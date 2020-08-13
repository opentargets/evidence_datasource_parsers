import pandas as pd
import numpy as np
import python_jsonschema_objects as pjo
import math
from datetime import datetime
from multiprocessing import Pool
import argparse
import gzip
import requests
import logging
import functools as fn

# Importing settings:
from settings import Config

class genetics_portal_evidence_generator(Config):
    """
    Based on the provided schema this class generates evidence string for the Genetics portal
    derived associations.

    Parameters inherited from Config class
    """

    # Source name:
    source_id = 'ot_genetics_portal'

    def __init__(self, schema_json, schema_version):
        """
        The init function loads the json schema into a builder object and a namespace:
        """

        # Initialize json builder based on the schema:
        self.builder = pjo.ObjectBuilder(schema_json)
        self.evidence_builder = self.builder.build_classes()
        self.schema_version = schema_version

    def get_evidence_string(self,row):

        # Variant level information:
        self.chromosome = row['chrom']
        self.position = row['pos']
        self.references_allele = row['ref']
        self.alternative_allele = row['alt']
        self.rs_id = row['rsid']
        self.consequence_code = row['most_severe_gene_csq']
        self.consequence_link = row['consequence_link']

        # Association level data:
        self.sample_size = row['sample_size']
        self.odds_ratio = row['odds_ratio']
        self.pval_mantissa = row['pval_mantissa']
        self.pval_exponent = row['pval_exponent']
        self.oddsr_ci_lower = row['oddsr_ci_lower']
        self.oddsr_ci_upper = row['oddsr_ci_upper']

        # Study level information:
        self.study_id = row['study_id']
        self.date_added = row['pub_date']
        self.author = row['pub_author']

        # Target level information:
        self.ensembl_gene_id = row['gene_id']

        # Disease level information:
        self.reported_trait = row['trait_reported']
        self.efo_id = row['efo']

        # Locus to gene score:
        self.l2g_score = row['y_proba_full_model']

        ##
        ## Applying a few checks before proceeding:
        ##

        # If the association has no EFO attached:
        if self.efo_id is None or self.efo_id == '':
            logging.warning('No EFO id for association row: {}'.format(row.name))
            return None

        # If because of some reason no l2g score is given:
        if self.l2g_score is None or self.l2g_score == '':
            logging.warning('No l2g score is available for association row: {}'.format(row.name))
            return None

        ##
        ## Process input:
        ##

        # Test if odds ratio is present:
        if math.isnan(row['odds_ratio']):
            self.odds_ratio = None
        else:
            self.odds_ratio = row['odds_ratio']

        # Test if odds ratio confidence interval is present:
        if math.isnan(self.oddsr_ci_lower):
            self.confidence_interval = None
        else:
            self.confidence_interval = '{}-{}'.format(self.oddsr_ci_lower,self.oddsr_ci_upper)

        # Test if Pubmed ID is present:
        try:
            self.pmid = row['pmid'].split(':')[1]
        except IndexError:
            self.pmid = row['pmid']
        except AttributeError:
            self.pmid = None

        date_stamp = datetime.strptime(self.date_added, '%Y-%m-%d')
        self.date_added = date_stamp.replace(microsecond=0).isoformat()
        self.year = date_stamp.year

        # Generating genetics portal variant id:
        self.genetics_portal_variant_id = '{}_{}_{}_{}'.format(self.chromosome,self.position,self.references_allele,self.alternative_allele)

        ##
        ## Prepare evidence:
        ##

        # Generate literature property:
        self.literature_prop = self.generate_literature_value()

        # Generate target property:
        target_prop = self.generate_target_value()

        # Generate unique association fields:
        unique_association_prop = self.generate_unique_fields()

        # Generate disease prop:
        disease_prop = self.generate_disease_value()

        # Generate features for the returned json:
        variant_prop = self.generate_variant_value()

        # Generating evidence property:
        evidence_value = self.generate_evidence_field()

        # Compiling properties into evidence:
        try:
            evidence = self.evidence_builder.Opentargets4(
                type='genetic_association',
                access_level='public',
                sourceID=self.source_id,
                variant=variant_prop,
                evidence=evidence_value,
                target=target_prop,
                disease=disease_prop,
                unique_association_fields=unique_association_prop,
                validated_against_schema_version=self.schema_version
            )
            return evidence.serialize()
        except:
            logging.warning('Evidence generation failed for row: {}'.format(row.name))
            raise

    def generate_literature_value(self):
        """
        This function generates values for publication link.

        return: dictionary
        """

        # Generate return value:
        if self.pmid:
            return {
                'references': [
                  {
                    'lit_id': '{}/{}'.format(self.LITERATURE_URL, self.pmid),
                    'author': self.author,
                    'year': self.year
                  }
                ]
            }
        else:
            return None

    def generate_target_value(self):
        """
        this function generates value for the target key of the gwas evidence

        return: dict
        """

        # Testing for type:
        if not isinstance(self.ensembl_gene_id, str):
            raise TypeError('[Error] Target could not be generated. The provided Ensembl ID ({}) is not a string.'.format(self.ensembl_gene_id))

        # Generate target object:
        return {
            'activity': '{}/predicted_damaging'.format(self.ACTIVITY_URL),
            'id': '{}/{}'.format(self.TARGET_URL, self.ensembl_gene_id),
            'target_type': '{}/gene_evidence'.format(self.TARGET_TYPE_URL)
        }

    def generate_disease_value(self):
        """
        This function generates value for the disease key of the gwas evidence

        return: dict
        """

        # Testing for efo_id type:
        if not isinstance(self.efo_id, str):
            raise TypeError('[Error] Disease could not be generated. The provided EFO ID ({}) is not a string.'.format(self.efo_id))

        # Testing type of reported trait:
        if self.reported_trait and not isinstance(self.reported_trait, str):
            raise TypeError('[Error] Disease could not be generated. The provided reported trait ({}) is not a string.'.format(self.reported_trait))

        # Generate
        return {
            'id': self.efo_id,
            'reported_trait': self.reported_trait
        }

    def generate_unique_fields(self):
        """
        Function to generate the unique terms. All fields are mandatory.

        Returns:
            dictionary with the following keys:
            - target: link to target on ensembl
            - disease_id: link to disease on efo
            - variant: link to variant on genetics portal
            - study: link to study page on genetics portal
        """

        return {
            'target': '{}/{}'.format(self.TARGET_URL, self.ensembl_gene_id),
            'disease_id': '{}/{}'.format(self.DISEASE_URL, self.efo_id),
            'variant': '{}/variant/{}'.format(self.GENETICS_PORTAL_URL, self.genetics_portal_variant_id),
            'study': '{}/study/{}'.format(self.GENETICS_PORTAL_URL, self.study_id)
        }

    @staticmethod
    def map_variant_type(ref_allele,alt_allele,consequence_base_url = 'http://purl.obolibrary.org/obo'):
        """
        This function maps variants to variant ontology
        based on the provided reference and alternative alleles.
        """

        # These are the available ontology terms to annotate variant types:
        term_mapper = {
            'SNP': 'SO_0000694',
            'deletion': 'SO_0000159',
            'insertion': 'SO_0000667'
        }

        # If any of the allele is missing, we cannot infer variant type:
        if (not ref_allele) or (not alt_allele):
            return OrderedDict({
                'type': None,
                'type_link': None
            })

        # Infer variant type:
        if (len(ref_allele) == 1) and (len(alt_allele) ==1):
            variant_type = 'SNP'
        elif (len(ref_allele) == 1) and (len(alt_allele) > 1):
            variant_type = 'deletion'
        elif (len(ref_allele) > 1) and (len(alt_allele) == 1):
            variant_type = 'insertion'
        else:
            variant_type = None

        if variant_type:
            return {
                'type': variant_type,
                'type_link': '{}/{}'.format(consequence_base_url,term_mapper[variant_type])
            }
        else:
            return {
                'type': None,
                'type_link': None
            }

    def generate_variant_value(self):
        return_data = {
            'id': self.genetics_portal_variant_id,
            'rs_id': self.rs_id,
            'source_link': '{}/variant/{}'.format(self.GENETICS_PORTAL_URL, self.genetics_portal_variant_id),
        }

        # Adding variant type annotation:
        return_data.update(self.map_variant_type(self.references_allele, self.alternative_allele, self.CONSEQUENCE_URL))

        return return_data

    def generate_evidence_field(self):
        """
        This function returns with the
        """

        return_value = {
            'variant2disease': self.generate_variant2disease(),
            'gene2variant': self.generate_gene2variant(),
        }

        return return_value

    def generate_variant2disease(self):

        # Generating p-value: if the exponent is too low, we apply a lower minimum.
        pval = float('{}e{}'.format(self.pval_mantissa,self.pval_exponent)) if self.pval_exponent > -300 else 1e-302

        # Generate resource score field:
        resource_score = {
            'type': 'pvalue',
            'method': {
                'description': 'pvalue for the snp to disease association'
            },
            'mantissa': int(round(self.pval_mantissa)),
            'exponent': int(self.pval_exponent),
            'value': pval
        }

        # Generate evidence code field:
        evidence_codes = [
            self.EVIDENCE_CODE_INFERENCE,
            self.EVIDENCE_CODE_EVIDENCE_TYPE
        ]

        return_value = dict(gwas_sample_size=self.sample_size, provenance_type=self.generate_provenance(),
                            is_associated=True,study_link='{}/study/{}'.format(self.GENETICS_PORTAL_URL, self.study_id),
                            resource_score=resource_score, evidence_codes=evidence_codes, date_asserted=self.date_added,
                            confidence_interval=self.confidence_interval, odds_ratio=self.odds_ratio)

        # Adding reported trait:
        return_value['reported_trait'] = self.reported_trait

        # The literature link is only added if the pubmed ID is present:
        return_value['unique_experiment_reference'] = '{}/{}'.format(self.LITERATURE_URL, self.pmid) if self.pmid else '{}/0000'.format(self.LITERATURE_URL)

        return return_value

    def generate_provenance(self):
        # Generate provenence type:
        provenence = {
            'expert': {
                'status': True,
                'statement': 'Primary submitter of the data'
            },
            'database': {
                'version': self.date_added,
                'id': self.source_id,
                'dbxref': {
                    'version': self.date_added,
                    'id': self.GENETICS_PORTAL_URL
                }
            }
        }

        if self.literature_prop:
            provenence['literature'] = self.literature_prop

        return provenence

    def generate_gene2variant(self):

        # Generate evidence code field:
        evidence_codes = [
            self.EVIDENCE_CODE_INFERENCE,
            self.EVIDENCE_CODE_SOURCE
        ]

        # Generate resource score field:
        resource_score = {
            'type': 'locus_to_gene_score',
            'method': {
                'description': 'Locus to gene score generated by OpenTargets Genetics portal',
            },
            'value': self.l2g_score
        }

        return_data = {
            'provenance_type': self.generate_provenance(),
            'is_associated': True,
            'date_asserted': self.date_added,
            'evidence_codes': evidence_codes,
            'resource_score': resource_score,
            'functional_consequence': self.consequence_link,
            'consequence_code': self.consequence_code
        }

        return return_data


def initialize_evidence_generation(df, schemaFile):
    """
    As it is not possible for multiprocessing.Pool to copy over states of objects,
    this function instantiates the evidence generator class and generates evidences.
    """

    # Fetching schema:
    json_schema = requests.get(schemaFile).json()

    # Initialize evidence builder object:
    evidence_builder = genetics_portal_evidence_generator(json_schema, '1.6.4')

    # Generate evidence for all rows of the dataframe:
    evidences = df.apply(evidence_builder.get_evidence_string, axis=1)

    return evidences


def parallelize_dataframe(df, schemaFile, n_cores=2):

    # Splitting dataframe as many chunks as many cores we have defined:
    df_split = np.array_split(df, n_cores)

    # Initialize pool object:
    pool = Pool(n_cores)

    # Create partial function:
    partial_function = fn.partial(initialize_evidence_generation, schemaFile = schemaFile)

    # Execute funciton on dataframe chunks and pool output:
    df = pd.concat(pool.map(partial_function, df_split))

    # Closing pool object:
    pool.close()
    pool.join()

    return df


def remove_duplicates(df):
    """
    This function reports and removes duplicated evidences to make sure
    all evidences are unique.
    """

    unique_fields = ['chrom', 'pos', 'ref', 'alt', 'gene_id', 'study_id', 'efo']

    # Are there any duplicated evidences:
    duplicates = df.loc[df.duplicated(unique_fields)]

    # Returning dataframe if no duplication was found:
    if len(duplicates) == 0:
        logging.info('Number of unique evidences: {}.'.format(len(df)))
        logging.info('No duplicated evidences were found.')
        return df 

    # Reporting duplication:
    logging.warning('Out of {} rows, the number of duplicated entries: {}'.format(len(df), len(duplicates)))

    # Report duplicates:
    logging.warning('The following rows are duplicated:')
    logging.warning(duplicates)
    logging.warning('Removing duplicates...')

    # Removing duplicates:
    df.drop_duplicates(unique_fields, keep=False, inplace=True)
    logging.info('Number of unique evidences: {}.'.format(len(df)))

    return df


def main():
    # Initialize logger:
    logging.basicConfig(
        filename='evidence_builder.log',
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    logging.info('Process started...')

    # Parsing input parameter:
    parser = argparse.ArgumentParser(description='This script generates Genetics portal sourced evidences.')

    # Database related input:
    parser.add_argument('--inputFile', help='Input parquet file with the table containing association details.', type=str, required=True)
    parser.add_argument('--schemaFile', help='OpenTargets JSON schema file (Draft-4 compatible!!).', type=str)
    parser.add_argument('--cores', help='Number of computing cores available for the evidence generation.', type=int, default=2)
    parser.add_argument('--outputFile', help='Name of the gzipped json output file.', type=str, default='output.json.gz')
    parser.add_argument('--sample', help='If provided this many randomly selected of evidences will be generated.', type=int, required=False)
    parser.add_argument('--threshold', help='If provided, evidences will be filtered for locus to gene score at the defined threshold.', type=float, required=False)
    args = parser.parse_args()

    # Parse input parameters:
    inputFile = args.inputFile
    cores = args.cores
    outputFile = args.outputFile
    sample = args.sample
    schemaFile = args.schemaFile
    threshold = args.threshold

    #  Opening input file parquet of tsv:
    if 'parquet' in inputFile:
        genetics_dataframe = pd.read_parquet(inputFile)
    elif '.tsv' in inputFile:
        genetics_dataframe = pd.read_csv(inputFile, sep='\t')

    logging.info('Loading data completed.')

    # Removing duplicates:
    genetics_dataframe = remove_duplicates(genetics_dataframe)

    # Applyting l2g score threshold if specified:
    if threshold:
        logging.info('Applying l2g score threshold {}'.format(threshold))
        rows = len(genetics_dataframe)
        genetics_dataframe = genetics_dataframe.loc[genetics_dataframe.y_proba_full_model > threshold]
        logging.info('Number of rows went from {} to {} after applying the threshold.'.format(rows, len(genetics_dataframe)))

    # If required, the dataframe is subset:
    if sample:
        genetics_dataframe = genetics_dataframe.sample(sample, random_state=12937)
        logging.info('Resample comleted. Number of rows: {}.'.format(len(genetics_dataframe)))

    # Evidences are generated in a parallel process spread across the defined number of cores:
    evidences = parallelize_dataframe(df=genetics_dataframe, schemaFile=schemaFile, n_cores=cores)
    evidences.dropna(inplace=True)

    logging.info('{} evidence strings have been successfully generated.'.format(len(evidences)))

    # Save gzipped json file:
    with gzip.open(outputFile, "wt") as f:
        evidences.apply(lambda x: f.write(str(x)+'\n'))

    logging.info('Evidence strings saved. Exiting.')


if __name__ == '__main__':
    
    main()
