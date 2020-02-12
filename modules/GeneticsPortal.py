import pandas as pd
import numpy as np
import json 
import python_jsonschema_objects as pjo
import math
from datetime import datetime
from multiprocessing import Pool
import argparse
import gzip
import requests
import logging

class genetics_portal_evidence_generator(object):
    """
    Based on the provided schema this class generates evidence string for the Genetics portal
    derived associations.
    """
    
    ## 
    ## The following values shouldx be read from a config file....
    ##

    # Defining base URLs for links:
    activity_base_url = 'http://identifiers.org/cttv.activity'
    target_id_base_url = 'http://identifiers.org/ensembl'
    target_type_base_url = 'http://identifiers.org/cttv.target'
    literature_base_url = 'http://europepmc.org/abstract/MED'
    genetics_portal_base_url = 'https://genetics.opentargets.org'
    disease_base_url = 'http://www.ebi.ac.uk/efo'
    consequence_base_url = 'http://purl.obolibrary.org/obo'
                           
    # Evidence code: computational inference:
    evidence_code = 'http://purl.obolibrary.org/obo/ECO_0000362'

    # Source name:
    source_id = 'ot_genetics_portal'

    def __init__(self, schema_json, schema_version):
        """
        The init function loads the json schema into a builder object and a namespace;
        """

        # Initialize json builder based on the schema:
        self.builder = pjo.ObjectBuilder(schema_json)
        self.evidence_builder = self.builder.build_classes()     
        self.schema_version = schema_version

    def get_evidence_string(self,row):

        # Variant level information:
        self.chromosome = row['chrom'] # OK
        self.position = row['pos'] # OK
        self.references_allele = row['ref'] # OK
        self.alternative_allele = row['alt'] # OK
        self.rs_id = row['rsid'] # OK
        self.consequence_code = row['most_severe_gene_csq']
        self.consequence_link = row['consequence_link']
        
        # Association level data:
        self.sample_size = row['sample_size'] # OK
        self.odds_ratio = row['odds_ratio'] # OK
        self.pval_mantissa = row['pval_mantissa']
        self.pval_exponent = row['pval_exponent']
        self.oddsr_ci_lower = row['oddsr_ci_lower']
        self.oddsr_ci_upper = row['oddsr_ci_upper']
        
        # Study level information:
        self.study_id = row['study_id'] # OK
        self.date_added = row['pub_date'] # OK
        self.author = row['pub_author'] # OK
        
        # Target level information:
        self.ensembl_gene_id = row['gene_id'] # OK

        # Disease level information:
        self.reported_trait = row['trait_reported'] # OK
        self.efo_id = row['efo']# OK
        
        # Locus to gene score:
        self.l2g_score = row['y_proba_full_model'] # OK

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
        
        # evidence_value = self.generate_evidence_field()
        evidence_value = self.generate_evidence_field()
        
        # Building return object:
        try:
            evidence = self.evidence_builder.Opentargets4(
                type = 'genetic_association',
                access_level = 'public',
                sourceID = self.source_id,
                variant = variant_prop,
                evidence = evidence_value,
                target = target_prop,
                disease = disease_prop,
                unique_association_fields = unique_association_prop,
                validated_against_schema_version = self.schema_version
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
                'references' : [
                  {
                    'lit_id' : '{}/{}'.format(self.literature_base_url, self.pmid),
                    'author' : self.author,
                    'year' : self.year
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
            'activity' : '{}/predicted_damaging'.format(self.activity_base_url),
            'id' : '{}/{}'.format(self.target_id_base_url, self.ensembl_gene_id),
            'target_type' : '{}/gene_evidence'.format(self.target_type_base_url)        
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
            'id' : self.efo_id,
            'reported_trait' : self.reported_trait
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
            'target' : '{}/{}'.format(self.target_id_base_url, self.ensembl_gene_id),
            'disease_id' : '{}/{}'.format(self.disease_base_url, self.efo_id),
            'variant' : '{}/variant/{}'.format(self.genetics_portal_base_url, self.genetics_portal_variant_id),
            'study' : '{}/study/{}'.format(self.genetics_portal_base_url, self.study_id)
        }

    @staticmethod
    def map_variant_type(ref_allele,alt_allele):
        """
        This function maps variants to variant ontology
        based on the provided reference and alternative alleles.
        """
        consequence_base_url = 'http://purl.obolibrary.org/obo'
        
        # If any of the allele is missing, we cannot infer variant type:
        if (not ref_allele) or (not alt_allele):
            return OrderedDict({
                'type' : None,
                'type_link' : None
            })
        
        # These are the available terms:
        term_mapper = {
            'SNP' : 'SO_0000694',
            'deletion' : 'SO_0000159',
            'insertion' : 'SO_0000667'
        }
        
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
                'type' : variant_type,
                'type_link' : '{}/{}'.format(consequence_base_url,term_mapper[variant_type])
            }
        else:
            return {
                'type' : None,
                'type_link' : None
            }
 
    def generate_variant_value(self):
        return_data = {
            'id' : self.genetics_portal_variant_id,
            'rs_id' : self.rs_id,
            'source_link' : '{}/variant/{}'.format(self.genetics_portal_base_url, self.genetics_portal_variant_id),
        }
        
        # Adding variant type annotation:
        return_data.update(self.map_variant_type(self.references_allele, self.alternative_allele))

        return return_data
    
           
    def generate_evidence_field(self):
        """
        This function returns with the 
        """
        
        return_value = {
            'variant2disease' : self.generate_variant2disease(),
            'gene2variant' : self.generate_gene2variant(),
        }
        
        return return_value
        
    def generate_variant2disease(self):
        
        # Generating p-value: if the exponent is too low, we apply a lower minimum.
        pval = float('{}e{}'.format(self.pval_mantissa,self.pval_exponent)) if self.pval_exponent > -300 else 1e-302
              
        # Generate resource score field:
        resource_score = {
            'type' : 'pvalue',
            'method' : {
                'description' : 'pvaluefor the snp to disease association'
            },
            'mantissa' : int(round(self.pval_mantissa)),
            'exponent' : int(self.pval_exponent),
            'value': pval
        }
        
        # Generate evidence code field:
        evidence_codes = [
            self.evidence_code,
            self.genetics_portal_base_url
        ]
        
        return_value = {
            'gwas_sample_size' : self.sample_size,
            'provenance_type' : self.generate_provenance(),
            'is_associated' : True,
            'study_link' : '{}/study/{}'.format(self.genetics_portal_base_url, self.study_id),
            'resource_score' : resource_score,
            'evidence_codes' : evidence_codes,
            'date_asserted' : self.date_added,
            'confidence_interval' : self.confidence_interval,
            'odds_ratio' : self.odds_ratio
        }
        
        # Adding reported trait:
        return_value['reported_trait'] = self.reported_trait
        
        # The literature link is only added if the pubmed ID is present:
        return_value['unique_experiment_reference'] = '{}/{}'.format(self.literature_base_url, self.pmid) if self.pmid else '{}/0000'.format(self.literature_base_url)

        return return_value
    
    def generate_provenance(self):
        # ("'['date_asserted', 'is_associated', 'provenance_type', 'resource_score']' are required attributes for evidence \nwhile setting 'evidence' in opentargets_3", 'occurred at index 0')
        # Generate provenence type:
        provenence = {
            'expert' : {
                'status' : True,
                'statement' : 'Primary submitter of the data'
            },
            'database' : {
                'version' : self.date_added,
                'id' : self.source_id,
                'dbxref' : {
                    'version' : self.date_added,
                    'id' : self.genetics_portal_base_url
                }
            }
        }
        
        if self.literature_prop:
            provenence['literature'] = self.literature_prop
        
        return provenence
    
    def generate_gene2variant(self):
        
        # Generate evidence code field:
        evidence_codes = [
            self.evidence_code,
            self.genetics_portal_base_url
        ]
        
        # Generate resource score field:
        resource_score = {
            'type' : 'locus_to_gene_score',
            'method' : {
                'description' : 'Locus to gene score generated by OpenTargets Genetics portal',
            },
            'value' : self.l2g_score
        }
        
        return_data = {
            'provenance_type' : self.generate_provenance(),
            'is_associated' : True,
            'date_asserted' : self.date_added,
            'evidence_codes' : evidence_codes,
            'resource_score' : resource_score,
            'functional_consequence' : self.consequence_link,
            'consequence_code' : self.consequence_code
        }
        
        return return_data
            

def initialize_evidence_generation(df):
    """
    As it is not possible for multiprocessing.Pool to copy over states of objects,
    this function instantiates the evidence generator class and generates evidences.
    """
    json_schema = requests.get('https://raw.githubusercontent.com/opentargets/json_schema/Draft-4_compatible/opentargets.json').json()
    evidence_builder = genetics_portal_evidence_generator(json_schema, '1.6.4')
    evidences = df.apply(evidence_builder.get_evidence_string, axis = 1)

    return evidences
    

def parallelize_dataframe(df, func, n_cores=2):
    """
    Splitting dataframe to distribute across cores:
    """
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)

    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
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
    parser.add_argument('--inputFile', help='Input parquet file with the table containing association details.', type=str, required = True)
    parser.add_argument('--schemaFile', help='OpenTargets JSON schema file (Draft-4 compatible!!).', type=str)
    parser.add_argument('--cores', help='tsv file to save output.', type = int, default = 2)
    parser.add_argument('--outputFile', help = 'Filename for the json file output.', type = str, default = 'output.json.gz')
    args = parser.parse_args()

    # Parse input parameters:
    json_schema = args.schemaFile
    inputFile = args.inputFile
    cores = args.cores
    outputFile = args.outputFile

    #  Opening input file:
    if 'parquet' in inputFile:
        genetics_dataframe = pd.read_parquet(inputFile)
    elif '.tsv' in inputFile:
        genetics_dataframe = pd.read_csv(inputFile, sep = '\t')

    logging.info('Loading data completed.')

    evidences = parallelize_dataframe(df = genetics_dataframe.head(), func = initialize_evidence_generation, n_cores = cores)

    logging.info('Evidence strings generated.')

    # Save file:
    with gzip.open(outputFile, "wt") as f:
        evidences.apply(lambda x: f.write(str(x)+'\n') if x is not None else 'cica')
        
    logging.info('Evidence strings saved. Exiting.')


if __name__ == '__main__':

    main()
