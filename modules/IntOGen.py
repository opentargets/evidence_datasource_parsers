import sys
import logging
import datetime
from collections import OrderedDict
import requests
import pandas as pd
import python_jsonschema_objects as pjo
from common.HGNCParser import GeneParser
from settings import Config, file_or_resource

__copyright__ ="Copyright 2014-2020, Open Targets"
__credits__   =["Gautier Koscielny", "David Tamborero", "ChuangKee Ong", "Michaela Spitzer", "Asier Gonzalez"]
__license__   ="Apache 2.0"
__version__   ="1.3.0"
__maintainer__="Open Targets Data Team"
__email__     =["data@opentargets.org"]
__status__    ="Production"

INTOGEN_ROLE_MAP={
    'Act': 'http://identifiers.org/cttv.activity/gain_of_function',
    'LoF': 'http://identifiers.org/cttv.activity/loss_of_function',
    'Ambiguous': 'http://identifiers.org/cttv.activity/unknown',
    'ambiguous': 'http://identifiers.org/cttv.activity/unknown',
    '': 'http://identifiers.org/cttv.activity/unknown'
}

INTOGEN_SYMBOL_MAPPING={
    'C15orf55': 'NUTM1',
    'CSDA': 'YBX3',
    'EIF2C3': 'AGO3',
    'ERBB2IP': 'ERBIN',
    'FAM123B': 'AMER1',
    'HNRPDL': 'HNRNPDL',
    'MLL': 'KMT2A',
    'MLL2': 'KMT2D',
    'MLL3': 'KMT2C',
    'RQCD1': 'CNOT9'
}

class IntOGen():

    def __init__(self, schema_json=Config.OT_JSON_SCHEMA, schema_version=Config.VALIDATED_AGAINST_SCHEMA_VERSION):

        self.evidence_strings=list()
        self.genes=None
        self.symbols={}
        self.cancer_acronym_to_efo_map = OrderedDict()
        self.cohort_info= OrderedDict()
        self.logger=logging.getLogger(__name__)

        # Initialize json builder based on the schema:
        json_schema = requests.get(schema_json).json()
        self.builder = pjo.ObjectBuilder(json_schema)
        self.evidence_builder = self.builder.build_classes()
        self.schema_version = schema_version

    def process_intogen(self, infile=Config.INTOGEN_DRIVER_GENES_FILENAME, outfile=Config.INTOGEN_EVIDENCE_FILENAME):

        gene_parser=GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes=gene_parser.genes
        self.get_cancer_acronym_to_efo_mappings()
        self.read_intogen_cohorts()
        self.read_intogen(driver_genes_filename=infile)
        self.write_evidence_strings(filename=outfile)

    def get_cancer_acronym_to_efo_mappings(self, cancer_mapping_filename=Config.INTOGEN_CANCER2EFO_MAPPING_FILENAME):
        """Parse intogen_cancer2EFO_mapping.tsv and create a dictionary to map IntOGen 3-letter cancer acronyms to EFO"""

        with open(cancer_mapping_filename, 'r') as cancer2efo:
            for line in cancer2efo:
                line = line.strip('\n')
                (intogen_cancer_full_name, cancer_acronym, efo_uri, efo_label) = line.split("\t")
                if cancer_acronym not in self.cancer_acronym_to_efo_map:
                    self.cancer_acronym_to_efo_map[cancer_acronym] = []
                self.cancer_acronym_to_efo_map[cancer_acronym].append({'efo_uri': efo_uri.strip(), 'efo_label': efo_label.strip(), 'intogen_full_name': intogen_cancer_full_name})

    def read_intogen_cohorts(self, cohorts_filename=Config.INTOGEN_COHORTS):
        """Parse IntOGen cohorts file and create a dictionary with necessary information"""

        cohorts_df = pd.read_csv(cohorts_filename, sep='\t', header=0)

        for index, row in cohorts_df.iterrows():
            if row['COHORT'] not in self.cohort_info:
                cohort_info_dict = {}
                cohort_info_dict['short_cohort_name'] = row['WEB_SHORT_COHORT_NAME']
                cohort_info_dict['cohort_description'] = row['WEB_LONG_COHORT_NAME']
                cohort_info_dict['number_samples'] = row['SAMPLES']
                self.cohort_info[row['COHORT']] = cohort_info_dict
            else:
                self.logger.error("Cohort %s appears multiple times in %s"%(row['COHORT'], cohorts_filename))
                sys.exit(1)


    def read_intogen(self, driver_genes_filename=Config.INTOGEN_DRIVER_GENES_FILENAME):

        # the database was created in 2014
        #now=datetime.datetime.now()
        now=datetime.datetime(2020, 2, 1, 0, 0)
        provenance_type =  {
            'database'  : {
                'id' : "IntOGen Cancer Drivers Database",
                'version' : '2020.02',
                'dbxref' : {
                    'url' : "https://www.intogen.org/search",
                    'id' : "IntOGen Cancer Drivers Database",
                    'version' : "2020.02"
                }
            },
           'literature' : {
               'references' : [
                   {
                       'lit_id': "http://europepmc.org/abstract/MED/25759023"
                   }
               ]
           }
        }

        intogen_driver_genes_df = pd.read_csv(driver_genes_filename, sep='\t', header=0, keep_default_na=False)

        n=0
        for index, line in intogen_driver_genes_df.iterrows():
            n += 1
            # One IntOGen cancer type may be mapped to multiple disease, e.g. PCPG - Pheochromocytoma and paraganglioma, generate evidence strings for all of them
            for disease in self.cancer_acronym_to_efo_map[line['CANCER_TYPE']]:

                validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION
                access_level = "public"
                type = "somatic_mutation"
                sourceID = "intogen"

                # *** Target information ***
                # get the Ensembl gene id from the symbol (mapping from 2014 won't work)
                # TODO: Check if this is still necessary with Nov 19 file
                if line['SYMBOL'] in INTOGEN_SYMBOL_MAPPING:
                    line['SYMBOL']=INTOGEN_SYMBOL_MAPPING[line['SYMBOL']]

                ensembl_gene_id = self.genes.get(line['SYMBOL'])
                if not ensembl_gene_id:
                    self.logger.error("%s is not found in Ensembl" % line['SYMBOL'])
                    continue

                target_info = {
                    'id' : "http://identifiers.org/ensembl/{0}".format(ensembl_gene_id),
                    'target_name' : line['SYMBOL'],
                    'activity' : INTOGEN_ROLE_MAP[line['ROLE']],
                    'target_type' : 'http://identifiers.org/cttv.target/gene_evidence'
                }
                # *** Disease information ***
                disease_info = {
                    'id' : disease['efo_uri'],
                    'name' : disease['efo_label'],
                    'source_name' : disease['intogen_full_name']
                }
                # *** Evidence ***
                linkout = [
                    {
                        'url' : 'https://www.intogen.org/search?gene=%s&cohort=%s' % (line['SYMBOL'], line['COHORT']),
                        'nice_name' : 'IntOGen -  %s gene cancer mutations in %s (%s)' % (line['SYMBOL'], disease['efo_label'], line['CANCER_TYPE'])
                    }
                ]

                # inheritance_pattern - 'gain_of_function' = 'Dominant', 'loss_of_function' = 'Recessive'
                inheritance_pattern='unknown'
                if line['ROLE'] == 'Act':
                    inheritance_pattern='dominant'
                elif line['ROLE'] == 'LoF':
                    inheritance_pattern='recessive'
                # known_mutations
                mutation = [
                    {
                        'functional_consequence': 'http://purl.obolibrary.org/obo/SO_0001564',
                        'preferred_name': 'gene_variant',
                        'inheritance_pattern': inheritance_pattern,
                        'number_samples_tested': self.cohort_info[line['COHORT']]['number_samples'],
                        'number_mutated_samples': line['SAMPLES']
                    }
                ]


                # resource_score
                resource_score = {
                    'type' : "pvalue",
                    'method' : {
                        'description': "IntOGen Driver identification methods as described in Rubio-Perez, C., Tamborero, D., Schroeder, MP., Antolin, AA., Deu-Pons,J., Perez-Llamas, C., Mestres, J., Gonzalez-Perez, A., Lopez-Bigas, N. In silico prescription of anticancer drugs to cohorts of 28 tumor types reveals novel targeting opportunities. Cancer Cell 27 (2015), pp. 382-396",
                        'reference': "http://europepmc.org/abstract/MED/25759023",
                        'url': "https://www.intogen.org/about"
                    },
                    'value' : line['QVALUE_COMBINATION']

                }

                # Extract driver gene methods and store in list
                significant_driver_methods = line['METHODS'].split(",")

                # Dictionary with cohort info
                cohort = {
                    'cohort_id': line['COHORT'],
                    'cohort_short_name': self.cohort_info[line['COHORT']]['short_cohort_name'],
                    'cohort_description': self.cohort_info[line['COHORT']]['cohort_description']
                }

                evidence_info = {
                    'date_asserted' : now.isoformat(),
                    'is_associated' : True,
                    'evidence_codes' : ["http://purl.obolibrary.org/obo/ECO_0000053"],
                    'provenance_type' : provenance_type,
                    'resource_score' : resource_score,
                    'urls' : linkout,
                    'known_mutations' : mutation,
                    'significant_driver_methods' : significant_driver_methods,
                    'cohort' : cohort
                }

                # *** unique_association_fields ***
                unique_association_fields = {
                    'target_id': target_info['id'],
                    'disease_id': disease_info['id'],
                    'cohort_short_name': self.cohort_info[line['COHORT']]['short_cohort_name'],
                }

                try:
                    evidence = self.evidence_builder.Opentargets(
                         type = type,
                         access_level = access_level,
                         sourceID = sourceID,
                         evidence = evidence_info,
                         target = target_info,
                         disease = disease_info,
                         unique_association_fields = unique_association_fields,
                         validated_against_schema_version = validated_against_schema_version
                    )
                    self.evidence_strings.append(evidence)
                except:
                    self.logger.warning('Evidence generation failed for row: {}'.format(n))
                    raise



        self.logger.info("%s evidence parsed"%(n-1))
        self.logger.info("%s evidence created"%len(self.evidence_strings))
        print("%s evidence created"%len(self.evidence_strings))


    def write_evidence_strings(self, filename=Config.INTOGEN_EVIDENCE_FILENAME):
        self.logger.info("Writing IntOGen evidence strings")
        with open(filename, 'w') as tp_file:
            n=0
            for evidence_string in self.evidence_strings:
                n+=1
                self.logger.info(evidence_string.disease['id'])
                tp_file.write(evidence_string.serialize()+"\n")
        tp_file.close()

def main():
    into=IntOGen()
    into.process_intogen()

if __name__ == "__main__":
    main()
