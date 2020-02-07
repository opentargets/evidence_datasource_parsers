import sys
import logging
import datetime
from collections import OrderedDict
import pandas as pd
from common.HGNCParser import GeneParser
from settings import Config, file_or_resource
import opentargets.model.core as opentargets
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.linkout as evidence_linkout
import opentargets.model.evidence.association_score as association_score
import opentargets.model.evidence.mutation as evidence_mutation

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

    def __init__(self):

        self.evidence_strings=list()
        self.genes=None
        self.symbols={}
        self.cancer_acronym_to_efo_map = OrderedDict()
        self.cohort_info= OrderedDict()
        self.logger=logging.getLogger(__name__)

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
                self.cancer_acronym_to_efo_map[cancer_acronym].append({'efo_uri': efo_uri, 'efo_label': efo_label, 'intogen_full_name': intogen_cancer_full_name})

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
        now=datetime.datetime(2019, 11, 12, 0, 0)
        provenance_type = evidence_core.BaseProvenance_Type(
            database=evidence_core.BaseDatabase(
                id="IntOGen Cancer Drivers Database",
                version='2019.11',
                dbxref=evidence_core.BaseDbxref(url="https://www.intogen.org/search", id="IntOGen Cancer Drivers Database", version="2019.11")),
            literature=evidence_core.BaseLiterature(
                references=[evidence_core.Single_Lit_Reference(lit_id="http://europepmc.org/abstract/MED/25759023")]
            )
        )
        error=provenance_type.validate(self.logger)
        if error > 0:
            self.logger.error(provenance_type.to_JSON(indentation=4))
            print(provenance_type.to_JSON(indentation=4))
            sys.exit(1)

        intogen_driver_genes_df = pd.read_csv(driver_genes_filename, sep='\t', header=0, keep_default_na=False)

        n=0
        for index, line in intogen_driver_genes_df.iterrows():
            n += 1
            if n > 1:
                # One IntOGen cancer type may be mapped to multiple disease, e.g. PCPG - Pheochromocytoma and paraganglioma, generate evidence strings for all of them
                for disease in self.cancer_acronym_to_efo_map[line['CANCER_TYPE']]:
                    # *** General properties ***
                    evidenceString=opentargets.Literature_Curated()
                    evidenceString.validated_against_schema_version=Config.VALIDATED_AGAINST_SCHEMA_VERSION
                    evidenceString.access_level="public"
                    evidenceString.type="somatic_mutation"
                    evidenceString.sourceID="intogen"

                    # *** Target information ***
                    # get the Ensembl gene id from the symbol (mapping from 2014 won't work)
                    # TODO: Check if this is still necessary with Nov 19 file
                    if line['SYMBOL'] in INTOGEN_SYMBOL_MAPPING:
                        line['SYMBOL']=INTOGEN_SYMBOL_MAPPING[line['SYMBOL']]

                    ensembl_gene_id = self.genes.get(line['SYMBOL'])
                    if not ensembl_gene_id:
                        self.logger.error("%s is not found in Ensembl" % line['SYMBOL'])
                        continue

                    evidenceString.target = bioentity.Target(
                        id="http://identifiers.org/ensembl/{0}".format(ensembl_gene_id),
                        target_name=line['SYMBOL'],
                        activity=INTOGEN_ROLE_MAP[line['ROLE']],
                        target_type='http://identifiers.org/cttv.target/gene_evidence'
                    )
                    # *** Disease information ***
                    evidenceString.disease = bioentity.Disease(
                        id=disease['efo_uri'],
                        name=disease['efo_label'],
                        source_name=disease['intogen_full_name']
                    )
                    # *** Evidence ***
                    linkout = evidence_linkout.Linkout(
                        url='https://www.intogen.org/search?gene=%s&cohort=%s'%(line['SYMBOL'], line['COHORT']),
                        nice_name='IntOGen -  %s gene cancer mutations in %s (%s)'%(line['SYMBOL'], disease['efo_label'], line['CANCER_TYPE'])
                    )
                    # inheritance_pattern - 'gain_of_function' = 'Dominant', 'loss_of_function' = 'Recessive'
                    inheritance_pattern='unknown'
                    if line['ROLE'] == 'Act':
                        inheritance_pattern='dominant'
                    elif line['ROLE'] == 'LoF':
                        inheritance_pattern='recessive'
                    # known_mutations
                    mutation = evidence_mutation.Mutation(
                        functional_consequence='http://purl.obolibrary.org/obo/SO_0001564',
                        preferred_name='gene_variant',
                        inheritance_pattern=inheritance_pattern,
                        number_samples_tested= self.cohort_info[line['COHORT']]['number_samples'],
                        number_mutated_samples= line['SAMPLES']
                    )
                    # resource_score
                    resource_score = association_score.Pvalue(
                        type="pvalue",
                        method=association_score.Method(
                            description="IntOGen Driver identification methods as described in Rubio-Perez, C., Tamborero, D., Schroeder, MP., Antolin, AA., Deu-Pons,J., Perez-Llamas, C., Mestres, J., Gonzalez-Perez, A., Lopez-Bigas, N. In silico prescription of anticancer drugs to cohorts of 28 tumor types reveals novel targeting opportunities. Cancer Cell 27 (2015), pp. 382-396",
                            reference="http://europepmc.org/abstract/MED/25759023",
                            url="https://www.intogen.org/about"),
                        value=line['QVALUE_COMBINATION']
                    )

                    evidenceString.evidence = evidence_core.Literature_Curated(
                        date_asserted=now.isoformat(),
                        is_associated=True,
                        evidence_codes=["http://purl.obolibrary.org/obo/ECO_0000053"],
                        provenance_type=provenance_type,
                        resource_score=resource_score,
                        urls=[linkout],
                        known_mutations=[mutation]
                    )

                    # *** unique_association_fields ***
                    # (target_id & disease_id are sufficient)
                    evidenceString.unique_association_fields = {
                        'target_id': evidenceString.target.id,
                        'disease_id': evidenceString.disease.id,
                        'cohort_id': line['COHORT'],
                        'cohort_short_name': self.cohort_info[line['COHORT']]['short_cohort_name'],
                        'cohort_description': self.cohort_info[line['COHORT']]['cohort_description'],
                        'methods': line['METHODS']
                    }

                    error=evidenceString.validate(logging)
                    if error > 0:
                        self.logger.error(evidenceString.to_JSON())
                        sys.exit(1)

                    self.evidence_strings.append(evidenceString)

            self.logger.info("%s evidence parsed"%(n-1))
            self.logger.info("%s evidence created"%len(self.evidence_strings))
            print("%s evidence created"%len(self.evidence_strings))


    def write_evidence_strings(self, filename=Config.INTOGEN_EVIDENCE_FILENAME):
        self.logger.info("Writing IntOGen evidence strings")
        with open(filename, 'w') as tp_file:
            n=0
            for evidence_string in self.evidence_strings:
                n+=1
                self.logger.info(evidence_string.disease.id[0])
                error=evidence_string.validate(logging)
                if error == 0:
                    tp_file.write(evidence_string.to_JSON(indentation=None)+"\n")
                else:
                    self.logger.error("REPORTING ERROR %i" % n)
                    self.logger.error(evidence_string.to_JSON(indentation=4))
        tp_file.close()

def main():
    into=IntOGen()
    into.process_intogen()

if __name__ == "__main__":
    main()
