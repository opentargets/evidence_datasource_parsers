from settings import Config
from common.HGNCParser import GeneParser
import sys
import logging
import datetime
import opentargets.model.core as opentargets
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.association_score as association_score

__copyright__ = "Copyright 2014-2018, Open Targets"
__credits__   = ["Michaela Spitzer"]
__license__   = "Apache 2.0"
__version__   = "1.2.8"
__maintainer__= "Michaela Spitzer"
__email__     = ["data@opentargets.org"]
__status__    = "Production"

SYSBIO_DATABASE_ID='SYSBIO'
SYSBIO_VERSION='2018.11'

class SYSBIO:
    def __init__(self):
        self.evidence_strings = list()
        self.symbols = {}
        self.logger = logging.getLogger(__name__)

    def process_sysbio(self):
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()

        self.symbols = gene_parser.genes
        self.build_evidence()
        self.write_evidence()

    def build_evidence(self, filename1=Config.SYSBIO_FILENAME1, filename2=Config.SYSBIO_FILENAME2):

        now = datetime.datetime.now()

        # *** Build evidence.provenance_type object ***
        provenance_type = evidence_core.BaseProvenance_Type(
            database=evidence_core.BaseDatabase(
                id=SYSBIO_DATABASE_ID,
                version=SYSBIO_VERSION
            )
        )
        error = provenance_type.validate(logging)

        if error > 0:
            self.logger.error(provenance_type.to_JSON(indentation=4))
            sys.exit(1)

        # Build dictionary with publication info (pmid, method description, score type)
        # TODO: Change to (pmid, gene set name, method description, score type) so that scores can be categorised per gene set
        SYSBIO_PMID_MAP = {}
        with open(filename2, 'r') as sysbio_publications:

            for line in sysbio_publications:
                # pmid	method	score_type
                (sysbio_pmid, sysbio_method, score_type) = tuple(line.rstrip().split('\t'))
                SYSBIO_PMID_MAP[sysbio_pmid]={'method':sysbio_method, 'score_type':score_type}

        with open(filename1, 'r') as sysbio_input:
            n = 0
            for line in sysbio_input:
                n +=1
                if n>1:
                    '''
                        pmid	gene_set_name	target_id	disease_id	disease_name	score
                        28892060	Intestine Key Driver Genes	CD53	EFO_0003767	Inflammatory bowel disease	0.083369979
                        28892060	Intestine Key Driver Genes	RHOH	EFO_0003767	Inflammatory bowel disease	0.11673251
                        28892060	Intestine Key Driver Genes	DOCK2	EFO_0003767	Inflammatory bowel disease	0.122212311
                        28892060	Intestine Key Driver Genes	FGR	EFO_0003767	Inflammatory bowel disease	0.13290268

                    '''
                    (pmid, gene_set_name, target_name, disease_id, disease_name, score) = tuple(line.rstrip().split('\t'))

                    # *** Build evidence.resource_score object ***
                    resource_score = association_score.Pvalue(
                        type="pvalue",
                        method=association_score.Method(
                            description = SYSBIO_PMID_MAP[pmid]['method'],
                            reference = "http://europepmc.org/abstract/MED/{0}".format(pmid)
                        ),
                        # TODO: add scores as p-value or a rank-based score (to be decided)
                        # Default score is 0.5 if no scores are provided
                        value=0.5
                    )

                    evidenceString = opentargets.Literature_Curated()
                    evidenceString.validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION
                    evidenceString.access_level = "public"
                    evidenceString.type = "affected_pathway"
                    evidenceString.sourceID = "sysbio"

                    # *** Build unique_association_field object ***
                    evidenceString.unique_association_fields = {}
                    evidenceString.unique_association_fields['pmid'] = pmid
                    evidenceString.unique_association_fields['gene_set'] = gene_set_name
                    evidenceString.unique_association_fields['gene_name'] = target_name

                    # *** Build target object ***
                    if target_name in self.symbols:
                        ensembl_gene_id = self.symbols[target_name]

                        evidenceString.target = bioentity.Target(
                            id="http://identifiers.org/ensembl/{0}".format(ensembl_gene_id),
                            target_name=target_name,
                            # TODO activity is a required field in target object, currently set as unknown
                            activity="http://identifiers.org/cttv.activity/unknown",
                            target_type='http://identifiers.org/cttv.target/gene_evidence'
                        )
                        # *** Build disease object ***
                        evidenceString.disease = bioentity.Disease(
                            id=disease_id,
                            name=disease_name
                        )
                        # *** Build evidence oject ***
                        evidenceString.evidence = evidence_core.Literature_Curated()
                        evidenceString.evidence.date_asserted = now.isoformat()
                        evidenceString.evidence.is_associated = True
                        evidenceString.evidence.evidence_codes = ["http://purl.obolibrary.org/obo/ECO_0000053"]
                        evidenceString.evidence.provenance_type = provenance_type
                        evidenceString.evidence.resource_score = resource_score


                        error = evidenceString.validate(logging)

                        #print(evidenceString.to_JSON(indentation=None))
                        ##TODO issue with append, take only last item of the gene
                        self.evidence_strings.append(evidenceString)

                    else:
                        self.logger.error("%s is not found in Ensembl" % target_name)

            self.logger.error("%s evidence parsed"%(n-1))
            self.logger.error("%s evidence created"%len(self.evidence_strings))

        sysbio_input.close()

    def write_evidence(self, filename=Config.SYSBIO_EVIDENCE_FILENAME):
        self.logger.info("Writing SYSBIO evidence strings")

        with open(filename, 'w') as sysbio_output:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                self.logger.info(evidence_string.disease.id[0])
                # get max_phase_for_all_diseases
                error = evidence_string.validate(logging)
                if error == 0:
                    sysbio_output.write(evidence_string.to_JSON(indentation=None)+"\n")
                else:
                    self.logger.error("REPORTING ERROR %i" %n)
                    self.logger.error(evidence_string.to_JSON(indentation=4))
            sysbio_output.close()


def main():
    SYSBIO().process_sysbio()

if __name__ == "__main__":
    main()

