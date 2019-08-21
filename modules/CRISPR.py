from settings import Config
from common.HGNCParser import GeneParser
import sys
import logging
import datetime
import opentargets.model.core as opentargets
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.association_score as association_score

__copyright__ = "Copyright 2014-2019, Open Targets"
__credits__   = ["Michaela Spitzer"]
__license__   = "Apache 2.0"
__version__   = "1.2.8"
__maintainer__= "Open Targets Data Team"
__email__     = ["data@opentargets.org"]
__status__    = "Production"

CRISPR_DATABASE_ID='CRISPR'
CRISPR_VERSION='2019.03'

class CRISPR:
    def __init__(self):
        self.evidence_strings = list()
        self.symbols = {}
        self.logger = logging.getLogger(__name__)

    def process_crispr(self):
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()

        self.symbols = gene_parser.genes
        self.build_evidence()
        self.write_evidence()

    def build_evidence(self, filename1=Config.CRISPR_FILENAME1, filename2=Config.CRISPR_FILENAME2):

        now = datetime.datetime.now()

        # *** Build evidence.provenance_type object ***
        provenance_type = evidence_core.BaseProvenance_Type(
            database=evidence_core.BaseDatabase(
                id=CRISPR_DATABASE_ID,
                version=CRISPR_VERSION
            )
        )
        error = provenance_type.validate(logging)

        if error > 0:
            self.logger.error(provenance_type.to_JSON(indentation=4))
            sys.exit(1)

        # Build dictionary with publication info (pmid, method description, score type, min_score, max_score)
        CRISPR_METHOD_MAP = {}
        with open(filename2, 'r') as crispr_descriptions:

            for line in crispr_descriptions:
                # pmid	method	score_type
                (efo_id, crispr_method) = tuple(line.rstrip().split('\t'))
                CRISPR_METHOD_MAP[efo_id]={'method':crispr_method}

        with open(filename1, 'r') as crispr_input:
            n = 0
            for line in crispr_input:
                n +=1
                if n>1:
                    # pmid	gene_set_name	target_id	disease_id	disease_name	score
                    # 28892060	Intestine Key Driver Genes	CD53	EFO_0003767	Inflammatory bowel disease	0.083369979
                    # 28892060	Intestine Key Driver Genes	RHOH	EFO_0003767	Inflammatory bowel disease	0.11673251
                    # 28892060	Intestine Key Driver Genes	DOCK2	EFO_0003767	Inflammatory bowel disease	0.122212311
                    # 28892060	Intestine Key Driver Genes	FGR	EFO_0003767	Inflammatory bowel disease	0.13290268
                    (pmid, gene_set_name, target_name, disease_id, disease_name, score) = tuple(line.rstrip().split('\t'))


                    # *** Build evidence.resource_score object ***
                    resource_score = association_score.Probability(
                        type="probability",
                        method=association_score.Method(
                            description = CRISPR_METHOD_MAP[disease_id]['method'],
                            reference = "http://europepmc.org/abstract/MED/{0}".format(pmid)
                        ),
                        value=float(score)/100
                    )

                    evidenceString = opentargets.Literature_Curated(
                        validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION,
                        access_level = "public",
                        type = "affected_pathway",
                        sourceID = "crispr"
                    )

                    # *** Build unique_association_field object ***
                    evidenceString.unique_association_fields = {
                        'pmid': pmid,
                        'gene_set': gene_set_name,
                        'gene_name': target_name,
                        'disease_id': disease_id
                    }

                    # *** Build target object ***
                    if target_name in self.symbols.values():
                        ensembl_gene_id = target_name

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
                        evidenceString.evidence = evidence_core.Literature_Curated(
                            date_asserted = now.isoformat(),
                            is_associated = True,
                            evidence_codes = ["http://purl.obolibrary.org/obo/ECO_0000053"],
                            provenance_type = provenance_type,
                            resource_score = resource_score
                        )

                        # Append current evidence string
                        self.evidence_strings.append(evidenceString)
                    else:
                        self.logger.error("%s is not found in Ensembl" % target_name)

            self.logger.error("%s evidence parsed"%(n-1))
            self.logger.error("%s evidence created"%len(self.evidence_strings))

    def write_evidence(self, filename=Config.CRISPR_EVIDENCE_FILENAME):
        self.logger.info("Writing CRISPR evidence strings")

        with open(filename, 'w') as crispr_output:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                self.logger.info(evidence_string.disease.id[0])
                error = evidence_string.validate(logging)
                if error == 0:
                    crispr_output.write(evidence_string.to_JSON(indentation=None)+"\n")
                else:
                    self.logger.error("REPORTING ERROR %i" %n)
                    self.logger.error(evidence_string.to_JSON(indentation=4))

def main():
    CRISPR().process_crispr()

if __name__ == "__main__":
    main()

