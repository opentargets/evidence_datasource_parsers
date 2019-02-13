from settings import Config
from common.HGNCParser import GeneParser
import sys
import logging
import datetime
import math
import opentargets.model.core as opentargets
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.association_score as association_score

# To check scores:
# cat sysbio-29-01-2019.json | jq -r '[.evidence.resource_score.value, .unique_association_fields.gene_set, .unique_association_fields.gene_name]| @csv'

__copyright__ = "Copyright 2014-2018, Open Targets"
__credits__   = ["Michaela Spitzer"]
__license__   = "Apache 2.0"
__version__   = "1.2.8"
__maintainer__= "Michaela Spitzer"
__email__     = ["data@opentargets.org"]
__status__    = "Production"

SYSBIO_DATABASE_ID='SYSBIO'
SYSBIO_VERSION='2018.11'

class SysBio:
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

        # Build dictionary with publication info (pmid, method description, score type, min_score, max_score)
        SYSBIO_PMID_MAP = {}
        with open(filename2, 'r') as sysbio_publications:

            for line in sysbio_publications:
                # pmid	method	score_type
                (sysbio_pmid, gene_set_name, sysbio_method, score_type, min_score, max_score) = tuple(line.rstrip().split('\t'))
                SYSBIO_PMID_MAP[gene_set_name]={'pmid':sysbio_pmid, 'method':sysbio_method, 'score_type':score_type, 'min_score':min_score, 'max_score':max_score}

        with open(filename1, 'r') as sysbio_input:
            n = 0
            for line in sysbio_input:
                n +=1
                if n>1:
                    # pmid	gene_set_name	target_id	disease_id	disease_name	score
                    # 28892060	Intestine Key Driver Genes	CD53	EFO_0003767	Inflammatory bowel disease	0.083369979
                    # 28892060	Intestine Key Driver Genes	RHOH	EFO_0003767	Inflammatory bowel disease	0.11673251
                    # 28892060	Intestine Key Driver Genes	DOCK2	EFO_0003767	Inflammatory bowel disease	0.122212311
                    # 28892060	Intestine Key Driver Genes	FGR	EFO_0003767	Inflammatory bowel disease	0.13290268
                    (pmid, gene_set_name, target_name, disease_id, disease_name, score) = tuple(line.rstrip().split('\t'))


                    # *** Build evidence.resource_score object ***
                    # Option 1: score is a 'p-value'   value = log10(pvalue) scaled to [0.5,1]
                    # Option 2: score is rank-based    value = score scaled to [0.5,1]
                    # Option 3: no score               value = 0.75
                    if SYSBIO_PMID_MAP[gene_set_name]['score_type'] == "p-value":
                        resource_score = association_score.Probability(
                            type="probability",
                            method=association_score.Method(
                                description = SYSBIO_PMID_MAP[gene_set_name]['method'],
                                reference = "http://europepmc.org/abstract/MED/{0}".format(pmid)
                            ),
                            value=SysBio.renormalize(math.log10(float(score)), [math.log10(float(SYSBIO_PMID_MAP[gene_set_name]['max_score'])), math.log10(float(SYSBIO_PMID_MAP[gene_set_name]['min_score']))], [0.5, 1])
                        )
                    elif SYSBIO_PMID_MAP[gene_set_name]['score_type'] == "rank":
                         resource_score = association_score.Probability(
                            type="probability",
                            method=association_score.Method(
                                description = SYSBIO_PMID_MAP[gene_set_name]['method'],
                                reference = "http://europepmc.org/abstract/MED/{0}".format(pmid)
                            ),
                            value=SysBio.renormalize(float(score), [float(SYSBIO_PMID_MAP[gene_set_name]['min_score']),float(SYSBIO_PMID_MAP[gene_set_name]['max_score'])], [0.5, 1])
                        )
                    else:
                        resource_score = association_score.Probability(
                            type="probability",
                            method=association_score.Method(
                                description=SYSBIO_PMID_MAP[gene_set_name]['method'],
                                reference="http://europepmc.org/abstract/MED/{0}".format(pmid)
                            ),
                            value=0.75
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
                        ##TODO issue with append, take only last item of the gene
                        self.evidence_strings.append(evidenceString)
                    else:
                        self.logger.error("%s is not found in Ensembl" % target_name)

            self.logger.error("%s evidence parsed"%(n-1))
            self.logger.error("%s evidence created"%len(self.evidence_strings))

    def write_evidence(self, filename=Config.SYSBIO_EVIDENCE_FILENAME):
        self.logger.info("Writing SysBio evidence strings")

        with open(filename, 'w') as sysbio_output:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                self.logger.info(evidence_string.disease.id[0])
                error = evidence_string.validate(logging)
                if error == 0:
                    sysbio_output.write(evidence_string.to_JSON(indentation=None)+"\n")
                else:
                    self.logger.error("REPORTING ERROR %i" %n)
                    self.logger.error(evidence_string.to_JSON(indentation=4))

    @staticmethod
    def renormalize(n, start_range, new_range):
        # apply the function f(x) to n using and old (start_range) and a new range
        # where f(x) = (dNewRange / dOldRange * (n - old_range_lower_bound)) + new_lower
        delta1 = start_range[1] - start_range[0]
        delta2 = new_range[1] - new_range[0]
        max_new_range = max(new_range)
        min_new_range = min(new_range)

        if delta1 or delta2:
            try:
                normalized = (delta2 * (n - start_range[0]) / delta1) + new_range[0]
            except ZeroDivisionError:
                normalized = new_range[0]
        else:
            normalized = n
        # The formula results in values slightly smaller and larger than the boundaries of the new range
        if normalized > max_new_range:
            return max_new_range
        elif normalized < min_new_range:
            return min_new_range

        return round(normalized,4)

def main():
    SysBio().process_sysbio()

if __name__ == "__main__":
    main()

