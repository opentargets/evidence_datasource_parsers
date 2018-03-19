from settings import Config
from common.HGNCParser import GeneParser
import sys
import logging
import datetime
import opentargets.model.core as opentargets
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.linkout as evidence_linkout
import opentargets.model.evidence.association_score as association_score

__copyright__ = "Copyright 2014-2018, Open Targets"
__credits__   = ["ChuangKee Ong", "Luz Garcia Alonso"]
__license__   = "Apache 2.0"
__version__   = "1.2.8"
__maintainer__= "ChuangKee Ong"
__email__     = ["data@opentargets.org"]
__status__    = "Production"

''' TCGA -> EFO mapping '''
TUMOR_TYPE_EFO_MAP = {
    'BLCA': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000292', 'label': 'bladder carcinoma'},
    'BRCA': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000305', 'label': 'breast carcinoma'},
    'HNSC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000181', 'label': 'head and neck squamous cell carcinoma'},
    'KIRC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000349', 'label': 'clear cell renal carcinoma'},
    'LIHC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000182', 'label': 'hepatocellular carcinoma'},
    'LUAD': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000571', 'label': 'lung adenocarcinoma'},
    'LUSC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000708', 'label': 'squamous cell lung carcinoma'},
    'PRAD': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000673', 'label': 'prostate adenocarcinoma'},
    'STAD': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000503', 'label': 'stomach adenocarcinoma'},
    'THCA': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0002892', 'label': 'thyroid carcinoma'},
    'UCEC': {'uri': 'http://www.ebi.ac.uk/efo/EFO_1000233', 'label': 'endometrial endometrioid adenocarcinoma'},
    ##TODO tumor uri & label need update
    'KICH': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000000', 'label': 'kidney chromophobe'},
    'KIRP': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000000', 'label': 'kidney renal papillary cell carcinoma'},
    'COREAD': {'uri': 'http://www.ebi.ac.uk/efo/EFO_0000000', 'label': 'colorectoral adenoma'}
}

''' TCGA acronyms '''
TUMOR_TYPE_MAP = {
    'BLCA': 'bladder carcinoma',
    'BRCA': 'breast carcinoma',
    'HNSC': 'head and neck squamous cell carcinoma',
    'KIRC': 'clear cell renal carcinoma',
    'LIHC': 'hepatocellular carcinoma',
    'LUAD': 'lung adenocarcinoma',
    'LUSC': 'squamous cell lung carcinoma',
    'PRAD': 'prostate adenocarcinoma',
    'STAD': 'stomach adenocarcinoma',
    'THCA': 'thyroid carcinoma',
    'UCEC': 'endometrial endometrioid adenocarcinoma',
    'KICH': 'kidney chromophobe',
    'KIRP': 'kidney renal papillary cell carcinoma',
    #TODO tumor description need update
    'COREAD': ''
}

''' Pathway -> Perturbed Targets
    https://drive.google.com/drive/folders/1L5Y_umEZiccWJnXiiaYMNKUKYTjnp3ZU
'''
PATHWAY_TARGET_MAP = {
    'Androgen' : ['AR'],
    'EGFR'     : ['EGFR'],
    'Estrogen' : ['ESR1'],
    'Hypoxia'  : ['HIF1A'],
    'JAK.STAT' : ['JAK1', 'JAK2', 'STAT1', 'STAT2'],
    'MAPK'     : ['MAPK2K1', 'MAP2K2', 'RAF1'],
    'NFkB'     : ['TLR4', 'NKFB1', 'RELA'],
    'PI3K'     : ['PIK3CA', 'PI3K(Class1)'],
    'TGFb'     : ['TGFBR1', 'TGFBR2'],
    'TNFa'     : ['TNFRSF1A'],
    'Trail'    : ['TNFSF10', 'BCL2', 'BCL-XL', 'BCL-W', 'MCL1'],
    'VEGF'     : ['VEGFR'],
    'WNT'      : ['WNT3A', 'GSK3A', 'GSK3B'],
    'p53'      : ['TP53']
}

''' Pathway -> Reactome Pathway ID '''
PATHWAY_REACTOME_MAP = {
    #TODO Pathway Reactome ID needs update
    'Androgen' : '',
    'EGFR'     : '',
    'Estrogen' : '',
    'Hypoxia'  : '',
    'JAK.STAT' : '',
    'MAPK'     : '',
    'NFkB'     : '',
    'PI3K'     : '',
    'TGFb'     : '',
    'TNFa'     : '',
    'Trail'    : '',
    'VEGF'     : '',
    'WNT'      : '',
    'p53'      : ''
}

''' These symbols are secondary/generic/typo that needs update '''
PROGENY_SYMBOL_MAPPING = {
    'NKFB1'   : 'NFKB1',
    'MAPK2K1' : 'PRKMK1',
    'PI3K(Class1)' : 'PIK3CA',
    'VEGFR'   : 'KDR',
    'BCL-W'   : 'BCL2L2',
    'BCL-XL'  : 'BCL2L1'
}

class PROGENY():
    def __init__(self):
        self.evidence_strings = list()
        self.symbols = {}
        self.logger = logging.getLogger(__name__)

    def process_progeny(self):
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()

        self.symbols = gene_parser.genes
        self.build_evidence()
        self.write_evidence()

    def build_evidence(self, filename=Config.PROGENY_FILENAME):

        now = datetime.datetime.now()
        '''
            build evidence.provenance_type object
        '''
        provenance_type = evidence_core.BaseProvenance_Type(
            database=evidence_core.BaseDatabase(
                id="PROGENY",
                version='2018.04',
                dbxref=evidence_core.BaseDbxref(url="https://saezlab.github.io/progeny/", id="PROGENY Pathway association analysis of TCGA tumor types", version="2018.04")),
            literature=evidence_core.BaseLiterature(
                references=[evidence_core.Single_Lit_Reference(lit_id="http://europepmc.org/abstract/MED/29295995")]
            )
        )
        error = provenance_type.validate(logging)

        if error > 0:
            self.logger.error(provenance_type.to_JSON(indentation=4))
            sys.exit(1)

        with open(filename, 'r') as progeny_input:
            n = 0

            for line in progeny_input:
                n +=1
                if n>1:
                    '''
                        Pathway Cancer_type     logFC   AveExpr t       P.Value adj.P.Val       B       Sample_size
                        VEGF    BRCA    82.5722740654623        -33.5309221366032       20.8841002773029        2.09462188741988e-55    2.93247064238783e-54    18.4517211174995        235
                        Hypoxia KIRC    625.176272105012        317.864789440655        23.2893493089142        1.4197622817373e-50     1.98766719443222e-49    11.4589892248308        144
                        p53     BRCA    -167.884630928143       11.9146505664384        -19.2135398711684       4.89209545526163e-50    3.42446681868314e-49    16.9318129144987        235
                    '''
                    (pathway_id, tumor_type, logfc, aveexpr, t, pval, fdr, b, sample) = tuple(line.rstrip().split('\t'))

                    # TODO pathway_id need to be mapped to Reactome pathway id:description via PATHWAY_REACTOME_MAP
                    # pathway = pathway_id.split(":")
                    # pathway_id = pathway[0].rstrip()
                    # pathway_desc = pathway[1].rstrip()
                    pathway_desc = 'Reactome pathway descriptions'

                    # TODO this is a temporary placeholder until gene_symbol is identify for each evidence
                    # gene_symbol = 'BRAF'

                    '''
                        build evidence.resource_score object
                    '''
                    resource_score = association_score.Pvalue(
                        type="pvalue",
                        method=association_score.Method(
                            description="PROGENY Pathway association analysis of TCGA tumor types as described in Schubert et al (2018)",
                            reference  ="http://europepmc.org/abstract/MED/29295995",
                            url="https://saezlab.github.io/progeny/"
                        ),
                        value=float(pval)
                    )

                    evidenceString = opentargets.Literature_Curated()
                    evidenceString.validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION
                    evidenceString.access_level = "public"
                    evidenceString.type = "affected_pathway"
                    evidenceString.sourceID = "progeny"

                    '''
                        build unique_association_field object
                    '''
                    evidenceString.unique_association_fields = {}
                    evidenceString.unique_association_fields['tumor_type_acronym'] = tumor_type
                    evidenceString.unique_association_fields['tumor_type'] = TUMOR_TYPE_MAP[tumor_type]
                    evidenceString.unique_association_fields['pathway_id'] = 'http://www.reactome.org/PathwayBrowser/#%s' % (pathway_id)
                    evidenceString.unique_association_fields['efo_id'] = TUMOR_TYPE_EFO_MAP[tumor_type]['uri']

                    target_type = 'http://identifiers.org/cttv.target/gene_evidence'
                    ensembl_gene_id = None

                    '''
                        loop through perturbed targets for each Pathway
                    '''
                    if pathway_id in PATHWAY_TARGET_MAP:
                        for gene_symbol in PATHWAY_TARGET_MAP[pathway_id]:

                            if gene_symbol in PROGENY_SYMBOL_MAPPING:
                                gene_symbol = PROGENY_SYMBOL_MAPPING[gene_symbol]
                            '''
                                build target object,
                            '''
                            if gene_symbol in self.symbols:

                                ensembl_gene_id = self.symbols[gene_symbol]

                                evidenceString.target = bioentity.Target(
                                    id="http://identifiers.org/ensembl/{0}".format(ensembl_gene_id),
                                    target_name=gene_symbol,
                                    # TODO activity is a required field in target object, currently set as unknown
                                    activity="http://identifiers.org/cttv.activity/unknown",
                                    target_type=target_type
                                )
                                '''
                                    build disease object
                                '''
                                evidenceString.disease = bioentity.Disease(
                                    id=TUMOR_TYPE_EFO_MAP[tumor_type]['uri'],
                                    name=TUMOR_TYPE_EFO_MAP[tumor_type]['label']
                                )
                                '''
                                    build evidence object
                                '''
                                evidenceString.evidence = evidence_core.Literature_Curated()
                                evidenceString.evidence.date_asserted = now.isoformat()
                                evidenceString.evidence.is_associated = True
                                # TODO check is this the correct evidence code "computational combinatorial evidence"
                                evidenceString.evidence.evidence_codes = ["http://purl.obolibrary.org/obo/ECO_0000053"]
                                evidenceString.evidence.provenance_type = provenance_type
                                evidenceString.evidence.resource_score = resource_score
                                '''
                                    build evidence.url object
                                '''
                                linkout = evidence_linkout.Linkout(
                                    url='http://www.reactome.org/PathwayBrowser/#%s' % (pathway_id),
                                    nice_name='%s' % (pathway_desc)
                                )

                                evidenceString.evidence.urls = [linkout]

                                error = evidenceString.validate(logging)

                                # if error > 0:
                                #    self.logger.error(evidenceString.to_JSON())
                                # sys.exit(1)

                                '''
                                    add gene_symbol in unique_association_field object
                                '''
                                evidenceString.unique_association_fields['symbol'] = gene_symbol

                                print(evidenceString.to_JSON(indentation=None))
                                ##TODO issue with append, take only last item of the gene
                                self.evidence_strings.append(evidenceString)

                            else:
                                self.logger.error("%s is not found in Ensembl" % gene_symbol)

            self.logger.error("%s evidence parsed"%(n-1))
            self.logger.error("%s evidence created"%len(self.evidence_strings))

        progeny_input.close()

    def write_evidence(self, filename=Config.PROGENY_EVIDENCE_FILENAME):
        self.logger.info("Writing PROGENY evidence strings")

        with open(filename, 'w') as progeny_output:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                self.logger.info(evidence_string.disease.id[0])
                # get max_phase_for_all_diseases
                error = evidence_string.validate(logging)
                if error == 0:
                    progeny_output.write(evidence_string.to_JSON(indentation=None)+"\n")
                else:
                    self.logger.error("REPORTING ERROR %i" %n)
                    self.logger.error(evidence_string.to_JSON(indentation=4))
            progeny_output.close()

'''
def main():
    import logging
    logger = logging.getLogger(__name__)

if __name__ == "__main__":
    main()
'''
