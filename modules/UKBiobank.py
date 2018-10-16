from settings import Config
from common.HGNCParser import GeneParser
import sys
import logging
import datetime
import opentargets.model.core as opentargets
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.genetics as evidence_genetics
import opentargets.model.evidence.linkout as evidence_linkout
import opentargets.model.evidence.association_score as association_score

__copyright__ = "Copyright 2014-2018, Open Targets"
__credits__   = ["ChuangKee Ong"]
__license__   = "Apache 2.0"
__version__   = "1.2.8"
__maintainer__= "ChuangKee Ong"
__email__     = ["data@opentargets.org"]
__status__    = "Production"

class UKBiobank():
    def __init__(self):
        self.evidence_strings = list()
        self.symbols = {}
        self.logger = logging.getLogger(__name__)

    def process_ukbiobank(self):
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()

        self.symbols = gene_parser.genes
        self.build_evidence()
        self.write_evidence()

    def build_evidence(self, filename=Config.UKBIOBANK_FILENAME):

        now = datetime.datetime.now()

        '''
            build evidence.provenance_type object
        '''
        provenance_type = evidence_core.BaseProvenance_Type(
            database=evidence_core.BaseDatabase(
                id="UKBiobank",
                version='2018.04',
                dbxref=evidence_core.BaseDbxref(url="http://www.ukbiobank.ac.uk/", id="UKBiobank", version="2018.04")),
            literature = evidence_core.BaseLiterature(
                references = [evidence_core.Single_Lit_Reference(lit_id="http://europepmc.org/abstract/MED/29518141")]
            )
        )
        error = provenance_type.validate(logging)

        if error > 0:
            self.logger.error(provenance_type.to_JSON(indentation=4))
            sys.exit(1)

        with open(filename, 'r') as ukbiobank_input:
            n = 0

            for line in ukbiobank_input:
                n +=1
                if n>1:
                    '''
                    ! Filter out lines where $11~/^ENSG/

                    rs1000007       4079    "Diastolic blood pressure, automated reading"   1.28E-07        317756  Neale et al. UK Biobank 2017-09-15      Common variants (MAF>1%)                -0.0144916 -ENSG00000198612 COPS8   nearest_gene_five_prime_end     -241901
                    rs10000093      23108   Impedance of leg (left) 2.41E-06        331296  Neale et al. UK Biobank 2017-09-15      Common variants (MAF>1%)                -0.0106638      -       ENSG00000145431     PDGFC   nearest_gene_five_prime_end     -489342
                    rs10000174      23116   Leg fat mass (left)     6.93E-06        331275  Neale et al. UK Biobank 2017-09-15      Common variants (MAF>1%)                -0.00866718     -       No nearest_gene_five_prime_end found!
                    '''
                    (variant_id, disease_id, trait, pval, sample_size, source, study, odd_ratio, beta_coef, beta_coef_dir, gene_id, gene_symbol, func_conseq, distance) = tuple(line.rstrip().split('\t'))

                    disease_id = 'http://www.ebi.ac.uk/efo/EFO_0000000' #+ disease_id
                    variant_id = 'http://identifiers.org/dbsnp/' + variant_id
                    #TODO study & func_conseq needs update
                    study      = 'STUDYID_' + study
                    func_conseq= 'http://purl.obolibrary.org/obo/SO_0001631'
                    '''
                        build evidence.resource_score object
                    '''
                    resource_score = association_score.Pvalue(
                        type="pvalue",
                        method=association_score.Method(
                            description="UKBiobank",
                            reference  ="http://europepmc.org/abstract/MED/29518141",
                            url="http://www.ukbiobank.ac.uk/"
                        ),
                        value=float(pval)
                    )

                    evidenceString = opentargets.Genetics()
                    evidenceString.validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION
                    evidenceString.access_level = "public"
                    evidenceString.type = "genetic_association"
                    evidenceString.sourceID = "ukbiobank"
                    '''
                        build unique_association_field object
                    '''
                    evidenceString.unique_association_fields = {}
                    evidenceString.unique_association_fields['symbol'] = gene_symbol
                    evidenceString.unique_association_fields['gene_id'] = gene_id
                    evidenceString.unique_association_fields['disease_id'] = disease_id
                    evidenceString.unique_association_fields['variant'] = variant_id
                    evidenceString.unique_association_fields['sample_size'] = sample_size
                    evidenceString.unique_association_fields['pvalue'] = pval
                    evidenceString.unique_association_fields['odd_ratio'] = odd_ratio
                    evidenceString.unique_association_fields['study'] = study

                    target_type = 'http://identifiers.org/cttv.target/gene_evidence'
                    ensembl_gene_id = None

                    if gene_symbol in self.symbols:
                        ensembl_gene_id = self.symbols[gene_symbol]

                        '''
                            build target object,
                        '''
                        evidenceString.target = bioentity.Target(
                            id="http://identifiers.org/ensembl/{0}".format(ensembl_gene_id),
                            target_name=gene_symbol,
                            #TODO activity is a required field in target object, currently set as unknown
                            activity="http://identifiers.org/cttv.activity/unknown",
                            target_type=target_type
                        )

                        '''
                            build disease object
                        '''
                        evidenceString.disease = bioentity.Disease(
                            id=disease_id,
                            #TODO update disease label
                            name=disease_id
                        )

                        '''
                            build variant object
                        '''
                        evidenceString.variant = bioentity.Variant(
                            id=variant_id,
                            type="snp single"
                        )

                        '''
                            build evidence variant2disease object
                        '''
                        #TODO variant2disease & gene2variant object missing
                        evidenceString.evidence = opentargets.GeneticsEvidence(
                            variant2disease=evidence_genetics.Variant2Disease(
                                unique_experiment_reference=study,
                                date_asserted=now.isoformat(),
                                resource_score=resource_score,#evidence_score.Pvalue(value=2.000000039082963e-25),
                                provenance_type=provenance_type,
                                evidence_codes=["http://purl.obolibrary.org/obo/ECO_0000205"],
                                gwas_sample_size=int(sample_size)
                            ),
                            gene2variant=evidence_genetics.Gene2Variant(
                                date_asserted=now.isoformat(),
                                provenance_type=provenance_type,
                                evidence_codes=["http://purl.obolibrary.org/obo/ECO_0000205"],
                                functional_consequence=func_conseq
                            ),
                        )
                        error = evidenceString.validate(logging)

                        if error > 0:
                            self.logger.error(evidenceString.to_JSON())
                            sys.exit(1)

                        self.evidence_strings.append(evidenceString)

                    else:
                        self.logger.error("%s is not found in Ensembl" % gene_symbol)

            self.logger.info("%s evidence parsed"%(n-1))
            self.logger.info("%s evidence created"%len(self.evidence_strings))

        ukbiobank_input.close()

    def write_evidence(self, filename=Config.UKBIOBANK_EVIDENCE_FILENAME):
        self.logger.info("Writing UKBiobank evidence strings")
        with open(filename, 'w') as ukbiobank_output:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                self.logger.info(evidence_string.disease.id[0])
                # get max_phase_for_all_diseases
                error = evidence_string.validate(logging)

                if error == 0:
                    print(evidence_string.to_JSON(indentation=None))
                    ukbiobank_output.write(evidence_string.to_JSON(indentation=None)+"\n")
                else:
                    self.logger.error("REPORTING ERROR %i" %n)
                    self.logger.error(evidence_string.to_JSON(indentation=4))
            ukbiobank_output.close()


def main():
    UKBiobank().process_ukbiobank()

if __name__ == "__main__":
    main()
