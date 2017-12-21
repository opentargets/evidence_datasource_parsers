from settings import Config, file_or_resource
from common.HGNCParser import GeneParser
from common.RareDiseasesUtils import RareDiseaseMapper
from tqdm import tqdm
import google.cloud.storage as gcs
import google.cloud.exceptions as gce

import opentargets.model.core as opentargets
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.phenotype as evidence_phenotype
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.linkout as evidence_linkout
import opentargets.model.evidence.association_score as association_score
import opentargets.model.evidence.mutation as evidence_mutation

import sys
import gzip
import logging
import urllib3
import csv

__copyright__  = "Copyright 2014-2017, Open Targets"
__credits__    = ["Gautier Koscielny", "ChuangKee Ong"]
__license__    = "Apache 2.0"
__version__    = "1.2.7"
__maintainer__ = "ChuangKee Ong"
__email__      = ["gautierk@targetvalidation.org", "ckong@ebi.ac.uk"]
__status__     = "Production"

G2P_FILENAME = Config.G2P_FILENAME
G2P_EVIDENCE_FILENAME = Config.G2P_EVIDENCE_FILENAME

class G2P(RareDiseaseMapper):
    def __init__(self, es=None):

        super().__init__()
        self.genes = None
        self.evidence_strings = list()
        self.gcs_client = gcs.Client()
        self.bucket = None
        self._logger = logging.getLogger(__name__)

        try:
            for bucket in self.gcs_client.list_buckets():
                print(bucket)
            self.bucket = self.gcs_client.get_bucket(Config.GOOGLE_BUCKET_EVIDENCE_INPUT)
        except gce.NotFound:
            self._logger.error('Please set export GOOGLE_APPLICATION_CREDENTIALS=/path/to/key.json')
            print('Sorry, that bucket does not exist!')

    def process_g2p(self):

        self.get_omim_to_efo_mappings()

        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

        self.generate_evidence_strings(G2P_FILENAME)
        self.write_evidence_strings(G2P_EVIDENCE_FILENAME)

    def generate_evidence_strings(self, filename):
        total_efo = 0
        blob = self.bucket.get_blob(filename)
        with open('/tmp/%s'%(G2P_FILENAME), 'wb') as file_obj:
            blob.download_to_file(file_obj)

        with gzip.open('/tmp/%s'%(G2P_FILENAME), mode='rt') as zf:
            reader = csv.reader(zf, delimiter=',', quotechar='"')
            c = 0
            for row in reader:

                c += 1
                if c > 1:
                    '''
                    "gene symbol","gene mim","disease name","disease mim","DDD category","allelic requirement","mutation consequence",phenotypes,"organ specificity list",pmids,panel,"prev symbols","hgnc id"
                    '''
                    (gene_symbol, gene_mim, disease_name, disease_mim, DDD_category, allelic_requirement, mutation_consequence, phenotypes, organ_specificity_list,pmids,panel, prev_symbols, hgnc_id) = row
                    gene_symbol.rstrip()

                    if gene_symbol in self.genes:
                        ''' map gene symbol to ensembl '''
                        target = self.genes[gene_symbol]
                        ensembl_iri = "http://identifiers.org/ensembl/" + target

                        ''' Map disease to EFO or Orphanet '''
                        if disease_mim in self.omim_to_efo_map:
                            total_efo +=1
                            diseases = self.omim_to_efo_map[disease_mim]

                            for disease in diseases:
                                self._logger.info("%s %s %s %s"%(gene_symbol, target, disease_name, disease['efo_uri']))

                                obj = opentargets.Literature_Curated(type='genetic_literature')
                                provenance_type = evidence_core.BaseProvenance_Type(
                                    database=evidence_core.BaseDatabase(
                                        id="Gene2Phenotype",
                                        version='v0.2',
                                        dbxref=evidence_core.BaseDbxref(
                                            url="http://www.ebi.ac.uk/gene2phenotype",
                                            id="Gene2Phenotype", version="v0.2")),
                                    literature=evidence_core.BaseLiterature(
                                        references=[evidence_core.Single_Lit_Reference(lit_id="http://europepmc.org/abstract/MED/25529582")]
                                    )
                                )

                            obj.access_level = "public"
                            obj.sourceID = "gene2phenotype"
                            obj.validated_against_schema_version = "1.2.7"
                            obj.unique_association_fields = {"target": ensembl_iri, "original_disease_label" : disease_name, "disease_uri": disease['efo_uri'], "source_id": "gene2phenotype"}
                            obj.target = bioentity.Target(id=ensembl_iri,
                                                          activity="http://identifiers.org/cttv.activity/unknown",
                                                          target_type='http://identifiers.org/cttv.target/gene_evidence',
                                                          target_name=gene_symbol)
                            # http://www.ontobee.org/ontology/ECO?iri=http://purl.obolibrary.org/obo/ECO_0000204 -- An evidence type that is based on an assertion by the author of a paper, which is read by a curator.
                            resource_score = association_score.Probability(
                                type="probability",
                                value=1)

                            obj.disease = bioentity.Disease(id=disease['efo_uri'], name=disease['efo_label'], source_name=disease_name)
                            obj.evidence = evidence_core.Literature_Curated()
                            obj.evidence.is_associated = True
                            obj.evidence.evidence_codes = ["http://purl.obolibrary.org/obo/ECO_0000204"]
                            obj.evidence.provenance_type = provenance_type
                            obj.evidence.date_asserted = '2017-06-14T00:00:00'
                            obj.evidence.provenance_type = provenance_type
                            obj.evidence.resource_score = resource_score
                            linkout = evidence_linkout.Linkout(
                                url='http://www.ebi.ac.uk/gene2phenotype/search?panel=ALL&search_term=%s' % (
                                gene_symbol,),
                                nice_name='Gene2Phenotype%s' % (gene_symbol))
                            obj.evidence.urls = [linkout]
                            error = obj.validate(logging)

                            if error > 0:
                                logging.error(obj.to_JSON())
                                sys.exit(1)
                            else:
                                self.evidence_strings.append(obj)
                    else:
                        self._logger.error("%s\t%s not mapped: please check manually"%(disease_name, disease_mim))

            print("%i %i" % (total_efo, c))

    def write_evidence_strings(self, filename):
        self._logger.info("Writing IntOGen evidence strings")
        with open('/tmp/' + filename, 'w') as tp_file:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                self._logger.info(evidence_string.disease.id[0])
                # get max_phase_for_all_diseases
                error = evidence_string.validate(logging)
                if error == 0:
                    tp_file.write(evidence_string.to_JSON(indentation=None) + "\n")
                else:
                    self._logger.error("REPORTING ERROR %i" % n)
                    self._logger.error(evidence_string.to_JSON(indentation=4))
                    # sys.exit(1)


        blob = self.bucket.blob(filename)
        with open('/tmp/' + filename, 'rb') as my_file:
            blob.upload_from_file(my_file)

def main():
    g2p = G2P()

if __name__ == "__main__":
    main()