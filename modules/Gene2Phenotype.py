from settings import Config
from common.HGNCParser import GeneParser
from common.RareDiseasesUtils import RareDiseaseMapper

import python_jsonschema_objects as pjo

import sys
import logging
import csv
import gzip
import requests
import datetime

__copyright__  = "Copyright 2014-2020, Open Targets"
__credits__    = ["Gautier Koscielny", "ChuangKee Ong", "Michaela Spitzer", "Asier Gonzalez" ]
__license__    = "Apache 2.0"
__version__    = "1.3.0"
__maintainer__ = "Open Targets Data Team"
__email__      = ["data@opentargets.org"]
__status__     = "Production"

class G2P(RareDiseaseMapper):
    def __init__(self, schema_json=Config.OT_JSON_SCHEMA, schema_version=Config.VALIDATED_AGAINST_SCHEMA_VERSION):
        super(G2P, self).__init__()
        self.genes = None
        self.evidence_strings = list()

        self._logger = logging.getLogger(__name__)

        # Initialize json builder based on the schema:
        json_schema = requests.get(schema_json).json()
        self.builder = pjo.ObjectBuilder(json_schema)
        self.evidence_builder = self.builder.build_classes()
        self.schema_version = schema_version

    def process_g2p(self):

        self.get_omim_to_efo_mappings()

        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

        # Parser DD file
        self.generate_evidence_strings(Config.G2P_DD_FILENAME)
        # Parser DD file
        self.generate_evidence_strings(Config.G2P_eye_FILENAME)
        # Parser DD file
        self.generate_evidence_strings(Config.G2P_skin_FILENAME)
        # Parser DD file
        self.generate_evidence_strings(Config.G2P_cancer_FILENAME)

        # Save results to file
        self.write_evidence_strings(Config.G2P_EVIDENCE_FILENAME)

    def generate_evidence_strings(self, filename):

        total_efo = 0

        with gzip.open(filename, mode='rt') as zf:

            reader = csv.DictReader(zf, delimiter=',', quotechar='"')
            c = 0
            for row in reader:
                c += 1
                if c > 1:
                    # Column names are:
                    # "gene symbol","gene mim","disease name","disease mim","DDD category","allelic requirement",
                    # "mutation consequence",phenotypes,"organ specificity list",pmids,panel,"prev symbols","hgnc id",
                    # "gene disease pair entry date"
                    gene_symbol = row["gene symbol"]
                    disease_name = row["disease name"]
                    disease_mim = row["disease mim"]
                    allelic_requirement = row["allelic requirement"]
                    mutation_consequence = row["mutation consequence"]
                    confidence = row["DDD category"]
                    panel = row["panel"]

                    date = row["gene disease pair entry date"]
                    # Handle missing dates ("No date" in file)
                    try:
                        date = datetime.datetime.strptime(date, "%Y-%m-%d %H:%M:%S").isoformat()
                    except ValueError:
                        date = "N/A"


                    gene_symbol.rstrip()

                    if gene_symbol in self.genes:
                        # Map gene symbol to ensembl
                        target = self.genes[gene_symbol]
                        ensembl_iri = "http://identifiers.org/ensembl/" + target

                        # Map disease to EFO or Orphanet
                        if disease_mim in self.omim_to_efo_map:
                            total_efo +=1
                            diseases = self.omim_to_efo_map[disease_mim]

                            for disease in diseases:
                                self._logger.info("%s %s %s %s"%(gene_symbol, target, disease_name, disease['efo_uri']))

                                type = "genetic_literature"

                                provenance_type = {
                                    'database' : {
                                        'id' : "Gene2Phenotype",
                                        'version' : '2020.04.02',
                                        'dbxref' : {
                                            'url': "http://www.ebi.ac.uk/gene2phenotype",
                                            'id' : "Gene2Phenotype",
                                            'version' : "2020.04.02"

                                        }
                                    },
                                    'literature' : {
                                        'references' : [
                                            {
                                                'lit_id' : "http://europepmc.org/abstract/MED/25529582"
                                            }
                                        ]
                                    }
                                }

                                # *** General properties ***
                                access_level = "public"
                                sourceID = "gene2phenotype"
                                validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION

                                # *** Target info ***
                                target = {
                                    'id' : ensembl_iri,
                                    'activity' : "http://identifiers.org/cttv.activity/unknown",
                                    'target_type' : "http://identifiers.org/cttv.target/gene_evidence",
                                    'target_name' : gene_symbol
                                }
                                # http://www.ontobee.org/ontology/ECO?iri=http://purl.obolibrary.org/obo/ECO_0000204 -- An evidence type that is based on an assertion by the author of a paper, which is read by a curator.

                                # *** Disease info ***
                                disease_info = {
                                    'id' : disease['efo_uri'],
                                    'name' : disease['efo_label'],
                                    'source_name' : disease_name
                                }
                                # *** Evidence info ***
                                # Score based on mutational consequence
                                if confidence in Config.G2P_confidence2score:
                                    score = Config.G2P_confidence2score[confidence]
                                else:
                                    self.logger.error('{} is not a recognised G2P confidence, assigning an score of 0'.format(confidence))
                                    score = 0
                                resource_score = {
                                    'type': "probability",
                                    'value': score
                                }

                                # Linkout
                                linkout = [
                                    {
                                        'url' : 'http://www.ebi.ac.uk/gene2phenotype/search?panel=ALL&search_term=%s' % (gene_symbol,),
                                        'nice_name' : 'Gene2Phenotype%s' % (gene_symbol)
                                    }
                                ]

                                evidence = {
                                    'is_associated' : True,
                                    'confidence' : confidence,
                                    'allelic_requirement' : allelic_requirement,
                                    'mutation_consequence' : mutation_consequence,
                                    'evidence_codes' : ["http://purl.obolibrary.org/obo/ECO_0000204"],
                                    'provenance_type' : provenance_type,
                                    'date_asserted' : date,
                                    'resource_score' : resource_score,
                                    'urls' : linkout
                                }
                                # *** unique_association_fields ***
                                unique_association_fields = {
                                    'target_id' : ensembl_iri,
                                    'original_disease_label' : disease_name,
                                    'disease_id' : disease['efo_uri'],
                                    'gene_panel': panel
                                }


                                try:
                                    evidence = self.evidence_builder.Opentargets(
                                        type = type,
                                        access_level = access_level,
                                        sourceID = sourceID,
                                        evidence = evidence,
                                        target = target,
                                        disease = disease_info,
                                        unique_association_fields = unique_association_fields,
                                        validated_against_schema_version = validated_against_schema_version
                                    )
                                    self.evidence_strings.append(evidence)
                                except:
                                    self.logger.warning('Evidence generation failed for row: {}'.format(c))
                                    raise
                    else:
                        self._logger.error("%s\t%s not mapped: please check manually"%(disease_name, disease_mim))

                print("%i %i" % (total_efo, c))

    def write_evidence_strings(self, filename):
        self._logger.info("Writing Gene2Phenotype evidence strings to %s", filename)
        with open(filename, 'w') as tp_file:
            n = 0
            for evidence_string in self.evidence_strings:
                n += 1
                self._logger.info(evidence_string['disease']['id'])
                tp_file.write(evidence_string.serialize() + "\n")
        tp_file.close()


def main():
    g2p = G2P()
    g2p.process_g2p()

if __name__ == "__main__":
    main()