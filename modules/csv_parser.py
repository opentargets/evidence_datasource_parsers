import csv
from common.EFOData import OBOParser
from common.HGNCParser import GeneParser
import json
import logging
from datetime import datetime

logging.basicConfig(filename='output_phewascatalog.log',
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                            level=logging.INFO)


class Phewas(object):

    def __init__(self, phewas_id, snp = None, ensg_id = None, efo_id = None, phenotype = None, gene_name = None, xref = None):
        self.ensg_id = ensg_id
        self.efo_id = efo_id
        self.snp = snp
        self.phenotype = phenotype
        self.gene_name = gene_name
        self.phewas_id = phewas_id
        self.xref = xref

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,
                          sort_keys=True, indent=4)


class CSVParser(object):



    def parse(self,filename):
        with open(filename, "rb") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                yield row


class PhewasProcessor(object):
    def __init__(self,schema_version ):
        self.genes = dict()
        self.efos = list()
        self.icd9 = dict()
        self.csv_parser = CSVParser()
        self.schema_version = schema_version




    def find_efo(self,phewas_phenotype, phecode):

        matched_efos = [efo_dict for efo_dict in self.efos if phewas_phenotype.lower() in efo_dict['synonyms'] ]
        if not matched_efos:


            # matched_efos = [efo_dict for efo_dict in self.efos if
            #                 any([True for e in self.icd9.get(phecode, []) if e in efo_dict['icd9']])]
            matched_efos = []
            icd9s = self.icd9.get(phecode, [])
            for efo_dict in self.efos:
                if icd9s:

                    for icd9 in icd9s:
                        if icd9 in efo_dict['icd9']:
                            matched_efos.append(efo_dict)
                            break



            # if not matched_efos:
            #     self.find_zooma_phenotype_mapping(phewas_phenotype)

        return matched_efos


    def setup(self):
        obo_parser = OBOParser('../resources/efo.obo')
        obo_parser.parse()
        self.efos = obo_parser.efos
        self.obsolete_efos = obo_parser.get_obsolete_efos()

        hp_obo_parser = OBOParser('../resources/hp.obo')
        hp_obo_parser.parse()
        self.efos.extend(hp_obo_parser.efos)

        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes


        for icd9_row in self.csv_parser.parse('../resources/phecode_icd9_rolled.csv'):
            pheCode = icd9_row['PheCode']
            try:

                icd9_code = float(icd9_row['ICD9'])
            except ValueError:
                icd9_code = icd9_row['ICD9']

            if self.icd9.get(pheCode):
                self.icd9[pheCode].append(icd9_code)
            else:
                self.icd9[pheCode] = list()
                self.icd9[pheCode].append(icd9_code)


    def convert_phewas_catalog_evidence_json(self):

        logging.info('Start processing ')

        missing_efo_fieldnames = ['phenotype', 'similar_efo']
        fieldnames = ['phenotype', 'efo_id', 'gene_name','ensg_id','cases','p-value','odds-ratio','snp']
        logging.info('Start the phewas catalog mapping')
        with open('missing_efo.csv', 'w') as out_missing_csv , open('phewas_efo_ensg.csv', 'w') as out_csv, open('../output/phewas_catalog.json', 'w') as out_json:
            writer = csv.DictWriter(out_csv, fieldnames)
            writer.writeheader()
            missing_efo_writer = csv.DictWriter(out_missing_csv, missing_efo_fieldnames)
            missing_efo_writer.writeheader()
            for phewas_row in self.csv_parser.parse('../resources/phewas-catalog.csv'):

                ensg_id = self.genes.get(phewas_row['gene_name'])
                if ensg_id:
                    matched_efos = self.find_efo(phewas_row['phewas phenotype'],phewas_row['phewas code'])
                    if matched_efos :
                        for efo in matched_efos:
                            inner_dict = dict(zip(fieldnames, [phewas_row['phewas phenotype'], efo['id'],phewas_row['gene_name'],ensg_id,phewas_row['cases'],phewas_row['p-value'],phewas_row['odds-ratio'],phewas_row['snp']]))
                            evidence = self.generate_evidence(phewas_row,efo['id'], ensg_id)
                            writer.writerow(inner_dict)
                            if evidence:
                                json.dump(evidence,out_json)
                                out_json.write('\n')
                    else:
                        inner_dict = dict(zip(missing_efo_fieldnames, [phewas_row['phewas phenotype'], '']))
                        missing_efo_writer.writerow(inner_dict)

        logging.info('Completed')

    def generate_evidence(self,phewas_dict, disease_id, target_id):
        phewas_evidence = dict()
        disease_id = disease_id.replace(':','_')
        if disease_id in self.obsolete_efos:
            new_efo = self.obsolete_efos[disease_id]
            if new_efo:
                print 'Obsolete efo- {} New EFO - {}'.format(disease_id, new_efo)
                disease_id = new_efo
            else:
                return None

        if disease_id.startswith('EFO'):
            phewas_evidence['disease'] = {'id': 'http://www.ebi.ac.uk/efo/'+disease_id}
        elif disease_id.startswith('HP') or disease_id.startswith('MP') :
            phewas_evidence['disease'] = {'id': 'http://purl.obolibrary.org/obo/' + disease_id}
        elif disease_id.startswith('Orphanet') :
            phewas_evidence['disease'] = {'id': 'http://www.orpha.net/ORDO/' + disease_id}
        elif disease_id.startswith('NCBITaxon') :
            phewas_evidence['disease'] = {'id': 'http://purl.obolibrary.org/obo/' + disease_id}

        phewas_evidence['target'] = {"activity": "http://identifiers.org/cttv.activity/predicted_damaging",
                    "id": "http://identifiers.org/ensembl/{}".format(target_id),
                    "target_type": "http://identifiers.org/cttv.target/gene_evidence"}
        if phewas_evidence.get('target') and phewas_evidence.get('disease'):
            phewas_evidence['validated_against_schema_version'] = self.schema_version
            phewas_evidence["access_level"] = "public"
            phewas_evidence["sourceID"] = "phewas_catalog"
            phewas_evidence['type'] = 'genetic_association'
            phewas_evidence["variant"]= {"type": "snp single", "id": "http://identifiers.org/dbsnp/{}".format(phewas_dict['snp'])}
            phewas_evidence['unique_association_fields'] = {'odds_ratio':phewas_dict['odds-ratio'], 'cases' : phewas_dict['cases'], 'phenotype' : phewas_dict['phewas phenotype']}

            #phewas_evidence['resource_score'] = {'type': 'pvalue', 'method': {"description":"pvalue for the phenotype to snp association."},"value":phewas_dict['p-value']}
            i = datetime.now()

            evidence = dict()
            evidence['variant2disease'] = {'unique_experiment_reference':'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3969265/',
                                           'provenance_type': {"literature":{"references":[{"lit_id":"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3969265/"}]},
                                                               "expert":{"status":True,"statement":"Primary submitter of data"},
                                                               "database":{"version":"2017-07-01T09:53:37+00:00","id":"PHEWAS Catalog",
                                                                           "dbxref":{"version":"2017-07-01T09:53:37+00:00","id":"http://identifiers.org/phewascatalog"}}},
                                           'is_associated': True,
                                           'resource_score':{'type': 'pvalue', 'method': {"description":"pvalue for the phenotype to snp association."},"value":float(phewas_dict['p-value'])},
                                           'date_asserted': i.strftime('%Y-%m-%dT%H:%M:%S+00:00'),
                                           'evidence_codes': ['http://identifiers.org/eco/GWAS','http://purl.obolibrary.org/obo/ECO_0000205'],
                                           }
            evidence['gene2variant'] = {'provenance_type': {"expert":{"status":True,"statement":"Primary submitter of data"},
                                                            "database":{"version":"2017-07-01T09:53:37+00:00","id":"PHEWAS Catalog","dbxref":{"version":"2017-07-01T09:53:37+00:00","id":"http://identifiers.org/phewascatalog"}}},
                                        'is_associated': True, 'date_asserted' : i.strftime('%Y-%m-%dT%H:%M:%S+00:00'),
                                        'evidence_codes':["http://identifiers.org/eco/cttv_mapping_pipeline", "http://purl.obolibrary.org/obo/ECO_0000205"],
                                        'functional_consequence':'http://purl.obolibrary.org/obo/SO_0001632'}
            phewas_evidence['evidence'] = evidence
        else:
            logging.info('Missing disease/target evidence : {}'.format(phewas_evidence))
            phewas_evidence = None

        return phewas_evidence



def remove_dup():


    with open('../missing_efo.csv', 'r') as in_file, open('../missing_efos.csv', 'w') as out_file:
        seen = set()  # set for fast O(1) amortized lookup
        reader = csv.reader(in_file, dialect=csv.excel_tab)
        for row in reader:
            if row[0] in seen: continue  # skip duplicate

            seen.add(row[0])
            out_file.write(row[0])
            out_file.write('\n')

def unique_phenotypes():
    with open('../phewas-catalog.csv', 'r') as in_file, open('../unique_efos.csv', 'w') as out_file:
        seen = set()  # set for fast O(1) amortized lookup
        reader = csv.reader(in_file)
        for row in reader:
            if row[2] in seen: continue  # skip duplicate

            seen.add(row[2])
            out_file.write(row[2])
            out_file.write('\n')

def main():
    phewas_processor = PhewasProcessor('1.2.7')
    phewas_processor.setup()
    phewas_processor.convert_phewas_catalog_evidence_json()
    #remove_dup()
    #unique_phenotypes()
    #phewas_processor.find_zooma_phenotype_mapping('test')


if __name__ == "__main__":
    main()
