from common.Redis import RedisQueueWorkerProcess,RedisQueue
from pymongo import MongoClient
from datetime import datetime
import multiprocessing
import uuid
from redislite import Redis
import os
import logging
import math
import json
import csv

from common.HGNCParser import GeneParser
from common.Settings import Config

NO_OF_WORKERS = 8
MAX_CHUNKS =100
UNIQUE_RUN_ID = str(uuid.uuid4()).replace('-', '')[:16]
TEMP_DIR = os.path.join(os.path.sep, 'tmp')
REDISLITE_DB_PATH = os.path.join(TEMP_DIR, 'phewas_redislite.rdb')

class MongoDataManager(object):

    def __init__(self):
        self.genes = dict()
        self.efos = dict()



    def setup(self):
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

        for phenotype_row in self.parse('../23andme phenotypes.csv'):
            phenotype_dict =  dict()
            phenotype_dict['phenotype'] = phenotype_row['phenotype']
            phenotype_dict['efo'] = phenotype_row['efo']
            phenotype_dict['no_of_cases'] = phenotype_row['case']
            phenotype_dict['notes'] = phenotype_row['notes']

            self.efos[phenotype_row['source']] = phenotype_dict



    def parse(self, filename):
            with open(filename, "rb") as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    yield row

    def process(self):
        logging.basicConfig(filename='output_mongo.log',
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                            level=logging.INFO)

        r_server = Redis(dbfilename=str(REDISLITE_DB_PATH), serverconfig={'save': [], 'maxclients': 10000})
        mongo_client = MongoClient(Config.MONGO_URL)
        db = mongo_client[Config.MONGO_DB]
        db.authenticate(Config.MONGO_USER,Config.MONGO_PWD)
        gwas23andme = db[Config.MONGO_TABLE]


        no_of_workers = NO_OF_WORKERS or multiprocessing.cpu_count()

        mongo_q = RedisQueue(queue_id=UNIQUE_RUN_ID + '|mongo_processor_q',
                                  max_size=MAX_CHUNKS * no_of_workers,
                                  job_timeout=120)
        transformer_q = RedisQueue(queue_id=UNIQUE_RUN_ID + '|data_transformer_q',
                             max_size=MAX_CHUNKS * no_of_workers,
                             job_timeout=120)
        writer_q = RedisQueue(queue_id=UNIQUE_RUN_ID + '|writer_q',
                             max_size=MAX_CHUNKS * no_of_workers,
                             job_timeout=120)


        workers = [MongoWorkerProcess(mongo_q,
                                      r_server.db,
                                      transformer_q,
                                               ) for i in range(no_of_workers)]

        for w in workers:
            w.start()

        transformers = [DataTransformerProcess(self.efos,
                                               self.genes,
                                               transformer_q,
                                      r_server.db,
                                      writer_q,
                                      ) for i in range(no_of_workers)]

        for t in transformers:
            t.start()

        # writers = [FileWriterProcess(writer_q,
        #                               r_server.db
        #                               )for i in range(no_of_workers)]
        #
        # for writer in writers:
        #     writer.start()

        writer = FileWriterProcess(writer_q,r_server.db)
        writer.start()


        logging.info('Getting distinct genes and phenotypes')
        distinct_phenotypes = list(gwas23andme.distinct("source"))
        #distinct_genes = list(gwas23andme.distinct("ingenes"))
        logging.info('Start the real deal here!!!! ')
        #i=1
        for phenotype in distinct_phenotypes:
            logging.info('Processing phenotype {} '.format(phenotype))
            #logging.info('Total docs for phenotype {} are {}'.format(phenotype, gwas23andme.count({"source" : phenotype})))
            # if not phenotype.startswith('PD'):
            mongo_q.put((phenotype, 'test'), r_server)
                # i = i + 1
                # if i>3:
                #     break
            # for gene in distinct_genes:
            #     mongo_q.put((phenotype,gene),r_server)


        mongo_q.set_submission_finished(r_server=r_server)

        for a in workers:
                a.join()
        for a in transformers:
                a.join()
        writer.join()
        # for w in writers:
        #     w.join()

class MongoWorkerProcess(RedisQueueWorkerProcess):

    def __init__(self, queue_in,
                 redis_path,
                 queue_out):
        super(MongoWorkerProcess, self).__init__(queue_in,redis_path,queue_out)
        self.gwas23andme = None
        self.val = '1e-5'




    def process(self, data):
        if not self.gwas23andme:
            mongo_client = MongoClient(Config.MONGO_URL)
            db = mongo_client[Config.MONGO_DB]
            db.authenticate(Config.MONGO_USER, Config.MONGO_PWD)
            gwas23andme = db[Config.MONGO_TABLE]


        phenotype, gene  = data

        # docs_cursor = self.gwas23andme.find({'source': phenotype, 'best_pvalue': {'$gt': 0.5}, 'ingenes' : gene},
        #                          {"source": 1, "adj_pvalue": 1, "im_num_1": 1, "im_num_0": 1, "alleles": 1,
        #                                 "effect": 1, "assay_name": 1, 'ingenes': 1, 'dtype': 1, "best_effect": 1 })
        docs_cursor = gwas23andme.find({"source": phenotype, 'best_pvalue': {'$exists': True, '$lte': float(self.val)}},
                                            {"source": 1, "adj_pvalue": 1, "im_num_1": 1, "im_num_0": 1, "alleles": 1,
                                         "effect": 1, "assay_name": 1, 'ingenes': 1, 'dtype': 1, "best_effect": 1 })

        gwas_coll = list(docs_cursor)
        logging.info('Found {} documents for phenotype {} meeting significant threshold'.format(len(gwas_coll),phenotype))
        # filename = 'output/23andme_' + phenotype + '.json'
        #return (filename, json_util.dumps(gwas_coll))

        return ( gwas_coll)

class DataTransformerProcess(RedisQueueWorkerProcess):

    def __init__(self,
                 efos,
                 genes,
                 queue_in,
                 redis_path,
                 queue_out):
        super(DataTransformerProcess, self).__init__(queue_in,redis_path,queue_out)
        self.efos = efos
        self.genes = genes


    def process(self, data):
        # filename,  other_info = data
        other_info = data
        evidences = list()
        for gwas_dict in other_info:
            if gwas_dict['ingenes']:
                if gwas_dict.get('dtype') == 'binary':
                    if gwas_dict.get('best_effect'):
                        try:
                            gwas_dict['odds_ratio'] = str(math.exp(gwas_dict.get('best_effect')))
                        except OverflowError:
                            logging.info('Overflow Error for document  {} with best_effect {} '.format(gwas_dict['_id'],gwas_dict.get('best_effect')))
                gwas_id_str = str(gwas_dict['_id'])
                gwas_dict['_id'] = gwas_id_str
                evidence_list = self.generate_evidence(gwas_dict)
                if evidence_list:
                    evidences.extend(evidence_list)

        return evidences
        #return (json.dumps(other_info))
        #return (filename, json_util.dumps(other_info))

    def generate_evidence(self, phewas_dict):
        ensg_ids = list()

        for gene_name in phewas_dict['ingenes']:

            ensg_ids.append(self.genes.get(gene_name))

        phenotype_dict = self.efos[phewas_dict['source']]
        disease_id = phenotype_dict['efo']

        if disease_id.startswith('EFO'):
            disease = {'id': 'http://www.ebi.ac.uk/efo/' + disease_id}
        elif disease_id.startswith('HP'):
            disease = {'id': 'http://purl.obolibrary.org/obo/' + disease_id}
        else:
            disease = None
        phewas_evidences = list()
        phewas_evidence = dict()
        for ensg_id in ensg_ids:

            if disease and phenotype_dict['no_of_cases']:

                phewas_evidence['disease'] = disease
                phewas_evidence['target'] = {"activity": "http://identifiers.org/cttv.activity/predicted_damaging",
                                             "id": "http://identifiers.org/ensembl/{}".format(ensg_id),
                                             "target_type": "http://identifiers.org/cttv.target/gene_evidence"}
                phewas_evidence['validated_against_schema_version'] = '1.2.6'
                phewas_evidence["access_level"] = "public"
                phewas_evidence["sourceID"] = "twentythreeandme"
                phewas_evidence['type'] = 'genetic_association'
                phewas_evidence["variant"] = {"type": "snp single",
                                              "id": "http://identifiers.org/dbsnp/{}".format(phewas_dict['assay_name'])}
                phewas_evidence['unique_association_fields'] = {'odds_ratio': phewas_dict.get('odds_ratio'),
                                                                'cases': phenotype_dict['no_of_cases'],
                                                                'phenotype': phewas_dict['source'],
                                                                '23andme_id': phewas_dict['_id'],
                                                                'notes': phenotype_dict['notes']}

                # phewas_evidence['resource_score'] = {'type': 'pvalue', 'method': {"description":"pvalue for the phenotype to snp association."},"value":phewas_dict['p-value']}
                i = datetime.now()

                evidence = dict()
                evidence['variant2disease'] = {'unique_experiment_reference': 'http://europepmc.org/abstract/MED/0',
                                               'provenance_type': {"literature": {"references": [{"lit_id": "http://europepmc.org/abstract/MED/1"}]},
                                                                   "expert": {"status": True,
                                                                              "statement": "Primary submitter of data"},
                                                                   "database": {"version": "2017-06-01T09:53:37+00:00",
                                                                                "id": "23andme",
                                                                                "dbxref": {
                                                                                    "version": "2017-06-01T09:53:37+00:00",
                                                                                    "id": "http://identifiers.org/23andme"}}},
                                               'is_associated': True,
                                               'resource_score': {'type': 'pvalue', 'method': {
                                                   "description": "pvalue for the phenotype to snp association."},
                                                                  "value": float(phewas_dict['adj_pvalue'])},
                                               'date_asserted': i.strftime('%Y-%m-%d %H:%M:%S'),
                                               'evidence_codes': ['http://identifiers.org/eco/GWAS',
                                                                  'http://purl.obolibrary.org/obo/ECO_0000205'],
                                               }
                evidence['gene2variant'] = {
                    'provenance_type': {"expert": {"status": True, "statement": "Primary submitter of data"},
                                        "database": {"version": "2017-06-01T09:53:37+00:00", "id": "23andme",
                                                     "dbxref": {"version": "2017-06-01T09:53:37+00:00",
                                                                "id": "http://identifiers.org/23andme"}}},
                    'is_associated': True, 'date_asserted': i.strftime('%Y-%m-%d %H:%M:%S'),
                    'evidence_codes': ["http://identifiers.org/eco/cttv_mapping_pipeline",
                                       "http://purl.obolibrary.org/obo/ECO_0000205"],
                    'functional_consequence': 'http://purl.obolibrary.org/obo/SO_0001632'}
                phewas_evidence['evidence'] = evidence
                phewas_evidences.append(phewas_evidence)

        return phewas_evidences


class FileWriterProcess(RedisQueueWorkerProcess):
    def __init__(self, queue_in,
                 redis_path):
        super(FileWriterProcess, self).__init__(queue_in, redis_path)
        self.output_file = open('output_bestpvalue/23andme_evidence.json','a')

    def process(self, data):
        evidences = data
        if evidences:
            for evidence in evidences:
                self.output_file.write(json.dumps(evidence))
                self.output_file.write('\n')




def main():
    mongo_processor = MongoDataManager()
    mongo_processor.setup()
    mongo_processor.process()


if __name__ == "__main__":
    main()
