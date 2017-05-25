from common.Redis import RedisQueueWorkerProcess,RedisQueue
from pymongo import MongoClient
from bson import json_util
import multiprocessing
import uuid
from redislite import Redis
import os
import logging
import math
import json

NO_OF_WORKERS = 8
MAX_CHUNKS =100
UNIQUE_RUN_ID = str(uuid.uuid4()).replace('-', '')[:16]
TEMP_DIR = os.path.join(os.path.sep, 'tmp')
REDISLITE_DB_PATH = os.path.join(TEMP_DIR, 'phewas_redislite.rdb')

class MongoDataManager(object):


    def process(self):
        logging.basicConfig(filename='output_mongo.log',
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                            level=logging.INFO)

        r_server = Redis(dbfilename=str(REDISLITE_DB_PATH), serverconfig={'save': [], 'maxclients': 10000})
        mongo_client = MongoClient('host:port')
        db = mongo_client['xxx']
        db.authenticate('uname','pwd')
        gwas23andme = db['xxxx']

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

        transformers = [DataTransformerProcess(transformer_q,
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
        distinct_genes = list(gwas23andme.distinct("ingenes"))
        logging.info('Start the real deal here!!!! ')
        for phenotype in distinct_phenotypes:
            logging.info('Processing phenotype {} '.format(phenotype))
            logging.info('Total docs for phenotype {} are {}'.format(phenotype, gwas23andme.count({"source" : phenotype})))

            mongo_q.put((phenotype, 'test'), r_server)
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
            mongo_client = MongoClient('host:port')
            db = mongo_client['xxx']
            db.authenticate('uname', 'pwd')
            gwas23andme = db['xxxx']

        phenotype, gene  = data

        # docs_cursor = self.gwas23andme.find({'source': phenotype, 'best_pvalue': {'$gt': 0.5}, 'ingenes' : gene},
        #                          {"source": 1, "adj_pvalue": 1, "im_num_1": 1, "im_num_0": 1, "alleles": 1,
        #                                 "effect": 1, "assay_name": 1, 'ingenes': 1, 'dtype': 1, "best_effect": 1 })
        docs_cursor = self.gwas23andme.find({"source": phenotype, 'best_pvalue': {'$exists': True, '$lte': float(self.val)}},
                                            {"source": 1, "adj_pvalue": 1, "im_num_1": 1, "im_num_0": 1, "alleles": 1,
                                         "effect": 1, "assay_name": 1, 'ingenes': 1, 'dtype': 1, "best_effect": 1 })

        gwas_coll = list(docs_cursor)
        logging.info('Found {} documents for phenotype {} meeting significant threshold'.format(len(gwas_coll),phenotype))
        # filename = 'output/23andme_' + phenotype + '.json'
        #return (filename, json_util.dumps(gwas_coll))

        return ( gwas_coll)

class DataTransformerProcess(RedisQueueWorkerProcess):

    def __init__(self, queue_in,
                 redis_path,
                 queue_out):
        super(DataTransformerProcess, self).__init__(queue_in,redis_path,queue_out)


    def process(self, data):
        # filename,  other_info = data
        other_info = data
        for gwas_dict in other_info:
            if gwas_dict.get('dtype') == 'binary':
                if gwas_dict.get('best_effect'):
                    try:
                        gwas_dict['odds_ratio'] = math.exp(gwas_dict.get('best_effect'))
                    except OverflowError:
                        logging.info('Overflow Error for document  {} with best_effect {} '.format(gwas_dict['_id'],gwas_dict.get('best_effect')))
            gwas_id_str = str(gwas_dict['_id'])
            gwas_dict['_id'] = gwas_id_str
        return (json.dumps(other_info))
        #return (filename, json_util.dumps(other_info))



class FileWriterProcess(RedisQueueWorkerProcess):
    def __init__(self, queue_in,
                 redis_path):
        super(FileWriterProcess, self).__init__(queue_in, redis_path)
        self.output_file = open('output_bestpvalue/23andme.json','a')

    def process(self, data):
        other_info = data
        if other_info:
            self.output_file.write(other_info)




def main():
    mongo_processor = MongoDataManager()
    mongo_processor.process()

if __name__ == "__main__":
    main()
