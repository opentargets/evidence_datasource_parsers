import json
import logging
import os

import uuid
from redislite import Redis
from datetime import datetime
import multiprocessing
from common.Redis import RedisQueueWorkerProcess,RedisQueue
from common.EFOData import OBOParser


TEMP_DIR = os.path.join(os.path.sep, 'tmp')
REDISLITE_DB_PATH = os.path.join(TEMP_DIR, 'nferx_redislite.rdb')
NO_OF_WORKERS = 8
MAX_CHUNKS =100
UNIQUE_RUN_ID = str(uuid.uuid4()).replace('-', '')[:16]


class NferxManager(object):

    def __init__(self,schema_version ):

        self.schema_version = schema_version

    def process(self, file):

        obo_parser = OBOParser('../resources/efo.obo')
        obsolete_efos = obo_parser.get_obsolete_efos()

        output_file = open(file, 'a')
        logging.basicConfig(filename='output_nferx.log',
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                            level=logging.INFO)

        r_server = Redis(dbfilename=str(REDISLITE_DB_PATH), serverconfig={'save': [], 'maxclients': 10000})


        no_of_workers = NO_OF_WORKERS or multiprocessing.cpu_count()


        transformer_q = RedisQueue(queue_id=UNIQUE_RUN_ID + '|data_transformer_q',
                             max_size=MAX_CHUNKS * no_of_workers,
                             job_timeout=120)
        writer_q = RedisQueue(queue_id=UNIQUE_RUN_ID + '|writer_q',
                             max_size=MAX_CHUNKS * no_of_workers,
                             job_timeout=120)




        transformers = [DataTransformerProcess(
                                      transformer_q,
                                      r_server.db,
                                      writer_q,
                                      obsolete_efos,
                                      self.schema_version
                                      ) for i in range(no_of_workers)]

        for t in transformers:
            t.start()

        writer = FileWriterProcess(output_file,writer_q,r_server.db)
        writer.start()

        logging.info('Getting distinct genes and phenotypes')

        json_file = '../resources/nferx_release5.json'
        with open(json_file) as input_file:
            for line in input_file:
                nferx_dict = json.loads(line)
                transformer_q.put(nferx_dict,r_server)

        transformer_q.set_submission_finished(r_server=r_server)


        for a in transformers:
                a.join()
        writer.join()
        output_file.flush()


class DataTransformerProcess(RedisQueueWorkerProcess):

    def __init__(self,

                 queue_in,
                 redis_path,
                 queue_out,
                 obsolete_efos,
                 schema_version):
        super(DataTransformerProcess, self).__init__(queue_in,redis_path,queue_out)
        self.obsolete_efos = obsolete_efos
        self.schema_version = schema_version

    def process(self, data):
        nferx_dict = data
        nferx_evidence = dict()

        cosine_dist = float(nferx_dict['nfer_data']['cosine_distance'])
        document_no = nferx_dict['nfer_data']['document_cooccurrence']
        if document_no == 0 or cosine_dist == 0 :
            return None
        efo_id = nferx_dict['EFO_ID']
        if efo_id in self.obsolete_efos:
            new_efo = self.obsolete_efos[efo_id]
            if new_efo:
                print 'Obsolete efo- {} New EFO - {}'.format(efo_id, new_efo)
                efo_id = new_efo
            else:
                return None
        if efo_id.startswith('EFO'):
            nferx_evidence['disease'] = {'id': 'http://www.ebi.ac.uk/efo/'+efo_id}
        elif efo_id.startswith('HP') or efo_id.startswith('MP') :
            nferx_evidence['disease'] = {'id': 'http://purl.obolibrary.org/obo/' + efo_id}
        elif efo_id.startswith('Orphanet') :
            nferx_evidence['disease'] = {'id': 'http://www.orpha.net/ORDO/' + efo_id}
        elif efo_id.startswith('NCBITaxon') :
            nferx_evidence['disease'] = {'id': 'http://purl.obolibrary.org/obo/' + efo_id}
        else:
            return None
        # evidence_stats['document_count']= str(document_no)
        # evidence_stats['cosine_distance'] = str(cosine_dist)
        #
        #

        nferx_evidence['target'] = {"activity": "http://identifiers.org/cttv.activity/predicted_damaging",
                                    "id": "http://identifiers.org/ensembl/{}".format(nferx_dict['ENSGID']),
                                    "target_type":"http:/identifiers.org/cttv.target/gene_evidence"}
        nferx_evidence['validated_against_schema_version'] = self.schema_version
        nferx_evidence["access_level"] = "public"
        nferx_evidence["sourceID"] = "nferx"
        nferx_evidence['type'] = 'literature'

        nferx_evidence['unique_association_fields'] = {'disease_count': nferx_dict['nfer_data']['disease_count'],
                                                       'gene_count': nferx_dict['nfer_data']['gene_count'],
                                                       'serial_id': str(nferx_dict['nfer_data']['serial_id']),
                                                       'document_count': str(document_no),
                                                       'link_url': nferx_dict['nfer_data']['link_url']}
        i = datetime.now()
        evidence = {'unique_experiment_reference': 'http://europepmc.org/abstract/MED/0',
                    'provenance_type': {"database": {"version": "2017-07-01T09:53:37+00:00", "id": "Nferx"}},
                    'is_associated': True,
                    'resource_score': {'type': 'pvalue',
                                       'method': {"description": "custom scoring - cosine distance between vectors."},
                                       "value": cosine_dist},
                    'date_asserted': i.strftime('%Y-%m-%dT%H:%M:%S+00:00'),
                    'evidence_codes': ['http://www.targetvalidation.org/evidence/literature_mining',
                                       'http://purl.obolibrary.org/obo/ECO_0000213'],
                    'literature_ref': {'lit_id': 'http://europepmc.org/abstract/MED/1'},
                    }
        nferx_evidence['evidence'] = evidence
        return nferx_evidence

class FileWriterProcess(RedisQueueWorkerProcess):
    def __init__(self,
                 output_file,
                 queue_in,
                 redis_path):
        super(FileWriterProcess, self).__init__(queue_in, redis_path)
        self.output_file = output_file


    def process(self, data):
        nferx_evidence = data
        if nferx_evidence:
            json.dump(nferx_evidence, self.output_file)
            self.output_file.write('\n')




def main():
    nferx_manager = NferxManager('1.2.7')
    nferx_manager.process('../output/nferx.json')


if __name__ == "__main__":
    main()
