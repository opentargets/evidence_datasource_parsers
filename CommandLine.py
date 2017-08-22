import argparse
import sys

from  modules.csv_parser import PhewasProcessor
# from modules.mongo_multiprocessing import MongoDataManager
# from modules.nferx_parser import NferxManager
def main():

    parser = argparse.ArgumentParser(description='Open Targets evidence generator')

    parser.add_argument("--phewas", dest='phewas',
                        help="process phewas data and generate evidences for open targets pipeline",
                        action="append_const",const=str)
    parser.add_argument("--23andme", dest='biogen_23andme',
                        help="process 23andme data and generate evidences for open targets pipeline",
                        action="append_const", const=str)
    parser.add_argument("--nferx", dest='biogen_nferx',
                        help="process nferx data and generate evidences for open targets pipeline",
                        action="append_const", const=str)
    parser.add_argument("--schema-version", dest='schema_version',
                        help="set the schema version",
                        action='store', default='1.2.6')
    args = parser.parse_args()

    if args.phewas :
        phewas_processor = PhewasProcessor(schema_version = args.schema_version)
        phewas_processor.setup()
        phewas_processor.convert_phewas_catalog_evidence_json()
    # if args.biogen_23andme :
    #     mongo_processor = MongoDataManager()
    #     mongo_processor.setup()
    #     mongo_processor.process()
    # if args.biogen_nferx :
    #     nferx_manager = NferxManager()
    #     nferx_manager.process('../nferx/nferx.json')


if __name__ == '__main__':
    sys.exit(main())
