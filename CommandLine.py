import argparse
import sys

from  modules.PheWAS import PhewasProcessor
from modules.Gene2Phenotype import G2P
from modules.GenomicsEnglandPanelApp import GE
from modules.MouseModels import Phenodigm
from modules.IntOGen import IntOGen
from modules.SLAPEnrich import SLAPEnrich


from settings import Config

def main():

    parser = argparse.ArgumentParser(description='Open Targets evidence generator')

    parser.add_argument("--phewas", dest='phewas',
                        help="process phewas data and generate evidences for open targets pipeline",
                        action="append_const",const=str)
    parser.add_argument("--genomicsengland", dest='genomicsengland',
                        help="process genomics england data and generate evidences for open targets pipeline",
                        action="append_const",const=str)
    parser.add_argument("--intogen", dest='intogen',
                        help="process IntoGen data and generate evidences for open targets pipeline",
                        action="append_const", const=str)
    parser.add_argument("--gene2phenotype", dest='gene2phenotype',
                        help="process phewas data and generate evidences for open targets pipeline",
                        action="append_const", const=str)
    parser.add_argument("--phenodigm", dest='phenodigm',
                        help="process phenodigm data and generate evidences for open targets pipeline",
                        action="append_const", const=str)
    parser.add_argument("--slapenrich", dest='slapenrich',
                        help="process slapenrich data and generate evidences for open targets pipeline",
                        action="append_const", const=str)
    parser.add_argument("--schema-version", dest='schema_version',
                        help="set the schema version",
                        action='store', default=Config.VALIDATED_AGAINST_SCHEMA_VERSION)
    args = parser.parse_args()

    if args.phewas:
        phewas_processor = PhewasProcessor(schema_version = args.schema_version)
        phewas_processor.setup()
        phewas_processor.convert_phewas_catalog_evidence_json()
    elif args.genomicsengland:
        GE().process_all()
    elif args.intogen:
        IntOGen().process_intogen()
    elif args.gene2phenotype:
        G2P().process_g2p()
    elif args.phenodigm:
        Phenodigm().generate_evidence()
    elif args.slapenrich:
        SLAPEnrich().process_slapenrich()



if __name__ == '__main__':
    sys.exit(main())
