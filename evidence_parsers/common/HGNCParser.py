import ujson as json
import requests
from tqdm import tqdm
from settings import Config
import logging

logger = logging.getLogger(__name__)

class Gene(object):
    def __init__(self, id=None):

        self.id = id
        self.hgnc_id = None
        self.approved_symbol = ""
        self.approved_name = ""
        self.status = ""
        self.locus_group = ""
        self.previous_symbols = []
        self.previous_names = []
        self.symbol_synonyms = []
        self.name_synonyms = []
        self.chromosome = ""

        self.ensembl_gene_id = ""

        self.gene_family_tag = ""
        self.gene_family_description = ""

        self.ensembl_external_name = ""
        self.ensembl_gene_version = None

        self.ensembl_release = None

    def load_hgnc_data_from_json(self, data):

        if 'ensembl_gene_id' in data:
            self.ensembl_gene_id = data['ensembl_gene_id']
            if not self.ensembl_gene_id:
                self.ensembl_gene_id = data['ensembl_id_supplied_by_ensembl']
            if 'hgnc_id' in data:
                self.hgnc_id = data['hgnc_id']
            if 'symbol' in data:
                self.approved_symbol = data['symbol']
            if 'name' in data:
                self.approved_name = data['name']
            if 'status' in data:
                self.status = data['status']
            if 'locus_group' in data:
                self.locus_group = data['locus_group']
            if 'prev_symbols' in data:
                self.previous_symbols = data['prev_symbols']
            if 'prev_names' in data:
                self.previous_names = data['prev_names']
            if 'alias_symbol' in data:
                self.symbol_synonyms.extend( data['alias_symbol'])
            if 'alias_name' in data:
                self.name_synonyms = data['alias_name']

            if 'gene_family_tag' in data:
                self.gene_family_tag = data['gene_family_tag']
            if 'gene_family_description' in data:
                self.gene_family_description = data['gene_family_description']


            if 'pubmed_id' in data:
                self.pubmed_ids = data['pubmed_id']


class GeneParser(object):
    ''' Parses the HGNC data and allows to lookup ENSGIDs from
    standard gene symbol.

    >>> gene_parser = GeneParser()
    >>> gene_parser._get_hgnc_data_from_json()
    >>> gene_parser.genes['BRAF']
    'ENSG00000157764'
    >>> gene_parser.genes['TOMM40']
    'ENSG00000130204'
    '''
    def __init__(self):
        self.genes = dict()

    def _get_hgnc_data_from_json(self):

        r = requests.get(Config.GENES_HGNC)
        data = r.json()

        for row in tqdm(data['response']['docs'],
                desc='Downloading HGNC genes from json response',
                unit='genes'):

            ensembl_gene_id = ''

            if 'ensembl_gene_id' in row:
                ensembl_gene_id = row['ensembl_gene_id']
                if not ensembl_gene_id:
                    ensembl_gene_id = data['ensembl_id_supplied_by_ensembl']
            self.genes[row['symbol']] = ensembl_gene_id

            if 'prev_symbol' in row:
                # to handle obsolete gene symbols like EFCAB4B
                for prev_symbol in row['prev_symbol']:
                    self.genes[prev_symbol] = ensembl_gene_id
        logger.info('All HGNC genes parsed')