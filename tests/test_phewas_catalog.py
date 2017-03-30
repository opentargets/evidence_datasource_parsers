
import logging
import unittest
from modules.csv_parser import PhewasProcessor


logger = logging.getLogger(__name__)


class LiteratureNLPTestCase(unittest.TestCase):


    def test_efo_match(self):
        phewas_processor = PhewasProcessor()
        phewas_processor.setup()


        '''direct match'''
        matched_efos = phewas_processor.find_efo('Dementias', '290.1')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn('EFO:0003862',matched_efo_ids)

        '''synonym match'''
        matched_efos = phewas_processor.find_efo('Prostate cancer', '185')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn('EFO:0001663',matched_efo_ids )

        '''icd9 match'''
        matched_efos = phewas_processor.find_efo('Other dermatoses', '702')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn('EFO:0002496', matched_efo_ids )

        matched_efos = phewas_processor.find_efo('Alzheimer\'s disease', '290.11')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn( 'EFO:0000249',matched_efo_ids)

        '''hpo ontology match'''
        matched_efos = phewas_processor.find_efo('Jaundice', '573.5')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn( 'HP:0000952',matched_efo_ids )

        '''no match'''
        matched_efos = phewas_processor.find_efo('Iron metabolism disorder', '275.1')

        self.assertEqual(len(matched_efos),0 )


    def test_gene_match(self):
        phewas_processor = PhewasProcessor()
        phewas_processor.setup()
        matched_ensembl_id = phewas_processor.genes.get('TOMM40')
        self.assertEqual(matched_ensembl_id, 'ENSG00000130204' )

    def test_phewas_processor(self):
        phewas_processor = PhewasProcessor()
        phewas_processor.setup()
        phewas_processor.convert_phewas_catalog_evidence_json()






