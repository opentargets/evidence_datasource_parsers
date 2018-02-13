
import logging
import unittest
from modules.PheWAS import PhewasProcessor

logger = logging.getLogger(__name__)


class PhewasTestCase(unittest.TestCase):


    def test_efo_match(self):
        phewas_processor = PhewasProcessor()
        phewas_processor.setup()


        '''direct match'''
        matched_efos = phewas_processor.find_efo(b'Dementias', b'290.1')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn(b'EFO:0003862',matched_efo_ids)

        '''synonym match'''
        matched_efos = phewas_processor.find_efo(b'Prostate cancer', b'185')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn(b'EFO:0001663',matched_efo_ids )

        '''icd9 match'''
        matched_efos = phewas_processor.find_efo(b'Other dermatoses', b'702')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn(b'EFO:0002496', matched_efo_ids )

        matched_efos = phewas_processor.find_efo(b'Alzheimer\'s disease', b'290.11')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn( b'EFO:0000249',matched_efo_ids)

        '''hpo ontology match'''
        matched_efos = phewas_processor.find_efo(b'Jaundice', b'573.5')
        matched_efo_ids = [efo['id'] for efo in matched_efos]
        self.assertIn( b'HP:0000952',matched_efo_ids )

        '''no match'''
        matched_efos = phewas_processor.find_efo(b'Iron metabolism disorder', b'275.1')

        self.assertEqual(len(matched_efos),0 )


    def test_gene_match(self):
        phewas_processor = PhewasProcessor()
        phewas_processor.setup()
        matched_ensembl_id = phewas_processor.genes.get('TOMM40')
        self.assertEqual(matched_ensembl_id, 'ENSG00000130204' )

    #NOTE: this is not a test, but actually the way to run the Phewasprocessor
    # def test_phewas_processor(self):
    #     phewas_processor = PhewasProcessor()
    #     phewas_processor.setup()
    #     phewas_processor.convert_phewas_catalog_evidence_json()






