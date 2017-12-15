import logging
import unittest
from modules.GenomicsEnglandPanelApp import GE

logger = logging.getLogger(__name__)


class GETestCase(unittest.TestCase):


    def test_get_omim_mappings(self):
        ge = GE()
        lc = ge.get_omim_to_efo_mappings()
        print(lc)
        self.assertTrue(lc > 0)

    def test_get_opentargets_zooma_to_efo_mappings(self):
        ge = GE()
        lc = ge.get_opentargets_zooma_to_efo_mappings()
        print(lc)
        self.assertTrue(lc >0)

if __name__ == '__main__':
    unittest.main()