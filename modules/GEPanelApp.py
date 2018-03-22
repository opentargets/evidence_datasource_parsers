import re
import csv
import requests
from ontoma import OnToma
from settings import Config
from collections import OrderedDict
from common.HGNCParser import GeneParser
from common.GCSUtils import GCSBucketManager

import logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright 2014-2018, Open Targets"
__credits__   = ["Gautier Koscielny", "ChuangKee Ong"]
__license__   = "Apache 2.0"
__version__   = "1.2.8"
__maintainer__= "ChuangKee Ong"
__email__     = ["data@opentargets.org"]
__status__    = "Production"

##TODO: 1) implement logging


'''
* Get 'Phenotypes' list for each 'HighEvidence' gene from each GE Panel.
* For each phenotype 1st check that it has an OMIM id [6 digits].
* If OMIM id found, try map to EFO using mapping file
  OMIM_EFO_MAP => https://github.com/opentargets/OnToma/blob/master/ontoma/constants.py
  * If found EFO mapping create a 'self.map_strings' with 'panel_name, gene_symbol, levelOfConfidence, phenotypes, efo_id'
    i.e Neutropenia, severe congenital 3, autosomal recessive, 610738
        OMIM  130650
        #120330:Papillorenal syndrome
  * If NO EFO mapping found
    i.e Verheij syndrome, 615583
* If NO OMIM id found .....

'''

#class GEPanelApp(GCSBucketManager):
class GEPanelApp():
    #ge_object = GE()
    #ge_object.execute_ge_request()
    #ge_object.use_zooma()
    #ge_object.process_panel_app_file()
    #ge_object.write_evidence_strings(Config.GE_EVIDENCE_FILENAME)

    def __init__(self):
        super(GEPanelApp, self).__init__()
        self.panel_app_id_map = OrderedDict()
        self.panel_app_info = list()
        self.map_strings = OrderedDict()

    def process_ge(self):

        self.ontoma = OnToma()
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()

        self.omim_to_efo_map = self.ontoma._omim_to_efo
        self.symbol_to_ensembl_map = gene_parser.genes

        print("Parsed %s OMIM to EFO mapping " % len(self.omim_to_efo_map))
        print("Parsed %s ZOOMA to EFO mapping " % len(self.ontoma._zooma_to_efo_map))
        print("Parsed %s Name to EFO mapping " % len(self.ontoma.name_to_efo))
        print("Parsed %s Name to HP mapping " % len(self.ontoma.name_to_hp))
        self.get_panel_code_mapping()
        self.process_ge_panel()

    def get_panel_code_mapping(self, filename=Config.GE_PANEL_MAPPING_FILENAME):

        print('Mapping old GE panel id to new panel id, using : ' + filename)

        #with open(self.get_gcs_filename(filename), mode='rt') as fh:
        with open(filename, mode='rt') as fh:
            reader = csv.reader(fh, delimiter=',', quotechar='"')
            next(reader, None) # skip header

            for row in reader:
                (old_panel_id, new_panel_id) = row
                print ("panel id :%s => is now %s"%(old_panel_id, new_panel_id))
                self.panel_app_id_map[old_panel_id] = new_panel_id

    def process_ge_panel(self):
        '''
         Create panel app info list and phenotype set
         :return: Unique phenotype list
        '''
        print("Retrieving GE PanelApp data...")
        url = 'https://panelapp.genomicsengland.co.uk/WebServices/search_genes/all/'
        phenotype_list = []
        c = 0
        d = 0
        e = 0
        f = 0

        for panel_name, panel_id, panel_version, panel_diseasegroup, panel_diseasesubgroup in self.get_panel_list():
            #c += 1
            print("Reading panel : %s %s %s %s %s" % (panel_name, panel_id, panel_version, panel_diseasegroup, panel_diseasesubgroup))
            r = requests.get(url, params={"panel_name": panel_name}, timeout=30)
            results = r.json()

            ''' process genes in a panel '''
            for row in results['results']:
                ensembl_iri = None

                ensembl_gene_id = self.symbol_to_ensembl_map.get(row['GeneSymbol'])

                if not ensembl_gene_id:
                    print("%s is not found in Ensembl" % row['GeneSymbol'])
                    continue

                ensembl_iri = "http://identifiers.org/ensembl/" + ensembl_gene_id

                if ensembl_iri and row['EnsembleGeneIds'] and row['Phenotypes'] \
                        and row['LevelOfConfidence'] == 'HighEvidence':

                    for item in row['Phenotypes']:
                        e +=1
                        item = item.rstrip().lstrip().rstrip("?")
                        if len(item) > 0:
                            print("High confidence Gene %s : %s has phenotypes: %s" % (row['GeneSymbol'], ensembl_iri, item))
                            '''
                             "search_omim"
                                Neutropenia, severe congenital 3, autosomal recessive, 610738
                                OMIM  130650
                             "search_omim2"
                                #120330:Papillorenal syndrome
                             "search_omim3"
                                Leukemia, acute myeloid, 601626(1)
                                Heimler Syndrome 2, 616617 (includes amelogenesis imperfecta)
                            '''
                            #search_omim  = re.search('(\d{6})$', item)
                            #search_omim2 = re.search('^#(\d{6})', item)
                            #search_omim3 = re.search('(\d{6})', item)
                            search_omim = re.search('^#(\d{6})|(\d{6})|(\d{6})$',item)

                            if search_omim:
                                omim_id = search_omim.group(1)

                            #if search_omim or search_omim2 or search_omim3:
                            #    if search_omim:
                            #        omim_id = search_omim.group(1)
                            #    elif search_omim2:
                            #        omim_id = search_omim2.group(1)
                            #    elif search_omim3:
                            #        omim_id = search_omim3.group(1)

                                #if self.ontoma.find_efo(omim_id, code='OMIM'):
                                if omim_id in self.omim_to_efo_map:
                                    for efo_id in self.omim_to_efo_map[omim_id]:
                                        self.panel_app_info.append([panel_name,
                                                                    panel_id,
                                                                    panel_version,
                                                                    panel_diseasegroup,
                                                                    panel_diseasesubgroup,
                                                                    row['GeneSymbol'],
                                                                    row['EnsembleGeneIds'][0],
                                                                    ensembl_iri,
                                                                    row['LevelOfConfidence'],
                                                                    item,
                                                                    omim_id,
                                                                    row['Publications'],
                                                                    row['Evidences'],
                                                                    [omim_id],
                                                                    "OMIM MAPPING"
                                                                    ])
                                        c += 1
                                        print("\tOMIM id found: " + omim_id + " mapped to EFO id: " + efo_id)
                                        self.map_strings = "%s\t%s\t%s\t%s\t%s" % (panel_name, row['GeneSymbol'], row['LevelOfConfidence'], item, efo_id)
                                else:
                                    '''
                                     Verheij syndrome, 615583
                                     Sifrim-Hitz-Weiss syndrome  617159
                                    '''
                                    f += 1
                                    print("\tOMIM id found: %s with NO EFO mapping" % omim_id)

                            else:
                                d += 1
                                print ("NO OMIM id found %s" % item)##


#        for item in self.panel_app_info:
#            print(item)

        print("total OMIM_id found with EFO mapping " + str(c))
        print("total OMIM_id found withOUT EFO mapping " + str(f))
        print("total phenotypes without OMIM_id " + str(d))
        print("total phenotypes " + str(e))

        return

    @staticmethod
    def get_panel_list():
        '''
         Get list of GE panels
         :return: tuple of list of panel name and panel id's
        '''
        try:
            r = requests.get(Config.GE_PANEL_LIST, params={})
            results = r.json()

            for row in results['result']:
                yield (row['Name'], row['Panel_Id'], row["CurrentVersion"],
                       row["DiseaseGroup"], row["DiseaseSubGroup"])
        except requests.exceptions.HTTPError as e:
            return "Error: " + str(e)





