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
  ** If EFO mapping found
    create a 'self.map_strings' with 'panel_name, gene_symbol, levelOfConfidence, phenotypes, efo_id'
    i.e Neutropenia, severe congenital 3, autosomal recessive, 610738
        OMIM  130650
        #120330:Papillorenal syndrome
  ** If NO EFO mapping found try OntoMA
    i.e Verheij syndrome,  615583
        => Ontoma: http://www.ebi.ac.uk/efo/EFO_0000699,
* If NO OMIM id found, try OntoMa
    i.e. Diffuse keratoderma with knuckle pads
        => Ontoma: 'http://www.orpha.net/ORDO/Orphanet_2698'

* The OntoMa mapping process
  ** https://github.com/opentargets/OnToma/blob/master/ontoma/interface.py#L365
  ** Searches for a matching EFO code for a given phenotype/disease string
     operations roughly ordered from least expensive to most expensive
     and also from most authorative to least authorative
        1. EFO OBO lookup
        2. Zooma mappings lookup
        3. Zooma API high confidence lookup
        4. OLS API EFO lookup - exact match
        --- below this line we might not have a term in the platform ---
        5. HP OBO lookup
        6. OLS API HP lookup - exact match
        7. OLS API EFO lookup - not exact
        (8. ?Zooma medium)

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

        self.ontoma = OnToma(exclude=['zooma'])
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

                            ''' These should have OMIM Id
                              610219
                              Orofacial cleft 6, 608864
                              Nephrotic syndrome 14	617575
                              [Hair morphology 1, hair thickness], 612630 -3
                            '''

                            '''
                              Need to be strictly 6 digits &
                              without "Orpha|ORPHA|HP|PMID" within
                            '''
                            # 7702/241/8925/14972
                            search_omim  = re.search('[^0-9]*(\d{6})[^0-9]*', item) #7928/212/9363/15575
                            search_omim2 = re.search('(Orpha|ORPHA|HP|PMID)', item)

                            if search_omim and search_omim2 is None:
                                omim_id = search_omim.group(1)
#                                >> > t.find_term('615877', code='OMIM')
#                                >> > t.find_term('notadisease') is None
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
                                        print("\tRESULTS: %s OMIM id found: %s => Curated EFO id: %s " %(item, omim_id, efo_id))
                                        self.map_strings = "%s\t%s\t%s\t%s\t%s" % (panel_name, row['GeneSymbol'], row['LevelOfConfidence'], item, efo_id)
                                else:
                                    '''
                                     Verheij syndrome, 615583
                                     Sifrim-Hitz-Weiss syndrome  617159
                                     Stromme syndrome 	243605
                                    '''
                                    f += 1
                                    #print("\tRESULTS: %s OMIM id found: %s => NO EFO id " % (item, omim_id))

                                    if self.ontoma.find_term(item):
                                        phenotype_label = self.ontoma.find_term(item, suggest=True, verbose=True)
                                        mapping_type = "OMIM automatic mapping via Ontoma"
                                        print("\tRESULTS: OMIM id found: %s => Ontoma: %s, No Curated EFO id, Action %s" %(item, phenotype_label['term'], phenotype_label['action']))

                                    else:
                                        print("%s => ZOOMA" %(item))
                            else:
                                ''' These should NOT have OMIM Id
                                 PMID: 26204423
                                 Peeling skin HP:0040189
                                 OrphaNet ORPHA284984
                                 ORPHA169147
                                 ORPHA:220493  Joubert syndrome  with ocular defect
                                 Orphanet:231178
                                 Movement disorder, autonomic dysfunction, developmental delay, behavioural difficulties
                                 Brain Dopamine Serotonin Vesicular Transport Disease (Other disorders of neurotransmitter metabolism)
                                 Vörner type palmoplantar keratoderma
                                '''
                                d += 1
                                print("\tRESULTS: %s NO OMIM id found" % (item))

                                #if self.ontoma.find_efo(item):
                                #    phenotype_label = self.ontoma.find_efo(item)
                                #    mapping_type = "Automatic mapping via Ontoma"
                                #    print("\tRESULTS: NO OMIM id found: %s => Ontoma: %s" % (item, phenotype_label))
                                #else:
                                #    print("%s => ZOOMA" % (item))


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
            r = requests.get('https://panelapp.genomicsengland.co.uk/WebServices/list_panels', params={})
            results = r.json()

            for row in results['result']:
                yield (row['Name'], row['Panel_Id'], row["CurrentVersion"],
                       row["DiseaseGroup"], row["DiseaseSubGroup"])
        except requests.exceptions.HTTPError as e:
            return "Error: " + str(e)





