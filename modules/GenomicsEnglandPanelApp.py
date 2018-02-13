from collections import OrderedDict
import opentargets.model.core as opentargets
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.linkout as evidence_linkout
import opentargets.model.evidence.association_score as association_score
from ontologyutils.rdf_utils import OntologyClassReader
from common.RareDiseasesUtils import RareDiseaseMapper
from common.GCSUtils import GCSBucketManager
import logging
import datetime
import csv
import re
import requests
from requests.packages.urllib3.poolmanager import PoolManager
import urllib.request
import json
import hashlib
from common.HGNCParser import GeneParser
import urllib.parse
from settings import Config
import ssl

logger = logging.getLogger(__name__)

TEST_SAMPLE=True

class MyAdapter(requests.adapters.HTTPAdapter):
    def init_poolmanager(self, connections, maxsize, block=False):
        self.poolmanager = PoolManager(
            num_pools=connections,
            maxsize=maxsize,
            block=block,
            ssl_version=ssl.PROTOCOL_TLSv1_1,
        )
        print("poolmanager set")

class GE(RareDiseaseMapper, GCSBucketManager):

    def __init__(self):

        super(GE, self).__init__()
        self.lookup_data = None
        self.hashkeys = OrderedDict()
        self.efo = OntologyClassReader()
        self.hpo = OntologyClassReader()
        self.hpo_labels = OrderedDict()
        self.efo_labels = OrderedDict()
        self.ols_synonyms = OrderedDict()
        self.not_ols_synonyms = set()
        self.panel_app_info = list()
        self.high_confidence_mappings = OrderedDict()
        self.genes = None
        self.other_zooma_mappings = OrderedDict()
        self.phenotype_set = set()
        self.evidence_strings = list()
        self.map_omim = OrderedDict()
        self.fh_zooma_high = None
        self.fh_zooma_low = None
        self._logger = logging.getLogger(__name__)
        self._logger.warning("GE init")
        self.map_strings = OrderedDict()

    def process_all(self):
        self._logger.warning("Process all")

        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

        # get EFO and HPO
        self._logger.warning("Get EFO")
        self.efo.load_efo_classes()
        for k, v in self.efo.current_classes.items():
            self.efo_labels[v.lower()] = k

        self._logger.warning("Get HPO")
        self.hpo.load_hpo_classes()
        for k, v in self.hpo.current_classes.items():
            self.hpo_labels[v.lower()] = k

        self.get_omim_to_efo_mappings()
        self.get_opentargets_zooma_to_efo_mappings()

        self.execute_ge_request(sample=TEST_SAMPLE)
        self.use_zooma()
        self.process_panel_app_file(sample=TEST_SAMPLE)
        self.write_evidence_strings(Config.GE_EVIDENCE_FILENAME)



    @staticmethod
    def request_to_panel_app():
        '''
        Makes a request to panel app to get the list of all panels
        :return: tuple of list of panel name and panel id's
        '''
        #requests_cache.install_cache('GE_results_cache_Feb', backend='sqlite', expire_after=3000000)
        #s = requests.session()
        #s.get('https://panelapp.genomicsengland.co.uk/', MyAdapter())
        r = requests.get('https://panelapp.genomicsengland.co.uk/WebServices/list_panels', params={})
        results = r.json()

        for item in results['result']:
            yield (item['Name'], item['Panel_Id'])

    def execute_ge_request(self, sample=TEST_SAMPLE):
        '''
        Create panel app info list and phenotype set
        :return: Unique phenotype list
        '''
        self._logger.warning("execute_ge_request...")
        phenotype_list = []
        nb_panels = 0
        for panel_name, panel_id in self.request_to_panel_app():
            nb_panels +=1
            self._logger.warning("reading panel %s %s" % (panel_name, panel_id))
            url = 'https://panelapp.genomicsengland.co.uk/WebServices/search_genes/all/'
            r = requests.get(url, params={"panel_name": panel_name}, timeout=30)
            results = r.json()
            for item in results['results']:

                ensembl_iri = None

                ensembl_gene_id = self.genes.get(item['GeneSymbol'])
                if not ensembl_gene_id:
                    self._logger.error("%s is not found in Ensembl" % item['GeneSymbol'])
                    continue

                ensembl_iri = "http://identifiers.org/ensembl/" + ensembl_gene_id
                self._logger.warning("Gene %s %s" % (item['GeneSymbol'], ensembl_iri))

                if ensembl_iri and item['EnsembleGeneIds'] and item['Phenotypes'] and item['LevelOfConfidence'] == 'HighEvidence':
                    for element in item['Phenotypes']:

                        element = element.rstrip().lstrip().rstrip("?")
                        if len(element) > 0:

                            '''
                            First check whether it's an OMIM identifier
                            '''
                            match_omim = re.match('^(\d{6,})$', element)
                            match_omim2 = re.match('^\s*OMIM:?[#\s]*(\d{6,})$', element)
                            if match_omim or match_omim2:
                                if match_omim:
                                    omim_id = match_omim.groups()[0]
                                elif match_omim2:
                                    omim_id = match_omim2.groups()[0]
                                self._logger.info("Found OMIM ID: %s" % (omim_id))
                                if omim_id in self.omim_to_efo_map:
                                    self._logger.info("Maps to EFO")
                                    for mapping in self.omim_to_efo_map[omim_id]:
                                        disease_label = mapping['efo_label']
                                        disease_uri = mapping['efo_uri']
                                        self.panel_app_info.append([panel_name,
                                                                    panel_id,
                                                                    item['GeneSymbol'],
                                                                    item['EnsembleGeneIds'][0],
                                                                    ensembl_iri,
                                                                    item['LevelOfConfidence'],
                                                                    element,
                                                                    omim_id,
                                                                    item['Publications'],
                                                                    item['Evidences'],
                                                                    [omim_id],
                                                                    disease_uri,
                                                                    disease_label,
                                                                    "OMIM MAPPING"
                                                                    ])
                                self.map_strings = "%s\t%s\t%s\t%s\t%s\t%s"%(panel_name, item['GeneSymbol'], item['LevelOfConfidence'], element, disease_uri, disease_label)
                            else:

                                '''
                                if there is already an OMIM xref to EFO/Orphanet, no need to map
                                '''
                                disease_uri = None
                                disease_label = None
                                is_hpo = False
                                is_efo = False

                                omim_ids = []
                                phenotype_label = None
                                mapping_type = "EFO MAPPING"

                                match_hpo_id = re.match('^(.+)\s+(HP:\d+)$', element)
                                match_curly_brackets_omim = re.match('^{([^\}]+)},\s+(\d+)', element)
                                match_no_curly_brackets_omim = re.match('^(.+),\s+(\d{6,})$', element)
                                match_no_curly_brackets_omim_continued = re.match('^(.+),\s+(\d{6,})\s+.*$', element)
                                # Myopathy, early-onset, with fatal cardiomyopathy 611705
                                match_no_curly_brackets_no_comma_omim = re.match('^(.+)\s+(\d{6,})\s*$', element)
                                if element.lower() in self.efo_labels:
                                    disease_uri = self.efo_labels[element.lower()]
                                    disease_label = element
                                    phenotype_label = disease_label
                                    is_efo = True
                                    mapping_type = "EFO MAPPING"
                                elif element.lower() in self.hpo_labels:
                                    disease_uri = self.hpo_labels[element.lower()]
                                    disease_label = self.hpo.current_classes[disease_uri]
                                    phenotype_label = disease_label
                                    is_hpo = True
                                    mapping_type = "HPO MAPPING (in HPO labels)"
                                elif match_hpo_id:
                                    # Ichthyosis HP:64
                                    disease_label = match_hpo_id.groups()[0]
                                    phenotype_label = disease_label
                                    phenotype_id = match_hpo_id.groups()[1]
                                    disease_uri = "http://purl.obolibrary.org/obo/" + phenotype_id.replace(":", "_")
                                    if disease_uri in self.hpo.current_classes:
                                        disease_label = self.hpo.current_classes[disease_uri]
                                        is_hpo = True
                                    mapping_type = "HPO MAPPING (match HPO id)"
                                elif self.search_request_to_ols(query=element) is True:
                                    disease_uri = self.ols_synonyms[element]["iri"]
                                    disease_label = self.ols_synonyms[element]["label"]
                                    phenotype_label = disease_label
                                    is_hpo = True
                                    mapping_type = "HPO MAPPING (request to OLS)"
                                elif match_curly_brackets_omim:
                                    #[{Pancreatitis, idiopathic}, 167800]
                                    phenotype_label = match_curly_brackets_omim.groups()[0]
                                    omim_ids.append(match_curly_brackets_omim.groups()[1])
                                    mapping_type = "OMIM MAPPING"
                                elif match_no_curly_brackets_omim:
                                    #[{Pancreatitis, idiopathic}, 167800]
                                    phenotype_label = match_no_curly_brackets_omim.groups()[0]
                                    omim_ids.append(match_no_curly_brackets_omim.groups()[1])
                                    mapping_type = "OMIM MAPPING"
                                elif match_no_curly_brackets_omim_continued:
                                    #[{Pancreatitis, idiopathic}, 167800]
                                    phenotype_label = match_no_curly_brackets_omim_continued.groups()[0]
                                    omim_ids.append(match_no_curly_brackets_omim_continued.groups()[1])
                                    mapping_type = "OMIM MAPPING"
                                elif match_no_curly_brackets_no_comma_omim:
                                    #[{Pancreatitis, idiopathic}, 167800]
                                    phenotype_label = match_no_curly_brackets_no_comma_omim.groups()[0]
                                    omim_ids.append(match_no_curly_brackets_no_comma_omim.groups()[1])
                                    mapping_type = "OMIM MAPPING"
                                else:
                                    phenotype_label = None
                                    #if isinstance(element, str):
                                    phenotype_label = element.strip()
                                    #else:
                                    #    phenotype_label = element.decode('iso-8859-1').encode('utf-8').strip()
                                    phenotype_label = re.sub(r"\#", "", phenotype_label)
                                    phenotype_label = re.sub(r"\t", "", phenotype_label)
                                    omim_ids = re.findall(r"\d{5}",phenotype_label)
                                    phenotype_label = re.sub(r"\d{5}", "", phenotype_label)
                                    phenotype_label = re.sub(r"\{", "", phenotype_label)
                                    phenotype_label = re.sub(r"\}", "", phenotype_label)
                                    mapping_type = "GENERAL MAPPING"

                                self._logger.info("[%s] => [%s]" % (element, phenotype_label))



                                if omim_ids is None:
                                    omim_ids = []

                                self.map_omim[phenotype_label] = omim_ids



                                if not is_hpo and not is_efo and all(l not in self.omim_to_efo_map for l in omim_ids) and phenotype_label.lower() not in self.zooma_to_efo_map:
                                    self._logger.info("Unknown term '%s' with unknown OMIM ID(s): %s"%(phenotype_label, ";".join(omim_ids)))
                                    phenotype_list.append(phenotype_label)

                                    self.panel_app_info.append([panel_name,
                                                                panel_id,
                                                                item['GeneSymbol'],
                                                                item['EnsembleGeneIds'][0],
                                                                ensembl_iri,
                                                                item['LevelOfConfidence'],
                                                                element,
                                                                phenotype_label,
                                                                item['Publications'],
                                                                item['Evidences'],
                                                                omim_ids,
                                                                disease_uri,
                                                                disease_label,
                                                                mapping_type])

                                else:
                                    self._logger.info("THERE IS A MATCH")

                                    if is_hpo or is_efo:
                                        self.panel_app_info.append([panel_name,
                                                                    panel_id,
                                                                    item['GeneSymbol'],
                                                                    item['EnsembleGeneIds'][0],
                                                                    ensembl_iri,
                                                                    item['LevelOfConfidence'],
                                                                    element,
                                                                    phenotype_label,
                                                                    item['Publications'],
                                                                    item['Evidences'],
                                                                    omim_ids,
                                                                    disease_uri,
                                                                    disease_label,
                                                                    mapping_type])

                                    elif omim_ids and any(l  in self.omim_to_efo_map for l in omim_ids):
                                        for omim_id in omim_ids:
                                            if omim_id in self.omim_to_efo_map:
                                                for mapping in self.omim_to_efo_map[omim_id]:
                                                    disease_label = mapping['efo_label']
                                                    disease_uri = mapping['efo_uri']

                                                    self.panel_app_info.append([panel_name,
                                                                                panel_id,
                                                                                item['GeneSymbol'],
                                                                                item['EnsembleGeneIds'][0],
                                                                                ensembl_iri,
                                                                                item['LevelOfConfidence'],
                                                                                element,
                                                                                phenotype_label,
                                                                                item['Publications'],
                                                                                item['Evidences'],
                                                                                omim_ids,
                                                                                disease_uri,
                                                                                disease_label,
                                                                                mapping_type])

                                    elif phenotype_label.lower() in self.zooma_to_efo_map:
                                        for mapping in self.zooma_to_efo_map[phenotype_label.lower()]:

                                            disease_label = mapping['efo_label']
                                            disease_uri = mapping['efo_uri']

                                            self.panel_app_info.append([panel_name,
                                                                        panel_id,
                                                                        item['GeneSymbol'],
                                                                        item['EnsembleGeneIds'][0],
                                                                        ensembl_iri,
                                                                        item['LevelOfConfidence'],
                                                                        element,
                                                                        phenotype_label,
                                                                        item['Publications'],
                                                                        item['Evidences'],
                                                                        omim_ids,
                                                                        disease_uri,
                                                                        disease_label,
                                                                        mapping_type])




                                self.map_strings = "%s\t%s\t%s\t%s\t%s\t%s" % (
                                panel_name, item['GeneSymbol'], item['LevelOfConfidence'], element, disease_uri,
                                disease_label)

            self.phenotype_set = set(phenotype_list)

            if sample == True:
                print(json.dumps(self.panel_app_info, indent=2))
                break

        return self.phenotype_set


    def search_request_to_ols(self, query=None, ontology='efo', type='class'):
        in_ols = False
        if query in self.ols_synonyms:
            return True
        elif query in self.not_ols_synonyms:
            return False
        else:
            # curl 'http://www.ebi.ac.uk/ols/api/search?q=Myofascial%20Pain%20Syndromes&queryFields=synonym&exact=true&ontology=efo,ordo' -i -H 'Accept: application/json'
            url = 'http://www.ebi.ac.uk/ols/api/search?q=%s&queryFields=synonym&exact=true&ontology=efo,ordo,hpo'%(query)
            self._logger.info("Requesting '%s' as synonym to OLS"%(query))
            print("Requesting '%s' as synonym to OLS" %(query))
            r = requests.get(
                'http://www.ebi.ac.uk/ols/api/search',
                params={'q': query, 'queryFields': 'synonym', 'exact': 'true', 'ontology': 'efo,ordo,hpo'})
            print(r.text)
            results = r.json()
            if results["response"]["numFound"] > 0:
                in_ols = True
                docs = results["response"]["docs"]
                for doc in docs:
                    if doc["ontology_name"] == 'efo':
                        self.ols_synonyms[query] = { "iri" : doc["iri"], "label" : doc["label"], "ontology_name" : doc["ontology_name"] }
                        return True
                for doc in docs:
                    if doc["ontology_name"] == 'hpo':
                        self.ols_synonyms[query] = { "iri" : doc["iri"], "label" : doc["label"], "ontology_name" : doc["ontology_name"] }
                        return True
                for doc in docs:
                    if doc["ontology_name"] == 'ordo':
                        self.ols_synonyms[query] = { "iri" : doc["iri"], "label" : doc["label"], "ontology_name" : doc["ontology_name"] }
                        return True
            else:
                self.not_ols_synonyms.add(query)
        return in_ols

    def request_to_zooma(self, property_value=None):
        '''
        Make a request to Zooma to get correct phenotype mapping and disease label
        :param property_value: Phenotype name from Genomics England
        :return: High confidence mappings .  Writes the output in the input file
        see docs: http://www.ebi.ac.uk/spot/zooma/docs/api.html
        '''
        #requests_cache.install_cache('zooma_results_cache_jan', backend='sqlite', expire_after=3000000)
        self._logger.info("Requesting")
        r = requests.get('http://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate',
                             params={'propertyValue': property_value, 'propertyType': 'phenotype'})
        results = r.json()
        for item in results:
            if item['confidence'] == "HIGH":
                self.high_confidence_mappings[property_value] = {
                    'uri': item['_links']['olslinks'][0]['semanticTag'],
                    'label': item['derivedFrom']['annotatedProperty']['propertyValue'],
                    'omim_id': self.map_omim[property_value]
                }

            else:
                self.other_zooma_mappings[property_value] = {
                    'uri': item['_links']['olslinks'][0]['semanticTag'],
                    'label': item['derivedFrom']['annotatedProperty']['propertyValue'],
                    'omim_id': self.map_omim[property_value]
                }

        return self.high_confidence_mappings

    def use_zooma(self):
        '''
        Call request to Zooma function
        :return: None.
        '''


        logger.info("use Zooma")
        for phenotype in self.phenotype_set:
            if phenotype:
                self._logger.info("Mapping '%s' with zooma..."%(phenotype))
                self.request_to_zooma(phenotype)

        with open(Config.GE_ZOOMA_DISEASE_MAPPING, 'w') as outfile:
            tsv_writer = csv.writer(outfile, delimiter='\t')
            for phenotype, value in self.high_confidence_mappings.items():
                tsv_writer.writerow([phenotype, value['uri'], value['label'], value['omim_id']])

        with open(Config.GE_ZOOMA_DISEASE_MAPPING_NOT_HIGH_CONFIDENT, 'w') as map_file:
            csv_writer = csv.writer(map_file, delimiter='\t')
            for phenotype, value in self.other_zooma_mappings.items():
                csv_writer.writerow([phenotype, value['uri'], value['label'], value['omim_id']])

    def process_panel_app_file(self, sample=TEST_SAMPLE):
        '''
        Core method to create Evidence strings
        :return: None
        '''
        logger.info("Process panel app file")
        now = datetime.datetime.now()


        with open(Config.GE_ZOOMA_DISEASE_MAPPING, 'w') as outfile:
            tsv_writer = csv.writer(outfile, delimiter='\t')
            for phenotype, value in self.high_confidence_mappings.items():
                tsv_writer.writerow([phenotype, value['uri'], value['label'], value['omim_id']])

        mapping_fh = open("/tmp/genomics_england_mapped.txt", 'w')
        mapping_tsv_writer = csv.writer(mapping_fh, delimiter='\t')
        unmapped_fh = open("/tmp/genomics_england_unmapped.txt", 'w')
        unmapped_tsv_writer = csv.writer(unmapped_fh, delimiter='\t')

        for row in self.panel_app_info:
            # encode to UTF8
            panel_name, panel_id, gene_symbol, ensemble_gene_ids, ensembl_iri, level_of_confidence, original_label, phenotype, publications, evidences, omim_ids, disease_uri, disease_label, mapping_method = row
            #panel_name = panel_name.encode("UTF-8")
            #panel_id = panel_id.encode("UTF-8")
            #original_label = original_label.encode("UTF-8")

            if len(omim_ids) > 0 and disease_uri:
                self.generate_single_evidence(panel_name, panel_id, gene_symbol, ensemble_gene_ids, ensembl_iri, level_of_confidence,
                                              original_label, phenotype, publications, evidences, omim_ids, disease_uri, disease_label, now)
                mapping_tsv_writer.writerow([panel_name, panel_id, gene_symbol, phenotype, disease_uri, disease_label])
            elif phenotype in self.high_confidence_mappings:
                disease_label = self.high_confidence_mappings[phenotype]['label']
                disease_uri = self.high_confidence_mappings[phenotype]['uri']
                self.generate_single_evidence(panel_name, panel_id, gene_symbol, ensemble_gene_ids, ensembl_iri, level_of_confidence, original_label, phenotype, publications, evidences, omim_ids, disease_uri , disease_label, now)
                mapping_tsv_writer.writerow([panel_name, panel_id, gene_symbol, phenotype, disease_uri, disease_label])
            elif phenotype.lower() in self.zooma_to_efo_map:
                for item in self.zooma_to_efo_map[phenotype.lower()]:
                    disease_uri = item['efo_uri']
                    disease_label = "N/A"
                    self.generate_single_evidence(panel_name, panel_id, gene_symbol, ensemble_gene_ids, ensembl_iri, level_of_confidence, original_label, phenotype, publications, evidences, omim_ids, disease_uri , disease_label, now)
                    mapping_tsv_writer.writerow(
                        [panel_name, panel_id, gene_symbol, phenotype, disease_uri, disease_label])
            elif disease_uri:
                self.generate_single_evidence(panel_name, panel_id, gene_symbol, ensemble_gene_ids, ensembl_iri,
                                              level_of_confidence, original_label, phenotype, publications, evidences, omim_ids,
                                              disease_uri, disease_label, now)
                mapping_tsv_writer.writerow(
                    [panel_name, panel_id, gene_symbol, original_label, disease_uri, disease_label])
            else:
                unmapped_tsv_writer.writerow(
                        [panel_name, panel_id, gene_symbol, original_label])

            if sample == True:
                break

        mapping_fh.close()
        unmapped_fh.close()

    def generate_single_evidence(self, panel_name, panel_id, gene_symbol, ensembl_gene_id, ensembl_iri, level_of_confidence, original_label, phenotype, publications, evidences, omim_ids, disease_uri, disease_label, now):

        '''
        Generate evidence string using the Literature Curated schema
        :param panel_name:
        :param panel_id:
        :param gene_symbol:
        :param ensembl_gene_id:
        :param ensembl_iri:
        :param level_of_confidence:
        :param original_label:
        :param phenotype:
        :param publications:
        :param evidences:
        :param omim_ids:
        :param disease_uri:
        :param disease_label:
        :param now:
        :return:
        '''

        '''
        Create a list of publications
        '''
        single_lit_ref_list = list()
        if publications is not None and len(str(publications).strip()) > 0:
            publications = re.findall(r"\'(.+?)\'", str(publications))
            if len(publications) > 0:
                pset = set()
                for paper_set in publications:
                    print(paper_set)
                    paper_set = re.findall(r"\d{7,12}", paper_set)
                    for paper in paper_set:
                        print(paper)
                        if paper not in pset:
                            pset.add(paper)
                            lit_url = "http://europepmc.org/abstract/MED/" + paper
                            single_lit_ref_list.append(evidence_core.Single_Lit_Reference(lit_id=lit_url))

        if len(single_lit_ref_list) == 0:
            print("no publication found for %s %s %s" % (panel_name, panel_id, publications))
            lit_url = "NA"
            single_lit_ref_list.append(evidence_core.Single_Lit_Reference(lit_id=lit_url))

        obj = opentargets.Literature_Curated(type='genetic_literature')
        target = "http://identifiers.org/ensembl/" + ensembl_gene_id
        provenance_type = evidence_core.BaseProvenance_Type(
            database=evidence_core.BaseDatabase(
                id="Genomics England PanelApp",
                version='v4.1',
                dbxref=evidence_core.BaseDbxref(
                    url="https://panelapp.genomicsengland.co.uk/",
                    id="Genomics England PanelApp", version="v5.7")),
            literature=evidence_core.BaseLiterature(
                references=single_lit_ref_list
            )
        )
        obj.access_level = "public"
        obj.sourceID = "genomics_england"
        obj.validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION

        obj.unique_association_fields = dict(
            panel_name=panel_name,
            original_disease_name=disease_label,
            panel_id=panel_id,
            target_id=ensembl_iri,
            disease_iri=disease_uri)

        json.dumps(obj.unique_association_fields)
        provenance_type.to_JSON()

        # it's unicode
        as_str = json.dumps(obj.unique_association_fields)
        md5 = hashlib.md5(as_str.encode("UTF-8"))
        hashkey = md5.hexdigest()
        if hashkey in self.hashkeys:

            self._logger.warning(
                "Doc {0} - Duplicated evidence string for {1} to disease {2} URI: {3}".format(panel_name,
                                                                                              panel_id, disease_label,
                                                                                              disease_uri))
        else:

            self.hashkeys[hashkey] = obj

            obj.target = bioentity.Target(id=ensembl_iri,
                                          activity="http://identifiers.org/cttv.activity/predicted_damaging",
                                          target_type='http://identifiers.org/cttv.target/gene_evidence',
                                          target_name=gene_symbol)
            if level_of_confidence == 'LowEvidence':
                resource_score = association_score.Probability(
                    type="probability",
                    method=association_score.Method(
                        description="Further details in the Genomics England PanelApp.",
                        #reference="NA",
                        url="https://panelapp.genomicsengland.co.uk"),
                value = 0.25)
            elif level_of_confidence == 'HighEvidence':
                resource_score = association_score.Probability(
                    type="probability",
                    method=association_score.Method(
                        description="Further details in the Genomics England PanelApp.",
                        #reference="NA",
                        url="https://panelapp.genomicsengland.co.uk"),
                value = 1)
            else :
                resource_score = association_score.Probability(
                    type="probability",
                    method=association_score.Method(
                        description="Further details in the Genomics England PanelApp.",
                        #reference="NA",
                        url="https://panelapp.genomicsengland.co.uk"),
                value = 0.5)

            obj.disease = bioentity.Disease(id=disease_uri, source_name=phenotype, name=disease_label)
            obj.evidence = evidence_core.Literature_Curated()
            obj.evidence.is_associated = True
            obj.evidence.evidence_codes = ["http://purl.obolibrary.org/obo/ECO_0000205"]
            obj.evidence.provenance_type = provenance_type
            obj.evidence.date_asserted = now.isoformat()
            obj.evidence.provenance_type = provenance_type
            obj.evidence.resource_score = resource_score
            #specific
            linkout = evidence_linkout.Linkout(
                url=urllib.parse.urljoin(Config.GE_LINKOUT_URL, panel_id, gene_symbol),
                nice_name='Further details in the Genomics England PanelApp')
            obj.evidence.urls = [linkout]
            linkout.to_JSON()

    def write_evidence_strings(self, filename):
        '''
        Validate the evidence string file internally and write it to a file
        :param filename: name of the empty input file
        :return: None .  Writes the output in the input file
        '''
        logger.info("Writing Genomics England evidence strings")
        with open(filename, 'w') as tp_file:
            for hashkey, evidence_string in self.hashkeys.items():
                error = evidence_string.validate(logger)
                if error > 1:
                    logger.error(evidence_string.to_JSON(indentation=4))
                tp_file.write(evidence_string.to_JSON(indentation=None) + "\n")

        blob = self.bucket.blob(filename)
        with open(filename, 'rb') as my_file:
            blob.upload_from_file(my_file)

def main():
    ge_object = GE()
    ge_object.execute_ge_request()
    ge_object.use_zooma()
    ge_object.process_panel_app_file()
    ge_object.write_evidence_strings(Config.GE_EVIDENCE_FILENAME)


if __name__ == "__main__":
    main()
