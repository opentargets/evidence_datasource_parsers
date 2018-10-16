import re
import csv
import json
import requests
import datetime
import hashlib
import logging
import urllib.request
import opentargets.model.core as opentargets
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.linkout as evidence_linkout
import opentargets.model.evidence.association_score as association_score

from ontoma import OnToma
from settings import Config
from collections import OrderedDict
from common.HGNCParser import GeneParser

__copyright__ = "Copyright 2014-2018, Open Targets"
__credits__   = ["Gautier Koscielny", "ChuangKee Ong"]
__license__   = "Apache 2.0"
__version__   = "1.2.8"
__maintainer__= "ChuangKee Ong"
__email__     = ["data@opentargets.org"]
__status__    = "Production"

'''
* Get 'Phenotypes' list for each 'HighEvidence' gene for all GE Panel.
* For each phenotype check that if it has an OMIM id.

* If omim_id - Y , try map to EFO using mapping file (https://github.com/opentargets/OnToma/blob/master/ontoma/constants.py)
  ** omim_id - Y, curated EFO - Y
  ** omim_id - Y, curated EFO - N, ontoma - Y (using omim_id)
  ** omim_id - Y, curated EFO - N, ontoma - N
* If omim_id - N , try OnToma (using phenotype string)
  ** omim_id - N , ontoma - Y
  ** omim_id - N , ontoma - N

notes:
* Configure such that ZOOMA API is not use in OntoMA
* ZOOMA API is then used for both category below to get high and low confidence mapping
    ** omim_id - Y, curated EFO - N, ontoma - N
    ** omim_id - N , ontoma - N
'''
class GEPanelApp():
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.hashkeys = OrderedDict()
        self.panel_app_info = list()
        self.panel_app_id_map = OrderedDict()
        self.confidence_zooma_mappings = OrderedDict()
        self.other_zooma_mappings = OrderedDict()

    def process_ge(self):
        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()

        self.symbol_to_gene = gene_parser.genes
        self.ontoma = OnToma(exclude=['zooma'])

        self.logger.info("Parsed %s OMIM to EFO mapping " % len(self.ontoma._omim_to_efo))
        self.logger.info("Parsed %s ZOOMA to EFO mapping " % len(self.ontoma._zooma_to_efo_map))
        self.logger.info("Parsed %s Name to EFO mapping " % len(self.ontoma.name_to_efo))
        self.logger.info("Parsed %s Name to HP mapping " % len(self.ontoma.name_to_hp))

        self.process_ge_panel()
        self.map_with_zooma()
        self.get_panel_id_mapping()
        self.build_evidence()
        self.write_evidence(Config.GE_EVIDENCE_FILENAME)

    def process_ge_panel(self):
        '''
         Create panel app info list and phenotype set
         :return: Unique phenotype list with no mappings
        '''
        self.logger.info("Retrieving GE PanelApp data")
        url = 'https://panelapp.genomicsengland.co.uk/WebServices/search_genes/all/'

        phenotype_list = []

        c = 0
        d = 0
        e = 0
        f = 0
        g = 0

        for panel_name, panel_id, panel_version, panel_diseasegroup, panel_diseasesubgroup in self.get_panel_list():
            self.logger.info("Reading panel : %s %s %s %s %s"
                             % (panel_name, panel_id, panel_version, panel_diseasegroup, panel_diseasesubgroup))
            r = requests.get(url, params={"panel_name": panel_name}, timeout=30)
            results = r.json()

            for row in results['results']:
                ensembl_iri = None
                ensembl_gene_id = self.symbol_to_gene.get(row['GeneSymbol'])

                if not ensembl_gene_id:
                    self.logger.warning("%s is not found in Ensembl" % row['GeneSymbol'])
                    continue

                ensembl_iri = "http://identifiers.org/ensembl/" + ensembl_gene_id

                if ensembl_iri and row['EnsembleGeneIds'] and row['Phenotypes'] \
                        and row['LevelOfConfidence'] == 'HighEvidence':
                    for item in row['Phenotypes']:
                        item = item.rstrip().lstrip().rstrip("?")
                        if len(item) > 0:
                            self.logger.info("HighEvidence Gene %s : %s has phenotypes: %s"
                                             % (row['GeneSymbol'], ensembl_iri, item))
                            '''
                              These phenotypes have OMIM Id
                                610219
                                Orofacial cleft 6, 608864
                                Nephrotic syndrome 14	617575
                                [Hair morphology 1, hair thickness], 612630 -3

                                note: need to be strictly 6 digits & without "Orpha|ORPHA|HP|PMID" within
                            '''
                            # MAY 20
                            # omim - Y curated EFO - Y => 8035
                            # omim - Y curated EFO - N ontoma - Y => 87
                            # omim - Y curated EFO - N ontoma - N => 132 --> ZOOMA
                            # omim - N ontoma - Y => 9415
                            # omim - N ontoma - N => 369 ---> ZOOMA
                            search_omim  = re.search('[^0-9]*(\d{6})[^0-9]*', item)
                            search_omim2 = re.search('(Orpha|ORPHA|HP|PMID)', item)

                            if search_omim and search_omim2 is None:
                                omim_id = search_omim.group(1)
                                if omim_id in self.ontoma._omim_to_efo:
                                    for efo_id in self.ontoma._omim_to_efo[omim_id]:
                                        c += 1
                                        mapping_type = 'OMIM to EFO curated mapping'
                                        disease_uri = efo_id['iri']
                                        disease_label = efo_id['label']

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
                                                                    disease_uri,
                                                                    disease_label,
                                                                    mapping_type])

                                        self.logger.info("[%s] => [%s], omim_id - Y, curated EFO - Y" % (item, efo_id))
                                else:
                                    if self.ontoma.find_term(omim_id, suggest=True, verbose=True):
                                        d += 1
                                        '''
                                          omim_id - Y, curated EFO - N, ontoma - Y

                                          Camptodactyly-arthropathy-coxa vara-pericarditis syndrome   208250
                                          Chondrocalcinosis 2 118600
                                          Neutropenia, nonimmune chronic idiopathic, of adults, 607847
                                        '''
                                        phenotype_label = self.ontoma.find_term(omim_id, suggest=True, verbose=True)
                                        mapping_type = 'OMIM to EFO mapping via OnToma'
                                        disease_uri = phenotype_label['term']
                                        disease_label = phenotype_label['label']

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
                                                                    disease_uri,
                                                                    disease_label,
                                                                    mapping_type])

                                        self.logger.info("[%s] => [%s], omim_id - Y, curated EFO - N, ontoma - Y action [%s]"
                                            % (item, phenotype_label['term'], phenotype_label['action']))
                                    else:
                                        '''
                                         omim_id - Y, curated EFO - N, ontoma - N

                                         Alkuraya-Kucinskas syndrome	617822
                                         Glycosylphosphatidylinositol biosynthesis defect 15, 617810
                                         Polycystic kidney disease, adult type I, 173900
                                        '''
                                        e += 1
                                        phenotype_list.append(item)
                                        self.logger.info("[%s] => [None], omim_id - Y, curated EFO - N, ontoma - N " % (item))
                            else:
                                if self.ontoma.find_term(item, suggest=True, verbose=True):
                                    f += 1
                                    '''
                                     omim_id - N, curated EFO - N, ontoma - Y

                                     Peeling skin HP:0040189
                                     ORPHA:220493  Joubert syndrome with ocular defect
                                     Brain Dopamine Serotonin Vesicular Transport Disease (Other disorders of neurotransmitter metabolism)
                                    '''
                                    phenotype_label = self.ontoma.find_term(item, suggest=True, verbose=True)
                                    mapping_type = 'PHENOTYPE string to EFO mapping via OnToma'
                                    disease_uri = phenotype_label['term']
                                    disease_label = phenotype_label['label']

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
                                                                disease_uri,
                                                                disease_label,
                                                                mapping_type])

                                    self.logger.info("[%s] => [%s], omim_id - N, curated EFO - N, ontoma - Y action [%s] " % (item, phenotype_label['term'], phenotype_label['action']))
                                else:
                                    '''
                                     omim_id - N, curated EFO - N, ontoma - N

                                     hypogammaglobulinaemia
                                     Oral and GI squamous cell carcinoma
                                     CEREBELLAR ATAXIA MENTAL RETARDATION AND DYSEQUILIBRIUM SYNDROME TYPE 1 (CMARQ1)
                                    '''
                                    g += 1
                                    phenotype_list.append(item)
                                    self.logger.info("[%s] => [None], omim_id - N, curated EFO - N, ontoma - N " % (item))

            self.phenotype_set = set(phenotype_list)

        self.logger.info("total OMIM_id found with EFO mapping " + str(c))
        self.logger.info("total OMIM_id found with NO EFO mapping but map via Ontoma " + str(d))
        self.logger.info("total OMIM_id found with NO EFO mapping & Ontoma " + str(e))
        self.logger.info("NO OMIM_id found, with NO EFO mapping but map via Ontoma " + str(f))
        self.logger.info("NO OMIM_id found, with NO EFO mapping & Ontoma " + str(g))

        return self.phenotype_set

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
                yield (row['Name'], row['Panel_Id'], row["CurrentVersion"], row["DiseaseGroup"], row["DiseaseSubGroup"])

        except requests.exceptions.HTTPError as e:
            return "Error: " + str(e)

    def map_with_zooma(self):
        '''
         Call ZOOMA API to map terms NOT map via omim_id & OnToma
        :return: None.
        '''
        self.logger.info("Mapping using ZOOMA API for terms NOT map through omim_id and OnToma")

        for phenotype in self.phenotype_set:
            if phenotype:
                self.logger.info("Mapping term: '%s' with ZOOMA API" % (phenotype))
                self.request_to_zooma(phenotype)

        with open(Config.GE_ZOOMA_DISEASE_MAPPING, 'w') as outfile:
            tsv_writer = csv.writer(outfile, delimiter='\t')
            for phenotype, value in self.confidence_zooma_mappings.items():
                tsv_writer.writerow([phenotype, value['uri'], value['label']])

        with open(Config.GE_ZOOMA_DISEASE_MAPPING_NOT_HIGH_CONFIDENT, 'w') as map_file:
            csv_writer = csv.writer(map_file, delimiter='\t')
            for phenotype, value in self.other_zooma_mappings.items():
                csv_writer.writerow([phenotype, value['uri'], value['label']])

    def request_to_zooma(self, property_value=None):
        '''
         Request ZOOMA API to get correct phenotype mapping and disease label

         :param property_value: Phenotype name from Genomics England
         :return: High confidence mappings .
         see docs: http://www.ebi.ac.uk/spot/zooma/docs/api.html
        '''
        r = requests.get('http://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate',
                         params={'propertyValue': property_value, 'propertyType': 'phenotype'})
        results = r.json()

        for item in results:
            if item['confidence'] == "HIGH":
                self.confidence_zooma_mappings[property_value] = {
                  'uri': item['_links']['olslinks'][0]['semanticTag'],
                  'label': item['derivedFrom']['annotatedProperty']['propertyValue'],
                }
            else:
                self.other_zooma_mappings[property_value] = {
                    'uri': item['_links']['olslinks'][0]['semanticTag'],
                    'label': item['derivedFrom']['annotatedProperty']['propertyValue'],
                }

        return self.confidence_zooma_mappings

    def get_panel_id_mapping(self, filename=Config.GE_PANEL_MAPPING_FILENAME):

        self.logger.info('Mapping old GE panel id to new panel id, using : ' + filename)

        with open(filename, mode='rt') as fh:
            reader = csv.reader(fh, delimiter=',', quotechar='"')
            next(reader, None) # skip header

            for row in reader:
                (old_panel_id, new_panel_id) = row
                self.logger.info("panel id :%s => is now %s" % (old_panel_id, new_panel_id))
                self.panel_app_id_map[old_panel_id] = new_panel_id

    def build_evidence(self):
        '''
         Core method to create Evidence strings
         :return: None
        '''
        now = datetime.datetime.now()

        for row in self.panel_app_info:
            panel_name, panel_id, panel_version, panel_diseasegroup, panel_diseasesubgroup, gene_symbol, ensembl_gene_ids, ensembl_iri, level_of_confidence, original_label, phenotype, publications, evidences, omim_ids, disease_uri, disease_label, mapping_method = row
            # encode to UTF8 panel_name = panel_name.encode("UTF-8")

            self.generate_single_evidence(panel_name, panel_id, panel_version, panel_diseasesubgroup, panel_diseasesubgroup,
                                          gene_symbol, ensembl_gene_ids, ensembl_iri, level_of_confidence,
                                          original_label, phenotype, publications, evidences, omim_ids,
                                          disease_uri, disease_label, now)

    def generate_single_evidence(self, panel_name, panel_id, panel_version, panel_diseasegroup, panel_diseasesubgroup,
                                 gene_symbol, ensembl_gene_id, ensembl_iri, level_of_confidence,
                                 original_label, phenotype, publications, evidences, omim_ids,
                                 disease_uri, disease_label, now):
        '''
        Generate evidence string using the Literature Curated schema
            :param panel_name:
            :param panel_id:
            :param panel_version:
            :param panel_diseasegroup:
            :param panel_diseasesubgroup:
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
        single_lit_ref_list = list()

        if publications is not None and len(str(publications).strip()) > 0:
            publications = re.findall(r"\'(.+?)\'", str(publications))
            if len(publications) > 0:
                pset = set()
                for paper_set in publications:
                    paper_set = re.findall(r"\d{7,12}", paper_set)
                    for paper in paper_set:
                        if paper not in pset:
                            pset.add(paper)
                            lit_url = "http://europepmc.org/abstract/MED/" + paper
                            single_lit_ref_list.append(evidence_core.Single_Lit_Reference(lit_id=lit_url))

        if len(single_lit_ref_list) == 0:
            self.logger.warning("no publication found for %s %s %s" % (panel_name, panel_id, publications))
            lit_url = "NA"
            single_lit_ref_list.append(evidence_core.Single_Lit_Reference(lit_id=lit_url))

        obj = opentargets.Literature_Curated(type='genetic_literature')
        target = "http://identifiers.org/ensembl/" + ensembl_gene_id

        provenance_type = evidence_core.BaseProvenance_Type(
            database=evidence_core.BaseDatabase(
                id="Genomics England PanelApp",
                version=Config.GE_PANEL_VERSION,
                dbxref=evidence_core.BaseDbxref(
                    url="https://panelapp.genomicsengland.co.uk/",
                    id="Genomics England PanelApp", version=Config.GE_PANEL_VERSION)),
            literature=evidence_core.BaseLiterature(
                references=single_lit_ref_list
            )
        )
        obj.access_level = "public"
        obj.sourceID = "genomics_england"
        obj.validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION

        if panel_id in self.panel_app_id_map:
            new_panel_id = self.panel_app_id_map[panel_id]
        else:
            new_panel_id = panel_id

        obj.unique_association_fields = dict(
            panel_name=panel_name,
            panel_id=new_panel_id,
            panel_version=panel_version,
            panel_diseasegroup=panel_diseasegroup,
            panel_diseasesubgroup=panel_diseasesubgroup,
            previous_code=panel_id,
            original_disease_name=disease_label,
            disease_iri=disease_uri,
            target_id=ensembl_iri)

        json.dumps(obj.unique_association_fields)
        provenance_type.to_JSON()

        # it's unicode
        as_str = json.dumps(obj.unique_association_fields)
        md5 = hashlib.md5(as_str.encode("UTF-8"))
        hashkey = md5.hexdigest()

        if hashkey in self.hashkeys:
            self.logger.warning("Panel Name: {0} - Duplicated evidence string for {1} to disease {2} URI: {3}"
                                .format(panel_name, panel_id, disease_label, disease_uri))
        else:
            self.hashkeys[hashkey] = obj

            obj.target = bioentity.Target(id=ensembl_iri,
                                          activity="http://identifiers.org/cttv.activity/predicted_damaging",
                                          target_type='http://identifiers.org/cttv.target/gene_evidence',
                                          target_name=gene_symbol)

            #TODO : should we have only High confidence genes from #L105
            if level_of_confidence == 'LowEvidence':
                resource_score = association_score.Probability(
                    type="probability",
                    method=association_score.Method(
                        description="Further details in the Genomics England PanelApp.",
                        url="https://panelapp.genomicsengland.co.uk"),
                    value=0.25)
            elif level_of_confidence == 'HighEvidence':
                resource_score = association_score.Probability(
                    type="probability",
                    method=association_score.Method(
                        description="Further details in the Genomics England PanelApp.",
                        url="https://panelapp.genomicsengland.co.uk"),
                    value=1)
            else:
                resource_score = association_score.Probability(
                    type="probability",
                    method=association_score.Method(
                        description="Further details in the Genomics England PanelApp.",
                        url="https://panelapp.genomicsengland.co.uk"),
                    value=0.5)

            obj.disease = bioentity.Disease(id=disease_uri, source_name=phenotype, name=disease_label)
            obj.evidence = evidence_core.Literature_Curated()
            obj.evidence.is_associated = True
            obj.evidence.date_asserted = now.isoformat()
            obj.evidence.evidence_codes = ["http://purl.obolibrary.org/obo/ECO_0000205"]
            obj.evidence.provenance_type = provenance_type
            obj.evidence.resource_score = resource_score

            linkout = evidence_linkout.Linkout(
                url=urllib.parse.urljoin(Config.GE_LINKOUT_URL, "%s/%s" % (new_panel_id, gene_symbol)),
                nice_name='Further details in the Genomics England PanelApp for panel %s and gene %s ' % (new_panel_id, gene_symbol))
            obj.evidence.urls = [linkout]
            linkout.to_JSON()


    def write_evidence(self, filename):
        '''
         Validate evidence string file and write it to a file
         :param filename: name of the empty input file
         :return: None . Writes the output in the input file
        '''
        self.logger.info("Writing Genomics England evidence strings")

        with open(filename, 'w') as ge_output:
            c = 0
            for hashkey, evidence_string in self.hashkeys.items():
                c +=1

                ge_output.write(evidence_string.to_JSON(indentation=None) + "\n")

                error = evidence_string.validate(self.logger)

                if error== 0:
                    ge_output.write(evidence_string.to_JSON(indentation=None) + "\n")
                else:
                    self.logger.error("Reporting validation error in line %s" % c)
                    self.logger.error(evidence_string.to_JSON(indentation=4))
        ge_output.close()