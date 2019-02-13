from collections import OrderedDict
import logging
import urllib
from settings import Config

class RareDiseaseMapper(object):

    def __init__(self):
        print("init RareDiseaseMapper")
        super(RareDiseaseMapper, self).__init__()
        self._logger = logging.getLogger(__name__+".RareDiseaseMapper")
        self.omim_to_efo_map = OrderedDict()
        self.zooma_to_efo_map = OrderedDict()

    def get_omim_to_efo_mappings(self):
        self._logger.info("OMIM to EFO parsing - requesting from URL %s" % Config.OMIM_TO_EFO_MAP_URL)
        response = urllib.request.urlopen(Config.OMIM_TO_EFO_MAP_URL)
        self._logger.info("OMIM to EFO parsing - response code %s" % response.status)
        line_count = 0
        for line in response.readlines():
            line = line.decode('utf8').strip()
            #if its an empty line after stripping, skip it
            if len(line) == 0:
                continue
            self._logger.debug("Parsing line: %s",line)
            '''
            omim	efo_uri	efo_label	source	status
            '''
            line_count += 1
<<<<<<< HEAD
            if len(line.decode('utf8').strip()) > 0:
                (omim, efo_uri, efo_label, source, status) = line.decode('utf8').strip().split("\t")
                if omim not in self.omim_to_efo_map:
                    self.omim_to_efo_map[omim] = []
                self.omim_to_efo_map[omim].append({'efo_uri': efo_uri, 'efo_label': efo_label})
=======
            (omim, efo_uri, efo_label, source, status) = line.split("\t")
            if omim not in self.omim_to_efo_map:
                self.omim_to_efo_map[omim] = []
            self.omim_to_efo_map[omim].append({'efo_uri': efo_uri, 'efo_label': efo_label})
>>>>>>> c24d58009581c131d6b58d29642f3d8f78c2e620
        return line_count

    def get_opentargets_zooma_to_efo_mappings(self):
        self._logger.info("ZOOMA to EFO parsing - requesting from URL %s" % Config.ZOOMA_TO_EFO_MAP_URL)
        response = urllib.request.urlopen(Config.ZOOMA_TO_EFO_MAP_URL)
        self._logger.info("ZOOMA to EFO parsing - response code %s" % response.status)
        n = 0
        for line in response.readlines():
            line = line.decode('utf8').strip()
            #if its an empty line after stripping, skip it
            if len(line) == 0:
                continue
            self._logger.debug("Parsing line: %s",line)
            '''
            STUDY	BIOENTITY	PROPERTY_TYPE	PROPERTY_VALUE	SEMANTIC_TAG	ANNOTATOR	ANNOTATION_DATE
            disease	Amyotrophic lateral sclerosis 1	http://www.ebi.ac.uk/efo/EFO_0000253
            '''
            n +=1
            if n > 1:
                #self._logger.info("[%s]"%line)
                (study, bioentity, property_type, property_value, semantic_tag, annotator, annotation_date) = line.split("\t")
                if property_value.lower() not in self.omim_to_efo_map:
                    self.zooma_to_efo_map[property_value.lower()] = []
                self.zooma_to_efo_map[property_value.lower()].append({'efo_uri': semantic_tag, 'efo_label': semantic_tag})
        return n
