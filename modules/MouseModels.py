#!/usr/bin/env python3

import argparse
import collections
import datetime
import gzip
import hashlib
import json
import logging
import os
import re
import tempfile
import urllib.request

import pyspark
import requests
from retry import retry

from common.RareDiseasesUtils import RareDiseaseMapper
from ontologyutils.rdf_utils import OntologyClassReader
import opentargets.model.core as cttv
import opentargets.model.bioentity as bioentity
import opentargets.model.evidence.phenotype as evidence_phenotype
import opentargets.model.evidence.core as evidence_core
import opentargets.model.evidence.association_score as association_score
from settings import Config

HGNC_DATASET_URI = 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt'
MGI_DATASET_URI = 'http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt'


class Phenodigm(RareDiseaseMapper):

    def __init__(self, logging):
        super(Phenodigm, self).__init__()
        self.spark = pyspark.sql.SparkSession.builder.appName('phenodigm_parser').getOrCreate()
        self.efo = OntologyClassReader()
        self.hpo = OntologyClassReader()
        self.mp = OntologyClassReader()
        self.cache = collections.OrderedDict()
        self.counter = 0
        self.mmGenes = None
        self.hsGenes = None
        self.hgnc2mgis = collections.OrderedDict()
        self.mgi2symbols = collections.OrderedDict()
        self.symbol2hgncids = collections.OrderedDict()
        self.mgi2mouse_models = collections.OrderedDict()
        self.mouse_model2diseases = collections.OrderedDict()
        self.disease_gene_locus = collections.OrderedDict()
        self.ontology = collections.OrderedDict()
        self.ontology_ontology = collections.OrderedDict()
        self.mouse_models = collections.OrderedDict()
        self.diseases = collections.OrderedDict()
        self.hashkeys = collections.OrderedDict()
        self._logger = logging

    def fetch_human_ensembl_mappings(self):
        """Loads HGNC ID → Ensembl gene ID mappings from the HGNC database."""
        # TODO: possibly unify fetching via temp file
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        urllib.request.urlretrieve(HGNC_DATASET_URI, tmp_file.name)
        human_genes = (
            self.spark.read.csv(tmp_file.name, sep='\t', header=True)
                .select('hgnc_id', 'ensembl_gene_id')
        )
        # TODO: use the Spark dataframe directly for joins
        self.hsGenes = {r.hgnc_id: r.ensembl_gene_id for r in human_genes.collect()}
        os.remove(tmp_file.name)

    def fetch_mouse_ensembl_mappings(self):
        """Loads MGI ID → Ensembl gene ID mappings from the MGI database."""
        # TODO: possibly unify fetching via temp file
        tmp_file = tempfile.NamedTemporaryFile(delete=False)
        urllib.request.urlretrieve(MGI_DATASET_URI, tmp_file.name)
        mouse_genes = (
            self.spark.read.csv(tmp_file.name, sep='\t', header=True)
                .withColumnRenamed('1. MGI accession id', 'mgi_id')
                .withColumnRenamed('11. Ensembl gene id', 'ensembl_gene_id')
                .select('mgi_id', 'ensembl_gene_id')
        )
        # TODO: use the Spark dataframe directly for joins
        self.mmGenes = {r.mgi_id: r.ensembl_gene_id for r in mouse_genes.collect()}
        os.remove(tmp_file.name)

    def access_solr(self):

        start = 0
        rows = 20000
        counter = 0
        total = 0
        numFound = 1

        while (total < numFound):
            counter += 1
            self._logger.info(start)

            rsp = self.query_solr(counter, start, rows)

            numFound = rsp['response']['numFound']
            self._logger.info("NumFound %i"%(numFound))
            nbItems = len(rsp['response']['docs'])
            total += nbItems
            self._logger.info("Number of items: %i" % (nbItems))

            self.parse_phenodigm(docs=rsp['response']['docs'])

            start += nbItems

    # Use @retry decorator to ensure that errors like the query failing because server was overloaded, are handled
    # correctly and the request is retried
    @retry(tries=3, delay=5, backoff=1.2, jitter=(1, 3))
    def query_solr(self, counter, start, rows):
        '''
        Queries IMPC solr API

        :param url: Query
        :return: Dictionary with response
        '''

        # SOLR base URL
        url = Config.MOUSEMODELS_PHENODIGM_SOLR + '/solr/phenodigm/select'

        # Full URI for printing
        uri = url + '?q=*:*&wt=json&indent=true&start=%i&rows=%i' % (start, rows)

        # Params for querying
        params = dict(q="*", wt="json", indent="true", start="%i" % start, rows="%i" % rows)
        self._logger.info("REQUEST {0}. {1}".format(counter, uri))

        # Query
        r = requests.get(url, params=params, timeout=30)
        self._logger.info("REQUEST %s" % (r.url))

        # Check for erroneous HTTP response statuses
        r.raise_for_status()
        rsp = r.json()
        return rsp

    def parse_phenodigm(self, docs=None):

        if docs:
            for doc in docs:
                '''
                1. Load all gene documents
                2. Load all mouse models of diseases.
                . Load all ontology mapping
                '''
                try:
                    if doc['type'] == 'gene':
                        # this is a human gene
                        if 'hgnc_gene_id' in doc:
                            hgnc_gene_id = doc['hgnc_gene_id']
                            hgnc_gene_symbol = doc['hgnc_gene_symbol']
                            self.symbol2hgncids[hgnc_gene_symbol] = hgnc_gene_id
                        elif 'gene_id' in doc:
                            gene_id = doc['gene_id']
                            gene_symbol = doc['gene_symbol']
                            self.mgi2symbols[gene_id] = gene_symbol
                            # this is a mouse gene
                    elif doc['type'] == 'gene_gene':
                        # gene_id # hgnc_gene_id
                        hgnc_gene_id = doc['hgnc_gene_id']
                        gene_id = doc['gene_id']
                        if hgnc_gene_id and not hgnc_gene_id in self.hgnc2mgis:
                            self.hgnc2mgis[hgnc_gene_id] = []
                        self.hgnc2mgis[hgnc_gene_id].append(gene_id)
                    elif doc['type'] == 'mouse_model':
                        marker_symbol = doc['marker_symbol']
                        marker_id = doc['marker_id']
                        model_id = doc['model_id']

                        # if there is a mouse model then add the marker id
                        if marker_id not in self.mgi2symbols:
                            self.mgi2symbols[marker_id] = marker_symbol

                        if not marker_symbol in self.mgi2mouse_models:
                            self.mgi2mouse_models[marker_symbol] = []
                        self.mgi2mouse_models[marker_symbol].append(model_id)

                        if not model_id in self.mouse_models:
                            self.mouse_models[model_id] = doc

                            model_phenotypes = []
                            for raw_mp in doc['model_phenotypes']:
                                mt = re.match("^(MP\:\d+)\s+", raw_mp)
                                if mt:
                                    mp_id = mt.groups()[0]
                                    model_phenotypes.append(mp_id)
                            self.mouse_models[model_id]['model_phenotypes'] = model_phenotypes
                    elif doc['type'] == 'disease_model_summary':
                        model_id = doc['model_id']
                        if not model_id in self.mouse_model2diseases:
                            self.mouse_model2diseases[model_id] = []
                        self.mouse_model2diseases[model_id].append(doc)
                    elif doc['type'] == 'disease_gene_summary' and 'marker_symbol' in doc and 'hgnc_gene_locus' in doc and 'disease_id' in doc:
                        hgnc_gene_id = doc['hgnc_gene_id']
                        hgnc_gene_symbol = doc['hgnc_gene_symbol']
                        marker_symbol = doc['marker_symbol']
                        try:
                            disease_id = doc['disease_id']
                            #if 'disease_id' in doc:
                            #    raise KeyError()
                        except Exception as error:

                            if isinstance(error, KeyError):
                                self._logger.error("Error checking disease in document: %s" % (str(error)))
                                self._logger.error(json.dumps(doc, indent=4))
                                raise Exception()
                        if not disease_id in self.disease_gene_locus:
                            self.disease_gene_locus[disease_id] = { hgnc_gene_id: [ marker_symbol ] }
                        elif not hgnc_gene_id in self.disease_gene_locus[disease_id]:
                            self.disease_gene_locus[disease_id][hgnc_gene_id] = [ marker_symbol ]
                        else:
                            self.disease_gene_locus[disease_id][hgnc_gene_id].append(marker_symbol)
                    elif doc['type'] == 'disease' or (doc['type'] == 'disease_gene_summary' and 'disease_id' in doc and not doc['disease_id'] in self.diseases):
                        '''and doc['disease_id'].startswith('ORPHANET')):'''
                        disease_phenotypes = []
                        if 'disease_phenotypes' in doc:
                            raw_disease_phenotypes = doc['disease_phenotypes']
                            for raw_phenotype in raw_disease_phenotypes:
                                mt = re.match("^(HP\:\d+)\s+", raw_phenotype)
                                hp_id = mt.groups()[0]
                                disease_phenotypes.append(hp_id)
                        disease_id = doc['disease_id']
                        self.diseases[disease_id] = doc
                        self.diseases[disease_id]['disease_phenotypes'] = disease_phenotypes

                        #matchOMIM = re.match("^OMIM:(.+)$", disease_id)
                        #if matchOMIM:
                        #    terms = efo.getTermsByDbXref(disease_id)
                        #    if terms == None:
                        #        terms = OMIMmap[disease_id]
                        #    if terms == None:
                        #        self._logger.error("{0} '{1}' not in EFO".format(disease_id, doc['disease_term']))
                    #else:
                    #    self._logger.error("Never heard of this '%s' type of documents in PhenoDigm. Exiting..."%(doc['type']))
                    #    sys.exit(1)
                    elif doc['type'] == 'ontology_ontology':
                        if doc['mp_id'] not in self.ontology_ontology:
                            self.ontology_ontology[doc['mp_id']] = []
                        self.ontology_ontology[doc['mp_id']].append(doc)
                    elif doc['type'] == 'ontology':
                        self.ontology[doc['phenotype_id']] = doc
                except KeyError as ke:
                    self._logger.error("KeyError \n", ke)
                    self._logger.error(json.dumps(doc))
                    break

    def generate_phenodigm_evidence_strings(self, upper_limit=0):
        '''
         Once you have retrieved all the genes,and mouse models
         Create an evidence string for every gene to disease relationship
        '''
        now = datetime.datetime.now()
        efoMapping = collections.OrderedDict()
        index_g = 0
        for hs_symbol in self.hsGenes:
            index_g += 1
            self._logger.debug('Now processing human gene: {}'.format(hs_symbol))
            hgnc_gene_id = self.symbol2hgncids[hs_symbol]
            hs_ensembl_gene_id = self.hsGenes[hs_symbol]

            '''
            Useful if you want to test on a subset of the data
            '''
            if upper_limit > 0 and len(self.hashkeys) > upper_limit:
                break

            if hgnc_gene_id and hs_ensembl_gene_id and re.match('^ENSG.*', hs_ensembl_gene_id) and hgnc_gene_id in self.hgnc2mgis:

                # self._logger.info("processing human gene {0} {1} {2} {3}".format(index_g, hs_symbol, hs_ensembl_gene_id, hgnc_gene_id))
                #print("processing human gene {0} {1} {2} {3}".format(index_g, hs_symbol, hs_ensembl_gene_id, hgnc_gene_id ))

                '''
                 Retrieve mouse models
                '''
                for mouse_gene_id in self.hgnc2mgis[hgnc_gene_id]:

                    #if not marker_symbol == "Il13":
                    #    continue;
                    if mouse_gene_id not in self.mgi2symbols:
                        continue

                    marker_symbol = self.mgi2symbols[mouse_gene_id]
                    #print(marker_symbol)
                    # self._logger.info("\tProcessing mouse gene symbol %s" % (marker_symbol))

                    '''
                    Some mouse symbol are not mapped in Ensembl
                    We will check that once we get the symbol
                    '''

                    if marker_symbol in self.mmGenes:
                        mmEnsemblId = self.mmGenes[marker_symbol]
                        '''
                         Some genes don't have any mouse models...
                        '''
                        if marker_symbol in self.mgi2mouse_models:

                            self._logger.debug("\tFound %i mouse models for this marker" % (len(self.mgi2mouse_models[marker_symbol])))

                            '''
                             if the marker is associated to mouse models
                            '''
                            for model_id in self.mgi2mouse_models[marker_symbol]:
                                '''
                                 loop through every mouse model
                                '''
                                mouse_model = self.mouse_models[model_id]

                                '''
                                List all mouse phenotypes
                                '''
                                model_mouse_phenotypes = mouse_model['model_phenotypes']
                                # get all potential human phenotypes
                                model_human_phenotypes = set()
                                for mp_id in model_mouse_phenotypes:
                                    if mp_id in self.ontology_ontology:
                                        for doc in self.ontology_ontology[mp_id]:
                                            model_human_phenotypes.add(doc['hp_id'])

                                marker_id = mouse_model['marker_id']
                                model_description = mouse_model['model_description']
                                self._logger.debug("\t\tMouse model {0} for {1}".format(model_id, marker_id))
                                '''
                                 Check the model_id is in the dictionary containing all the models
                                '''
                                if model_id in self.mouse_model2diseases:

                                    self._logger.debug("\t\tMouse model {0} is in the dictionary containing all the models".format(model_id))

                                    for mouse_model2disease in self.mouse_model2diseases[model_id]:

                                        '''
                                         get the disease identifier
                                        '''
                                        disease_id = mouse_model2disease['disease_id']

                                        '''
                                         Retrieve the disease document
                                        '''
                                        disease = None
                                        if disease_id in self.diseases:
                                            disease = self.diseases[disease_id]
                                            self._logger.debug("\t\t\tdisease: %s %s"%(disease_id, disease['disease_term']))
                                        else:
                                            self._logger.debug("\t\t\tdisease: %s"%(disease_id))

                                        '''
                                            Get all phenotypes and check the ones that match the mouse ones 
                                            we don't want to show all the phenotypes
                                        '''
                                        disease_phenotypes = self.diseases[disease_id]['disease_phenotypes']
                                        model_disease_human_phenotypes = list(filter(lambda x: x in model_human_phenotypes, disease_phenotypes))


                                        '''
                                        Map the disease ID to EFO
                                        Can be a one to many mapping
                                        '''
                                        disease_terms = list()

                                        diseaseName = None
                                        if not disease_id in efoMapping:
                                            # find corresponding EFO disease term
                                            matchOMIM = re.match("^OMIM:(.+)$", disease_id)
                                            matchORPHANET = re.match("^ORPHANET:(.+)$", disease_id)
                                            self._logger.debug('disease is: {}'.format(disease_id))
                                            if matchOMIM:
                                                if ';PMID' in disease_id:
                                                    (source, omim_id) = disease_id.split(';')[0].split(':')
                                                else:
                                                    (source, omim_id) = disease_id.split(':')
                                                self._logger.debug("\t\t\tCheck OMIM id = %s is in efo"%omim_id)
                                                if omim_id in self.omim_to_efo_map:
                                                    efoMapping[disease_id] = self.omim_to_efo_map[omim_id]
                                                if len(disease_terms) > 0:
                                                    #for disease_term in disease_terms:
                                                    #    if disease_term['efo_uri'] in self.efo.current_classes:
                                                    #        self._logger.info("{0} => {1} {2}".format(disease_id, disease_term['efo_uri'], self.efo.current_classes[disease_term['efo_uri']]))
                                                    #    else:
                                                    #        self._logger.info("{0} => {1} (no EFO mapping)".format(disease_id, disease_term['efo_uri']))
                                                    disease_terms = efoMapping[disease_id]

                                            elif matchORPHANET:
                                                    suffix = matchORPHANET.groups()[0]
                                                    orphanetId = "Orphanet:{0}".format(suffix)
                                                    orphanet_uri = "http://www.orpha.net/ORDO/Orphanet_{0}".format(suffix)
                                                    if orphanet_uri in self.efo.current_classes:
                                                        efoMapping[disease_id] = [ {'efo_uri': orphanet_uri, 'efo_label': self.diseases[disease_id]['disease_term'] } ]
                                                        disease_terms = efoMapping[disease_id]
                                        else:
                                            disease_terms = efoMapping[disease_id]

                                        '''
                                        OK, we have a disease mapped to EFO
                                        we can proceed to the next stage
                                        we don't filter on the score anymore.
                                        this will be adjusted in Open Targets
                                        If the score >= 0.5 or in_locus for the same disease
                                        and mouse_model2disease['model_to_disease_score'] >= 50
                                        '''
                                        if disease_terms is not None and (
                                                ('disease_model_max_norm' in mouse_model2disease and
                                                 mouse_model2disease['disease_model_max_norm']>=50) or
                                            (disease_id in self.disease_gene_locus and
                                             hgnc_gene_id in self.disease_gene_locus[disease_id] and
                                             marker_symbol in self.disease_gene_locus[disease_id][hgnc_gene_id])):

                                            for disease_term in disease_terms:

                                                '''
                                                Create a new evidence string
                                                '''



                                                # 1.2.6 create an Animal_Models class
                                                evidenceString = cttv.Animal_Models()
                                                evidenceString.validated_against_schema_version = Config.VALIDATED_AGAINST_SCHEMA_VERSION
                                                evidenceString.access_level = "public"
                                                evidenceString.type = "animal_model"
                                                evidenceString.sourceID = "phenodigm"
                                                evidenceString.unique_association_fields = collections.OrderedDict()
                                                #evidenceString.unique_association_fields['projectName'] = 'otar_external_mousemodels'
                                                evidenceString.evidence = cttv.Animal_ModelsEvidence()
                                                evidenceString.evidence.date_asserted = now.isoformat()
                                                evidenceString.evidence.is_associated = True

                                                '''
                                                Target
                                                '''
                                                evidenceString.target = bioentity.Target(
                                                    id="http://identifiers.org/ensembl/{0}".format(hs_ensembl_gene_id),
                                                    activity="http://identifiers.org/cttv.activity/predicted_damaging",
                                                    target_type="http://identifiers.org/cttv.target/gene_evidence"
                                                    )

                                                '''
                                                Disease
                                                '''
                                                disease_name = 'N/A'
                                                if disease_term['efo_uri'] in self.efo.current_classes:
                                                    disease_name = self.efo.current_classes[disease_term['efo_uri']]
                                                self._logger.debug("Disease id is %s"%(disease_term['efo_uri']))
                                                evidenceString.disease = bioentity.Disease(
                                                    id=disease_term['efo_uri'],
                                                    name=disease_name,
                                                    source_name=self.diseases[disease_id]['disease_term']
                                                    )

                                                human_gene_id = "http://identifiers.org/ensembl/%s"%(hs_ensembl_gene_id)
                                                model_gene_id = "http://identifiers.org/ensembl/%s"%(mmEnsemblId)

                                                '''
                                                Evidence Codes
                                                '''
                                                if mouse_model2disease['model_source'] == 'MGI':
                                                    #evidenceString.evidence.evidence_codes.append("http://identifiers.org/eco/ECO:0000057")
                                                    evidenceString.unique_association_fields['predictionModel'] = 'mgi_predicted'
                                                else:
                                                    #evidenceString.evidence.evidence_codes.append("http://identifiers.org/eco/ECO:0000057")
                                                    evidenceString.unique_association_fields['predictionModel'] = 'impc_predicted'

                                                #evidenceString.unique_association_fields['disease_phenodigm_name'] = self.diseases[disease_id]['disease_term']
                                                evidenceString.unique_association_fields['disease_phenodigm_id'] = disease_id
                                                evidenceString.unique_association_fields['disease_id'] = disease_term['efo_uri']
                                                evidenceString.unique_association_fields['target_id'] = human_gene_id
                                                evidenceString.unique_association_fields['model_gene_id'] = model_gene_id
                                                #evidenceString.unique_association_fields['species'] = "Mus musculus"
                                                evidenceString.unique_association_fields['model_genetic_background'] = mouse_model['model_genetic_background']
                                                evidenceString.unique_association_fields['model_description'] = mouse_model['model_description']
                                                '''
                                                Orthologs
                                                Human gene => Mouse marker
                                                '''
                                                evidenceString.evidence.orthologs = evidence_phenotype.Orthologs(
                                                    evidence_codes = ["http://identifiers.org/eco/ECO:0000265"],
                                                    provenance_type= evidence_core.BaseProvenance_Type(database=evidence_core.BaseDatabase(id="MGI", version="2016")),
                                                    resource_score= association_score.Probability(type="probability", method= association_score.Method(description ="orthology from MGI"), value=1.0),
                                                    date_asserted= now.isoformat(),
                                                    human_gene_id = "http://identifiers.org/ensembl/{0}".format(hs_ensembl_gene_id),
                                                    model_gene_id = "http://identifiers.org/ensembl/{0}".format(mmEnsemblId),
                                                    species = "mouse"
                                                    )

                                                '''
                                                Biological Model
                                                'allele_ids':'MGI:1862010|MGI:1862010',
                                                'allelic_composition':'Pmp22<Tr-1H>/Pmp22<+>',
                                                'genetic_background':'involves: BALB/cAnNCrl * C3H/HeN',
                                                '''
                                                zygosity = 'oth'
                                                allelic_composition = mouse_model['model_description']
                                                self._logger.debug(mouse_model['model_description'])
                                                if mouse_model['model_description'].endswith("hom") or mouse_model['model_description'].endswith("het") or mouse_model['model_description'].endswith("hem"):
                                                    zygosity = mouse_model['model_description'][-3:]
                                                    allelic_composition = mouse_model['model_description'][:-4]
                                                evidenceString.evidence.biological_model = evidence_phenotype.Biological_Model(
                                                    evidence_codes = ["http://identifiers.org/eco/ECO:0000179"],
                                                    resource_score= association_score.Probability(
                                                        type="probability",
                                                        method= association_score.Method(description =""),
                                                        value=1.0),
                                                    date_asserted= now.isoformat(),
                                                    model_id = mouse_model['model_id'],
                                                    zygosity = zygosity,
                                                    genetic_background = mouse_model['model_genetic_background'],
                                                    allelic_composition = allelic_composition,
                                                    model_gene_id = "http://identifiers.org/ensembl/{0}".format(mmEnsemblId),
                                                    species = "mouse"
                                                    )

                                                evidenceString.evidence.biological_model.provenance_type= evidence_core.BaseProvenance_Type(database=evidence_core.BaseDatabase(id=mouse_model2disease['model_source'], version="2017"))
                                                evidenceString.evidence.biological_model.allele_ids = ""
                                                ''' add all mouse phenotypes '''
                                                evidenceString.evidence.biological_model.phenotypes = []

                                                for mp_id in model_mouse_phenotypes:
                                                    evidenceString.evidence.biological_model.phenotypes.append(
                                                        bioentity.Phenotype(
                                                            id = mp_id,
                                                            term_id = "http://purl.obolibrary.org/obo/" + mp_id.replace(":", "_"),
                                                            label = self.ontology[mp_id]['phenotype_term'],
                                                            species="mouse"
                                                            )
                                                        )

                                                ''' get all human phenotypes '''
                                                human_phenotypes = []
                                                if len(model_disease_human_phenotypes) > 0:
                                                    for hp_id in model_disease_human_phenotypes:
                                                        term_id = "http://purl.obolibrary.org/obo/" + hp_id.replace(":","_")
                                                        if term_id in self.hpo.current_classes:
                                                            term_name = self.hpo.current_classes[term_id]
                                                        else:
                                                            term_name = self.ontology[hp_id]['phenotype_term']
                                                        # self._logger.debug("HPO term is {0} {1}".format(hp, hp in hpo.terms))
                                                        # self._logger.debug("HPO term retrieved is {0} {1}".format( termId, termName ))
                                                        human_phenotypes.append(
                                                            bioentity.Phenotype(
                                                                id = hp_id,
                                                                term_id = term_id,
                                                                label = term_name,
                                                                species="human"
                                                                )
                                                            )
                                                '''
                                                get all matched mouse phenotypes
                                                Format is:
                                                '''
                                                '''
                                                mouse_phenotypes = []
                                                if mp_matched_ids:
                                                    for mp in mp_matched_ids:
                                                        term_id = "http://purl.obolibrary.org/obo/" + mp.replace(":", "_")
                                                        term_name = None
                                                        if term_id in self.mp.current_classes:
                                                            term_name = self.mp.current_classes[term_id]
                                                        elif term_id in self.mp.obsolete_classes:
                                                            new_id = self.mp.get_new_from_obsolete_uri(term_id)
                                                            if new_id and new_id in self.mp.current_classes:
                                                                term_name = self.mp.current_classes[new_id]
                                                            else:
                                                                term_name = self.mp.obsoletes[term_id]['label']

                                                        #termId = mp
                                                        #termName = "TO BE DERTERMINED"
                                                        #self._logger.info("MP term is {0}".format(mp))
                                                        #term = mpo.getTermById(mp)
                                                        #termId = term['tags']['id'][0]
                                                        #termName = term['tags']['name'][0]
                                                        mouse_phenotypes.append(
                                                            bioentity.Phenotype(
                                                                id = mp,
                                                                term_id = term_id,
                                                                label = term_name,
                                                                species='mouse'
                                                                )
                                                            )
                                                '''
                                                mouse_phenotypes = evidenceString.evidence.biological_model.phenotypes

                                                '''
                                                Disease model association
                                                '''
                                                score = 1.0
                                                method = 'Unknown method'
                                                if 'disease_model_max_norm' in mouse_model2disease:
                                                    score = (mouse_model2disease['disease_model_max_norm'])/100
                                                    method = 'phenodigm_model_to_disease_score'
                                                self._logger.debug("score: {0}\n".format(score))
                                                evidenceString.unique_association_fields['score'] = "%.15f"%(score)

                                                evidenceString.evidence.disease_model_association = evidence_phenotype.Disease_Model_Association(
                                                    disease_id = disease_term['efo_uri'],
                                                    model_id = "{0}".format(mouse_model['model_id']),
                                                    provenance_type= evidence_core.BaseProvenance_Type(database=evidence_core.BaseDatabase(id="PhenoDigm", version="June 2016")),
                                                    evidence_codes = ["http://identifiers.org/eco/ECO:0000057"],
                                                    resource_score= association_score.Summed_Total(type = "summed_total", method = association_score.Method(description = method), value = score),
                                                    date_asserted= now.isoformat(),
                                                    model_phenotypes = mouse_phenotypes,
                                                    human_phenotypes = human_phenotypes
                                                )

                                                '''
                                                Make sure we don't create duplicates
                                                If the model is associated already to a disease
                                                take the best score
                                                '''
                                                # it's unicode
                                                as_str = json.dumps(evidenceString.unique_association_fields)
                                                md5 = hashlib.md5(as_str.encode("UTF-8"))
                                                hashkey = md5.hexdigest()
                                                if not hashkey in self.hashkeys:
                                                    self.hashkeys[hashkey] = evidenceString
                                                else:
                                                    self._logger.warn("Duplicated mouse model {0} to disease {1} URI: {2}".format(model_id, disease_id, disease_term))
                                                    if self.hashkeys[hashkey].unique_association_fields['score'] > evidenceString.unique_association_fields['score']:
                                                        self.hashkeys[hashkey] = evidenceString
                                        else:

                                            self._logger.error("Unable to incorporate this strain for this disease: {0}".format(disease_id))
                                            self._logger.error("No disease terms {0}".format(disease_terms == None))
                                            self._logger.error("disease_model_max_norm in mouse_model2disease: {0}".format('disease_model_max_norm' in mouse_model2disease))
                                            self._logger.error("disease_id in disease_gene_locus: {0}".format(disease_id in self.disease_gene_locus))
                                            #self._logger.error("hs_symbol in disease_gene_locus[disease_id]: {0}".format(not disease_term_uris == None and disease_id in self.disease_gene_locus and hgnc_gene_id in self.disease_gene_locus[disease_id]))
                                            #self._logger.error("marker_symbol in disease_gene_locus[disease_id][hgnc_gene_id]): {0}".format(disease_term_uris is not None and disease_id in self.disease_gene_locus and marker_symbol in self.disease_gene_locus[disease_id][hgnc_gene_id]))

    def convert_evidence(self, data):
        '''
        This function maps the old evidence structure into the new json schema.

        Args:
            evidenceString (dict): old PhenoDigm evidence string
        Returns:
            dict: dictionary for the new evidence if all fields could be successfully mapped.
        '''

        # Converting evidence string from this weird object to a more managable dictionary:
        data = json.loads(data.to_JSON())

        try:
            return collections.OrderedDict({
            'targetFromSourceId' : data['target']['id'].split('/')[-1],
            'diseaseFromSourceMappedId': data['disease']['id'].split('/')[-1],
            'resourceScore': float(data['unique_association_fields']['score']),
            'datasourceId': 'phenodigm',
            'datatypeId': data['type'],
            'diseaseFromSource': data['disease']['source_name'],
            'targetInModel': data['unique_association_fields']['model_gene_id'].split('/')[-1],
            'diseaseFromSourceId': data['unique_association_fields']['disease_phenodigm_id'],
            'biologicalModelAllelicComposition': data['unique_association_fields']['model_description'],
            'biologicalModelGeneticBackground': data['unique_association_fields']['model_genetic_background'],
            'diseaseModelAssociatedModelPhenotypes': [{'id': p['id'], 'label': p['label']} for p in data['evidence']['disease_model_association']['model_phenotypes']],
            'diseaseModelAssociatedHumanPhenotypes': [{'id': p['id'], 'label': p['label']} for p in data['evidence']['disease_model_association']['human_phenotypes']],
            'biologicalModelId': data['evidence']['disease_model_association']['model_id'].split('#')[0]
            })
        except:
            self._logger.warning(f'The following evidence could not be mapped to the new schema:\n{data}')
            return None

    def write_evidence_strings(self, filename):

        countExported = 0
        self._logger.info("Writing Phenodigm evidence strings")
        with gzip.open(filename, 'wt') as tp_file:
            self._logger.info("Processing %i records" % (len(self.hashkeys)))
            for hashkey in self.hashkeys:
                self._logger.info("Processing key %s"%(hashkey))
                evidenceString = self.hashkeys[hashkey]
                error = evidenceString.validate(self._logger)
                score = evidenceString.evidence.disease_model_association.resource_score.value
                if error == 0 and score >= 0.9:
                    new_evidence = self.convert_evidence(evidenceString)
                    if new_evidence:
                        tp_file.write(json.dumps(new_evidence) + "\n")
                        countExported+=1

        self._logger.info("Exported %i evidence" % (countExported))

    def process_ontologies(self):

        self._logger.info("Load MP classes")
        self.mp.load_mp_classes()

        self._logger.info("Load HP classes")
        self.hpo.load_hpo_classes()

        self._logger.info("Load EFO classes")
        self.efo.load_efo_classes()

    def get_ontology_mappings(self):

        self._logger.info("Get ontology mappings")
        self.get_omim_to_efo_mappings()
        self._logger.info("Get all Zooma mapping for Open Targets")
        self.get_opentargets_zooma_to_efo_mappings()
        #self.omim_to_efo_map["OMIM:191390"] = ["http://www.ebi.ac.uk/efo/EFO_0003767"]
        #self.omim_to_efo_map["OMIM:266600"] = ["http://www.ebi.ac.uk/efo/EFO_0003767"]
        #self.omim_to_efo_map["OM:612278"] = ["http://www.ebi.ac.uk/efo/EFO_0003767"]
        #self.omim_to_efo_map["OMIM:608049"] = ["http://www.ebi.ac.uk/efo/EFO_0003756"]
        #self.omim_to_efo_map["OMIM:300494"] = ["http://www.ebi.ac.uk/efo/EFO_0003757"]

    def process_all(self):
        self._logger.info('1. Fetch Ensembl gene ID mappings for human.')
        self.fetch_human_ensembl_mappings()

        self._logger.info('2. Fetch Ensembl gene ID mappings for mouse.')
        self.fetch_mouse_ensembl_mappings()

        self._logger.info('3. Process ontologies.')
        self.process_ontologies()
        self.get_ontology_mappings()

        self._logger.info('4. Fetch information from IMPC SOLR.')
        self.access_solr()

        self._logger.info('5. Build evidence')
        self.generate_phenodigm_evidence_strings()

        self._logger.info('6. Write evidence strings')
        self.write_evidence_strings(Config.MOUSEMODELS_EVIDENCE_FILENAME)


def main(log_file):
    # Initialize logger based on the provided log file. If no log file is specified, logs are written to STDERR.
    logging_config = {
        'level': logging.INFO,
        'format': '%(asctime)s %(levelname)s %(module)s - %(funcName)s: %(message)s',
        'datefmt': '%Y-%m-%d %H:%M:%S',
    }
    if log_file:
        logging_config['filename'] = log_file
    logging.basicConfig(**logging_config)

    # Process the data.
    Phenodigm(logging).process_all()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evidence parser for animal models sources from PhenoDigm.')
    parser.add_argument('-l', '--log-file', help='Optional filename to redirect the logs into.', required=False)
    args = parser.parse_args()
    main(args.log_file)
