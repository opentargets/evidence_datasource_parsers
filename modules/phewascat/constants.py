PHEWAS_CATALOG_URL = 'https://storage.googleapis.com/otar000-evidence_input/PheWAScatalog/phewas-catalog.csv'

PHEWAS_PHECODE_MAP_URL = 'https://phewascatalog.org/files/phecode_icd9_map_unrolled.csv.zip'

PHEWAS_SCHEMA = 'https://raw.githubusercontent.com/opentargets/json_schema/ep-fixrelative/src/genetics.json'

'''define once where evidence is coming from
'''
PROVENANCE = {'literature': {
            "references":[{"lit_id":"http://europepmc.org/articles/PMC3969265"}]
            },
             "database":{
                 "version":"2013-12-31T09:53:37+00:00",
                "id":"PHEWAS Catalog"
                }
}

TOTAL_NUM_PHEWAS = 4269549
