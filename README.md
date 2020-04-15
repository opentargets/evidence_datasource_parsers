# OT evidence generators

Each folder in module corresponds corresponds to a datasource.

In each folder we have one or more standalone python scripts.

Generally these scripts:
1. map the disease terms (if any) to our ontology, sometimes using [OnToma](https://ontoma.readthedocs.io)
2. save the mappings in https://github.com/opentargets/mappings
3. Read the **github mappings** to generate evidence objects (JSON strings) according to our JSON schema

Code used by more than one script (that does not live in a python package)
is stored in the `common` folder and imported as follows:

```python
from common.<module> import <function>
```



### Install
Install (requires python 3):

```sh
virtualenv -p python3 venv
source venv/bin/activate
pip3 install -r requirements.txt
export PYTHONPATH=.
```
### Usage

Each script is a standalone python script.
Common dependencies are stored in the `common` folder.

Hence to run each parser, simply run the standalone script with your python
interpreter:
```sh
(venv)$ python3 modules/<parser you want>.py
```

### Contributor guidelines

Further development of this repository should follow the next premises:

1. Data QC
1. Configuration
1. Documentation
1. Logging
1. Standardize disease mappings
1. Reliability (e.g. network retries)
1. JSON schema validation
1. Uniqueness of unique_association_fields
1. Reproducibility (e.g. saving API results as files)
1. Containerization - Dockerize


### PheWAS catalog

The `PheWAS.py` script parses the PheWAS Catalog CSV file specified as `PHEWAS_CATALOG_FILENAME` in `settings.py` and that should be located either in the working directory or the `resources` folder. The mappings between the Phecodes and EFO are read from the [phewascat.mappings.tsv](https://raw.githubusercontent.com/opentargets/mappings/master/phewascat.mappings.tsv) file in the `mappings` repository.

```sh
(venv)$ python3 modules/PheWAS.py
```

### PhenoDigm - What version should I run?
The PhenoDigm parser used to generate the 19.04 release data is the [solr_phenodigm_1904](https://github.com/opentargets/evidence_datasource_parsers/tree/solr_phenodigm_1904) branch.

The version in _master_ is an older version (October 2018) that **DOES NOT** call the IMPC SOLR API and it is **unlikely to work** but it has not been tested.

[solr_phenodigm](https://github.com/opentargets/evidence_datasource_parsers/tree/solr_phenodigm) is the version that Gautier handed over to the OT data team in February 2019. It *DOES* call the IMPC SOLR API but it has a number of bugs and **DOES NOT WORK**.

**TODO**
- [ ] map `intergenic` rsIDs to genes (~900k evidences)
- [ ] improve mappings with manual curation

