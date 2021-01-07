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

### ClinGen
The ClinGen parser processes the _Gene Validity Curations_ table that can be downloaded from https://search.clinicalgenome.org/kb/gene-validity. At the time the parser was written (August 2020), the downloadable file contained six header lines that look as follows:
```tsv
CLINGEN GENE VALIDITY CURATIONS
FILE CREATED: 2020-08-04
WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity
+++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++
GENE SYMBOL,GENE ID (HGNC),DISEASE LABEL,DISEASE ID (MONDO),MOI,SOP,CLASSIFICATION,ONLINE REPORT,CLASSIFICATION DATE
+++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++

```

The mapping of the diseases is done on the fly as a three step process:
1. The MONDO id provided by ClinGen is used if it exists in EFO.
2. EFO is searched for xrefs to the MONDO id and all the EFO hits are used.
3. OnToma is used to search for perfect matches of the ClinGen disease name.

The unmapped diseases are saved to a file called `unmapped_diseases.tsv` so that they can be reported to [EFO](https://github.com/EBISPOT/efo/issues/).

The parser requires three parameters:
- `-i`, `--input_file`: Name of csv file downloaded from https://search.clinicalgenome.org/kb/gene-validity
- `-o`, `--output_file`: Name of evidence JSON file
- `-s`, `--schema_version`: JSON schema version to use, e.g. 1.6.8. It must be branch or a tag available in https://github.com/opentargets/json_schema

There is also an optional parameter to save the unmapped diseases to a file:
- `-u`, `--unmapped_diseases_file`: If specified, the diseases not mapped to EFO will be stored in this file'

To use the parser configure the python environment and run it as follows:
```bash
(venv)$ python3 modules/ClinGen.py -i ClinGen-Gene-Disease-Summary-2020-08-04.csv -o clingen_2020-08-04.json -s 1.6.9 -u unmapped_diseases_clingen.tsv
```

### Gene2Phenotype
The Gene2Phenotype parser processes the four gene panels (Developmental Disorders - DD, eye disorders, skin disorders and cancer) that can be downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/.

The mapping of the diseases, i.e. the "disease name" column, is done on the fly using [OnToma](https://pypi.org/project/ontoma/):
- Exact matches to EFO are used directly.
- If the match comes from HP or ORDO the disease is searched in MONDO and if there is an exact match it is used. If not the HP or ORDO term is used.
- If OnToma returns a fuzzy match it is ignore and MONDO is searched for exact matches.
- When no exact matches are found the disease is considered unmapped and it's saved to a file (see the `-u`/ `--unmapped_disease_file` option below).

The parser requires two parameters to run:
- `-s`, `--schema_version`: JSON schema version to use, e.g. 1.6.8. It must be branch or a tag available in https://github.com/opentargets/json_schema.
- `-v`, `--g2p_version`: Version of the Gene2Phenotype data used. If not available please use the date in which the data was downloaded in YYYY-MM-DD format.

There are also a number of optional parameters to specify the name of the input and out files:
- `-d`, `--dd_panel`: Name of Developmental Disorders (DD) panel file. It uses the value of G2P_DD_FILENAME in setting.py if not specified.
- `-e`, `--eye_panel`: Name of eye panel file. It uses the value of G2P_eye_FILENAME in setting.py if not specified.
- `-k`, `--skin_panel`: Name of skin panel. It uses the value of G2P_skin_FILENAME in setting.py if not specified.
- `-c`, `--cancer_panel`: Name of cancer panel file. It uses the value of G2P_cancer_FILENAME in setting.py if not specified.
- `-o`, `--output_file`: Name of output evidence file. It uses the value of G2P_EVIDENCE_FILENAME in setting.py if not specified.
- `-u`, `--unmapped_diseases_file`: If specified, the diseases not mapped to EFO will be stored in this file.

Note that when using the default file names, the input files have to exist in the working directory or in the _resources_ directory:

To use the parser configure the python environment and run it as follows:
```bash
(venv)$ python3 modules/Gene2Phenotype.py -s 1.7.1 -v 2020-08-19 -d DDG2P_19_8_2020.csv.gz -e EyeG2P_19_8_2020.csv.gz -k SkinG2P_19_8_2020.csv.gz -c CancerG2P_19_8_2020.csv.gz -o gene2phenotype-19-08-2020.json -u gene2phenotype-19-08-2020_unmapped_diseases.txt 
```

### IntOGen

The intOGen parser generates evidence strings from three files that need to be in the working directory or in the _resources_ directory:

- _intogen_cancer2EFO_mapping.tsv_: Mappings of intOGen cancer type acronyms to EFO terms. Currently is a manually curated file given that in intOGen they do not use an specific ontology to model cancers.
- _intogen_cohorts.tsv_: It contains information about the analysed cohorts and it can be downloaded from the [intOGen website](https://www.intogen.org/download). In the current implementation, the total number of samples included in each cohort is used to calculate the percentage of samples that carry a relevant mutation in a driver gene.
- _intogen_Compendium_Cancer_Genes.tsv_: It contains the genes that have been identified as _drivers_ in at least one cohort and information related to  the methods that yield significant results, the q-value and additional data about mutations. It can be downloaded from the same place as the cohorts file.

### PheWAS catalog

The `PheWAS.py` script parses the PheWAS Catalog CSV file specified as `PHEWAS_CATALOG_FILENAME` in `settings.py` and that should be located either in the working directory or the `resources` folder. The mappings between the Phecodes and EFO are read from the [phewascat.mappings.tsv](https://raw.githubusercontent.com/opentargets/mappings/master/phewascat.mappings.tsv) file in the `mappings` repository.

```sh
(venv)$ python3 modules/PheWAS.py
```

### CRISPR

The CRISPR parser processes three files available in the `resources` directory:

- _crispr_evidence.tsv_: Main file that contains the prioritisation score as well as target, disease and publication information. It was adapted from supplementary table 6 in the [Behan et al. 2019](https://www.nature.com/articles/s41586-019-1103-9) paper.
- _crispr_descriptions.tsv_: File used to extract the number of targets prioritised per cancer types as well as to retrieve the cell lines from the next file.
- _crispr_cell_lines.tsv_: It contains three columns with the cell line, tissue and the cancer type information. It is used to extract the cell lines in which the target has been identified as essential. The file has been adapted from the supplementary table 1 in the [Behan et al. 2019](https://www.nature.com/articles/s41586-019-1103-9) paper.

The parser has five parameters but only two of them are compulsory, the other three have defaults specified in the `settings.py`:

- `-i`, `--input_file`: Name or full path for _crispr_evidence.tsv_. The default value is taken from the `CRISPR_FILENAME1` variable in `settings.py`.
- `-d`, `--description_file`: Name or full path for the _crispr_descriptions.tsv_ file. The default value is taken from the `CRISPR_FILENAME2` variable in `settings.py`.
- `-c`, `--cell_types_file`: Name or full path for _crispr_cell_lines.tsv_. There is no default value, it **must** be specified.
- `-o`, `--output_file`:  Name of the evidence JSON file. The default value is taken from the `CRISPR_EVIDENCE_FILENAME` variable in `settings.py`.
- `-s`, `--schema_version`:  JSON schema version to use, e.g. 1.7.5. It must be branch or a tag available in https://github.com/opentargets/json_schema. There is no default value, it **must** be specified.

This is how it is run only specifying the required parameters:

```sh
(venv)$ python3 modules/CRISPR.py -s 1.7.5 -c crispr_cell_lines.tsv -o crispr-18-11-2020.json
```

### PhenoDigm

Generates target-disease evidence querying the IMPC SOLR API.

#### System requirements and running time
The PhenoDigm parser has been run both in a MacBook Pro and a Google Cloud machine. The current version is very memory greedy, using up to 60 GB. 

* MacBook Pro: It takes around 8 hours to run in a 16 GB laptop.
* Google Cloud Machine: It runs in around 2 hours in a 60 GB machine (_n1-standard-16_). There used to be such a machine set up but now one called _ag-ubuntu-20-04-lts_ in _open-targets-eu-dev_ can be used instead (see below). This machine has specs higher than needed (32 vCPU and 208 GB memory).

#### Installation
Before the parser can be run there is one dependency that needs to be installed manually:

```sh
# Create and activate virtual environment (requires python3)
virtualenv -p python3 phenodigm_venv
source phenodigm_venv/bin/activate

# Install dependencies
pip3 install -r requirements.txt

# Manually download and install ontology-utils as it requires specific tag
wget https://github.com/opentargets/ontology-utils/archive/bff0f189a4c6e8613e99a5d47e9ad4ceb6a375fc.zip
pip3 install bff0f189a4c6e8613e99a5d47e9ad4ceb6a375fc.zip

# Set environment variables
export PYTHONPATH=.
```
Additionally, the parser needs access to Google Cloud, which requires downloading a key JSON file and setting the `GOOGLE_APPLICATION_CREDENTIALS` variable. Ask someone from the back-end team about it.

#### Running MouseModels in a Google Cloud virtual machine

Start the machine and connect to it via SSH. When you are ready set-up the environment and run the parser:

```sh
# Move to installation directory
cd opentargets/evidence_datasource_parsers/

# Activate virtual environment
. phenodigm_venv/bin/activate

# Set up PYTHONPATH AND GOOGLE_APPLICATION_CREDENTIALS environment variables
. phenodigm_venv/bin/set_env_variables.sh

# Update cache
# NOTE: This is only required once per release, i.e. is not needed if PhenoDigm is
# rerun multiple times in a short period of time
python3 modules/MouseModels.py --update-cache

# Run PhenoDigm parser
python3 modules/MouseModels.py 

# Upload the evidence file to Google Cloud Storage
python3 modules/MouseModels.py -w
```
**NOTE:** Remember to stop the machine once you are done as it costs money to have it on! 

### Open Targets Genetics Portal

The evidence generation is done in two steps:
1. Pooling together data from Genetics portal buckets: `GeneticsPortal_prepare_data.py`
2. The resulting table is formatted into compressed set of JSON lines: `GeneticsPortal.py`

The first step is run on Google Cloud Dataproc cluster to take advantage of the proximity of the Genetics Portal data stored in Google buckets. Once you have the proper cluster set up, submit the job:

```bash
gcloud dataproc jobs submit pyspark \
    --cluster=${clusterName} \
    --project=${projectName} \
    --region=${region} \
    ${repo_directory}/modules/GeneticsPortal_prepare_data.py -- \
    --locus2gene=gs://genetics-portal-data/l2g/200127 \
    --toploci=gs://genetics-portal-data/v2d/200207/toploci.parquet \
    --study=gs://genetics-portal-data/v2d/200207/studies.parquet \
    --variantIndex=gs://genetics-portal-data/variant-annotation/190129/variant-annotation.parquet \
    --ecoCodes=gs://genetics-portal-data/lut/vep_consequences.tsv \
    --output=gs://genetics-portal-analysis/l2g-platform-export/data/l2g_joined.2020.03.02_exploded.parquet
```

The output is saved in a parquet file in one of the buckets. This then needs to be further processed:

```bash
python ${repo_directory}/modules/GeneticsPortal.py \
    --inputFile l2g_joined.2020.01.30_exploded.parquet \
    --schemaFile https://raw.githubusercontent.com/opentargets/json_schema/Draft-4_compatible/opentargets.json \
    --cores 32 --sample 1000 --outputFile output.json.gz
```

**Important**: to ensure the resulting json schema is valid, we are using the [python_jsonschema_objects](https://pypi.org/project/python-jsonschema-objects/0.0.13/) library, which enforces the proper structure. The only caveat is that the this library uses draft-4 JSON schema, while our JSON schema is written on draft-7. To resolve this discrepancy, our JSON schema repository has a parallel draft-4 compatible branch that we are using for evidence generation.
