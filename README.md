# OT evidence generators

Each folder in module corresponds to a datasource.

In each folder we have one or more standalone python scripts.

Generally these scripts:
1. map the disease terms (if any) to our ontology in various ways:
      - by using [OnToma](https://ontoma.readthedocs.io)
      - by using the [RareDiseasesUtils](https://github.com/opentargets/evidence_datasource_parsers/blob/master/common/RareDiseasesUtils.py) script
      - by using [Ontology Utils](https://github.com/opentargets/ontology-utils)
      - by importing manually curated files. Some of these are stored in the [mappings repo](https://github.com/opentargets/mappings)
2. Once the mapping is handled, evidence objects are generated in the form of JSON strings according to our JSON schema

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
The ClinGen parser processes the _Gene Validity Curations_ table that can be downloaded from https://search.clinicalgenome.org/kb/gene-validity. As of January 2021 the downloadable file contains six header lines that look as follows:
```tsv
CLINGEN GENE VALIDITY CURATIONS
FILE CREATED: 2021-01-18
WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity
++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++
GENE SYMBOL,GENE ID (HGNC),DISEASE LABEL,DISEASE ID (MONDO),MOI,SOP,CLASSIFICATION,ONLINE REPORT,CLASSIFICATION DATE,GCEP
+++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++
```

The mapping of the diseases is done on the fly as a three step process:
1. The MONDO id provided by ClinGen is used if it exists in EFO.
2. EFO is searched for xrefs to the MONDO id and all the EFO hits are used.
3. OnToma is used to search for perfect matches of the ClinGen disease name.

The unmapped diseases are saved to a file called `unmapped_diseases.tsv` so that they can be reported to [EFO](https://github.com/EBISPOT/efo/issues/).

The parser requires two parameters:
- `-i`, `--input_file`: Name of csv file downloaded from https://search.clinicalgenome.org/kb/gene-validity
- `-o`, `--output_file`: Name of evidence JSON file

There is also an optional parameter to save the unmapped diseases to a file:
- `-u`, `--unmapped_diseases_file`: If specified, the diseases not mapped to EFO will be stored in this file'

To use the parser configure the python environment and run it as follows:
```bash
(venv)$ python3 modules/ClinGen.py -i ClinGen-Gene-Disease-Summary-2020-08-04.csv -o clingen_2020-08-04.json -u unmapped_diseases_clingen.tsv
```

### Gene2Phenotype

The Gene2Phenotype parser processes the four gene panels (Developmental Disorders - DD, eye disorders, skin disorders and cancer) that can be downloaded from https://www.ebi.ac.uk/gene2phenotype/downloads/.

The mapping of the diseases, i.e. the "disease name" column, is done on the fly using [OnToma](https://pypi.org/project/ontoma/):
- Exact matches to EFO are used directly.
- If the match comes from HP or ORDO the disease is searched in MONDO and if there is an exact match it is used. If not the HP or ORDO term is used.
- If OnToma returns a fuzzy match it is ignore and MONDO is searched for exact matches.
- When no exact matches are found the disease is considered unmapped and it's saved to a file (see the `-u`/ `--unmapped_disease_file` option below).


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
(venv)$ python3 modules/Gene2Phenotype.py -d DDG2P_19_8_2020.csv.gz -e EyeG2P_19_8_2020.csv.gz -k SkinG2P_19_8_2020.csv.gz -c CancerG2P_19_8_2020.csv.gz -o gene2phenotype-19-08-2020.json -u gene2phenotype-19-08-2020_unmapped_diseases.txt 
```

### Genomics England Panel App

The Genomics England parser processes the associations between genes and diseases described in the _Gene Panels Data_ table. This data is provided by Genomics England and can be downloaded [here](https://storage.googleapis.com/otar000-evidence_input/PanelApp/20.11/All_genes_20200928-1959.tsv) from the _otar000-evidence_input_ bucket.

The source table is then formatted into a compressed set of JSON lines following the schema of the version to be used.

The mapping of the diseases is done on the fly using [OnToma](https://pypi.org/project/ontoma/):
1. Exact matches to an EFO term are used directly.
2. Sometimes an OMIM code can be present in the disease string. OnToma is then queried for both the OMIM code and the respective disease term. If OnToma returns a fuzzy match for both, it is checked whether they both point to the same EFO term. Being this the case, the term is considered as an exact match.

By default the result of the diseases and codes mappings are stored locally as of _disease_queries.json_ and _disease_queries.json_ respectively. This is intended for analysys purposes and to ease up a potential rerun of the parser.

The parser requires three parameters:
- `-i`, `--input_file`: Name of tsv file located in the [Panel App bucket](https://storage.googleapis.com/otar000-evidence_input/PanelApp/20.11/All_genes_20200928-1959.tsv).
- `-o`, `--output_file`: Name of evidence JSON file containing the evidence strings.
- `-s`, `--schema_version`: JSON schema version to use, e.g. 1.7.5. It must be branch or a tag available in https://github.com/opentargets/json_schema.
  
There is also an optional parameter to load a dictionary containing the results of querying OnToma with the disease terms:
- `-d`, `--dictionary`: If specified, the diseases mappings will be imported from this JSON file.'

To use the parser configure the python environment and run it as follows:
```bash
(venv)$ python3 modules/GenomicsEnglandPanelApp.py -i All_genes_20200928-1959.tsv -o genomics_england-2021-01-05.json -s 1.7.5 -d disease_queries.json
```

### IntOGen

The intOGen parser generates evidence strings from three files that need to be in the working directory or in the _resources_ directory:

- _intogen_cancer2EFO_mapping.tsv_: Mappings of intOGen cancer type acronyms to EFO terms. Currently is a manually curated file given that in intOGen they do not use an specific ontology to model cancers.
- _intogen_cohorts.tsv_: It contains information about the analysed cohorts and it can be downloaded from the [intOGen website](https://www.intogen.org/download). In the current implementation, the total number of samples included in each cohort is used to calculate the percentage of samples that carry a relevant mutation in a driver gene.
- _intogen_Compendium_Cancer_Genes.tsv_: It contains the genes that have been identified as _drivers_ in at least one cohort and information related to  the methods that yield significant results, the q-value and additional data about mutations. It can be downloaded from the same place as the cohorts file.

The parser uses the following parameters:
- `-g`, `--inputGenes`: Input source .tsv file listing the driver genes across the analyzed cohorts.
- `-c`, `--inputCohorts`: Input source .tsv file with information about the analyzed cohorts.
- `-d`, `--diseaseMapping`: optional; input look-up table containing the cancer type mappings to an EFO ID.
- `-o`, `--outputFile`: Gzipped JSON file containing the evidence strings.
- `-s`, `--skipMapping`: optional; state whether to skip the disease to EFO term mapping step. If used this step is not performed.
- `-l`, `--logFile`: optional; if not specified, logs are written to standard error.

To use the parser configure the python environment and run it as follows:
```bash
python modules/IntOGen.py \
    --inputGenes Compendium_Cancer_Genes.tsv \
    --inputCohorts cohorts.tsv \
    --diseaseMapping resources/cancer2EFO_mappings.tsv \
    --outputFile intogen-2021-03-08.json.gz
```

### PheWAS catalog

The PheWAS parser processes three files:
- phewas-catalog-19-10-2018.csv: Main tsv file containing information from phenome-wide association studies. The file is located in the [PheWAS bucket](https://storage.googleapis.com/otar000-evidence_input/PheWAS/data_files/phewas-catalog-19-10-2018.csv).
- phewas_w_consequences.csv: File containing PheWAS data enriched with data from the Genetics Portal on the variant and its functional consequence. This file is located in the [PheWAS bucket](https://storage.googleapis.com/otar000-evidence_input/PheWAS/data_files/phewas_w_consequences.csv).
- phewascat.mappings.tsv: File containing the mappings between the phenotypes provided by PheWAS and its respective disease listed in EFO. This file is located in [our mappings repo](https://raw.githubusercontent.com/opentargets/mappings/master/phewascat.mappings.tsv).

The source table is then formatted into a compressed set of JSON lines following the schema of the version to be used.

The parser uses the following parameters:
- `-i`, `--inputFile`: Main tsv file coming from PheWAS.
- `-c`, `--consequencesFile`: Input look-up table containing the variant data and consequences coming from the Variant Index.
- `-d`, `--diseaseMapping`: optional; input look-up table containing the PheWAS phenotypes mappings to an EFO IDs.
- `-s`, `--skipMapping`: optional; state whether to skip the disease to EFO term mapping step. If used this step is not performed.
- `-o`, `--outputFile`: Gzipped JSON file containing the evidence strings.
- `-l`, `--logFile`: optional; if not specified, logs are written to standard error.

To use the parser configure the python environment and run it as follows:
```bash
python modules/PheWAS.py \
    --inputFile phewas-catalog-19-10-2018.csv \
    --consequencesFile https://storage.googleapis.com/otar000-evidence_input/PheWAS/data_files/phewas_w_consequences.csv) \
    --diseaseMapping https://raw.githubusercontent.com/opentargets/mappings/master/phewascat.mappings.tsv \
    --outputFile phewas-2021-02-02.json.gz
```

### CRISPR

The CRISPR parser processes three files available in the `resources` directory:

- _crispr_evidence.tsv_: Main file that contains the prioritisation score as well as target, disease and publication information. It was adapted from supplementary table 6 in the [Behan et al. 2019](https://www.nature.com/articles/s41586-019-1103-9) paper.
- _crispr_descriptions.tsv_: File used to extract the number of targets prioritised per cancer types as well as to retrieve the cell lines from the next file.
- _crispr_cell_lines.tsv_: It contains three columns with the cell line, tissue and the cancer type information. It is used to extract the cell lines in which the target has been identified as essential. The file has been adapted from the supplementary table 1 in the [Behan et al. 2019](https://www.nature.com/articles/s41586-019-1103-9) paper.

*Usage:*

```sh
(venv)$ python3 modules/CRISPR.py -e <evidence_file> -d <description_file> -c <cell_line_file> -o <output_file> -l <log_file>
```

- `-e`, `--evidence_file`: name or full path for _crispr_evidence.tsv_.
- `-d`, `--descriptions_file`: name or full path for the _crispr_descriptions.tsv_ file.
- `-c`, `--cell_types_file`: name or full path for _crispr_cell_lines.tsv_.
- `-o`, `--output_file`: output is a gzipped JSON.
- `-l`, `--log_file`:  optional parameter. If not provided, logs are written to the standard error.

### PROGENy

The PROGENy parser processes three files:

- progeny_normalVStumor_opentargets.txt: Main file that contains a systematic comparison of pathway activities between normal and primary TCGA samples in 14 tumor types using PROGENy. This file can be downloaded [here](https://storage.googleapis.com/otar000-evidence_input/PROGENy/data_files/progeny_normalVStumor_opentargets.txt) from the _otar000-evidence_input_ bucket.
- cancer2EFO_mappings.tsv: File containing the mappings between the acronym of the type of cancer and its respective disease listed in EFO. This file can be found in the `resources` directory.
- pathway2Reactome_mappings.tsv: File containing the mappings between the analysed pathway and their respective target and ID in Reactome. This file can be found in the `resources` directory.

The source table is then formatted into a compressed set of JSON lines following the schema of the version to be used.

The parser uses the following parameters:
- `-i`, `--inputFile`: Main tsv file.
- `-d`, `--diseaseMapping`: optional; input look-up table containing the cancer type mappings to an EFO ID.
- `-s`, `--skipMapping`: optional; state whether to skip the disease to EFO term mapping step. If used this step is not performed.
- `-p`, `--pathwayMapping`: input look-up table containing the pathway mappings to a respective target and ID in Reactome.
- `-o`, `--outputFile`: gzipped JSON file containing the evidence strings.


To use the parser configure the python environment and run it as follows:
```bash
python modules/PROGENY.py \
    --inputFile progeny_normalVStumor_opentargets.txt \
    --diseaseMapping resources/cancer2EFO_mappings.tsv \
    --pathwayMapping resources/pathway2Reactome_mappings.tsv \
    --outputFile progeny-2021-01-18.json.gz
```

### SystemsBiology

This parser processes two files stored in Google cloud bucket: `gs://otar000-evidence_input/SysBio/data_files`. The `sysbio_evidence-31-01-2019.tsv` file contains evidence data, while `sysbio_publication_info_nov2018.tsv` contains study level information.

**Usage:**

```bash
python modules/SystemsBiology.py \
    --evidenceFile sysbio_evidence-31-01-2019.tsv \
    --studyFile    sysbio_publication_info_nov2018.tsv \
    --outputFile   sysbio_evidence.json.gz \
    --logFile      sysbio.parser.log
```

### SLAPenrich

The SLAPenrich parser processes twofiles:

- `slapenrich_opentargets.tsv`: Main file that contains a systematic comparison of on somatic mutations from TCGA across 25 different cancer types and a collection of pathway gene sets from Reactome. This file can be downloaded [here](https://storage.googleapis.com/otar000-evidence_input/SLAPEnrich/data_file/slapenrich_opentargets-21-12-2017.tsv) from the _otar000-evidence_input_ bucket.
- `cancer2EFO_mappings.tsv`: File containing the mappings between the acronym of the type of cancer and its respective disease listed in EFO. This file can be found in the `resources` directory.

The source table is then formatted into a compressed set of JSON lines following the schema of the version to be used.

The parser uses the following parameters:

- `-i`, `--inputFile`: Name of tsv file located in the [SLAPEnrich bucket](https://storage.googleapis.com/otar000-evidence_input/SLAPEnrich/data_file/slapenrich_opentargets-21-12-2017.tsv).
- `-d`, `--diseaseMapping`: optional; input look-up table containing the cancer type mappings to an EFO ID.
- `-s`, `--skipMapping`: optional; state whether to skip the disease to EFO term mapping step. If used this step is not performed.
- `-o`, `--outputFile`: gzipped JSON file containing the evidence strings.
- `-l`, `--logFile`: optional; if not specified, logs are written to standard error.

To use the parser configure the python environment and run it as follows:
```bash
python modules/SLAPEnrich.py \
    --inputFile slapenrich_opentargets.tsv \
    --diseaseMapping resources/cancer2EFO_mappings.tsv \
    --outputFile slapenrich-2021-01-18.json.gz
```

### PhenoDigm

Generates target-disease evidence querying the IMPC SOLR API.

The PhenoDigm parser has been run both in a MacBook Pro and a Google Cloud machine. The current version is very memory greedy, using up to 60 GB. 

* MacBook Pro: It takes around 8 hours to run in a 16 GB laptop.
* Google Cloud Machine: It runs in around 2 hours in a 60 GB machine (_n1-standard-16_). There used to be such a machine set up but now one called _ag-ubuntu-20-04-lts_ in _open-targets-eu-dev_ can be used instead (see below). This machine has specs higher than needed (32 vCPU and 208 GB memory).

```sh
# Prepare the virtual environment
python3 -m venv phenodigm_venv
source phenodigm_venv/bin/activate
pip3 install -r requirements.txt
export PYTHONPATH=.

# Manually download and install ontology-utils as it requires specific tag
wget https://github.com/opentargets/ontology-utils/archive/bff0f189a4c6e8613e99a5d47e9ad4ceb6a375fc.zip
pip3 install bff0f189a4c6e8613e99a5d47e9ad4ceb6a375fc.zip

# Run
python3 modules/MouseModels.py -l log.txt
```

If `-l` is unspecified, logs will be printed to STDERR.

### Open Targets Genetics Portal

Genetics portal evidence generation is done by calling this Python script:

```bash
python modules/GeneticsPortal.py \
    --locus2gene gs://genetics-portal-data/l2g/200127 \
    --toploci gs://genetics-portal-data/v2d/200207/toploci.parquet \
    --study gs://genetics-portal-data/v2d/200207/studies.parquet \
    --variantIndex gs://genetics-portal-data/variant-annotation/190129/variant-annotation.parquet \
    --ecoCodes gs://genetics-portal-data/lut/vep_consequences.tsv \
    --outputFile gentics-portal-evidence.json.gz \
    --logFile gentics-portal-evidence.log \
    --threshold 0.05
```

**Where:**

* `--threshold` is required. It provides a lower locus to gene score cutoff.
* `--logFile` is optional. If not specified, logs are written to standard error.



