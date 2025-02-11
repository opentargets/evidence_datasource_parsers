# Open Targets internal evidence generation

This repository contains a collection of modules which generate evidence for several internal data sources (“internal” meaning that the code is maintained by the data team; the data itself usually comes from sources outside Open Targets).

## How to set up and update the environment

This project uses `uv` for fast and reliable dependency management. The Python version is pinned to **3.11**, and dependencies are managed via `pyproject.toml`.

By executing the `setup.sh` script, you will install the required Python version, set up a virtual environment, and install the project dependencies:

```bash
bash setup.sh
```


1. **Install `uv` (if not already installed):**

   Run the following command to install `uv`:

```bash
# Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.${SHELL##*/}rc
```
## How to generate the evidence

This will create a Google Cloud instance, SSH into it, install the necessary dependencies, generate, validate, and upload the evidence. Tweak the commands as necessary.

To run this, conditions related to the service accounts need to be satisfied:

1. The service account must have a Storage Admin role for two buckets, `_otar000-evidence_input_` and `_otar001-core_`.
2. The service account must have a Compute Admin and Service Account User roles in the _open-targets-eu-dev_ project.
3. The user running the code must have access to use the service account.

The schema version which the evidence is validated against can be tweaked in [`configuration.yaml`](configuration.yaml) → global → schema.

```bash
# Set parameters.
export INSTANCE_NAME=evidence-generation
export INSTANCE_PROJECT=open-targets-eu-dev
export INSTANCE_ZONE=europe-west1-d
# Create the instance and SSH.
gcloud compute instances create \
  ${INSTANCE_NAME} \
  --project=${INSTANCE_PROJECT} \
  --zone=${INSTANCE_ZONE} \
  --machine-type=n1-highmem-32 \
  --service-account=426265110888-compute@developer.gserviceaccount.com \
  --scopes=https://www.googleapis.com/auth/cloud-platform \
  --create-disk=auto-delete=yes,boot=yes,device-name=${INSTANCE_NAME},image=projects/ubuntu-os-cloud/global/images/ubuntu-2004-focal-v20210927,mode=rw,size=2000,type=projects/open-targets-eu-dev/zones/europe-west1-d/diskTypes/pd-balanced
gcloud compute ssh --project ${INSTANCE_PROJECT} --zone ${INSTANCE_ZONE} ${INSTANCE_NAME}

screen

# Install the system dependencies.
sudo apt update
sudo apt install -y openjdk-8-jdk-headless python3-pip python3.8-venv r-base-core

# Activate the environment and install Python dependencies.
git clone https://github.com/opentargets/evidence_datasource_parsers
cd evidence_datasource_parsers
pip install uv
uv sync --all-groups --frozen
export PYTHONPATH="$PYTHONPATH:$(pwd)"

# Workaround for a potential OnToma race condition: pre-initialise cache directory.
# This prevents an issue where several OnToma instances are trying to do this at once and fail.
echo 'asthma' | ontoma --cache-dir cache_dir
```

At this point, we are ready to run the Snakemake pipeline. The following options are available:

* `snakemake --cores all`: Display help (the list of possible rules to be run) and do not run anything.
* `snakemake --cores all --until local`: Generate all files, but do not upload them to Google Cloud Storage. The files generated in this way do not have prefixes, e.g. `cancer_biomarkers.json.gz`. This is done intentionally, so that the pipeline can be re-run the next day without having to re-generate all the files.
  * It is also possible to locally run only a single rule by substituting its name instead of “local”.
* `snakemake --cores all --until all`: Generate all files and then upload them to Google Cloud Storage.

All individual parser rules are strictly local. The only rule which uploads files to Google Cloud Storage (all at once) is "all".

Some additional parameters which can be useful for debugging:

* `--keep-incomplete`: This will keep the output files of failed jobs. Must only be used with local runs. Note that Snakemake uses the presence of the output files to decide which jobs to run, so the incomplete files must be removed after investigation and before any re-runs of the workflow.
* `--dry-run`: Do not run the workflow, and only show the list of jobs to be run.

## Notes on how this repository is organised

Each module in [`src/modules/`](src/modules/) corresponds to one evidence generator.

Modules which are shared by multiple generators reside in [`src/common/`](src/common/).

## Historical notes on individual parsers

_Note that some information in this section may be obsolete. It is provided for historical purposes, and will be eventually migrated into the source code of individual parsers or into the Snakefile._

### ClinGen

The ClinGen parser processes the _Gene Validity Curations_ table that can be downloaded from <https://search.clinicalgenome.org/kb/gene-validity>. As of January 2021 the downloadable file contains six header lines that look as follows:

```tsv
CLINGEN GENE VALIDITY CURATIONS
FILE CREATED: 2021-01-18
WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity
++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++
GENE SYMBOL,GENE ID (HGNC),DISEASE LABEL,DISEASE ID (MONDO),MOI,SOP,CLASSIFICATION,ONLINE REPORT,CLASSIFICATION DATE,GCEP
+++++++++++,++++++++++++++,+++++++++++++,++++++++++++++++++,+++++++++,+++++++++,+++++++++,++++++++++++++,+++++++++++++,+++++++++++++++++++
```

### Gene2Phenotype

The Gene2Phenotype parser processes the four gene panels (Developmental Disorders - DD, eye disorders, skin disorders and cancer) that can be downloaded from <https://www.ebi.ac.uk/gene2phenotype/downloads/>.

### Genomics England Panel App

The Genomics England parser processes the associations between genes and diseases described in the _Gene Panels Data_ table. This data is provided by Genomics England and can be downloaded [here](https://storage.googleapis.com/otar000-evidence_input/PanelApp/20.11/All_genes_20200928-1959.tsv) from the _otar000-evidence_input_ bucket.

The source table is then formatted into a compressed set of JSON lines following the schema of the version to be used.

### IntOGen

The intOGen parser generates evidence strings from three files that need to be in the working directory or in the _resources_ directory:
* _intogen_cancer2EFO_mapping.tsv_: Mappings of intOGen cancer type acronyms to EFO terms. Currently is a manually curated file given that in intOGen they do not use an specific ontology to model cancers.
* _intogen_cohorts.tsv_: It contains information about the analysed cohorts and it can be downloaded from the [intOGen website](https://www.intogen.org/download). In the current implementation, the total number of samples included in each cohort is used to calculate the percentage of samples that carry a relevant mutation in a driver gene.
* _intogen_Compendium_Cancer_Genes.tsv_: It contains the genes that have been identified as _drivers_ in at least one cohort and information related to  the methods that yield significant results, the q-value and additional data about mutations. It can be downloaded from the same place as the cohorts file.

### CRISPR

The CRISPR parser processes three files available in the `resources` directory:
* _crispr_evidence.tsv_: Main file that contains the prioritisation score as well as target, disease and publication information. It was adapted from supplementary table 6 in the [Behan et al. 2019](https://www.nature.com/articles/s41586-019-1103-9) paper.
* _crispr_descriptions.tsv_: File used to extract the number of targets prioritised per cancer types as well as to retrieve the cell lines from the next file.
* _crispr_cell_lines.tsv_: It contains three columns with the cell line, tissue and the cancer type information. It is used to extract the cell lines in which the target has been identified as essential. The file has been adapted from the supplementary table 1 in the [Behan et al. 2019](https://www.nature.com/articles/s41586-019-1103-9) paper.

### PROGENy

The PROGENy parser processes three files:
* progeny_normalVStumor_opentargets.txt: Main file that contains a systematic comparison of pathway activities between normal and primary TCGA samples in 14 tumor types using PROGENy. This file can be downloaded [here](https://storage.googleapis.com/otar000-evidence_input/PROGENy/data_files/progeny_normalVStumor_opentargets.txt) from the _otar000-evidence_input_ bucket.
* cancer2EFO_mappings.tsv: File containing the mappings between the acronym of the type of cancer and its respective disease listed in EFO. This file can be found in the `resources` directory.
* pathway2Reactome_mappings.tsv: File containing the mappings between the analysed pathway and their respective target and ID in Reactome. This file can be found in the `resources` directory.

### SystemsBiology

This parser processes two files stored in Google cloud bucket: `gs://otar000-evidence_input/SysBio/data_files`. The `sysbio_evidence-31-01-2019.tsv` file contains evidence data, while `sysbio_publication_info_nov2018.tsv` contains study level information.

### SLAPenrich

The SLAPenrich parser processes two files:
* `slapenrich_opentargets.tsv`: Main file that contains a systematic comparison of on somatic mutations from TCGA across 25 different cancer types and a collection of pathway gene sets from Reactome. This file can be downloaded [here](https://storage.googleapis.com/otar000-evidence_input/SLAPEnrich/data_file/slapenrich_opentargets-21-12-2017.tsv) from the _otar000-evidence_input_ bucket.
* `cancer2EFO_mappings.tsv`: File containing the mappings between the acronym of the type of cancer and its respective disease listed in EFO. This file can be found in the `resources` directory.

### IMPC

Generates the mouse model target-disease evidence by querying the IMPC SOLR API.

The base of the evidence is the `disease_model_summary` table, which is unique on the combination of (`model_id`, `disease_id`). When target information is added, an original row may explode into multiple evidence strings. As a result, the final output is unique on the combination of (`biologicalModelId`, `targetFromSourceId`, `targetInModelId`, `diseaseFromSourceId`).



## How to set up and update the environment

The file `requirements.txt` contains the **direct** dependencies only with their exact versions pinned. Only this file should be edited directly.

Additionally, the file `requirements-frozen.txt` contains all **direct and transitive** dependencies with their exact versions pinned. It is generated automatically from the former file, and is intended for reproducing the exact working environment (as much as possible) as of any particular commit version. This file should not be directly edited.

These are the steps to update an environment:

* Add, delete or update packages as necessary in `requirements.txt`
* Install requirements with `uv pip install -r requirements.txt`
* Make sure everything works
* Update the frozen file with `uv pip freeze > requirements-frozen.txt`
* Include both files, `requirements.txt` and `requirements-frozen.txt`, with your PR

Make sure to always work in a clean virtual environment to avoid any surprises.