# PhenoDigm evidence generators

**NOTE:** This branch contains the working version of the MouseModels module. For all other parsers use master.

### Install
Clone repository, create a virtual environment (requires python 3) and install dependencies.
```sh
# Clone solr_phenodigm_1904 branch
git clone -b solr_phenodigm_1904 https://github.com/opentargets/evidence_datasource_parsers.git

# Create and activate virtual environment
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
Additionally, the parser needs access to Google Cloud, which requires downloading a key JSON file and setting the `GOOGLE_APPLICATION_CREDENTIALS` variable.
### Usage

Read the [wiki page](https://github.com/opentargets/evidence_datasource_parsers/wiki/Running-PhenoDigm-evidence-parser) for a detailed explanation on how to run it on Google Cloud.
