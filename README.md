# PhenoDigm evidence generator

**NOTE:** This branch contains the working version of the `MouseModels` module. For all other parsers use master.

### System requirements and running time
The PhenoDigm parser has been run both in a MacBook Pro and a Google Cloud machine. The current version is very memory greedy, using up to 60 GB. 

* MacBook Pro: It takes around 8 hours to run in a 16 GB laptop.
* Google Cloud Machine: It runs in less than 2 hours in a 60 GB machine (_n1-standard-16_). _ag-phenodigm-1904_ in _open-targets-eu-dev_ is set up to run it (see below).

### Installation
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

### Running MouseModels in _ag-phenodigm-1904_ Google Cloud machine
Start the machine and connect to it via SSH.

Before you start you may want to change the name of the output JSON by editing `MOUSEMODELS_EVIDENCE_FILENAME` variable in [settings.py](https://github.com/opentargets/evidence_datasource_parsers/blob/solr_phenodigm_1904/settings.py):
```python3
MOUSEMODELS_EVIDENCE_FILENAME = 'phenodigm-4-4-2019.json'
```
When you are ready follow these steps:
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
```
**NOTE:** Remember to stop the machine once you are done as it costs money to have it on! 
