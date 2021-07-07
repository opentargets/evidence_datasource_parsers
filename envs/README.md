# Conda + PIP environment for running the parsers

This repository follows advice outlined in [“Reproducible and upgradable Conda environments”.](https://pythonspeed.com/articles/conda-dependency-management/).

## To reproduce the existing environment
This recreates the exact environment which was present during the time of a particular commit in the repository. Make sure to substitute the appropriate lock file for your platform (Linux/OSX).
```bash
export LOCKFILE=envs/lock-conda-linux-64.txt
conda env remove -n evidence_datasource_parsers_lock
conda create --name evidence_datasource_parsers_lock --file $LOCKFILE
conda activate evidence_datasource_parsers_lock
pip install -r envs/lock-requirements.txt
```

## To update the environment
Add or update any direct dependencies to [`environment.yml`](environment.yml) (Conda) and `requirements.txt` (PIP), as needed.
* Do not add any transitive dependencies.
* Do not forget to remove any dependencies which are not required anymore.
* When a package is available through both Conda and PIP, it's recommended to use the Conda version.

Recreate and activate the environment:
```bash
conda env remove -n evidence_datasource_parsers
conda env create -f envs/environment.yml
conda activate evidence_datasource_parsers
pip install -r envs/requirements.txt
```
* Resolve the version conflicts if any arise.
* Make sure that the new environment works not only for the functionality which you are adding, but for the unchanged modules as well.
* Re-run the test suite.
* Ideally, run some or all of the modules using the old and the new environment and compare the results.
* Never install any packages directly to this environment without also adding them to `environment.yml` or `requirements.txt`.

Once the testing is complete, create a suite of lock files:
```bash
conda-lock lock -p linux-64 -p osx-64 -f envs/environment.yml --filename-template 'envs/lock-conda-{platform}.txt'
pip freeze | grep -v 'file://' > envs/lock-requirements.txt
```
* These can be used to exactly reproduce the environment you had at the time of creating them (see the section above).
* Two lock files are generated for Linux and for  to ensure cross-platform compatibility.

Commit the changes to `environment.yml`, `requirements.txt`, and all lock files alongside with your pull request.
