# phewascatalog_parser
script to convert phewascatalog database to JSONs, with ENSGIDs and EFO ids.

* Download the full phewas catalog using url - https://phewascatalog.org/phewas
* Clone this repo :  `git clone https://github.com/opentargets/phewascatalog_parser.git`

```sh
cd phewascatalog_parser
```

* Copy phewas-catalog.csv downloaded in first step to phewascatalog_parser/resources
* Create a virtual environment and activate it
```sh
virtualenv venv
source venv/bin/activate
```

* Install packages
```sh
pip install -r requirements.txt
```

* Run the phewas parser - with default schema-version - 1.2.6

```sh
python setup.py install
phewascatalog_parser -h
phewascatalog_parser --phewas
```

* Run the phewas parser - with updated schema-version
```sh
phewascatalog_parser --phewas --schema-version=1.2.7
```

* Output file - phewascatalog.json
