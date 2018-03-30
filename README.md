# OT evidence generator
Install (requires python 3):

```sh
virtualenv -p python3 venv
source venv/bin/activate
pip install -r requirements.txt
```

Run tests:
```sh
(venv)$ python -m pytest
```

Each script is a standalone python script.
Common dependencies are stored in the `common` folder.

Hence to run each parser, simply run the standalone script with your python
interpreter:
```sh
(venv)$ python modules/<parser you want>.py
```
