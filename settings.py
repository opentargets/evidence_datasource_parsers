import os
import pkg_resources as res

def file_or_resource(fname=None):
    '''get filename and check if in getcwd then get from
    the package resources folder
    '''
    filename = os.path.expanduser(fname)

    resource_package = __name__
    resource_path = '/'.join(('resources', filename))

    if filename is not None:
        abs_filename = os.path.join(os.path.abspath(os.getcwd()), filename) \
                       if not os.path.isabs(filename) else filename

        return abs_filename if os.path.isfile(abs_filename) \
            else res.resource_filename(resource_package, resource_path)

class Config:
    PHEWAS_CATALOG_LOCN = 'https://storage.googleapis.com/phewas-catalog'
    GENES_HGNC =  'http://ftp.ebi.ac.uk/pub/databases/genenames/new/json/hgnc_complete_set.json'
    PHEWAS_CATALOG_JSON = '../phewas.json'

    INTOGEN_FILENAME = file_or_resource('intogen_opentargets.tsv')
    INTOGEN_EVIDENCE_FILENAME = '/tmp/intogen.json'

