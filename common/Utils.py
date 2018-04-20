import logging
import requests
import tqdm

class TqdmLoggingHandler (logging.Handler):
    def __init__ (self, level = logging.NOTSET):
        super (self.__class__, self).__init__ (level)

    def emit (self, record):
        try:
            msg = self.format (record)
            tqdm.tqdm.write (msg)
            self.flush ()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

def mapping_on_github(modulename):
    '''
    send a HEAD request to github's mapping repo to see if we have a mapping
    file
    '''
    base = 'https://raw.githubusercontent.com/opentargets/mappings/master/{}.mappings.tsv'
    r = requests.head(base.format(modulename))
    return r.status_code == 200