from envparse import env, ConfigurationError
import logging


logger = logging.getLogger(__name__)



def read_option(option, cast=None,
                **kwargs):

    try:
        default_value = kwargs.pop('default')
    except KeyError:
        default_value = None

    try:
        # reading the environment variable with envparse
        return env(option, cast=cast, **kwargs)
    except ConfigurationError:
       return default_value



class Config():
    MONGO_URL = read_option('MONGO_URL', cast=str,
                                      default='')
    MONGO_DB = read_option('MONGO_DB', cast=str,
                            default='')
    MONGO_TABLE = read_option('MONGO_TABLE', cast=str,
                            default='')
    MONGO_USER = read_option('MONGO_USER', cast=str,
                            default='')
    MONGO_PWD = read_option('MONGO_PWD', cast=str,
                             default='')

