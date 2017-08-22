from setuptools import setup

setup(
    name='phewascatalog_parser',
    version='1.0',
    packages=['','common', 'modules'],
    url='',
    license='',
    author='priyankaw',
    author_email='',
    description='',
    entry_points={
          'console_scripts': ['phewascatalog_parser=CommandLine:main'],
      },
    data_files=['resources/efo.obo','resources/hp.obo','resources/phecode_icd9_rolled.csv','resources/phewas-catalog.csv']

)
