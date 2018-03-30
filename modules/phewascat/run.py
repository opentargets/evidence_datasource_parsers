import os
import click
from evidence_parsers.cli import pass_context
from settings import Config
from common.HGNCParser import GeneParser

@click.command('phewas', short_help='phewascatalog.org => JSON evidence strings')
@pass_context
def cli(ctx):
    """Transform phewascatalog.org data into JSON evidence strings"""
    ctx.log('Start processing phewascatalog.org data...')
    ctx.vlog('debug info')

def main():
    print('hello {}'.format(os.path.dirname(__file__)))

if __name__ == '__main__':
    main()