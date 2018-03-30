import click
from evidence_parsers.cli import pass_context

@click.command('phewas', short_help='phewascatalog.org => JSON evidence strings')
@pass_context
def main(ctx):
    """Transform phewascatalog.org data into JSON evidence strings"""
    ctx.log('Start processing phewascatalog.org data...')
    ctx.vlog('debug info')