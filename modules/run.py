import click
from evidence_parsers.cli import pass_context


@click.command('bigtest', short_help='Shows bigtest changes.')
@pass_context
def main(ctx):
    """Shows bigtes changes in the current working directory."""
    ctx.log('phewas stuff: none')
    ctx.vlog('bla bla bla, debug info')