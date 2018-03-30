import click
from evidence_parsers.cli import pass_context, EvParserCLI, Context


@click.command('status', short_help='Shows file changes.')
@pass_context
def cli(ctx):
    """Shows file changes in the current working directory."""
    ctx.log('Changed files: one')
    ctx.vlog('bla bla bla, debug info')

@click.command('main', short_help='Shows file changes.')
@click.option('--anoption')
@pass_context
def main(ctx, anoption):
    """Shows file changes in the current working directory."""
    ctx.log('Changed files: two')
    ctx.vlog('bla bla bla, debug info')


if __name__ == '__main__':
    main(Context(),'test')