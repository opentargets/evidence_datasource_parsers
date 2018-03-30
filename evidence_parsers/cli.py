import os
import sys
from glob import glob
import click

CONTEXT_SETTINGS = dict(auto_envvar_prefix='EVPARSER')


class Context(object):

    def __init__(self):
        self.verbose = False
        self.home = os.getcwd()

    def log(self, msg, *args):
        """Logs a message to stderr."""
        if args:
            msg %= args
        click.echo(msg, file=sys.stderr)

    def vlog(self, msg, *args):
        """Logs a message to stderr only if verbose is enabled."""
        if self.verbose:
            self.log(msg, *args)


pass_context = click.make_pass_decorator(Context, ensure=True)
cmd_folders = [d for d in glob(os.path.abspath(os.path.dirname(__file__)) + '/*/')
            if 'common' not in d]

class EvParserCLI(click.MultiCommand):

    def list_commands(self, ctx):
        rv = []
        for directory in cmd_folders:
            rv.append(os.path.basename(os.path.normpath(directory)))
        rv.sort()
        print(rv)
        return rv

    def get_command(self, ctx, name):
        try:
            if sys.version_info[0] == 2:
                name = name.encode('ascii', 'replace')
            mod = __import__('evidence_parsers.{}.run'.format(name),
                             None, None, ['main'])
        except ImportError:
            return
        return mod.main

@click.command(cls=EvParserCLI, context_settings=CONTEXT_SETTINGS)
@click.option('--output', type=click.Path(exists=True, file_okay=False,
                                        resolve_path=True),
              help='Changes the output folder')
@click.option('-v', '--verbose', is_flag=True,
              help='Enables verbose mode.')
@pass_context
def evparser(ctx, verbose, output):
    """A complex command line interface."""
    ctx.verbose = verbose
    if output is not None:
        ctx.output = output