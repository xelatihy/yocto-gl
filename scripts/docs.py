#! /usr/bin/env python3 -B

import click, glob, os

@click.group()
def cli():
    pass

@cli.command()
def serve():
    os.system('mkdocs serve')

@cli.command()
def deploy():
    os.system('mkdocs gh-deploy')

cli()
