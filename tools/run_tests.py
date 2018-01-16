#! /usr/bin/env python3 -b

import click, glob, os

@click.command()
@click.option('--tool','-t',help='which tool to use',default="yview")
def run(tool='yview'):
    if tool == "yview":
        for filename in glob.glob('tests/*/*.obj'):
            cmd = './bin/yview ' + filename
            print(cmd)
            os.system(cmd)

run()
