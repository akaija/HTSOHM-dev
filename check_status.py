#!/usr/bin/env python3

import click
from sqlalchemy import func

from htsohm.db import Material, session

@click.group()
def data():
    pass

@data.command()
def library_sizes():
    result = session.query(Material.run_id, func.count(Material.id)).group_by(Material.run_id)
    print("\nrun-id\t\t\t\t|\tlibrary size")
    print("{}+{}".format("-" * 32, "-" * 32))
    for row in result:
        print("{}\t|\t{}".format(row[0], row[1]))

if __name__ == '__main__':
    data()
