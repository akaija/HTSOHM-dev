#!/usr/bin/env python3

import os

import click
from sqlalchemy import or_

import htsohm
from htsohm.db import session, Material
from htsohm.hypothetical_MOFs.phi_utilities import load_data, rad_dist_func

@click.group()
def add_phi():
    pass

@add_phi.command()
@click.argument('run_path',type=click.Path())
def start(run_path):
    s, m = session, Material
    ids = [e[0] for e in s.query(m.id).filter(m.run_id==run_path, m.generation_index<100, or_(m.retest_passed==None, m.retest_passed==True)).all()]
    print(len(ids))

if __name__ == '__main__':
    add_phi()
