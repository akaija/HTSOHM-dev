#!/usr/bin/env python

import subprocess
from random import choice

import numpy as np
from sqlalchemy import text
from scipy.spatial.distance import euclidean

from man003.processing import post_processed_engine as engine
from man003.processing import run_id

class Limits:
    def __init__(self, column, min_, max_):
        self.column = column
        self.min = min_
        self.max = max_

class Quadrant:
    def __init__(self, tag, x_args, y_args):
        self.tag = tag
        self.x = Limits(*x_args)
        self.y = Limits(*y_args)

x = "co2_selectivity"
y = "co2_working_capacity"

quadrants = [Quadrant("A",   [x, 0,   130],  [y, 0,  10]),
             Quadrant("B",   [x, 260, 390],  [y, 10, 20]),
             Quadrant("C",   [x, 260, 390],  [y, 50, 60]),
             Quadrant("D",   [x, 260, 390],  [y, 90, 100]),
             Quadrant("E",   [x, 520, 650],  [y, 80, 90])]

def add_where_clause(q, query_string):
    query_string += """where """
    for a in ["x", "y"]:
        axis = getattr(q, a)
        query_string += """{0}>{1} and {0}<{2} and """.format(axis.column, axis.min, axis.max)
    query_string += """run_id='{}'""".format(run_id)
    return query_string

for q in quadrants:
    properties = ["average_vdw_radius", "average_well_depth", "number_density",
                  "helium_void_fraction", "volumetric_surface_area"]
    query_string = """select avg({}), avg({}), avg({}), avg({}), avg({}) from processed_data """.format(*properties)
    query_string = add_where_clause(q, query_string)
    rows = engine.connect().execute(text(query_string))
    for row in rows:
        center = np.array([*row])

    query_string = """select seed, {}, {}, {}, {}, {} from processed_data """.format(*properties)
    query_string = add_where_clause(q, query_string)

    rows = engine.connect().execute(text(query_string))
    seeds = []
    closest_distance = 10 ** 6
    for row in rows:
        point = np.array([*row[1:]])
        distance = euclidean(center, point)
        if distance < closest_distance:
            closest_distance = distance
            closest_seed = row[0]

    for color_by in ["epsilon", "charge"]:
        args = "~/builds/blender/blender-2.79/blender-2.79b-linux-glibc219-x86_64/blender empty_scene.blend --background --python render_wrapper.py {} {} {} side".format(
                closest_seed, color_by, q.tag)
        subprocess.call(args, shell=True)
