import os

import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from sqlalchemy.sql import func
from sqlalchemy import or_

import htsohm
from htsohm.db.__init__ import session, Material
from htsohm.files import load_config_file

def get_z(x ,y):
    if [x, y] == ['vf', 'sa']:
        return 'ga'
    elif [x, y] == ['vf', 'ga']:
        return 'sa'
    elif [x, y] == ['sa', 'ga']:
        return 'vf'

def query_all_counts_3D(config, run_id, generation, x_axis, y_axis, z):
    filters = [Material.run_id == run_id,
            or_(Material.retest_passed == None, Material.retest_passed == True),
            Material.generation_index < config['children_per_generation'],
            Material.generation <= generation]
    
    bin_group = [Material.gas_adsorption_bin, Material.surface_area_bin,
            Material.void_fraction_bin]

    if [x_axis, y_axis] == ['vf', 'sa']:
        bin_queries = [Material.void_fraction_bin,
                Material.surface_area_bin]
        filters.append(Material.gas_adsorption_bin == z)
    elif [x_axis, y_axis] == ['vf', 'ga']:
        bin_queries = [Material.void_fraction_bin,
                Material.gas_adsorption_bin]
        filters.append(Material.surface_area_bin == z)
    elif [x_axis, y_axis] == ['sa', 'ga']:
        bin_queries = [Material.surface_area_bin,
                Material.gas_adsorption_bin]
        filters.append(Material.void_fraction_bin == z)
    else:
        print('Unexpected structure-property flags.')

    tuples = session.query(*bin_queries, func.count(Material.uuid)) \
            .filter(*filters).group_by(*bin_group).all()
    data = [[ [e[0], e[1]], e[2] ] for e in tuples]

    count_tuples = session.query(func.count(Material.uuid)) \
            .filter(*filters).group_by(*bin_group).all()
    counts = [e[0] for e in count_tuples]

    if len(counts) == 0:
        max_counts = 0
    else:
        max_counts = max(counts)

    return data, max_counts

def get_label(x):
    if x == 'vf':
        return 'Void fraction bin'
    elif x == 'sa':
        return 'Surface area bin'
    elif x == 'ga':
        return 'Gas loading bin'
    else:
        print('Unexpected flag.')
        return None

def plot_3D(run_id, generations):
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    config_path = os.path.join(htsohm_dir, run_id, 'config.yaml')
    config = load_config_file(config_path)
    properties = config['material_properties']

    data_types = []
    if 'helium_void_fraction' in properties:
        data_types.append('vf')
    if 'surface_area' in properties:
        data_types.append('sa')
    if 'gas_adsorption_0' in properties:
        data_types.append('ga')
    data_combinations = []
    for i in range(len(data_types)):
        for j in range(i + 1, len(data_types)):
            combo = [data_types[i], data_types[j]]
            data_combinations.append(combo)

    for generation in generations:
        fig = plt.figure(figsize=(4 * len(data_combinations), 
            4 * config['number_of_convergence_bins']))
        G = gridspec.GridSpec(4 * config['number_of_convergence_bins'],
            4 * len(data_combinations))
        row_counter = 0
        print('Plot {} / {} ...'.format(generations.index(generation) + 1, len(generations)))

        for z in range(config['number_of_convergence_bins']):
            col_counter = 0
            for [x_axis, y_axis] in data_combinations:
                ax = plt.subplot2grid(
                    (4 * config['number_of_convergence_bins'],
                    4 * len(data_combinations)), (row_counter, col_counter),
                    rowspan = 4, colspan = 4)
                data, max_count = query_all_counts_3D(config, run_id, generation, x_axis,
                    y_axis, z)

                if [row_counter, col_counter] == [0, 0]:
                    ax.text(0.1 * config['number_of_convergence_bins'],
                            0.85 * config['number_of_convergence_bins'],
                            'Max count {}'.format(max_count))
                plt.xlim(0, config['number_of_convergence_bins'])
                plt.ylim(0, config['number_of_convergence_bins'])
                
                if row_counter == 0:
                    plt.xlabel(get_label(x_axis))
                    plt.ylabel(get_label(y_axis))

                if max_count != 0:
                    for d in data:
                        c = d[1] / max_count
                        [x, y] = [*d[0]]
                        ax.add_patch(patches.Rectangle((x, y), 1, 1,
                            facecolor = cm.Reds(c), edgecolor=None))
                col_counter += 4
            row_counter += 4
        plt.tight_layout()
        plt.savefig('{}_Gen{}_BinCounts_3D.png'.format(run_id, generation))
        plt.close()
