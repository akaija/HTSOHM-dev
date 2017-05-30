import os

import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

import htsohm
from htsohm.db.__init__ import session, MutationStrength
from htsohm.files import load_config_file

def get_prior(config, run_id, generation,
        gas_adsorption_bin, surface_area_bin, void_fraction_bin):
    ms = session.query(MutationStrength) \
                .filter(
                    MutationStrength.run_id == run_id,
                    MutationStrength.gas_adsorption_bin == gas_adsorption_bin,
                    MutationStrength.surface_area_bin == surface_area_bin,
                    MutationStrength.void_fraction_bin == void_fraction_bin,
                    MutationStrength.generation <= generation) \
                .order_by(MutationStrength.generation.desc()) \
                .first()

    if ms:
        return ms.strength
    else:
        return config['initial_mutation_strength'] 

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
                plt.xlim(0, config['number_of_convergence_bins'])
                plt.ylim(0, config['number_of_convergence_bins'])

                for x in range(config['number_of_convergence_bins']):
                    for y in range(config['number_of_convergence_bins']):
                        ms_key = [run_id, generation]
                        if [x_axis, y_axis] == ['vf', 'sa']:
                            ms_key = ms_key + [x, y, z]
                        elif [x_axis, y_axis] == ['vf', 'ga']:
                            ms_key = ms_key + [x, z, y]
                        elif [x_axis, y_axis] == ['sa', 'ga']:
                            ms_key = ms_key + [z, x, y]
                        s = get_prior(config, *ms_key)
                        ax.add_patch(patches.Rectangle((x, y), 1, 1,
                            facecolor=cm.Reds(s), edgecolor=None))
                col_counter += 4
            row_counter += 4
        plt.tight_layout()
        plt.savefig('{}_Gen{}_MutationStrengths_3D.png'.format(run_id, generation))
        plt.close()
