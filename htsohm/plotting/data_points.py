import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from sqlalchemy import or_

import htsohm
from htsohm.db.__init__ import session, Material
from htsohm.files import load_config_file

def query_data_points(run_id, generations):
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    config_path = os.path.join(htsohm_dir, run_id, 'config.yaml')
    config = load_config_file(config_path)
    properties = config['material_properties']

    filters = [Material.run_id == run_id, or_(Material.retest_passed == True,
        Material.retest_passed == None),
        Material.generation_index < config['children_per_generation']]

    data = {}
    for generation in generations:
        data['gen_{}'.format(generation)] = {}
        d = data['gen_{}'.format(generation)]

        d['current'] = {}
        c = d['current']
        d['previous'] = {}
        p = d['previous']
      
        if 'gas_adsorption_0' in properties:
            if 'gas_adsorption_1' not in properties:
                ga_query = Material.ga0_absolute_volumetric_loading
            elif 'gas_adsorption_1' in properties:
                ga_query = (Material.ga0_absolute_volumetric_loading -
                        Material.ga1_absolute_volumetric_loading)
            else:
                print('Unexpected gas adsorption property.')

            ga_c_tuples = session.query(ga_query) \
                    .filter(*filters, Material.generation == generation).all()
            c['ga'] = [e[0] for e in ga_c_tuples]

            ga_p_tuples = session.query(ga_query) \
                    .filter(*filters, Material.generation < generation).all()
            p['ga'] = [e[0] for e in ga_p_tuples]

        if 'surface_area' in properties:
            sa_c_tuples = session.query(Material.sa_volumetric_surface_area) \
                    .filter(*filters, Material.generation == generation).all()
            c['sa'] = [e[0] for e in sa_c_tuples]

            sa_p_tuples = session.query(Material.sa_volumetric_surface_area) \
                    .filter(*filters, Material.generation < generation).all()
            p['sa'] = [e[0] for e in sa_p_tuples]

        if 'helium_void_fraction' in properties:
            vf_c_tuples = session.query(Material.vf_helium_void_fraction) \
                    .filter(*filters, Material.generation == generation).all()
            c['vf'] = [e[0] for e in vf_c_tuples]

            vf_p_tuples = session.query(Material.vf_helium_void_fraction) \
                    .filter(*filters, Material.generation < generation).all()
            p['vf'] = [e[0] for e in vf_p_tuples]

    return data

def get_limits(config, x):
    if x == 'ga':
        return config['gas_adsorption_0']['limits']
    elif x == 'sa':
        return config['surface_area']['limits']
    elif x == 'vf':
        return config['helium_void_fraction']['limits']
    else:
        print('Unexpected flag.')

def get_label(x):
    if x == 'ga':
        return 'Volumetric gas adsorption ({0}/{0})'.format(r'$cm^3$')
    elif x == 'sa':
        return 'Volumetric surface area ({}/{})'.format(r'$m^2$', r'$cm^3$')
    elif x == 'vf':
        return 'Helium void fraction ({})'.format(r'$dim.$')
    else:
        print('Unexpected flag.')

def plot_2D(run_id, generations):
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

    data = query_data_points(run_id, generations)

    fig = plt.figure(figsize=(4 * len(generations), 4 * len(data_combinations)))
    G = gridspec.GridSpec(4 * len(data_combinations), 4 * len(generations))

    col_counter = 0
    for generation in generations:
        d = data['gen_{}'.format(generation)]
        c, p = d['current'], d['previous']

        row_counter = 0
        for combo in data_combinations:
            ax = plt.subplot2grid(
                    (4 * len(data_combinations), 4 * len(generations)),
                    (row_counter, col_counter), rowspan = 4, colspan = 4)
            plt.scatter(p[combo[0]], p[combo[1]], edgecolor = 'none',
                    facecolor = 'k', alpha = 0.7, s = 5)
            plt.scatter(c[combo[0]], c[combo[1]], edgecolor = 'none',
                    facecolor = 'r', alpha = 0.7, s = 5)

            [x_min, x_max] = [*get_limits(config, combo[0])]
            [y_min, y_max] = [*get_limits(config, combo[1])]
            plt.xlim(x_min, x_max)
            plt.ylim(y_min, y_max)

            if col_counter == 0:
                plt.xticks(rotation=45) 
                plt.xlabel(get_label(combo[0]))
                plt.ylabel(get_label(combo[1]))
            else:
                plt.tick_params(labelbottom = 'off', labelleft = 'off')

            if row_counter == 0:
                ax.text(0.1 * x_max, 0.85 * y_max, 'Gen. {}'.format(generation))

            ax.set_xticks(np.arange(x_min, x_max,
                (x_max - x_min) / config['number_of_convergence_bins']))
            ax.set_yticks(np.arange(y_min, y_max,
                (y_max - y_min) / config['number_of_convergence_bins']))
            plt.grid()

            row_counter += 4
        col_counter += 4

    plt.tight_layout()
    plt.savefig('{}_MaxGen{}_DataPoints_2D.png'.format(run_id, max(generations)))
