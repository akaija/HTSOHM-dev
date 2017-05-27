import os

import matplotlib.pyplot as plt
from sqlalchemy.sql import func
from sqlalchemy import or_

import htsohm
from htsohm.db.__init__ import session, Material
from htsohm.files import load_config_file

def query_number_of_empty_bins(run_id, generation):
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    config_path = os.path.join(htsohm_dir, run_id, 'config.yaml')
    config = load_config_file(config_path)
    properties = config['material_properties']

    filters = [Material.run_id == run_id, or_(Material.retest_passed == True,
        Material.retest_passed == None),
        Material.generation_index < config['children_per_generation']]

    bin_dimensions = []
    if 'gas_adsorption_0' in properties:
        bin_dimensions.append(Material.gas_adsorption_bin)
    if 'surface_area' in properties:
        bin_dimensions.append(Material.surface_area_bin)
    if 'helium_void_fraction' in properties:
        bin_dimensions.append(Material.void_fraction_bin)

    total_number_of_bins = config['number_of_convergence_bins'] ** len(bin_dimensions)
    all_accessed_bins = session.query(*bin_dimensions) \
            .filter(*filters, Material.generation <= generation) \
            .distinct().count()

    return total_number_of_bins - all_accessed_bins

def plot_empty_bin_convergence(run_id, max_generation = None):
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    config_path = os.path.join(htsohm_dir, run_id, 'config.yaml')
    config = load_config_file(config_path)

    filters = [Material.run_id == run_id, or_(Material.retest_passed == True,
        Material.retest_passed == None),
        Material.generation_index < config['children_per_generation']]

    if max_generation == None:
        max_generation = session.query(func.max(Material.generation)) \
                .filter(*filters).one()[0]

    data = []
    for generation in range(max_generation):
        data.append(query_number_of_empty_bins(run_id, generation))

    fig = plt.figure(figsize = (12, 8))
    plt.xlabel('Generation')
    plt.ylabel('Number of empty bins')
    plt.plot(range(len(data)), data)
    plt.savefig('{}_MaxGen{}_EmptyBins.png'.format(run_id, max_generation))


