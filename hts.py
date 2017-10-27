#!/usr/bin/env python3
from datetime import datetime
import os
from distutils.dir_util import copy_tree

import click
import yaml
import RASPA2

import htsohm
from htsohm.files import load_config_file
from htsohm.htsohm import worker_run_loop, calc_bin, run_all_one_off, load_pseudo_material, calculate_properties
from htsohm.db import session, Material

@click.group()
def hts():
    pass

@hts.command()
@click.argument('config_path',type=click.Path())
def start(config_path):
    """Create a new run.
    
    Args:
        config_path (str): path to config-file (ex: setting/htsohm.sample.yaml)

    Prints run_id, creates run-folder with config-file.

    """
    config = load_config_file(config_path)
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    run_id = datetime.now().isoformat()
    config['run_id'] = run_id
    config['raspa2_dir'] = os.path.dirname(RASPA2.__file__)
    config['htsohm_dir'] = htsohm_dir
    
    run_dir = os.path.join(htsohm_dir, run_id)
    os.makedirs(run_dir, exist_ok=True)
    file_name = os.path.join(run_dir, 'config.yaml')
    with open(file_name, 'w') as config_file:
        yaml.dump(config, config_file, default_flow_style=False)
    print('Run created with id: %s' % run_id)

@hts.command()
@click.argument('old_run_id')
@click.argument('generation')
@click.argument('config_path',type=click.Path())
def append(old_run_id, generation, config_path):
    """Create a new run, appending to previous dataset. Method creates a new
    run-directory, copying all structure-files from the previous run into this
    new directory. Then the method creates a new table, copying all records
    from the previous run, but updating the run_id and assigned bin for each
    material based on the number of bins specified by the user.
    
    Args:
        old_run_id (str): identification string for run to which new data will
            be appended.
        generation (int): last generation to include from the old run.
        config_path (str): path to run configuration file.

    """
    # load config, update with new run_id
    config = load_config_file(config_path)
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    run_id = datetime.now().isoformat()
    config['run_id'] = run_id
    config['raspa2_dir'] = os.path.dirname(RASPA2.__file__)
    config['htsohm_dir'] = htsohm_dir
    
    # create run-directory, dump config
    run_dir = os.path.join(htsohm_dir, run_id)
    os.makedirs(run_dir, exist_ok=True)
    file_name = os.path.join(run_dir, 'config.yaml')
    with open(file_name, 'w') as config_file:
        yaml.dump(config, config_file, default_flow_style=False)
    print('Run created with id: %s' % run_id)

    print('Copying data from run : '.format(old_run_id))

    query_filter = [Material.run_id==old_run_id, Material.generation<=generation]
    id_tuples = session.query(Material.id).filter(*query_filter).all()
    ids = [e[0] for e in id_tuples]

    # re-calculate structure-property bins for all records, copy to new table
    properties = config['material_properties']

    for some_id in ids:
        old_material = session.query(Material).get(some_id)
        new_material = old_material.clone()
        new_material.id = None
        new_material.run_id = run_id

        if 'gas_adsorption_0' in properties and 'gas_adsorption_1' not in properties:
            new_material.gas_adsorption_bin = calc_bin(
                    new_material.ga0_absolute_volumetric_loading,
                    *config['gas_adsorption_0']['limits'],
                    config['number_of_convergence_bins'])

        if 'gas_adsorption_0' in properties and 'gas_adsorption_1' in properties:
            new_material.gas_adsorption_bin = calc_bin(
                    (new_material.ga0_absolute_volumetric_loading - new_material.ga1_absolute_volumetric_loading),
                    *config['gas_adsorption_1']['limits'],
                    config['number_of_convergence_bins'])

        if 'helium_void_fraction' in properties:
            new_material.void_fraction_bin = calc_bin(
                    new_material.vf_helium_void_fraction,
                    *config['helium_void_fraction']['limits'],
                    config['number_of_convergence_bins'])
        
        if 'surface_area' in properties:
            new_material.surface_area_bin = calc_bin(
                    new_material.sa_volumetric_surface_area,
                    *config['surface_area']['limits'],
                    config['number_of_convergence_bins'])
        
        session.add(new_material)
    session.commit()

    # copy structure-data files, updating with new run_id
    old_pseudo_materials_dir = os.path.join(old_run_id, 'pseudo_materials')

    uuid_tuples = session.query(Material.uuid).filter(*query_filter).all()
    uuids = [e[0] for e in uuid_tuples]
    for uuid in uuids:
        file_path = os.path.join(old_pseudo_materials_dir, '{}.yaml'.format(uuid))
        with open(file_path) as structure_file:
            pseudo_material = yaml.load(structure_file)
        pseudo_material.run_id = run_id
        pseudo_material.dump()

    print('...copy complete!')

@hts.command()
@click.argument('run_id')
def launch_worker(run_id):
    """Start process to manage run.

    Args:
        run_id (str): identification string for run.

    Runs HTSOHM-method in one process.

    """
    htsohm._init(run_id)
    worker_run_loop(run_id)

@hts.command()
@click.argument('uuid')
def one_off(uuid):
    material_id = session.query(Material.id).filter(Material.uuid==uuid).one()[0]
    old_material = session.query(Material).get(material_id)

    material = Material('one_off')
    material.uuid = uuid
    run_id = session.query(Material.run_id).filter(Material.uuid==uuid).one()[0]
    htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
    config = load_config_file(os.path.join(htsohm_dir, run_id, 'config.yaml'))

    pseudo_material = load_pseudo_material(run_id, material)
    
    print('SCREENING : {}'.format(uuid))
    run_all_one_off(material, pseudo_material, config)
    print('...done!\n')

    print('Pseudomaterial UUID :\t\t{}'.format(uuid))

    print('\n\tOld values :')
    print('CH4 v/v 35bar :\t\t{}'.format(old_material.ga0_absolute_volumetric_loading))
    print('He void fraction :\t{}'.format(old_material.vf_helium_void_fraction))
    print('Vol. surface area :\t{}\n'.format(old_material.sa_volumetric_surface_area))

    print('\n\tNew values :')
    print('CH4 v/v 35bar :\t\t{}'.format(material.ga0_absolute_volumetric_loading))
    print('He void fraction :\t{}'.format(material.vf_helium_void_fraction))
    print('Vol. surface area :\t{}\n'.format(material.sa_volumetric_surface_area))

    calculate_properties(uuid)

if __name__ == '__main__':
    hts()
