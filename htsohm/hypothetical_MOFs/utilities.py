import os
from random import uniform
from collections import Counter
from math import pi, cos, sqrt

import numpy as np
import yaml

from htsohm.db import Material
from htsohm.pseudo_material import PseudoMaterial

def load_pseudomaterial(path):
    with open(path) as material_file:
        structure = yaml.load(material_file)
        
    material = Material(structure.uuid)
    material.generation = 0
    
    return material, structure

def normalise_fraction(value, min0, max0, min1, max1):
    return (((value - min0) *(max1 - min1)) / (max0 - min0)) + min1

def shift_frac(frac):
    if frac < -2:
        frac += 2
    elif frac < -1:
        frac += 1
    elif frac > 1:
        frac -= 1
    elif frac > 2:
        frac -= 2
    return frac

def calculate_average_sigma_epsilon(atom_sites, atom_types):
    atom_ids = []
    for atom_site in atom_sites:
        atom_ids.append(atom_site['chemical-id'])
    counts = Counter(atom_ids)

    sigma_products, epsilon_products = [], []
    for i in counts:
        count = counts[i]
        for atom_type in atom_types:
            if atom_type['chemical-id'] == i:
                sigma = atom_type['sigma']
                epsilon = atom_type['epsilon']
        sigma_product = count * sigma
        epsilon_product = count * epsilon
        sigma_products.append(sigma_product)
        epsilon_products.append(epsilon_product)
        
    return sum(sigma_products) / float(len(atom_sites)), sum(epsilon_products) / float(len(atom_sites))

def strip_numbers_and_symbols(string):
    valid_letters = 'QqWwEeRrTtYyUuIiOoPpAaSsDdFfGgHhJjKkLlZzXxCcVvBbNnMm'
    new_string = ''
    for char in string:
        if char in valid_letters:
            new_string += char
    return new_string

def load_MOF(path, config):
    for file in os.listdir(path):
        if file.endswith(".cif"):
            cif_path = os.path.join(path, file)
            name = file[:-4]
    try:
        print('.cif-file path : {}'.format(cif_path))
    except:
        print('.cif-file not found.')
    
    material = Material(name)
    material.generation = 0
    
    mixing_rules_path = os.path.join(path, 'force_field_mixing_rules.def')
    force_field_path = os.path.join(path, 'force_field.def')
    pseudo_atoms_path = os.path.join(path, 'pseudo_atoms.def')
    
    if 'hypothetical' in cif_path:
        lattice_constants = np.genfromtxt(cif_path, usecols=1, skip_header=8, max_rows=3)
        lattice_angles = np.genfromtxt(cif_path, usecols=1, skip_header=11, max_rows=3)
        cif_atom_types = np.genfromtxt(cif_path, usecols=0, skip_header=20, dtype=str)
        xfrac = np.genfromtxt(cif_path, usecols=2, skip_header=20)
        yfrac = np.genfromtxt(cif_path, usecols=3, skip_header=20)
        zfrac = np.genfromtxt(cif_path, usecols=4, skip_header=20)
    else:
        lattice_constants = np.genfromtxt(cif_path, usecols=1, skip_header=9, max_rows=3)
        lattice_angles = np.genfromtxt(cif_path, usecols=1, skip_header=12, max_rows=3)
        cif_atom_types = np.genfromtxt(cif_path, usecols=0, skip_header=24, dtype=str)
        xfrac = np.genfromtxt(cif_path, usecols=2, skip_header=24)
        yfrac = np.genfromtxt(cif_path, usecols=3, skip_header=24)
        zfrac = np.genfromtxt(cif_path, usecols=4, skip_header=24)
    
    mixing_atom_types = np.genfromtxt(mixing_rules_path, usecols=0, skip_header=7, skip_footer=4, dtype=str)
    sigma = np.genfromtxt(mixing_rules_path, usecols=3, skip_header=7, skip_footer=4)
    epsilon = np.genfromtxt(mixing_rules_path, usecols=2, skip_header=7, skip_footer=4)

    # shift *frac values that exceed abs(*frac) > 1...
    for i in range(len(xfrac)):
        xfrac[i] = shift_frac(xfrac[i])
        yfrac[i] = shift_frac(yfrac[i])
        zfrac[i] = shift_frac(zfrac[i])
    
    [lc_min, lc_max] = [*config['lattice_constant_limits']]
    if (lattice_constants[0] < lc_min or lattice_constants[0] > lc_max or
        lattice_constants[1] < lc_min or lattice_constants[1] > lc_max or
        lattice_constants[2] < lc_min or lattice_constants[2] > lc_max):
        print('\nWARNING : LATTICE CONSTANT OUT OF BOUNDS.\n')
        # double unit cell, if necessary
        if (lattice_constants[0] < lc_min and 2 * lattice_constants[0] <= lc_max and
            lattice_constants[1] < lc_min and 2 * lattice_constants[1] <= lc_max and
            lattice_constants[2] < lc_min and 2 * lattice_constants[2] <= lc_max):
            print('\nDOUBLING UNIT CELL.\n')
            a0 = lattice_constants[0]
            a1 = lattice_constants[1]
            a2 = lattice_constants[2]
            lattice_constants[0] += a0
            lattice_constants[1] += a1
            lattice_constants[2] += a2
        
            x0, y0, z0, at0 = xfrac, yfrac, zfrac, cif_atom_types
            xfrac, yfrac, zfrac, cif_atom_types = [], [], [], []
            for i in range(len(x0)):
                # 1st section
                cif_atom_types.append(at0[i])
                xfrac.append(x0[i] / 2 + 0.25)
                yfrac.append(y0[i] / 2 + 0.25)
                zfrac.append(z0[i] / 2 + 0.25)
    
                # 2nd section
                cif_atom_types.append(at0[i])
                xfrac.append(x0[i] / 2 - 0.25)
                yfrac.append(y0[i] / 2 + 0.25)
                zfrac.append(z0[i] / 2 + 0.25)
               
                # 3rd section
                cif_atom_types.append(at0[i])
                xfrac.append(x0[i] / 2 - 0.25)
                yfrac.append(y0[i] / 2 - 0.25)
                zfrac.append(z0[i] / 2 + 0.25)
    
                # 4th section
                cif_atom_types.append(at0[i])
                xfrac.append(x0[i] / 2 - 0.25)
                yfrac.append(y0[i] / 2 - 0.25)
                zfrac.append(z0[i] / 2 - 0.25)
    
                # 5th section
                cif_atom_types.append(at0[i])
                xfrac.append(x0[i] / 2 + 0.25)
                yfrac.append(y0[i] / 2 - 0.25)
                zfrac.append(z0[i] / 2 + 0.25)
    
                # 6th section
                cif_atom_types.append(at0[i])
                xfrac.append(x0[i] / 2 + 0.25)
                yfrac.append(y0[i] / 2 - 0.25)
                zfrac.append(z0[i] / 2 - 0.25)
               
                # 7th section
                cif_atom_types.append(at0[i])
                xfrac.append(x0[i] / 2 + 0.25)
                yfrac.append(y0[i] / 2 + 0.25)
                zfrac.append(z0[i] / 2 - 0.25)
                
                # 8th section
                cif_atom_types.append(at0[i])
                xfrac.append(x0[i] / 2 - 0.25)
                yfrac.append(y0[i] / 2 + 0.25)
                zfrac.append(z0[i] / 2 - 0.25)

    # remove numbers and symbols from atom-type names...
    cif_atom_types = [strip_numbers_and_symbols(e) for e in cif_atom_types]
    mixing_atom_types = [strip_numbers_and_symbols(e) for e in mixing_atom_types]

    structure = PseudoMaterial(material.uuid)
    structure.run_id = name
    
    new_atom_types = []
    structure.atom_types = []
    for i in range(len(mixing_atom_types)):
        new_atom_type = 'A_{}'.format(i)
        new_atom_types.append(new_atom_type)
        structure.atom_types.append({
            'chemical-id' : new_atom_type,
            'charge'      : 0.,
            'epsilon'     : epsilon[i],
            'sigma'       : sigma[i]
        })
    
    # update cif atom-type names...
    new_cif_atom_types = []
    for i in cif_atom_types:
       new_cif_atom_types.append(new_atom_types[mixing_atom_types.index(i)])

    structure.lattice_constants = {}
    structure.lattice_constants['a'] = lattice_constants[0]
    structure.lattice_constants['b'] = lattice_constants[1]
    structure.lattice_constants['c'] = lattice_constants[2]
    
    structure.lattice_angles = {}
    structure.lattice_angles['alpha'] = lattice_angles[0]
    structure.lattice_angles['beta'] = lattice_angles[1]
    structure.lattice_angles['gamma'] = lattice_angles[2]
 
    atom_sites = []
    atom_types_used = []
    for i in range(len(new_cif_atom_types)):
        if new_cif_atom_types[i] not in atom_types_used:
            atom_types_used.append(new_cif_atom_types[i])
        atom_sites.append({
            'chemical-id' : new_cif_atom_types[i],
            'x-frac'      : xfrac[i],
            'y-frac'      : yfrac[i],
            'z-frac'      : zfrac[i]
        })
    [a, b, c] = [*lattice_constants]
    [A, B, C] = [e * pi / 180. for e in lattice_angles]
    vol = a * b * c * sqrt( 1 + 2 * cos(A) * cos(B) * cos(C) - cos(A) ** 2 - cos(B) ** 2 - cos(C) ** 2)
    material.unit_cell_volume = vol
    nden = len(xfrac) / vol
    material.number_density = nden
    print('\nNUMBER DENSITY, original : {}'.format(nden))
    [nden_min, nden_max] = [*config['number_density_limits']]
    if nden > nden_max or nden < nden_min:
        new_nden = normalise_fraction(nden_min, 0.1852, nden_min, nden_max)
        new_n = new_nden * vol
        print('NUMBER DENSITY, new : {}\n'.format(nden_new))
        if new_n < len(atom_sites):
            selected_atom_sites = np.random.choice(atom_sites, new_n, replace=False)
            structure.atom_sites = selected_atom_sites
        else:
            structure.atom_sites = atom_sites
    else:
        print('NUMBER DENSITY NOT MODIFIED.')
        structure.atom_sites = atom_sites
#    structure.atom_sites = atom_sites
  
    # calculate average sigma, epsilon
    avg_sig, avg_ep = calculate_average_sigma_epsilon(structure.atom_sites, structure.atom_types)
    material.average_sigma = avg_sig
    material.average_epsilon = avg_ep

    return material, structure

import matplotlib.pyplot as plt
from sqlalchemy import create_engine, or_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

import htsohm
from htsohm.db import session, Material
from non_pseudo.db import Material as np_m

np_cs = "postgresql://htsohm:dentate-canst-freely@127.0.01/non_pseudo3"
np_engine = create_engine(np_cs)
np_s = sessionmaker(bind=np_engine)()

def make_query(querying, run_id, generation):
    return [e[0] for e in session.query(getattr(Material, querying)) \
            .filter(Material.run_id==run_id, Material.generation==generation).order_by(Material.id).all()]

def load_data(run_id, generation, ordering):

    ga_tuples = session.query(Material.ga0_absolute_volumetric_loading) \
        .filter(Material.run_id==run_id, Material.generation==generation).all()
    ga = [e[0] for e in ga_tuples]
    sa_tuples = session.query(Material.sa_volumetric_surface_area) \
        .filter(Material.run_id==run_id, Material.generation==generation).all()
    sa = [e[0] for e in sa_tuples]
    vf_tuples = session.query(Material.vf_helium_void_fraction) \
        .filter(Material.run_id==run_id, Material.generation==generation).all()
    vf = [e[0] for e in vf_tuples]
    return ga, sa, vf

def get_index(x_ax):
    if x_ax == 'ga0_absolute_volumetric_loading':
        return 'gas_adsorption_0'
    elif x_ax == 'sa_volumetric_surface_area':
        return 'surface_area'
    elif x_ax == 'vf_helium_void_fraction':
        return 'helium_void_fraction'

def plot_parent(x_ax, y_ax, run_id, config):
    r35 = '2017-08-01T14:13:50.410006'
    retest_filter = [or_(Material.retest_passed==True, Material.retest_passed==None)]
    filters = [*retest_filter, Material.generation<=3000, Material.run_id==r35,
              Material.generation_index<config['children_per_generation']]
    x = [e[0] for e in session.query(getattr(Material, x_ax)).filter(*filters).all()]
    y = [e[0] for e in session.query(getattr(Material, y_ax)).filter(*filters).all()]
    
    fig = plt.figure(figsize=(12,12))
    ax = plt.subplot(111)
    plt.scatter(x, y, alpha=0.2, edgecolor='none', facecolor='0.5', s=5)
    [xmin, xmax] = [*config[get_index(x_ax)]['limits']]
    plt.xlim(xmin, xmax)
    x_ticks = np.arange(xmin, xmax, (xmax - xmin) / config['number_of_convergence_bins'])
    ax.set_xticks(x_ticks)
    [ymin, ymax] = [*config[get_index(y_ax)]['limits']]
    plt.ylim(ymin, ymax)
    y_ticks = np.arange(ymin, ymax, (ymax - ymin) / config['number_of_convergence_bins']) 
    ax.set_yticks(y_ticks)
    plt.grid()
    plt.xlabel(x_ax)
    plt.ylabel(y_ax)
    plt.xticks(rotation=90)
#
    return fig, ax

def plot_parent_and_children(run_id, config):
    r35 = '2017-08-01T14:13:50.410006'
    retest_filter = [or_(Material.retest_passed==True, Material.retest_passed==None)]
    filters = [*retest_filter, Material.generation<=3000, Material.run_id==r35,
              Material.generation_index<config['children_per_generation']]
    ga_tuples = session \
        .query(Material.ga0_absolute_volumetric_loading) \
        .filter(*filters).all()
    ga = [e[0] for e in ga_tuples]
    sa_tuples = session.query(Material.sa_volumetric_surface_area).filter(*filters).all()
    sa = [e[0] for e in sa_tuples]
    vf_tuples = session.query(Material.vf_helium_void_fraction).filter(*filters).all()
    vf = [e[0] for e in vf_tuples]
    
#    ga0, sa0, vf0 = load_data(run_id, 0)
#    ga1, sa1, vf1 = load_data(run_id, 1)
    
    ga_o = np_s.query(np_m.ga1_absolute_volumetric_loading).filter(np_m.name==run_id).one()[0]
    sa_o = np_s.query(np_m.sa_volumetric_surface_area).filter(np_m.name==run_id).one()[0]
    vf_o = np_s.query(np_m.vf_helium_void_fraction).filter(np_m.name==run_id).one()[0]
    
    print('Void fraction v. volumetric surface area')
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(111)
    plt.scatter(vf, sa, alpha=1., edgecolor='none', facecolor='k', s=5)
    plt.scatter(vf0, sa0, alpha=1., edgecolor='none', facecolor='b', s=20)
    plt.scatter(vf1, sa1, alpha=1., edgecolor='none', facecolor='r', s=20)
    plt.scatter(vf_o, sa_o, alpha=1., edgecolor='none', facecolor='g', s=20)
    [xmin, xmax] = [*config['helium_void_fraction']['limits']]
    plt.xlim(xmin, xmax)
    x_ticks = np.arange(xmin, xmax, (xmax - xmin) / config['number_of_convergence_bins'])
    ax.set_xticks(x_ticks)
    [ymin, ymax] = [*config['surface_area']['limits']]
    plt.ylim(ymin, ymax)
    y_ticks = np.arange(ymin, ymax, (ymax - ymin) / config['number_of_convergence_bins']) 
    ax.set_yticks(y_ticks)
    plt.grid()
#    plt.tick_params(labelbottom='off', labelleft='off')
    plt.xlabel('He void fraction')
    plt.ylabel('Vol. surface area')
    plt.xticks(rotation=90)
#    plt.show()

    print('Void fraction v. methane deliverable capacity')
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(111)
    plt.scatter(vf, ga, alpha=1., edgecolor='none', facecolor='k', s=5)
    plt.scatter(vf0, ga0, alpha=1., edgecolor='none', facecolor='b', s=20)
    plt.scatter(vf1, ga1, alpha=1., edgecolor='none', facecolor='r', s=20)
    plt.scatter(vf_o, ga_o, alpha=1., edgecolor='none', facecolor='g', s=20)
    [xmin, xmax] = [*config['helium_void_fraction']['limits']]
    plt.xlim(xmin, xmax)
    x_ticks = np.arange(xmin, xmax, (xmax - xmin) / config['number_of_convergence_bins'])
    ax.set_xticks(x_ticks)
    [ymin, ymax] = [*config['gas_adsorption_0']['limits']]
    plt.ylim(ymin, ymax)
    y_ticks = np.arange(ymin, ymax, (ymax - ymin) / config['number_of_convergence_bins']) 
    ax.set_yticks(y_ticks)
    plt.grid()
#    plt.tick_params(labelbottom='off', labelleft='off')
    plt.xlabel('He void fraction')
    plt.ylabel('CH4 v/v STP')
    plt.xticks(rotation=90)
    plt.show()
    
    print('Surface area v. methane deliverable capacity')
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(111)
    plt.scatter(sa, ga, alpha=1., edgecolor='none', facecolor='k', s=5)
    plt.scatter(sa0, ga0, alpha=1., edgecolor='none', facecolor='b', s=20)
    plt.scatter(sa1, ga1, alpha=1., edgecolor='none', facecolor='r', s=20)
    plt.scatter(sa_o, ga_o, alpha=1., edgecolor='none', facecolor='g', s=20)
    [xmin, xmax] = [*config['surface_area']['limits']]
    plt.xlim(xmin, xmax)
    x_ticks = np.arange(xmin, xmax, (xmax - xmin) / config['number_of_convergence_bins'])
    ax.set_xticks(x_ticks)
    [ymin, ymax] = [*config['gas_adsorption_0']['limits']]
    plt.ylim(ymin, ymax)
    y_ticks = np.arange(ymin, ymax, (ymax - ymin) / config['number_of_convergence_bins']) 
    ax.set_yticks(y_ticks)
    plt.grid()
#    plt.tick_params(labelbottom='off', labelleft='off')
    plt.xticks(rotation=90)
    plt.xlabel('Vol. surface area')
    plt.ylabel('CH4 v/v STP')
    plt.show()

from numpy.random import choice

from htsohm.db import session
from htsohm.files import load_config_file
from htsohm.htsohm import retest2, run_all_one_off
#from htsohm.material_files import mutate_pseudo_material2

def simulate(pm, ps, config):
    run_all_one_off(pm, ps, config)
    session.add(pm)
    session.commit()

def mutate_hMOF(parent_material, parent_pseudo_material, mutation_strength, config):
    generation = 1

    ########################################################################
    # load boundaries from config-file
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]

    ########################################################################
    # create material object
    child_material = Material(parent_material.run_id)
    child_material.parent_id = parent_material.id
    child_material.generation = generation
    child_pseudo_material = PseudoMaterial(child_material.uuid)

    child_pseudo_material.run_id = parent_pseudo_material.run_id

    ########################################################################
    # perturb LJ-parameters
    child_pseudo_material.atom_types = parent_pseudo_material.atom_types.copy()
    for atom_type in child_pseudo_material.atom_types:    
        for x in ['epsilon', 'sigma']:
            old_x = atom_type[x]
            random_x = uniform(*config["{0}_limits".format(x)])
            atom_type[x] += mutation_strength * (random_x - old_x)

    ########################################################################
    # calculate new lattice constants
    child_pseudo_material.lattice_constants = parent_pseudo_material.\
            lattice_constants.copy()
    for i in ['a', 'b', 'c']:
        old_x = parent_pseudo_material.lattice_constants[i]
        random_x = uniform(*lattice_limits)
        child_pseudo_material.lattice_constants[i] += mutation_strength * (
                random_x - old_x)

    child_pseudo_material.lattice_angles = parent_pseudo_material.lattice_angles.copy()

    ########################################################################
    #perturb number density, calculate number of atoms
    child_ND = parent_pseudo_material.number_density()
    random_ND = uniform(*number_density_limits)
    child_ND += mutation_strength * (random_ND - child_ND)

    print('\nNUMBER DENSITY : {}\n'.format(child_ND))
    child_material.number_density = child_ND

    child_number_of_atoms = (
            int(child_ND * child_pseudo_material.volume()))

    ########################################################################
    # remove excess atom-sites, if any
    child_pseudo_material.atom_sites = np.random.choice(
            parent_pseudo_material.atom_sites,
            min(child_number_of_atoms,
                len(parent_pseudo_material.atom_sites)),
            replace=False).tolist()
    
    ########################################################################
    # perturb atom-site positions
    for atom_site in child_pseudo_material.atom_sites:
        for i in ['x-frac', 'y-frac', 'z-frac']:
            random_frac = uniform(-1, 1)
            atom_site[i] += mutation_strength * (random_frac - atom_site[i])

    ########################################################################
    # add atom-sites, if needed
    if child_number_of_atoms > len(child_pseudo_material.atom_sites):
        for new_sites in range(child_number_of_atoms -
                len(child_pseudo_material.atom_sites)):
            new_atom_site = {'chemical-id' : choice(
                child_pseudo_material.atom_types)['chemical-id']}
            for i in ['x-frac', 'y-frac', 'z-frac']:
                new_atom_site[i] = uniform(-1, 1)
            child_pseudo_material.atom_sites.append(new_atom_site)
    
    lcs = child_pseudo_material.lattice_constants
    child_material.unit_cell_volume = lcs['a'] * lcs['b'] * lcs['c']
    avg_sig, avg_ep = calculate_average_sigma_epsilon(child_pseudo_material.atom_sites, child_pseudo_material.atom_types)
    child_material.average_sigma = avg_sig
    child_material.average_epsilon = avg_ep

    return child_material, child_pseudo_material

def mutate_hMOF_sites(parent_material, parent_pseudo_material, mutation_strength, config):
    generation = 1

    ########################################################################
    # load boundaries from config-file
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]

    ########################################################################
    # create material object
    child_material = Material(parent_material.run_id)
    child_material.parent_id = parent_material.id
    child_material.generation = generation
    child_pseudo_material = PseudoMaterial(child_material.uuid)

    child_pseudo_material.run_id = parent_pseudo_material.run_id

    ########################################################################
    # perturb LJ-parameters
    child_pseudo_material.atom_types = parent_pseudo_material.atom_types.copy()

    ########################################################################
    # calculate new lattice constants
    child_pseudo_material.lattice_constants = parent_pseudo_material.\
            lattice_constants.copy()

    child_pseudo_material.lattice_angles = parent_pseudo_material.lattice_angles.copy()

    ########################################################################
    # remove excess atom-sites, if any
    child_pseudo_material.atom_sites = parent_pseudo_material.atom_sites.copy()
    
    ########################################################################
    # perturb atom-site positions
    for atom_site in child_pseudo_material.atom_sites:
        for i in ['x-frac', 'y-frac', 'z-frac']:
            random_frac = uniform(-1, 1)
            atom_site[i] += mutation_strength * (random_frac - atom_site[i])

    child_ND = parent_pseudo_material.number_density()
    child_material.number_density = child_ND

    lcs = child_pseudo_material.lattice_constants
    child_material.unit_cell_volume = lcs['a'] * lcs['b'] * lcs['c']
    avg_sig, avg_ep = calculate_average_sigma_epsilon(child_pseudo_material.atom_sites, child_pseudo_material.atom_types)
    child_material.average_sigma = avg_sig
    child_material.average_epsilon = avg_ep

    return child_material, child_pseudo_material
