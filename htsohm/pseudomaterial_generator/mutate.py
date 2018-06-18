import math
from random import choice, random, randrange, uniform

import numpy as np
from sqlalchemy.sql import func, or_

from htsohm import config
from htsohm.db import session, Material, MutationStrength
from htsohm.db import Structure, LennardJones, AtomSites
from htsohm.simulation.run_all import run_all_simulations

def select_parent(run_id, max_generation, generation_limit):
    """Use bin-counts to preferentially select a list of 'rare' parents.

    Args:
        run_id (str): identification string for run.
        max_generation (int): latest generation to include when counting number
            of materials ub each bin.
        generation_limit (int): number of materials to query in each generation
            (as materials are added to database they are assigned an index
            within the generation to bound the number of materials in each
            generation).

    Returns:
        The material id(int) corresponding to some parent-material selected
        from database with a bias favoring materials in bins with the lowest
        counts.

    """
    # Determine which dimensions are being binned
    simulations = config['simulations']
    queries = []
    if 'gas_adsorption' in simulations or 'artificial_gas_adsorption' in simulations:
        queries.append( getattr(Material, 'gas_adsorption_bin') )
    if 'surface_area' in simulations or 'artificial_surface_area' in simulations:
        queries.append( getattr(Material, 'surface_area_bin') )
    if 'helium_void_fraction' in simulations or 'artificial_void_fraction' in simulations:
        queries.append( getattr(Material, 'void_fraction_bin') )

    # Count the number of materials in each bin
    bins_and_counts = session \
        .query(
            func.count(Material.id),
            *queries
        ) \
        .filter(
            Material.run_id == run_id,
            or_(Material.retest_passed == True, Material.retest_passed == None),
            Material.generation <= max_generation,
            Material.generation_index < generation_limit,
        ) \
        .group_by(*queries).all()[1:]

    # Print bin-counts in interactive mode
    if config['selection_mode'] == 'manual':
        print('\n  BIN\t\t|  COUNT')
        print('----------------+-----------')
        for bin_and_count in bins_and_counts:
            print('  {}\t|  {}'.format(bin_and_count[1:len(queries) + 2], bin_and_count[0]))
        print('\n')
    
        input('Press Enter to calculate normalized weights...')

    # Extract only bin coordinate from bins_and_counts query
    bins = []
    for i in bins_and_counts:
        some_bin = {}
        for j in range(len(queries)):
            some_bin[queries[j]] = i[j + 1]
        bins.append(some_bin)

    # Assign weight to each bin
    total = sum([i[0] for i in bins_and_counts])
    weights = [ total / float(i[0]) for i in bins_and_counts ]
    normalized_weights = [ weight / sum(weights) for weight in weights ]

    # Print bins and weights in interactive mode
    if config['selection_mode'] == 'manual':
        print('\n  BIN\t\t|  COUNT\t|  WEIGHT\t|  NORMALIZED-WEIGHT')
        print('----------------+---------------+---------------+---------------------')
        for i in range(len(bins_and_counts)):
            print('  {}\t|  {}\t\t|  {}\t\t|  {}'.format(
                bins_and_counts[i][1:len(queries) + 1], bins_and_counts[i][0],
                weights[i], normalized_weights[i]))
        print('\n')
    
    # Have user select a parent-bin
        ga_bin = input('Press select a gas adsorption bin :\t')
        sa_bin = input('Press select a surface area bin :\t')
        vf_bin = input('Press select a void fraction bin :\t')
    
    # Query parent-ids corresponding to selected bin
        parent_query = session \
            .query(Material.uuid) \
            .filter(
                Material.run_id == run_id,
                or_(Material.retest_passed == True, Material.retest_passed == None),
                Material.gas_adsorption_bin == ga_bin,
                Material.surface_area_bin == sa_bin,
                Material.void_fraction_bin == vf_bin,
                Material.generation <= max_generation,
                Material.generation_index < generation_limit,
            ).all()
        potential_parents = [i[0] for i in parent_query]
    
        print('Potential parents :')
        for i in potential_parents:
            print('\t{}'.format(i))
    
    # Have user select parent-pseudomaterial
        parent_uuid = input('Please select a parent UUID :\t')
        parent_id = session.query(Material.id).filter(Material.uuid==parent_uuid).one()

    # Weighted-selection for parent bin (non-interactive mode)
    else:
        parent_bin = np.random.choice(bins, p = normalized_weights)
        parent_queries = [i == parent_bin[i] for i in queries]
        parent_query = session \
            .query(Material.id) \
            .filter(
                Material.run_id == run_id,
                or_(Material.retest_passed == True, Material.retest_passed == None),
                *parent_queries,
                Material.generation <= max_generation,
                Material.generation_index < generation_limit,
            ).all()
        potential_parents = [i[0] for i in parent_query]
    # Randomly-select parent-pseudomaterial (non-interactive mode)
        parent_id = int(np.random.choice(potential_parents))

    return parent_id

def retest(m_orig, retests, tolerance):
    """Reproduce simulations  to prevent statistical errors.

    Args:
        m_orig (sqlalchemy.orm.query.Query): material to retest.
        retests (int): number of times to reproduce each simulation.
        tolerance (float): acceptance criteria as percent deviation from
            originally calculated value.

    Queries database to determine if there are remained retests to  be run.
    Updates row in database with total number of retests and results.

    """
    m = m_orig.clone()

    run_all_simulations(m)
    print('Retest {} / {}'.format(m_orig.retest_num + 1, retests))

    simulations = config['simulations']

    # requery row from database, in case someone else has changed it, and lock it
    # if the row is presently locked, this method blocks until the row lock is released
    session.refresh(m_orig, lockmode='update')

    if config['interactive_mode'] == 'on':
        x0 = m_orig.retest_gas_adsorption_0_sum
        x1 = m_orig.retest_gas_adsorption_1_sum
        y = m_orig.retest_surface_area_sum
        z = m_orig.retest_void_fraction_sum
            
    if m_orig.retest_num < retests:
        if 'gas_adsorption' in simulations or 'artificial_gas_adsorption' in simulations:
            m_orig.retest_gas_adsorption_0_sum += m.ga0_absolute_volumetric_loading
            if 'gas_adsorption' in simulations:
                if isinstance(config['simulations']['gas_adsorption']['external_pressure'], list):
                    m_orig.retest_gas_adsorption_1_sum += m.ga1_absolute_volumetric_loading
        if 'surface_area' in simulations or 'artificial_surface_area' in simulations:
            m_orig.retest_surface_area_sum += m.sa_volumetric_surface_area
        if 'helium_void_fraction' in simulations or 'artificial_void_fraction' in simulations:
            m_orig.retest_void_fraction_sum += m.vf_helium_void_fraction
        m_orig.retest_num += 1
        session.commit()        

    if config['interactive_mode'] == 'on':
        print('\n  PROPERTY\t\t|  OLD SUM\t|  NEW VALUE\t|  NEW SUM')
        print('------------------------+---------------+---------------+-----------')
        if 'gas_adsorption' in simulations:
            print('  gas-0 adsorption\t|  {}\t|  {}\t|  {}'.format(
                round(x0, 4),
                round(m.ga0_absolute_volumetric_loading, 4),
                round(m_orig.retest_gas_adsorption_0_sum, 4)))
            if isinstance(config['simulations']['gas_adsorption']['external_pressure'], list):
                print('  gas-1 adsorption\t|  {}\t|  {}\t|  {}'.format(
                    round(x1, 4),
                    round(m.ga1_absolute_volumetric_loading, 4),
                    round(m_orig.retest_gas_adsorption_1_sum, 4)))
        if 'surface_area' in simulations:
            print('  surface area\t\t|  {}\t|  {}\t|  {}'.format(
                round(y, 4),
                round(m.sa_volumetric_surface_area, 4),
                round(m_orig.retest_surface_area_sum, 4)))
        if 'helium_void_fraction' in simulations:
            print('  void fraction\t\t|  {}\t|  {}\t|  {}'.format(
                round(z, 4),
                round(m.vf_helium_void_fraction, 4),
                round(m_orig.retest_void_fraction_sum, 4)))


    if m_orig.retest_num >= retests:
        if config['interactive_mode'] != 'on':
            try:
                m_orig.retest_passed = m.calculate_retest_result(tolerance)
                print('\nRETEST_PASSED :\t%s' % m_orig.retest_passed)
                session.commit()
            except ZeroDivisionError as e:
                print('WARNING: ZeroDivisionError - material.calculate_retest_result(tolerance)')
        elif config['interactive_mode'] == 'on':
            print('\tPlease check, using the values below, that the difference\n' +
                    '\tbetween the original value and averaged retest values is\n' +
                    '\tnot greater than the tolerance times bin-width.')
            print('==============================================================================')
            print('  TOLERANCE\t|  {}'.format(config['retests']['tolerance']))
            print('==============================================================================')
            print('  NUM. TESTS\t|  {}'.format(config['retests']['number']))
            print('==============================================================================')
            print('  PROPERTY\t\t| ORIGINAL VALUE\t|  RETEST SUM\t|  BIN-WIDTH')
            print('------------------------+-----------------------+---------------+-------------')

            s = config['simulations']
            if 'gas_adsorption_0' in simulations:
                print('  gas-0 adsorption\t|  {}\t\t|  {}\t|  {}'.format(
                    round(m_orig.ga0_absolute_volumetric_loading, 4),
                    round(m_orig.retest_gas_adsorption_0_sum, 4),
                    (s['gas_adsorption']['limits'][1] - s['gas_adsorption']['limits'][0]) / config['number_of_convergence_bins']))
                if isinstance(config['simulations']['gas_adsorption']['external_pressure'], list):
                    print('  gas-1 adsorption\t|  {}\t\t|  {}\t|  {}'.format(
                        round(m_orig.ga1_absolute_volumetric_loading, 4),
                        round(m_orig.retest_gas_adsorption_1_sum, 4),
                    (s['gas_adsorption']['limits'][1] - s['gas_adsorption']['limits'][0]) / config['number_of_convergence_bins']))
            if 'surface_area' in simulations:
                print('  surface area\t\t|  {}\t\t|  {}\t|  {}'.format(
                    round(m_orig.sa_volumetric_surface_area, 4),
                    round(m_orig.retest_surface_area_sum, 4),
                    (s['surface_area']['limits'][1] - s['surface_area']['limits'][0]) / config['number_of_convergence_bins']))
            if 'helium_void_fraction' in simulations:
                print('  void fraction\t\t|  {}\t\t|  {}\t|  {}'.format(
                    round(m_orig.vf_helium_void_fraction, 4),
                    round(m_orig.retest_void_fraction_sum, 4),
                    (s['helium_void_fraction']['limits'][1] - s['helium_void_fraction']['limits'][0]) / config['number_of_convergence_bins']))

            m_orig.retest_passed = set_variable(
                    'Does {} pass the retest? (True/False) :\t'.format(m_orig.uuid),
                    'm_orig.retest_passed')
            session.commit()

def retest_loop(parent_material):
    # run retests until we've run enough
    while parent_material.retest_passed is None:
        print("running retest...")
        retest(parent_material, config['retests']['number'],
                config['retests']['tolerance'])
        session.refresh(parent_material)

def determine_mutation_strength(run_id, gen, parent_material): 
    if config['mutation_strength_method'] == 'flat':
        mutation_strength = config['initial_mutation_strength']
    else:
        mutation_strength_key = [run_id, gen] + parent_material.bin
        mutation_strength = MutationStrength \
                .get_prior(*mutation_strength_key).clone().strength
    return mutation_strength

def closest_distance(x, y):
    """Finds closest distance between two points across periodic boundaries.

    Args:
        x (float): value in range [0,1].
        y (float): value in range [0,1].

    Returns:
        D (float): closest distance measured between `x` and `y`, with
            periodic boundaries at 0 and 1. For example, the disance between
            x = 0.2 and y = 0.8 would be 0.4, and not 0.6.

    """
    a = 1 - y + x
    b = abs(y - x)
    c = 1 - x + y
    return min(a, b, c)

def random_position(x_o, x_r, strength):
    """Produces a point along the closest path between two points.

    Args:
        x_o (float): value in range [0,1].
        x_r (float): value in range [0,1].
        strength (float): refers to mutation strength. Value determines
            fractional distance to `x_r` from `x_o`.

    Returns:
        xfrac (float): a value between `x_o` and `x_r` across the closest
            distance accross periodic boundaries at 0 and 1.

    """
    dx = closest_distance(x_o, x_r)
    if (x_o > x_r
            and (x_o - x_r) > 0.5):
        xfrac = (x_o + strength * dx) % 1.
    if (x_o < x_r
            and (x_r - x_o) > 0.5):
        xfrac = (x_o - strength * dx) % 1.
    if (x_o >= x_r
            and (x_o - x_r) < 0.5):
        xfrac = x_o - strength * dx
    if (x_o < x_r
            and (x_r - x_o) < 0.5):
        xfrac = x_o + strength * dx
    return xfrac

def mutate_material(parent_material, mutation_strength, generation):
    """Modifies a "parent" pseudomaterial structure by perturbing each
    parameter by some factor, dictated by the `mutation_strength`.
    
    Args:
        parent_material : record for parent pseudomaterial in 'materials' table.
        mutation_strength (float): perturbation factor [0, 1].
        generation (int): iteration count for overall bin-mutate-simulate routine.

    Returns:

    Todo:
        * Add methods for assigning and mutating charges.

    """
    
    ########################################################################
    # load boundaries from config-file
    lattice_limits          = config["lattice_constant_limits"]
    number_density_limits   = config["number_density_limits"]

    child_material = Material(parent_material.run_id)
    child_material.generation = generation
    child_material.parent_id = parent_material.id

    print('Parent UUID :\t{}'.format(parent_material.uuid))
    print('Child UUID :\t{}'.format(child_material.uuid))

    # perturb lennard-jones parameters
    for atom_type in parent_material.structure.lennard_jones:
        child_material.structure.lennard_jones.append(LennardJones(
            chemical_id = atom_type.chemical_id,
            sigma = atom_type.sigma + mutation_strength * (uniform(*config['sigma_limits']) - atom_type.sigma),
            epsilon = atom_type.epsilon + mutation_strength * (uniform(*config['epsilon_limits']) - atom_type.epsilon)))

    if config['interactive_mode'] == 'on':
        print('=========================================================================================')
        print('LENNARD JONES PARAMETERS')
        print('=========================================================================================')
        print('  chemical-id\t|  parent sigma\t|  child sigma\t|  parent epsilon\t| child epsilon')
        print('----------------+---------------+---------------+-----------------------+----------------')
        for i in range(len(child_material.structure.lennard_jones)):
            c_chem = child_material.structure.lennard_jones[i]
            p_chem = parent_material.structure.lennard_jones[i]
            print('  {}\t\t|  {}\t|  {}\t|  {}\t\t|  {}'.format(c_chem.chemical_id,
                round(p_chem.sigma, 4), round(c_chem.sigma, 4),
                round(p_chem.epsilon, 4), round(c_chem.epsilon, 4)))

    # perturb lattice constants
    child_material.structure.lattice_constant_a = parent_material.structure.lattice_constant_a \
            + mutation_strength * (uniform(*lattice_limits) - parent_material.structure.lattice_constant_a)
    child_material.structure.lattice_constant_b = parent_material.structure.lattice_constant_b \
            + mutation_strength * (uniform(*lattice_limits) - parent_material.structure.lattice_constant_b)
    child_material.structure.lattice_constant_c = parent_material.structure.lattice_constant_c \
            + mutation_strength * (uniform(*lattice_limits) - parent_material.structure.lattice_constant_c)

    if config['interactive_mode'] == 'on':
        print('==========================================')
        print('LATTICE CONSTANTS')
        print('==========================================')
        print('  direction\t|  parent \t|  child')
        print('----------------+---------------+---------')
        for i in ['a', 'b', 'c']:
            print('  {}\t\t|  {}\t|  {}'.format(i,
                round(getattr(parent_material.structure, 'lattice_constant_{}'.format(i)), 4),
                round(getattr(child_material.structure, 'lattice_constant_{}'.format(i)), 4)))

    # perturb number density/number of atom-sites
    parent_ND = len(parent_material.structure.atom_sites) / parent_material.structure.volume
    child_ND = parent_ND + mutation_strength * (uniform(*number_density_limits) - parent_ND)
    number_of_atoms = int(child_ND * child_material.structure.volume)

    # remove atom-sites, if necessary
    child_material.structure.atom_sites = np.random.choice(
        parent_material.structure.atom_sites,
        min(number_of_atoms, len(parent_material.structure.atom_sites)),
        replace = False).tolist()

    # store original atom-site positions before perturbation in order to compare
    # parent and child values later
    if config['interactive_mode'] == 'on':
        p_x, p_y, p_z = [], [], []
        for atom_site in child_material.structure.atom_sites:
            p_x.append(atom_site.x_frac)
            p_y.append(atom_site.y_frac)
            p_z.append(atom_site.z_frac)

    # perturb atom-site positions
    for atom_site in child_material.structure.atom_sites:
        atom_site.x_frac = random_position(atom_site.x_frac, random(), mutation_strength)
        atom_site.y_frac = random_position(atom_site.y_frac, random(), mutation_strength)
        atom_site.z_frac = random_position(atom_site.z_frac, random(), mutation_strength)

    # add atom-sites, if needed
    if number_of_atoms > len(parent_material.structure.atom_sites):
        for new_sites in range(number_of_atoms - len(parent_material.structure.atom_sites)):
            child_material.structure.atom_sites.append(AtomSites(
                chemical_id = 'A_{}'.format(choice(
                    range(len(parent_material.structure.lennard_jones)))),
                x_frac = random(), y_frac = random(), z_frac = random()))

    ##########################
    # mutate charges
    parent_charges = np.asarray([e.charge for e in parent_material.structure.atom_sites])

    if number_of_atoms < len(parent_charges):
        child_charges = parent_charges[:number_of_atoms]
    else:
        child_charges = np.zeros(number_of_atoms)
        for i in range(len(parent_charges)):
            child_charges[i] = parent_charges[i]

    # adjust charges to account for new number of atom sites
    total_charge = sum(child_charges)
    remaining_charge = total_charge
    while not math.isclose(0., total_charge, abs_tol=1e-9):
        # chose a random atom site
        i = np.random.choice(range(number_of_atoms))
        negative = config['charge_limit'] + child_charges[i]
        positive = config['charge_limit'] - child_charges[i]
        if total_charge > 0:
            delta_charge = uniform(0., min(negative, remaining_charge))
            child_charges[i] -= delta_charge
        else:
            delta_charge = random.uniform(0., min(positive, remaining_charge))
            child_charges[i] += delta_charge
        total_charge = sum(child_charges)
        remaining_charge -= delta_charge
    
    # randomize any new atom sites
    if number_of_atoms > len(parent_charges):
        for i in range(number_of_atoms - len(parent_charges), number_of_atoms):
            q0 = config['charge_limit'] - abs(child_charges[i])
            j = np.random.choice(range(number_of_atoms))
            q1 = config['charge_limit'] - abs(child_charges[j])
            delta_charge = uniform(0, min([q0, q1]))
            child_charges[i] += delta_charge
            child_charges[j] -= delta_charge

    # perturb all charges
    for i in range(number_of_atoms):
        q0 = config['charge_limit'] - abs(child_charges[i])
        j = np.random.choice(range(number_of_atoms))
        q1 = config['charge_limit'] - abs(child_charges[j])
        delta_charge = mutation_strength * uniform(0, min([q0, q1]))
        child_charges[i] += delta_charge
        child_charges[j] -= delta_charge

    print('NET CHARGE : {}'.format(sum(child_charges)))

    for i in range(number_of_atoms):
        child_material.structure.atom_sites[i].charge = child_charges[i]

    if config['interactive_mode'] == 'on':
        print('===================================================================')
        print('FIRST 10 ATOM-SITES')
        print('===================================================================')
        print('  chemical-id\t|  parent position\t\t|  child position')
        print('----------------+-------------------------------+------------------')
        for i in range(min([10, len(p_x)])):
            c = child_material.structure.atom_sites[i]
            print('  {}\t\t|  {}\t|  {}'.format(c.chemical_id,
                (round(p_x[i], 4), round(p_y[i], 4), round(p_z[i], 4)),
                (round(c.x_frac, 4), round(c.y_frac, 4), round(c.z_frac, 4))))

        print('==========================')
        print('NUMBER DENSITY')
        print('==========================')
        print('  parent\t|  child')
        print('----------------+---------')
        print('  {}\t| {}'.format(
            round(len(parent_material.structure.atom_sites) / parent_material.structure.volume, 6),
            round(len(child_material.structure.atom_sites) / child_material.structure.volume, 6)))

    

    return child_material

def new_material(run_id, gen):
    retest_passed = False
    while retest_passed != True:
        # select parent 
        parent_id = select_parent(run_id, gen, config['children_per_generation'])
        parent_material = session.query(Material).get(parent_id)

        # retest results
        retest_passed = parent_material.retest_passed
        retest_loop(parent_material)
        retest_passed = parent_material.retest_passed

    # determine mutation strength
    mutation_strength = determine_mutation_strength(run_id, gen, parent_material)
    if config['pseudomaterial_generator'] == 'hybrid':
        mutation_strength = np.random.choice([1, mutation_strength])
    
    # mutate material
    material = mutate_material(parent_material, mutation_strength, gen)
    return material
