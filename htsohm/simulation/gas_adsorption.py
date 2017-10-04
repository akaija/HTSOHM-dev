import sys
import os
import subprocess
import shutil
from datetime import datetime
from uuid import uuid4
from string import Template

import htsohm
from htsohm import config
from htsohm.simulation.forcefield.files import write_cif_file, write_mixing_rules
from htsohm.simulation.forcefield.files import write_pseudo_atoms, write_force_field
from htsohm.simulation.files import load_and_subs_template
from htsohm.simulation.calculate_bin import calc_bin

def write_raspa_file(filename, material, simulation_config):
    """Writes RASPA input file for simulating gas adsorption.

    Args:
        filename (str): path to input file.
        material (Material): material record.

    Writes RASPA input-file.

    """
    # Set default void fraction value if non was calculated
    if material.vf_helium_void_fraction == None:
        void_fraction = 0.
    else:
        void_fraction = material.vf_helium_void_fraction

    # Load simulation parameters from config
    values = {
            'NumberOfCycles'                : simulation_config['simulation_cycles'],
            'NumberOfInitializationCycles'  : simulation_config['initialization_cycles'],
            'FrameworkName'                 : material.uuid,
            'HeliumVoidFraction'            : void_fraction,
            'ExternalTemperature'           : simulation_config['external_temperature'],
            'MoleculeName'                  : simulation_config['adsorbate']}

    # Adjust printing multiple pressures (isotherms)
    if isinstance(simulation_config['external_pressure'], list):
        values.update({'ExternalPressure'   : '{} {}'.format(*simulation_config['external_pressure'])})
    else:
        values.update({'ExternalPressure'   : simulation_config['external_pressure']})

    # Load template and replace values
    input_data = load_and_subs_template('gas_adsorption.input', values)

    # Write simulation input-file
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(input_data)

def find_output_file(output_dir, pressure):
    # Currently this method only works for pressure values with only two
    # significant digits. This is due to the nature of RASPA's output file
    # nomenclature: any value greater than or equal to 1000000 is printed
    # in scientific notation with only two significant digits.
    if pressure >= 10 ** 6:
        p_string = '{:.1e}'.format(pressure)
    else:
        p_string = str(pressure)

    output_subdir = os.path.join(output_dir, 'Output', 'System_0')
    for file_name in os.listdir(output_subdir):
                if p_string in file_name:
                    return os.path.join(output_subdir, file_name)

def parse_output(output_dir, simulation_config):
    """Parse output file for gas adsorption data.

    Args:
        output_file (str): path to simulation output file.

    Returns:
        results (dict): absolute and excess molar, gravimetric, and volumetric
            gas loadings, as well as energy of average, van der Waals, and
            Coulombic host-host, host-adsorbate, and adsorbate-adsorbate
            interactions.

    """
    results = {}

    external_pressure = simulation_config['external_pressure']
    if isinstance(external_pressure, list):
        ga0_pressure = max(external_pressure)
        ga1_pressure = min(external_pressure)
    else:
        ga0_pressure = external_pressure
        ga1_pressure = None

    # parse first output file then second, in the case of single-pressure
    # adsorption this will be the only output file generated
    f = 'ga0'  # first pressure flag
    for p in [ga0_pressure, ga1_pressure]:
        if p != None:
            output_file = find_output_file(output_dir, p)
            with open(output_file) as origin:
                line_counter = 1
                for line in origin:
                    if "absolute [mol/kg" in line:
                        results['{}_absolute_molar_loading'.format(f)] = float(line.split()[5])
                    elif "absolute [cm^3 (STP)/g" in line:
                        results['{}_absolute_gravimetric_loading'.format(f)] = float(line.split()[6])
                    elif "absolute [cm^3 (STP)/c" in line:
                        results['{}_absolute_volumetric_loading'.format(f)] = float(line.split()[6])
                    elif "excess [mol/kg" in line:
                        results['{}_excess_molar_loading'.format(f)] = float(line.split()[5])
                    elif "excess [cm^3 (STP)/g" in line:
                        results['{}_excess_gravimetric_loading'.format(f)] = float(line.split()[6])
                    elif "excess [cm^3 (STP)/c" in line:
                        results['{}_excess_volumetric_loading'.format(f)] = float(line.split()[6])
                    elif "Average Host-Host energy:" in line:
                        host_host_line = line_counter + 8
                    elif "Average Adsorbate-Adsorbate energy:" in line:
                        adsorbate_adsorbate_line = line_counter + 8
                    elif "Average Host-Adsorbate energy:" in line:
                        host_adsorbate_line = line_counter + 8
                    line_counter += 1

            with open(output_file) as origin:
                line_counter = 1
                for line in origin:
                    if line_counter == host_host_line:
                        results['{}_host_host_avg'.format(f)] = float(line.split()[1])
                        results['{}_host_host_vdw'.format(f)] = float(line.split()[5])
                        results['{}_host_host_cou'.format(f)] = float(line.split()[7])
                    elif line_counter == adsorbate_adsorbate_line:
                        results['{}_adsorbate_adsorbate_avg'.format(f)] = float(line.split()[1])
                        results['{}_adsorbate_adsorbate_vdw'.format(f)] = float(line.split()[5])
                        results['{}_adsorbate_adsorbate_cou'.format(f)] = float(line.split()[7])
                    elif line_counter == host_adsorbate_line:
                        results['{}_host_adsorbate_avg'.format(f)] = float(line.split()[1])
                        results['{}_host_adsorbate_vdw'.format(f)] = float(line.split()[5])
                        results['{}_host_adsorbate_cou'.format(f)] = float(line.split()[7])
                    line_counter += 1
        
            f = 'ga1'  # flag for second pressure gas adsorption simulations

    # calculate bin
    if isinstance(external_pressure, list):
        value = abs(results['ga0_absolute_volumetric_loading'] - results['ga1_absolute_volumetric_loading'])
    else:
        value = results['ga0_absolute_volumetric_loading']
    results['gas_adsorption_bin'] = calc_bin(value, *simulation_config['limits'], simulation_config['bins'])

    return results

def run(material, simulation_config):
    """Runs gas loading simulation.

    Args:
        material_id (Material): material record.

    Returns:
        results (dict): gas loading simulation results.

    """
    # Determine where to write simulation input/output files, create directory
    simulation_directory = config['simulation_directory']
    if simulation_directory == 'HTSOHM':
        htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
        path = os.path.join(htsohm_dir, material.run_id)
    elif simulation_directory == 'SCRATCH':
        path = os.environ['SCRATCH']
    else:
        print('OUTPUT DIRECTORY NOT FOUND.')
    output_dir = os.path.join(path, 'output_%s_%s' % (material.uuid, uuid4()))
    print('Output directory :\t%s' % output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Write simulation input-files
    adsorbate = simulation_config['adsorbate']
    # RASPA input-file
    filename = os.path.join(output_dir, '%s_loading.input' % adsorbate)
    write_raspa_file(filename, material, simulation_config)
    # Pseudomaterial cif-file
    write_cif_file(material, output_dir)
    # Lennard-Jones parameters, force_field_mixing_rules.def
    write_mixing_rules(material, output_dir)
    # Pseudoatom definitions, pseudo_atoms.def (placeholder values)
    write_pseudo_atoms(material, output_dir)
    # Overwritten interactions, force_field.def (none overwritten by default)
    write_force_field(output_dir)

    # Run simulations
    print("Date :\t%s" % datetime.now().date().isoformat())
    print("Time :\t%s" % datetime.now().time().isoformat())
    print("Simulating %s loading in %s..." % (adsorbate, material.uuid))
    while True:
        try:
            subprocess.run(
                ['simulate', './%s_loading.input' % adsorbate],
                check = True,
                cwd = output_dir
            )
        
            print('OUTPUT DIR:\t%s' % output_dir)

            # Parse output
            results = parse_output(output_dir, simulation_config)
            shutil.rmtree(output_dir, ignore_errors=True)
            sys.stdout.flush()
        except FileNotFoundError as err:
            print(err)
            print(err.args)
            continue
        break

    return results
