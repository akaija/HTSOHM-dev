import sys
import os
import subprocess
import shutil
from datetime import datetime
from uuid import uuid4
from string import Template

import htsohm
from htsohm.material_files import write_cif_file, write_mixing_rules
from htsohm.material_files import write_pseudo_atoms, write_force_field
from htsohm.simulation.calculate_bin import calc_bin

def write_raspa_file(filename, uuid, params):
    """Writes RASPA input file for calculating surface area.

    Args:
        filename (str): path to input file.
        run_id (str): identification string for run.
        material_id (str): uuid for material.

    Writes RASPA input-file.

    """
    # Load simulation parameters from config
    values = {
            'NumberOfCycles'                : params['simulation_cycles'],
            'FrameworkName'                 : uuid}

    # Load template and replace values
    htsohm_dir = os.path.dirname(htsohm.__file__)
    input_template_path = os.path.join(htsohm_dir, 'simulation', 'surface_area.input')
    input_file = open(input_template_path)
    template = Template( input_file.read() )
    input_data = template.substitute(values)

    print(input_data)

    # Write simulation input-file
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(input_data)

def parse_output(output_file, params):
    """Parse output file for void fraction data.

    Args:
        output_file (str): path to simulation output file.

    Returns:
        results (dict): total unit cell, gravimetric, and volumetric surface
            areas.

    """
    results = {}
    with open(output_file) as origin:
        count = 0
        for line in origin:
            if "Surface area" in line:
                if count == 0:
                    results['sa_unit_cell_surface_area'] = float(line.split()[2])
                    count = count + 1
                elif count == 1:
                    results['sa_gravimetric_surface_area'] = float(line.split()[2])
                    count = count + 1
                elif count == 2:
                    results['sa_volumetric_surface_area'] = float(line.split()[2])

    print(
        "\nSURFACE AREA\n" +
        "%s\tA^2\n"      % (results['sa_unit_cell_surface_area']) +
        "%s\tm^2/g\n"    % (results['sa_gravimetric_surface_area']) +
        "%s\tm^2/cm^3"   % (results['sa_volumetric_surface_area']))

    results['surface_area_bin'] = calc_bin(results['sa_volumetric_surface_area'],
            *params['limits'], params['bins'])

    return results

def run(material, simulation_config):
    """Runs surface area simulation.

    Args:
        material (Material): material record.

    Returns:
        results (dict): surface area simulation results.

    """
    # Determine where to write simulation input/output files, create directory
    simulation_directory  = simulation_config['directory']
    if simulation_directory == 'HTSOHM':
        htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
        path = os.path.join(htsohm_dir, material.run_id)
    elif simulation_directory == 'SCRATCH':
        path = os.environ['SCRATCH']
    else:
        print('OUTPUT DIRECTORY NOT FOUND.')
    output_dir = os.path.join(path, 'output_%s_%s' % (material.uuid, uuid4()))
    print("Output directory :\t%s" % output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Write simulation input-files
    params = simulation_config['surface_area']
    # RASPA input-file
    filename = os.path.join(output_dir, "SurfaceArea.input")
    write_raspa_file(filename, material.uuid, params)
    # Pseudomaterial cif-file
    write_cif_file(material, output_dir)
    # Lennard-Jones parameters, force_field_mixing_rules.def
    write_mixing_rules(material, output_dir)
    # Pseudoatom definitions, pseudo_atoms.def (placeholder values)
    write_pseudo_atoms(material, output_dir)
    # Overwritten interactions, force_field.def (none overwritten by default)
    write_force_field(output_dir)

    # Run simulations
    while True:
        try:
            print("Date :\t%s" % datetime.now().date().isoformat())
            print("Time :\t%s" % datetime.now().time().isoformat())
            print("Calculating surface area of %s..." % (material.uuid))
            subprocess.run(['simulate', './SurfaceArea.input'], check=True, cwd=output_dir)
            filename = "output_%s_2.2.2_298.000000_0.data" % (material.uuid)
            output_file = os.path.join(output_dir, 'Output', 'System_0', filename)

            # Parse output
            results = parse_output(output_file, params)
            shutil.rmtree(output_dir, ignore_errors=True)
            sys.stdout.flush()
        except (FileNotFoundError, KeyError) as err:
            print(err)
            print(err.args)
            continue
        break

    return results
