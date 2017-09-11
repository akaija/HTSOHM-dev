import sys
import os
import subprocess
import shutil
from datetime import datetime
from uuid import uuid4

import htsohm
from htsohm import config
from htsohm.material_files import write_cif_file, write_mixing_rules
from htsohm.material_files import write_pseudo_atoms, write_force_field

def write_raspa_file(filename, uuid, helium_void_fraction=None):
    """Writes RASPA input file for simulating gas adsorption.

    Args:
        filename (str): path to input file.
        run_id (str): identification string for run.
        material_id (str): uuid for material.

    Writes RASPA input-file.

    """
    simulation_cycles      = config['gas_adsorption_0']['simulation_cycles']
    initialization_cycles  = config['gas_adsorption_0']['initialization_cycles']
    external_temperature   = config['gas_adsorption_0']['external_temperature']
    external_pressure      = config['gas_adsorption_0']['external_pressure']
    adsorbate              = config['gas_adsorption_0']['adsorbate']
       
    with open(filename, "w") as raspa_input_file:
        raspa_input_file.write(
            "SimulationType                 MonteCarlo\n" +
            "NumberOfCycles                 %s\n" % (simulation_cycles) +                  # number of MonteCarlo cycles
            "NumberOfInitializationCycles   %s\n" % (initialization_cycles) +    # number of initialization cycles
            "PrintEvery                     10\n" +
            "RestartFile                    no\n" +
            "\n" +
            "Forcefield     GenericMOFs\n" +
            "CutOff         12.8\n" +                                            # electrostatic cut-off, Angstroms
            "\n" +
            "Framework              0\n" +
            "FrameworkName          %s\n" % (uuid) +
            "UnitCells              1 1 1\n"
        )
        if 'helium_void_fraction' != None:
            raspa_input_file.write("HeliumVoidFraction     %s\n" % (helium_void_fraction))
        raspa_input_file.write(
            "ExternalTemperature    %s\n" % (external_temperature) +               # External temperature, K 
            "ExternalPressure       %s\n" % (external_pressure) +                    # External pressure, Pa
            "\n" +
            "Component 0 MoleculeName               %s\n" % (adsorbate) +
            "            MoleculeDefinition         TraPPE\n" +
            "            TranslationProbability     1.0\n" +
            "            ReinsertionProbability     1.0\n" +
            "            SwapProbability            1.0\n" +
            "            CreateNumberOfMolecules    0\n"
        )

def parse_output(output_file):
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
    with open(output_file) as origin:
        line_counter = 1
        for line in origin:
            if "absolute [mol/kg" in line:
                results['ga0_absolute_molar_loading'] = float(line.split()[5])
            elif "absolute [cm^3 (STP)/g" in line:
                results['ga0_absolute_gravimetric_loading'] = float(line.split()[6])
            elif "absolute [cm^3 (STP)/c" in line:
                results['ga0_absolute_volumetric_loading'] = float(line.split()[6])
            elif "excess [mol/kg" in line:
                results['ga0_excess_molar_loading'] = float(line.split()[5])
            elif "excess [cm^3 (STP)/g" in line:
                results['ga0_excess_gravimetric_loading'] = float(line.split()[6])
            elif "excess [cm^3 (STP)/c" in line:
                results['ga0_excess_volumetric_loading'] = float(line.split()[6])
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
                results['ga0_host_host_avg'] = float(line.split()[1])
                results['ga0_host_host_vdw'] = float(line.split()[5])
                results['ga0_host_host_cou'] = float(line.split()[7])
            elif line_counter == adsorbate_adsorbate_line:
                results['ga0_adsorbate_adsorbate_avg'] = float(line.split()[1])
                results['ga0_adsorbate_adsorbate_vdw'] = float(line.split()[5])
                results['ga0_adsorbate_adsorbate_cou'] = float(line.split()[7])
            elif line_counter == host_adsorbate_line:
                results['ga0_host_adsorbate_avg'] = float(line.split()[1])
                results['ga0_host_adsorbate_vdw'] = float(line.split()[5])
                results['ga0_host_adsorbate_cou'] = float(line.split()[7])
            line_counter += 1

    adsorbate = config['gas_adsorption_0']['adsorbate']
    print(
        "\n%s ADSORPTION\tabsolute\texcess\n" % adsorbate +
        "mol/kg\t\t\t%s\t%s\n" % (results['ga0_absolute_molar_loading'], results['ga0_excess_molar_loading']) +
        "cc/g\t\t\t%s\t%s\n"   % (results['ga0_absolute_gravimetric_loading'], results['ga0_excess_gravimetric_loading']) +
        "cc/cc\t\t\t%s\t%s\n"  % (results['ga0_absolute_volumetric_loading'], results['ga0_excess_volumetric_loading']) +
        "\nENERGIES\thost-host\tadsorbate-adsorbate\thost-adsorbate\n" +
        "avg\t\t%s\t\t%s\t\t%s\n" % (results['ga0_host_host_avg'], results['ga0_adsorbate_adsorbate_avg'], results['ga0_host_adsorbate_avg']) +
        "vdw\t\t%s\t\t%s\t\t%s\n" % (results['ga0_host_host_vdw'], results['ga0_adsorbate_adsorbate_vdw'], results['ga0_host_adsorbate_vdw']) +
        "cou\t\t%s\t\t%s\t\t\t%s\n" % (results['ga0_host_host_cou'], results['ga0_adsorbate_adsorbate_cou'], results['ga0_host_adsorbate_cou'])
    )

    return results

def run(run_id, pseudo_material, helium_void_fraction=None):
    """Runs gas loading simulation.

    Args:
        run_id (str): identification string for run.
        material_id (str): unique identifier for material.
        helium_void_fraction (float): material's calculated void fraction.

    Returns:
        results (dict): gas loading simulation results.

    """
    adsorbate             = config['gas_adsorption_0']['adsorbate']
    simulation_directory  = config['simulations_directory']
    if simulation_directory == 'HTSOHM':
        htsohm_dir = os.path.dirname(os.path.dirname(htsohm.__file__))
        path = os.path.join(htsohm_dir, run_id)
    elif simulation_directory == 'SCRATCH':
        path = os.environ['LOCAL']
    else:
        print('OUTPUT DIRECTORY NOT FOUND.')
    output_dir = os.path.join(path, 'output_%s_%s' % (pseudo_material.uuid, uuid4()))
    print('Output directory :\t%s' % output_dir)
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, '%s_loading.input' % adsorbate)
    write_raspa_file(filename, pseudo_material.uuid, helium_void_fraction)
    write_cif_file(pseudo_material, output_dir)
    write_mixing_rules(pseudo_material, output_dir)
    write_pseudo_atoms(pseudo_material, output_dir)
    write_force_field(output_dir)
    print("Date :\t%s" % datetime.now().date().isoformat())
    print("Time :\t%s" % datetime.now().time().isoformat())
    print("Simulating %s loading in %s..." % (adsorbate, pseudo_material.uuid))
    while True:
        try:
            subprocess.run(
                ['simulate', './%s_loading.input' % adsorbate],
                check = True,
                cwd = output_dir
            )
        
            file_name_part = "output_%s" % (pseudo_material.uuid)
            output_subdir = os.path.join(output_dir, 'Output', 'System_0')
            for file in os.listdir(output_subdir):
                if file_name_part in file:
                    output_file = os.path.join(output_subdir, file)
            print('OUTPUT FILE:\t%s' % output_file)
            results = parse_output(output_file)
            shutil.rmtree(output_dir, ignore_errors=True)
            sys.stdout.flush()
        except FileNotFoundError as err:
            print(err)
            print(err.args)
            continue
        break

    return results
