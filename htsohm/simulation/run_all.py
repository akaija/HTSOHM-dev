from htsohm import config
from htsohm import simulation

def get_simulation(simulation_type):
    if simulation_type == 'gas_adsorption':
        return simulation.gas_adsorption
    elif simulation_type == 'surface_area':
        return simulation.surface_area
    elif simulation_type == 'helium_void_fraction':
        return simulation.helium_void_fraction
    elif simulation_type == 'artificial_gas_adsorption':
        return simulation.artificial_gas_adsorption
    elif simulation_type == 'artificial_surface_area':
        return simulation.artificial_surface_area
    elif simulation_type == 'artificial_void_fraction':
        return simulation.artificial_void_fraction
    else:
        raise Exception('Simulation-type not found!')

def run_all_simulations(material):
    """Simulate helium void fraction, gas loading, and surface area.

    Args:
        material (sqlalchemy.orm.query.Query): material to be analyzed.

    Depending on properties specified in config, adds simulated data for helium
    void fraction, gas loading, heat of adsorption, surface area, and
    corresponding bins to row in database corresponding to the input-material.
        
    """
    simulation_config = config['simulations']
    for simulation_type in simulation_config:
        results = get_simulation(simulation_type).run(material, simulation_config[simulation_type])
        material.update_from_dict(results)


