from htsohm import config
from htsohm import simulation

def get_simulation(simulation_type):
    if  simulation_type == 'gas_adsorption_0':
        return simulation.gas_adsorption_0
    elif  simulation_type == 'gas_adsorption_1':
        return simulation.gas_adsorption_1
    elif simulation_type == 'surface_area':
        return simulation.surface_area
    elif simulation_type == 'void_fraction':
        return simulation.void_fraction
    else:
        raise Exception('Simulation-type not found!')

def run_all_simulations(material, structure):
    """Simulate helium void fraction, gas loading, and surface area.

    Args:
        material (sqlalchemy.orm.query.Query): material to be analyzed.

    Depending on properties specified in config, adds simulated data for helium
    void fraction, gas loading, heat of adsorption, surface area, and
    corresponding bins to row in database corresponding to the input-material.
        
    """
    for simulation_number in config["simulations"]:
        simulation_config = config["simulations"][simulation_number]
        getattr(simulation, simulation_config["type"]).run(material, structure, simulation_config)
        #material.update_from_dict(results)


