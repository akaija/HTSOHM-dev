import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from sqlalchemy import func, text

from htsohm import config
from htsohm.db import engine, session, Material
from htsohm import pseudomaterial_generator
from htsohm.simulation.run_all import run_all_simulations
from htsohm.simulation.utilities import calc_bin

def count_number_of_materials(run_id):
    return session.query(func.count(Material.id)).filter(Material.run_id==run_id)[0][0]

def query_data(run_id, sim_config):
    selects = "select m.average_epsilon, m.average_sigma, m.number_density, m.unit_cell_volume, "
    joins = " "
    wheres = " where run_id='{}' and ".format(run_id)

    columns = ["average_epsilon", "average_sigma", "number_density", "unit_cell_volume"]
    for i in sim_config:
        params = sim_config[i]
        if params["type"] == "gas_loading":
            selects += "g{}.absolute_volumetric_loading, ".format(i)
            joins += "join gas_loadings g{0} on m.id=g{0}.material_id ".format(i)
            wheres += "g{0}.adsorbate='{1}' and g{0}.pressure={2} and g{0}.temperature={3} and ".format(
                    i, params["adsorbate"], params["pressure"], params["temperature"])
        elif params["type"] == "surface_area":
            selects += "s{}.volumetric_surface_area, ".format(i)
            joins += "join surface_areas s{0} on m.id=s{0}.material_id ".format(i)
            wheres += "s{}.adsorbate='{}' and ".format(i, params["adsorbate"])
        elif params["type"] == "void_fraction":
            selects += "v{}.void_fraction, ".format(i)
            joins += "join void_fractions v{0} on m.id=v{0}.material_id ".format(i)
            wheres += "v{0}.adsorbate='{1}' and v{0}.temperature={2} and ".format(i, params["adsorbate"],
                    params["temperature"])
        columns.append("{}_{}".format(params["type"], i))
    sql = text(selects[:-2] + " from materials m " + joins[:-1] + wheres[:-5])
    result = engine.execute(sql)
    data = pd.DataFrame(result.fetchall())
    data.columns = columns
    X_train = data.loc[:, columns[:4]]
    y_train = data.loc[:, columns[4:]]
    return X_train, y_train

def query_bins(run_id, sim_config):
    selects = "select "
    joins = " "
    wheres = " where run_id='{}' and ".format(run_id)

    columns = []
    for i in sim_config:
        params = sim_config[i]
        if params["type"] == "gas_loading":
            selects += "g{}.bin_value, ".format(i)
            joins += "join gas_loadings g{0} on m.id=g{0}.material_id ".format(i)
            wheres += "g{0}.adsorbate='{1}' and g{0}.pressure={2} and g{0}.temperature={3} and ".format(
                    i, params["adsorbate"], params["pressure"], params["temperature"])
        elif params["type"] == "surface_area":
            selects += "s{}.bin_value, ".format(i)
            joins += "join surface_areas s{0} on m.id=s{0}.material_id ".format(i)
            wheres += "s{}.adsorbate='{}' and ".format(i, params["adsorbate"])
        elif params["type"] == "void_fraction":
            selects += "v{}.bin_value, ".format(i)
            joins += "join void_fractions v{0} on m.id=v{0}.material_id ".format(i)
            wheres += "v{0}.adsorbate='{1}' and v{0}.temperature={2} and ".format(i, params["adsorbate"],
                    params["temperature"])
        columns.append("{}_{}_bin".format(params["type"], i))
    sql = text(selects[:-2] + " from materials m " + joins[:-1] + wheres[:-5])
    result = engine.execute(sql)
    bins = []
    for row in result:
        bin_ = []
        index = 0
        for i in sim_config:
            bin_.append(row[index])
            index += 1
        bins.append(bin_)
    return bins

def worker_run_loop(run_id):
    """
    Args:
        run_id (str): identification string for run.

    Manages overall routine for generating pseudomaterials and simulating their properties.
    Method runs until convergence cutoff is reached.

    """
    # run until library contains one million materials
    while count_number_of_materials(run_id) < 1000000:
        if count_number_of_materials(run_id) < config["general"]["seed_population"]:
            # generate pseudomaterial
            material, structure = pseudomaterial_generator.random.new_material(run_id,
                    config["structure_parameters"])
        else:
            # Load training data
            X_train, y_train = query_data(run_id, config["simulations"])

            # Train model
            model = MultiOutputRegressor(GradientBoostingRegressor())
            model.fit(X_train, y_train)

            # Query bins
            bins = query_bins(run_id, config["simulations"])

            is_in_new_bin = False
            while is_in_new_bin != True:
                # Generate pseudomaterial
                material, structure = pseudomaterial_generator.random.new_material(run_id,
                        config["structure_parameters"])
                # Gather inputs
                X_test = np.array([getattr(material, e) for e in X_train.columns.values]).reshape(1, -1)
                # Predict properties
                y_test = model.predict(X_test)
                # Get bin from prediction
                predicted_bin = []
                index = 0
                for i in config["simulations"]:
                    params = config["simulations"][i]
                    value = y_test[0][index]
                    predicted_bin.append(calc_bin(value, *params["limits"], config["general"]["bins"]))
                    index += 1
                if predicted_bin not in bins:
                    is_in_new_bin = True

        # simulate properties of interest
        run_all_simulations(material, structure)

        # add pseudomaterial to database
        session.add(material)
        session.commit()
        
        sys.stdout.flush()
