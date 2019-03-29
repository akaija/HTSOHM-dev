# processing data from PORUS

pre_processed_db = "man003"
post_processed_db = pre_processed_db + "_processed"
connection_string_base = "postgresql://porus@127.0.01/"

run_id = "2018-08-22T10:53:27.019962"

bins = 20

# schema for post-processed data

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy import Column, Float, Integer, String

post_processed_connection_string = connection_string_base + post_processed_db
post_processed_engine = create_engine(post_processed_connection_string)
post_processed_session = sessionmaker(bind=post_processed_engine)()

class Base(object):

    def clone(self):
        copy = self.__class__()
        for col in self.__table__.columns:
            val = getattr(self, col.name)
            setattr(copy, col.name, val)
        return copy

    def update_from_dict(self, d):
        for k, v in d.items():
            setattr(self, k, v)

Base = declarative_base(cls=Base)

class ProcessedData(Base):
    __tablename__ = "processed_data"

    id             = Column(Integer, primary_key=True)
    seed           = Column(Integer)
    run_id         = Column(String(50))

    # properties of interest
    average_vdw_radius            = Column(Float)
    average_well_depth            = Column(Float)
    co2_selectivity               = Column(Float)
    co2_adsorption_uptake         = Column(Float)
    co2_working_capacity          = Column(Float)
    coulombic_interaction         = Column(Float)
    heat_of_adsorption            = Column(Float)
    helium_void_fraction          = Column(Float)
    number_density                = Column(Float)
    regenerability                = Column(Float)
    sorbent_selection_parameter   = Column(Float)
    vdw_interaction               = Column(Float)
    volumetric_surface_area       = Column(Float)

    # bins
    average_vdw_radius_bin            = Column(Integer)
    average_well_depth_bin            = Column(Integer)
    co2_selectivity_bin               = Column(Integer)
    co2_adsorption_uptake_bin         = Column(Integer)
    co2_working_capacity_bin          = Column(Integer)
    coulombic_interaction_bin         = Column(Integer)
    heat_of_adsorption_bin            = Column(Integer)
    helium_void_fraction_bin          = Column(Integer)
    number_density_bin                = Column(Integer)
    regenerability_bin                = Column(Integer)
    sorbent_selection_parameter_bin   = Column(Integer)
    vdw_interaction_bin               = Column(Integer)
    volumetric_surface_area_bin       = Column(Integer)

    def __init__(self, seed=None, run_id=None, ):
        self.seed = seed
        self.run_id = run_id

Base.metadata.create_all(post_processed_engine)
Base.metadata.bind = post_processed_engine

# populating properties of interest

from sqlalchemy import text

from htsohm.db import engine as pre_processed_engine

def query_all_data():
    print("Starting query...")
    pre_processed_connection_string = connection_string_base + pre_processed_db
    pre_processed_engine = create_engine(pre_processed_connection_string)
    pre_processed_session = sessionmaker(bind=pre_processed_engine)()
    
    processed_seeds = post_processed_session.query(ProcessedData.seed).all()
    
    sql = text("""
    select
      co2_ads.absolute_molar_loading, co2_ads.absolute_volumetric_loading, co2_ads.host_adsorbate_cou,
      co2_ads.host_adsorbate_vdw, n2_ads.absolute_molar_loading, n2_ads.host_adsorbate_cou, n2_ads.host_adsorbate_vdw,
      co2_des.absolute_molar_loading, co2_des.absolute_volumetric_loading, n2_des.absolute_molar_loading, m.seed,
      m.average_sigma, m.average_epsilon, m.number_density, v.void_fraction, s.volumetric_surface_area
    from materials m
    join gas_loadings co2_ads on m.id=co2_ads.material_id
    join gas_loadings n2_ads on m.id=n2_ads.material_id
    join gas_loadings co2_des on m.id=co2_des.material_id
    join gas_loadings n2_des on m.id=n2_des.material_id
    join void_fractions v on m.id=v.material_id
    join surface_areas s on m.id=s.material_id
    where run_id=:run_id
      and co2_ads.adsorbate='CO2' and co2_ads.pressure=10000 and co2_ads.temperature=298
      and n2_ads.adsorbate='N2' and n2_ads.pressure=90000 and n2_ads.temperature=298
      and co2_des.adsorbate='CO2' and co2_des.pressure=1000 and co2_des.temperature=298
      and n2_des.adsorbate='N2' and n2_des.pressure=9000 and n2_des.temperature=298
      and v.adsorbate='helium' and v.temperature=298
      and s.adsorbate='N2'
    limit 20000
    """)
    rows = pre_processed_engine.connect().execute(sql, run_id=run_id).fetchall()
    
    #################################################
    # index | value
    #-------+----------------------------------------
    #     0 | co2_ads.absolute_molar_loading,
    #     1 | co2_ads.absolute_volumetric_loading,
    #     2 | co2_ads.host_adsorbate_cou,
    #     3 | co2_ads.host_adsorbate_vdw,
    #     4 | n2_ads.absolute_molar_loading,
    #     5 | n2_ads.host_adsorbate_cou,
    #     6 | n2_ads.host_adsorbate_vdw,
    #     7 | co2_des.absolute_molar_loading,
    #     8 | co2_des.absolute_volumetric_loading,
    #     9 | n2_des.absolute_molar_loading,
    #    10 | m.seed,
    #    11 | m.average_sigma,
    #    12 | m.average_epsilon,
    #    13 | m.number_density,
    #    14 | v.void_fraction,
    #    15 | s.volumetric_surface_area
    #################################################
    
    for row in rows:
        if row[10] not in processed_seeds:
            pd = ProcessedData(row[10], run_id)
            
            ########################################
            # copied values
            ########################################
            pd.average_vdw_radius       = row[11]
            pd.average_well_depth       = row[12]
            pd.co2_adsorption_uptake    = row[1]
            pd.helium_void_fraction     = row[14]
            pd.number_density           = row[13]
            pd.volumetric_surface_area  = row[15]
        
            ########################################
            # aggregated values
            ########################################
            # selectivity
            pd.co2_selectivity = row[0] / row[4] * (90000 / 10000)
            # working capacity
            pd.co2_working_capacity = row[1] - row[8]
            # coulombic heat contribution
            pd.coulombic_interaction = (row[2] + row[5]) * 10 ** -6
            # van der waals heat contribution
            pd.vdw_interaction = (row[3] + row[6]) * 10 ** -6
            # heat of adsorption
            pd.heat_of_adsorption = pd.coulombic_interaction + pd.vdw_interaction
            # regenerability
            pd.regenerability = (row[0] - row[7]) / row[0] * 100
            # sorbent selection parameter
            des_s = (row[7] / row[9] * (9000 / 1000))
            co2_wc = row[0] - row[7]
            n2_wc = row[4] - row[9]
            pd.sorbent_selection_parameter = pd.co2_selectivity / des_s * (co2_wc / n2_wc)
            
            post_processed_session.add(pd)
    post_processed_session.commit()
    print("...done!")

# binning properties of interest

from math import ceil

from sqlalchemy import func

#temp
import numpy as np

def custom_round(x, y):
    return int(y * round(float(x) / y))

def get_min_max(column):
    # custom limits
    custom = {"average_vdw_radius"            : [0.5, 5.5],
              "sorbent_selection_parameter"   : [0, 20],
              "helium_void_fraction"          : [0, 1],
              "number_density"                : [0, 0.02],
              "vdw_interaction"               : [-10, 0],
              "coulombic_interaction"         : [-20, 0],
              "heat_of_adsorption"            : [-30, 0],
              "co2_working_capacity"          : [0, 120],
              "volumetric_surface_area"       : [0, 4500],
              "co2_selectivity"               : [0, 1100],
              "regenerability"                : [0, 100]
    }
    if column not in custom:
        x = getattr(ProcessedData, column)
        min_, max_ = post_processed_session.query(func.min(x), func.max(x)).filter(ProcessedData.run_id==run_id)[0]
        if (max_ - min_) >= bins:
            min_, max_ = round(min_ / bins) * bins, round(max_ / bins) * bins
            step = custom_round((max_ - min_) / bins, 5)
            max_ = min_ + bins * step
    else:
        [min_, max_] = [*custom[column]]
    return min_, max_

def print_ticks(column):
    min_, max_ = get_min_max(column)
    step = (max_ - min_) / bins
    ticks = np.arange(min_, max_ + step, step)
    print(column)
    print("\t{}".format(ticks))

def print_all_ticks():
    columns = ProcessedData.__table__.columns.keys()
    for column in columns:
        if column not in ["id", "seed", "run_id"] and "bin" not in column:
            print_ticks(column)

def bin_property(column, bins, run_id):
    min_, max_ = get_min_max(column)
#    print(column)
#    print(min_, max_)
    step = (max_ - min_) / bins
    print
    sql = text("""update processed_data
    set {}=({}-:min_)/:step
    where run_id=:run_id""".format("{}_bin".format(column), column))
    post_processed_engine.connect().execute(sql, min_=min_, step=step, run_id=run_id)

def bin_all_data():
    print("Starting binning...")
    columns = ProcessedData.__table__.columns.keys()
    for column in columns:
        if column not in ["id", "seed", "run_id"] and "bin" not in column:
            bin_property(column, bins, run_id)
    print("...done!")

def process_data():
    query_all_data()
    bin_all_data()
    print("Processing complete.")
