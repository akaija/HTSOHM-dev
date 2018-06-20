from sqlalchemy import Column, ForeignKey, Integer, String, Float

from htsohm.db import Base

class GasLoading(Base):
    __tablename__ = "gas_loadings"

    id = Column(Integer, primary_key=True)

    # relationship with `materials`
    material_id = Column(Integer, ForeignKey("materials.id"))

    # simulation input
    adsorbate      = Column(String(16))
    pressure       = Column(Float)
    temperature    = Column(Float)

    # simulation output
    absolute_volumetric_loading      = Column(Float)
    absolute_gravimetric_loading     = Column(Float)
    absolute_molar_loading           = Column(Float)
    excess_volumetric_loading        = Column(Float)
    excess_gravimetric_loading       = Column(Float)
    excess_molar_loading             = Column(Float)
    host_host_avg                    = Column(Float)
    host_host_vdw                    = Column(Float)
    host_host_cou                    = Column(Float) 
    adsorbate_adsorbate_avg          = Column(Float)
    adsorbate_adsorbate_vdw          = Column(Float)
    adsorbate_adsorbate_cou          = Column(Float)
    host_adsorbate_avg               = Column(Float)
    host_adsorbate_vdw               = Column(Float)
    host_adsorbate_cou               = Column(Float)

    # bin
    bin_value = Column(Integer)
