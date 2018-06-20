from sqlalchemy import Column, ForeignKey, Integer, String, Float

from htsohm.db import Base

class SurfaceArea(Base):
    __tablename__ = "surface_areas"

    id = Column(Integer, primary_key=True)

    # relationship with `materials`
    material_id = Column(Integer, ForeignKey("materials.id"))

    # simulation input
    adsorbate      = Column(String(16))

    # simulation output
    unit_cell_surface_area       = Column(Float)
    volumetric_surface_area      = Column(Float)
    gravimetric_surface_area     = Column(Float)

    # bin
    bin_value = Column(Integer)
