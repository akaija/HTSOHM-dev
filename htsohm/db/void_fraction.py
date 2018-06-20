from sqlalchemy import Column, ForeignKey, Integer, String, Float

from htsohm.db import Base

class VoidFraction(Base):
    __tablename__ = "void_fractions"

    id = Column(Integer, primary_key=True)

    # relationship with `materials`
    material_id = Column(Integer, ForeignKey("materials.id"))

    # simulation input
    adsorbate      = Column(String(16))
    temperature    = Column(Float)

    # simulation output
    void_fraction = Column(Float)

    # bin
    bin_value = Column(Integer)
