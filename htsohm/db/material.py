from sqlalchemy import Column, Integer, String, Float
from sqlalchemy.orm import relationship

from htsohm.db import Base, GasLoading, SurfaceArea, VoidFraction

class Material(Base):
    """Declarative class mapping to table storing material/simulation data.

    Attributes:
        id (int): database table primary_key.
        run_id (str): identification string for run.
    """
    __tablename__ = 'materials'
    
    id         = Column(Integer, primary_key=True)
    run_id     = Column(String(50)) 
    seed       = Column(Float)

    # structure properties
    unit_cell_volume     = Column(Float)
    number_density       = Column(Float)
    average_epsilon      = Column(Float)
    average_sigma        = Column(Float)

    # relationships
    gas_loading       = relationship("GasLoading")
    surface_area      = relationship("SurfaceArea")
    void_fraction     = relationship("VoidFraction")

    def __init__(self, run_id=None, seed=None, ):
        """Init material-row.

        Args:
            self (class): row in material table.
            run_id : identification string for run (default = None).

        Initializes row in materials datatable.

        """
        self.seed = seed
        self.run_id = run_id

    def clone(self):
        copy = super(Material, self).clone()
        return copy
