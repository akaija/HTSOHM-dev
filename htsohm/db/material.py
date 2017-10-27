
import sys
import uuid

from sqlalchemy import Column, ForeignKey, Integer, String, Float, Boolean
from sqlalchemy.sql import text

from htsohm import config
from htsohm.db import Base, session, engine

class Material(Base):
    """Declarative class mapping to table storing material/simulation data.

    Attributes:
        id (int): database table primary_key.
        run_id (str): identification string for run.
        uuid (str): unique identification string for material.
        parent_id (int): uuid of parent mutated to create material.
        generation (int): iteration in overall bin-mutate-simulate routine.
        generation_index (int): order material was created in generation (used
            to determine when all materials appear in database for a particular
            generation).
        retest_num (int): iteration in re-test routine for statistical errors.
        retest_methane_loading_sum (float): sum of all absolute volumetric
            methane loadings calculated in re-test routine.
        retest_surface_area_sum (float): sum of all volumetric surface areas
            calculated in re-test routine.
        retest_void_fraction_sum (float): sum of all helium void fractions
            calculated in re-test routine.
        retest_passed (bool): true if the average of all re-test results is
            within the acceptable range of deviation.
        ga_absolute_volumetric_loading (float): absolute volumetric loading.
        ga_absolute_gravimetric_loading (float): absolute gravimetric loading.
        ga_absolute_molar_loading (float): absolute molar loading.
        ga_excess_volumetric_loading (float): excess volumetric loading.
        ga_excess_gravimetric_loading (float): excess gravimetric loading.
        ga_excess_molar_loading (float): excess molar loading.
        ga_host_host_avg (float): average energy of host-host interactions.
        ga_host_host_vdw (float): energy of host-host van der Waals
            interactions.
        ga_host_host_cou (float): energy of host-host electrostatic
            interactions.
        ga_adsorbate_adsorbate_avg (float): average energy of adsorbate-
            adsorbate interactions.
        ga_adsorbate_adsorbate_vdw (float): energy of adsorbate-adsorbate van
            der Waals interactions.
        ga_adsorbate_adsorbate_cou (float): energy of adsorbate-adsorbate
            electrostatic interactions.
        ga_host_adsorbate_avg (float): average energy of host-adsorbate
            interactions.
        ga_host_adsorbate_vdw (float): energy of host-adsorbate van der Waals
            interactions.
        ga_host_adsorbate_cou (float): energy of host-adsorbate electrostatic
            interactions.
        sa_unit_cell_surface_area (float): surface area of unit-cell.
        sa_volumetric_surface_area (float): surface area per unit volume.
        sa_gravimetric_surface_area (float): surface area per unit mass.
        vf_helium_void_fraction (float): void fraction measured with helium
            probe.
        methane_loading_bin (int): region of methane loading-space
            corresponding to the material's simulation results.
        surface_area_bin (int): region of surface area-space corresponding to
            the material's simulation results.
        void_fraction_bin (int): region of void fraction-space corresponding to
            the material's simulation results.

    """
    __tablename__ = 'materials'
    # COLUMN                                                 UNITS
    id = Column(Integer, primary_key=True)                 # dimm.
    run_id = Column(String(50))                            # dimm.
    uuid = Column(String(40))
    parent_id = Column(Integer)                            # dimm.
    generation = Column(Integer)                           # generation#
    generation_index = Column(Integer)                     # index order of row in generation

    # retest columns
    retest_num = Column(Integer, default=0)
    retest_gas_adsorption_0_sum = Column(Float, default=0)
    retest_gas_adsorption_1_sum = Column(Float, default=0)
    retest_surface_area_sum = Column(Float, default=0)
    retest_void_fraction_sum = Column(Float, default=0)
    retest_passed = Column(Boolean)                        # will be NULL if retest hasn't been run

    # data collected
    #   gas adsorption 0
    ga0_absolute_volumetric_loading = Column(Float)            # cm^3 / cm^3
    ga0_absolute_gravimetric_loading = Column(Float)           # cm^3 / g
    ga0_absolute_molar_loading = Column(Float)                 # mol / kg
    ga0_excess_volumetric_loading = Column(Float)              # cm^3 / cm^3
    ga0_excess_gravimetric_loading = Column(Float)             # cm^3 / g
    ga0_excess_molar_loading = Column(Float)                   # mol /kg
    ga0_host_host_avg = Column(Float)                          # K
    ga0_host_host_vdw = Column(Float)                          # K
    ga0_host_host_cou = Column(Float)                          # K
    ga0_adsorbate_adsorbate_avg = Column(Float)                # K
    ga0_adsorbate_adsorbate_vdw = Column(Float)                # K
    ga0_adsorbate_adsorbate_cou = Column(Float)                # K
    ga0_host_adsorbate_avg = Column(Float)                     # K
    ga0_host_adsorbate_vdw = Column(Float)                     # K
    ga0_host_adsorbate_cou = Column(Float)                     # K
    #   gas adsorption 1
    ga1_absolute_volumetric_loading = Column(Float)            # cm^3 / cm^3
    ga1_absolute_gravimetric_loading = Column(Float)           # cm^3 / g
    ga1_absolute_molar_loading = Column(Float)                 # mol / kg
    ga1_excess_volumetric_loading = Column(Float)              # cm^3 / cm^3
    ga1_excess_gravimetric_loading = Column(Float)             # cm^3 / g
    ga1_excess_molar_loading = Column(Float)                   # mol /kg
    ga1_host_host_avg = Column(Float)                          # K
    ga1_host_host_vdw = Column(Float)                          # K
    ga1_host_host_cou = Column(Float)                          # K
    ga1_adsorbate_adsorbate_avg = Column(Float)                # K
    ga1_adsorbate_adsorbate_vdw = Column(Float)                # K
    ga1_adsorbate_adsorbate_cou = Column(Float)                # K
    ga1_host_adsorbate_avg = Column(Float)                     # K
    ga1_host_adsorbate_vdw = Column(Float)                     # K
    ga1_host_adsorbate_cou = Column(Float)                     # K
    #   surface area
    sa_unit_cell_surface_area = Column(Float)                 # angstroms ^ 2
    sa_volumetric_surface_area = Column(Float)                # m^2 / cm^3
    sa_gravimetric_surface_area = Column(Float)               # m^2 / g
    #   void fraction
    vf_helium_void_fraction = Column(Float)                   # dimm.

    # bins
    gas_adsorption_bin = Column(Integer)                     # dimm.
    surface_area_bin = Column(Integer)                        # dimm.
    void_fraction_bin = Column(Integer)                       # dimm.
    average_epsilon = Column(Float)
    average_sigma = Column(Float)
    number_density = Column(Float)
    unit_cell_volume = Column(Float)

    def __init__(self, run_id=None, ):
        """Init material-row.

        Args:
            self (class): row in material table.
            run_id : identification string for run (default = None).

        Initializes row in materials datatable.

        """
        self.uuid = str(uuid.uuid4())
        self.run_id = run_id

    @property
    def bin(self):
        """Determine material's structure-property bin.

        Args:
            self (class): row in material table.

        Returns:
            The bin corresponding to a material's gas loading, void
            fraction, and surface area data and their postion in this three-
            dimension parameter-space.

        """
        return [self.gas_adsorption_bin, self.surface_area_bin, self.void_fraction_bin]

    def calculate_generation_index(self):
        """Determine material's generation-index.

        Args:
            self (class): row in material table.

        Returns:
            The generation-index is used to count the number of materials
            present in the database (that is to have all definition-files in
            the RASPA library and simulation data in the materials datatable).
            This attribute is used to determine when to stop adding new
            materials to one generation and start another.

        """
        return session.query(Material).filter(
                Material.run_id==self.run_id,
                Material.generation==self.generation,
                Material.id < self.id,
            ).count()

    def calculate_percent_children_in_bin(self):
        """Determine number of children in the same bin as their parent.

        Args:
            self (class): row in material table.

        Returns:
            Fraction of children in the same bin as parent (self).
        """
        sql = text("""
            select
                m.gas_adsorption_bin,
                m.surface_area_bin,
                m.void_fraction_bin,
                (
                    m.gas_adsorption_bin = p.gas_adsorption_bin and
                    m.surface_area_bin = p.surface_area_bin and
                    m.void_fraction_bin = p.void_fraction_bin
                ) as in_bin
            from materials m
            join materials p on (m.parent_id = p.id)
            where m.generation = :gen
              and m.run_id = :run_id
              and p.gas_adsorption_bin = :ga_bin
              and p.surface_area_bin = :sa_bin
              and p.void_fraction_bin = :vf_bin
        """)

        rows = engine.connect().execute(
            sql,
            gen=self.generation,
            run_id=self.run_id,
            ga_bin=self.gas_adsorption_bin,
            sa_bin=self.surface_area_bin,
            vf_bin=self.void_fraction_bin
        ).fetchall()

        return len([ r for r in rows if r.in_bin ]) / len(rows)

    def calculate_retest_result(self, tolerance):
        """Determine if material has passed re-testing routine.

        Args:
            self (class): row in material table.
            tolerance (float): acceptable deviation as percent of parameter
                bin-width.

        Returns:
            (bool) True if material has NOT failed any of all re-tests.

        """
        simulations = config['material_properties']
        number_of_bins = config['number_of_convergence_bins']

        if 'gas_adsorption_0' in simulations:
            ga0_o = self.ga0_absolute_volumetric_loading    # initally-calculated values
            ga0_mean = self.retest_gas_adsorption_0_sum / self.retest_num
            ga0_limits = config['gas_adsorption_0']['limits']
            ga0_width = (ga0_limits[1] - ga0_limits[0]) / number_of_bins
        else:
            ga0_o = 0
            ga0_mean = 0
            ga0_width = 0

        if 'gas_adsorption_1' in simulations:
            ga1_o = self.ga1_absolute_volumetric_loading    # initally-calculated values
            ga1_mean = self.retest_gas_adsorption_1_sum / self.retest_num
            ga1_limits = config['gas_adsorption_1']['limits']
            ga1_width = (ga1_limits[1] - ga1_limits[0]) / number_of_bins
        else:
            ga1_o = 0
            ga1_mean = 0
            ga1_width = 0

        if 'surface_area' in simulations:
            sa_o = self.sa_volumetric_surface_area
            sa_mean = self.retest_surface_area_sum / self.retest_num
            sa_limits = config['surface_area']['limits']
            sa_width = (sa_limits[1] - sa_limits[0]) / number_of_bins
        else:
            sa_o = 0
            sa_mean = 0
            sa_width = 0

        if 'helium_void_fraction' in simulations:
            vf_o = self.vf_helium_void_fraction
            vf_mean = self.retest_void_fraction_sum / self.retest_num
            vf_limits = config['surface_area']['limits']
            vf_width = (vf_limits[1] - vf_limits[0]) / number_of_bins
        else:
            vf_o = 0
            vf_mean = 0
            vf_width = 0

        retest_failed = (
            abs(ga0_mean - ga0_o) >= tolerance * ga0_width and ga0_o != 0 or
            abs(ga1_mean - ga1_o) >= tolerance * ga1_width and ga1_o != 0 or
            abs(sa_mean - sa_o) >= tolerance * sa_width and sa_o != 0 or
            abs(vf_mean - vf_o) >= tolerance * vf_width and vf_o != 0
        )

        return not retest_failed

    def calculate_retest_result2(self, tolerance, config):
        """Determine if material has passed re-testing routine.

        Args:
            self (class): row in material table.
            tolerance (float): acceptable deviation as percent of parameter
                bin-width.

        Returns:
            (bool) True if material has NOT failed any of all re-tests.

        """
        simulations = config['material_properties']
        number_of_bins = config['number_of_convergence_bins']

        if 'gas_adsorption_0' in simulations:
            ga0_o = self.ga0_absolute_volumetric_loading    # initally-calculated values
            ga0_mean = self.retest_gas_adsorption_0_sum / self.retest_num
            ga0_limits = config['gas_adsorption_0']['limits']
            ga0_width = (ga0_limits[1] - ga0_limits[0]) / number_of_bins
        else:
            ga0_o = 0
            ga0_mean = 0
            ga0_width = 0

        if 'gas_adsorption_1' in simulations:
            ga1_o = self.ga1_absolute_volumetric_loading    # initally-calculated values
            ga1_mean = self.retest_gas_adsorption_1_sum / self.retest_num
            ga1_limits = config['gas_adsorption_1']['limits']
            ga1_width = (ga1_limits[1] - ga1_limits[0]) / number_of_bins
        else:
            ga1_o = 0
            ga1_mean = 0
            ga1_width = 0

        if 'surface_area' in simulations:
            sa_o = self.sa_volumetric_surface_area
            sa_mean = self.retest_surface_area_sum / self.retest_num
            sa_limits = config['surface_area']['limits']
            sa_width = (sa_limits[1] - sa_limits[0]) / number_of_bins
        else:
            sa_o = 0
            sa_mean = 0
            sa_width = 0

        if 'helium_void_fraction' in simulations:
            vf_o = self.vf_helium_void_fraction
            vf_mean = self.retest_void_fraction_sum / self.retest_num
            vf_limits = config['surface_area']['limits']
            vf_width = (vf_limits[1] - vf_limits[0]) / number_of_bins
        else:
            vf_o = 0
            vf_mean = 0
            vf_width = 0

        retest_failed = (
            abs(ga0_mean - ga0_o) >= tolerance * ga0_width and ga0_o != 0 or
            abs(ga1_mean - ga1_o) >= tolerance * ga1_width and ga1_o != 0 or
            abs(sa_mean - sa_o) >= tolerance * sa_width and sa_o != 0 or
            abs(vf_mean - vf_o) >= tolerance * vf_width and vf_o != 0
        )

        return not retest_failed
