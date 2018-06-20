from random import choice, random, randrange, uniform
import pytest

from htsohm.pseudomaterial_generator.random import random_number_density
from htsohm.structure import Structure

@pytest.fixture
def LCs():
    structure = Structure()
    structure.lattice_constants.a = 1
    structure.lattice_constants.b = 1
    structure.lattice_constants.c = 1
    return structure

def test_number_density_bounds():
    assert random_number_density((100000, 100000), LCs(), random()) == 100000

def test_two_site_constraint():
    assert random_number_density((1, 2), LCs(), random()) == 2

