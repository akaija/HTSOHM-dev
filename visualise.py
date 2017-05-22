import os

import bpy
import numpy as np

from htsohm.db.__init__ import session, Material, Structure
from htsohm.files import load_config_file

scene       = bpy.context.scene
add_cube    = bpy.ops.mesh.primitive_cube_add
add_sphere  = bpy.ops.mesh.primitive_uv_sphere_add

def MakeMaterial(name, diffuse, specular, alpha):
    material = bpy.data.materials.new(name)
    material.diffuse_color = diffuse
    material.diffuse_shader = 'LAMBERT'
    material.diffuse_intensity = 1.0
    material.specular_color = specular
    material.specular_shader = 'COOKTORR'
    material.specular_intensity = 0.5
    material.alpha = alpha
    material.ambient = 1
    return material

def SetMaterial(bpy_object, material):
    something = bpy_object.data
    something.materials.append(material)

def create_unit_cell(a, b, c, uuid):
    dimensions = np.array([a, b, c]) / 2.
    add_cube(location = dimensions)
    bpy.ops.transform.resize(value = dimensions)
    bpy.data.objects['Cube'].select = True
    cell_name = 'UnitCell_%s' % uuid
    bpy.context.scene.objects[0].name = cell_name
    black = MakeMaterial('black', (0, 0, 0), (1, 1, 1), 1)
    SetMaterial(bpy.context.object, black)
    bpy.data.objects[cell_name].select = False

def add_atom_site(x, y, z, radius,  material, uuid, site_count):
    add_sphere(location = (x, y, z), size = radius)
    bpy.data.objects['Sphere'].select = True
    bpy.ops.object.shade_smooth()
    SetMaterial(bpy.context.object, material)
    site_name = 'AtomSite_%s' % site_count
    bpy.context.scene.objects[0].name = site_name
    apply_cut_off(site_name, uuid)
    bpy.data.objects[site_name].select = False
    site_count += 1
    return site_count

def apply_cut_off(name, uuid):
    unit_cell = bpy.data.objects['UnitCell_%s' % uuid]
    atom_site = bpy.data.objects[name]
    cut_off = atom_site.modifiers.new(type='BOOLEAN', name="cut_off")
    cut_off.object = unit_cell
    cut_off.operation = 'INTERSECT'
    try:
        bpy.context.scene.objects.active = atom_site
        bpy.ops.object.modifier_apply(apply_as='DATA', modifier="cut_off")
    except:
        print('...something may have gone wrong applying cut-off...')
        pass

def add_periodic_segments(atom_site, chemical_species, mesh_material, uuid, a, b, c, site_count):
    for species in chemical_species:
        if species['chemical'] == atom_site.chemical_id:
            radius = species['radius']
            material = species['material']
    x = atom_site.x_frac * a
    y = atom_site.y_frac * b
    z = atom_site.z_frac * c
    if x <= radius:
        site_count = add_atom_site(
            x + a, y, z, radius, mesh_material, uuid, site_count
        )
    if x + radius >= a:
        site_count = add_atom_site(
            x - a, y, z, radius, mesh_material, uuid, site_count
        )
    if y <= radius:
        site_count = add_atom_site(
            x, y + b, z, radius, mesh_material, uuid, site_count
        )
    if y + radius >= b:
        site_count = add_atom_site(
            x, y - b, z, radius, mesh_material, uuid, site_count
        )
    if z <= radius:
        site_count = add_atom_site(
            x, y, z + c, radius, mesh_material, uuid, site_count
        )
    if z + radius >= c:
        site_count = add_atom_site(
            x, y, z - c, radius, mesh_material, uuid, site_count
        )
    return site_count

def visualise_pseudo_material(uuid):
    print('\nvisualising:\t%s' % uuid)

    material_id = session.query(Material.id).filter(Material.uuid==uuid).one()[0]
    material = session.query(Material).get(material_id)
    structure = material.structure

    print('\ncreating unit cell...')
    a = structure.lattice_constant_a
    b = structure.lattice_constant_b
    c = structure.lattice_constant_c
    create_unit_cell(a, b, c, uuid)
    print('...done!')

    config_path = os.path.join(material.run_id, 'config.yaml')
    config = load_config_file(config_path)

    epsilon_limits = config['epsilon_limits']

    chemical_species = []
    print('creating blender-materials for all atom-types...')
    for atom_type in structure.lennard_jones:
        chem_type = {
                'chemical' : atom_type.chemical_id,
                'radius' : atom_type.sigma,
                'material' : MakeMaterial(
                    '{}_{}'.format(atom_type.chemical_id, uuid),
                    (1,
                        1 - (0.3 + atom_type.epsilon / epsilon_limits[1]),
                        1 - (0.3 + atom_type.epsilon / epsilon_limits[1])),
                    (1, 1, 1), 1)}
        chemical_species.append(chem_type)
    print('...done!')

    print('adding atom-sites...')
    site_count = 0
    for atom_site in structure.atom_sites:
        print('%s of at least %s' % (site_count, len(structure.atom_sites)))
        for species in chemical_species:
            if species['chemical'] == atom_site.chemical_id:
                radius = species['radius']
                mesh_material = species['material']
        x = atom_site.x_frac * a
        y = atom_site.y_frac * b
        z = atom_site.z_frac * c

        site_count = add_atom_site(
            x, y, z, radius, mesh_material, uuid, site_count) 
        site_count = add_periodic_segments(
            atom_site, chemical_species, mesh_material, uuid, a, b, c, site_count)
    print('...done!')

    print('applying wireframe-modifier to unit cell...')
    unit_cell_name = 'UnitCell_%s' % uuid
    bpy.data.objects[unit_cell_name].select = True
    bpy.context.scene.objects.active = bpy.data.objects[unit_cell_name]
    bpy.ops.object.modifier_add(type='WIREFRAME')
    print('...done!')

    pwd_dir = os.environ['PWD']
    blender_file = os.path.join(pwd_dir, '%s.blend' % uuid)
    bpy.ops.wm.save_mainfile(filepath = blender_file)

    print('Finished.')
