import os
from math import pi

import bpy
import numpy as np
import yaml
import matplotlib.cm as cm

import htsohm
from htsohm.db import session, Material
from htsohm.files import load_config_file
from htsohm.pseudomaterial_generator.random import generate_material

scene       = bpy.context.scene
add_cube    = bpy.ops.mesh.primitive_cube_add
add_sphere  = bpy.ops.mesh.primitive_uv_sphere_add

scene.render.resolution_x = 1080
scene.render.resolution_y = 1080
scene.world.light_settings.use_environment_light = True
scene.world.horizon_color = (1, 1, 1)

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

def create_unit_cell(a, b, c, seed):
    dim = np.array([a, b, c]) / 2.
    add_cube(location = (dim[0] + 8., dim[1], dim[2]))
    bpy.ops.transform.resize(value = dim)
    bpy.data.objects["Cube"].select = True
    cell_name = "UnitCell_{}".format(seed)
    bpy.context.scene.objects[0].name = cell_name
    black = MakeMaterial('black', (0, 0, 0), (1, 1, 1), 1)
    SetMaterial(bpy.context.object, black)
    bpy.data.objects[cell_name].select = False

def add_atom_site(x, y, z, radius,  material, seed, site_count):
    add_sphere(location = (x, y, z), size = radius)
    bpy.data.objects['Sphere'].select = True
    bpy.ops.object.shade_smooth()
    SetMaterial(bpy.context.object, material)
    site_name = 'AtomSite_%s' % site_count
    bpy.context.scene.objects[0].name = site_name
    bpy.data.objects[site_name].select = False
    site_count += 1
    return site_count

def render_from_seed(seed, color_by, tag):
    print("seed   :\t{}".format(seed))
    seed = int(seed)

    # query run_id
    run_id = session.query(Material.run_id).filter(Material.seed==seed).one()[0]
    print("run_id :\t{}".format(run_id))

    # load config
    config_path = os.path.join(run_id, "config.yaml")
    config = load_config_file(config_path)["structure_parameters"]

    # generate structure
    m, s = generate_material(run_id, seed, config)
    
    # create unit cell
    print('\ncreating unit cell...')
    a, b, c = s.lattice_constants.a, s.lattice_constants.b, s.lattice_constants.c
    create_unit_cell(a, b, c, seed)
    print('...done!')

    ep_lims = config['epsilon_limits']

    chemical_species = []
    print('creating blender-materials for all atom-types...')
    for at in s.atom_types:
        val = (at.epsilon - ep_lims[0]) / (ep_lims[1] - ep_lims[0])
        color = cm.seismic(val)[:3]
        chem_type = {
                'chemical' : at.chemical_id,
                'radius' : at.sigma,
                'material' : MakeMaterial(
                    '{}_{}'.format(at.chemical_id, seed),
                    color,
                    (1, 1, 1), 1)}
        chemical_species.append(chem_type)
    print('...done!')

    print('adding atom-sites...')
    site_count = 0
    for atom_site in s.atom_sites:
#        print("{} of at least {}".format(site_count, len(s.atom_sites)))
        for species in chemical_species:
            if species['chemical'] == atom_site.chemical_id:
                radius = species['radius']
                if color_by == "epsilon":
                    mesh_material = species['material']
                elif color_by == "charge":
                    val = (atom_site.q + 1) / 2
                    color = cm.seismic(val)[:3]
                    material = MakeMaterial("site_{}".format(site_count), color, (1,1,1) ,1)
                    mesh_material = material
                else:
                    print("`color_by` not found...")
        x = atom_site.x * a
        y = atom_site.y * b
        z = atom_site.z * c

        site_count = add_atom_site(x, y, z, radius, mesh_material, seed, site_count) 
#        site_count = add_atom_site(x, y + b, z, radius, mesh_material, seed, site_count) 
#        site_count = add_atom_site(x, y - b, z, radius, mesh_material, seed, site_count) 
#        site_count = add_atom_site(x, y + b, z + c, radius, mesh_material, seed, site_count) 
#        site_count = add_atom_site(x, y + b, z - c, radius, mesh_material, seed, site_count) 
#        site_count = add_atom_site(x, y - b, z + c, radius, mesh_material, seed, site_count) 
#        site_count = add_atom_site(x, y - b, z - c, radius, mesh_material, seed, site_count) 
#        site_count = add_atom_site(x, y, z + c, radius, mesh_material, seed, site_count) 
#        site_count = add_atom_site(x, y, z - c, radius, mesh_material, seed, site_count) 

    print('...done!')

    print('applying wireframe-modifier to unit cell...')
    unit_cell_name = "UnitCell_{}".format(seed)
    bpy.data.objects[unit_cell_name].select = True
    bpy.context.scene.objects.active = bpy.data.objects[unit_cell_name]
    bpy.ops.object.modifier_add(type='WIREFRAME')
    print('...done!')

    # add camera
    bpy.ops.object.camera_add(location=(a + 10., b / 2., c / 2.), rotation=(pi / 2, 0., pi / 2))
    scene.camera = bpy.data.objects["Camera"]
    bpy.data.cameras["Camera"].type = "ORTHO"
    bpy.data.cameras["Camera"].ortho_scale = 64.0
    bpy.data.cameras["Camera"].show_mist = True

    pwd_dir = os.environ['PWD']
#    blender_file = os.path.join(pwd_dir, "man003/renders/{}_{}.blend".format(seed, color_by))
#    bpy.ops.wm.save_mainfile(filepath = blender_file)

    png_file = os.path.join(pwd_dir, "man003/renders/{}_{}_{}.png".format(tag, seed, color_by))
    bpy.data.scenes['Scene'].render.filepath = png_file
    bpy.ops.render.render( write_still=True )

    print('Finished.')
