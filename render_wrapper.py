#! /usr/bin/env python

import sys
import os
import bpy

blend_dir = os.path.basename(bpy.data.filepath)
if blend_dir not in sys.path:
    sys.path.append(blend_dir)

venv_dir = os.environ['VIRTUAL_ENV']
lib_dir = os.path.join(venv_dir, 'lib', 'python3.7', 'site-packages')
if lib_dir not in sys.path:
    sys.path.append(lib_dir)
lib64_dir = os.path.join(venv_dir, 'lib64', 'python3.7', 'site-packages')
if lib64_dir not in sys.path:
    sys.path.append(lib64_dir)

pwd_dir = os.environ['PWD']
if pwd_dir not in sys.path:
    sys.path.append(pwd_dir)

import render
import imp

imp.reload(render)

if __name__ == "__main__":
    render.render_from_seed(str(sys.argv[5]), str(sys.argv[6]), str(sys.argv[7]), str(sys.argv[8]))
