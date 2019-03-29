#!/usr/bin/env python

import subprocess
import psycopg2

seed = 262441180
args = "~/builds/blender/blender-2.79/blender-2.79b-linux-glibc219-x86_64/blender empty_scene.blend --background --python render_wrapper.py {} epsilon SINGLE side".format(seed)
subprocess.call(args, shell=True)
