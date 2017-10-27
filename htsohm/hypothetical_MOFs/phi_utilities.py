import os

import numpy as np
from math import pi, sqrt, cos
import yaml

import non_pseudo
from non_pseudo.db import session, Material

def load_data(uuid):
    np_path = os.path.dirname(os.path.dirname(non_pseudo.__file__))
    cif_path = os.path.join('cif_files', '{}.cif'.format(uuid))
    
    if 'hypothetical' in cif_path:
        [a, b, c] = np.genfromtxt(cif_path, usecols=1, skip_header=8, max_rows=3)
        [alpha, beta, gamma] = np.genfromtxt(cif_path, usecols=1, skip_header=11, max_rows=3)
        xfrac = np.genfromtxt(cif_path, usecols=2, skip_header=20)
        yfrac = np.genfromtxt(cif_path, usecols=3, skip_header=20)
        zfrac = np.genfromtxt(cif_path, usecols=4, skip_header=20)
    else:
        [a, b, c] = np.genfromtxt(cif_path, usecols=1, skip_header=9, max_rows=3)
        [alpha, beta, gamma] = np.genfromtxt(cif_path, usecols=1, skip_header=12, max_rows=3)
        xfrac = np.genfromtxt(cif_path, usecols=2, skip_header=24)
        yfrac = np.genfromtxt(cif_path, usecols=3, skip_header=24)
        zfrac = np.genfromtxt(cif_path, usecols=4, skip_header=24)
    
    return xfrac, yfrac, zfrac, a, b, c, alpha, beta, gamma

def distance(x0, y0, z0, x1, y1, z1):
    return sqrt((x1 - x0)**2 + (y1 - y0)**2 + (z1 - z0)**2)

def rad_dist_func(xfrac, yfrac, zfrac, a, b, c, alpha, beta, gamma, n_bins, side):
    x = [e * a for e in xfrac]
    y = [e * b for e in yfrac]
    z = [e * c for e in zfrac]
    num_part = len(x)
    [A, B, C] = [e * pi / 180. for e in [alpha, beta, gamma]]
    vol = a*b*c*sqrt(1 + 2*cos(A)*cos(B)*cos(C)-cos(A)**2-cos(B)**2-cos(C)**2)
    nden = num_part / vol
    dr = side / n_bins
    rs = np.arange(dr, side + dr, dr)
    counts = np.zeros([len(rs)])
    index = 0
    for r in rs:
        for i in range(num_part):
            for j in range(num_part):
                if i != j:
                    d = distance(x[i], y[i], z[i], x[j], y[j], z[j])
                    if d >= r and d < r + dr:
                        counts[index] += 2
                        
        counts[index] /= num_part
        counts[index] /= 4 * pi * r**2 * dr
        counts[index] /= nden
        index += 1
        
    return (rs, counts)
