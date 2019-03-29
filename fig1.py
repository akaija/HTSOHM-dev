#!/usr/bin/env python

import os

import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

def add_image(ax, path, name):
    img = mpimg.imread(path)
    #ax.imshow(img, cmap=plt.cm.viridis)
    ax.imshow(img, cmap=plt.cm.seismic)
    if name != None:
        ax.set_title(name, loc="left")
    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.axis("off")

lt = os.path.join("man003", "renders", "NU-125_side.png")
lb = os.path.join("man003", "renders", "NU-125_angle.png")
rt = os.path.join("man003", "renders", "pseudomaterial_side.png")
rb = os.path.join("man003", "renders", "pseudomaterial_angle.png")

font = "arial"
mpl.rcParams.update({"font.family": font,
                     "font.size": 8,
                     "mathtext.fontset": "custom",
                     "mathtext.rm": font,
                     "mathtext.it": font})

fig = plt.figure(figsize=(3.5, 3.5))

gs = gridspec.GridSpec(2, 2, hspace=0, wspace=0, left=0, right=0.8)
cb = gridspec.GridSpec(1, 1, left=0.825, right=0.85)

ax = fig.add_subplot(gs[0,0])
add_image(ax, lt, "    A")
ax = fig.add_subplot(gs[1,0])
add_image(ax, lb, "    B")
ax = fig.add_subplot(gs[0,1])
add_image(ax, rt, "    C")
ax.annotate("10 $\AA$", xy=(0.5, 0.05), xycoords="axes fraction", xytext=(0.5 + 1/6, 0.05), textcoords="axes fraction",
        horizontalalignment="left", verticalalignment="center",
        arrowprops=dict(arrowstyle="-", facecolor="black"))
ax = fig.add_subplot(gs[1,1])
add_image(ax, rb, "    D")

norm = Normalize(vmin=0, vmax=50)
#fig.subplots_adjust(right=0.85, hspace=0.45)
#cb_ax = fig.add_axes([0.87, 0.25/2, 0.02, 0.75])
cb_ax = fig.add_subplot(cb[0,0])
#cbar = ColorbarBase(cb_ax, cmap=cm.viridis, norm=norm, orientation="vertical")
cbar = ColorbarBase(cb_ax, cmap=cm.Reds, norm=norm, orientation="vertical")
cbar.set_label("Potential well depth (K)")

plt.tight_layout()
plt.savefig("fig1.png", dpi=600)
