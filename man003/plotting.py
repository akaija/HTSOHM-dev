# plotting functions

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker
import numpy as np
from sqlalchemy import text

from man003.processing import run_id, bins, get_min_max, post_processed_engine

threshold = 50

#font_path = "/usr/share/fonts/msttcore/times.ttf"
#prop = mpl.font_manager.FontProperties(fname=font_path)
font_name = "Arial"
mpl.rcParams.update({"font.family": font_name, "font.size": 8,
                     "mathtext.fontset": "custom",
                     "mathtext.rm": font_name,
                     "mathtext.it": "{}:italic".format(font_name),
                     "mathtext.bf": "{}:bold".format(font_name)})

def get_axis_label(x):
    if x == "average_vdw_radius":
        return "Avg. $\sigma$ value ($\AA$)"
    elif x == "average_well_depth":
        return "Avg. $\epsilon$ value (K)"
    elif x == "co2_selectivity":
        return "CO$_{2}$ selectivity (dim.)"
    elif x == "co2_adsorption_uptake":
        return "CO$_{2}$ uptake$_{ads}$ ($v_{STP}/v$)"
    elif x == "co2_working_capacity":
        return "CO$_{2}$ working capacity ($v_{STP}/v$)"
    elif x == "coulombic_interaction":
        return "Coulombic heat of adsorption (GJ/mol)"
    elif x == "heat_of_adsorption":
        return "Avg. heat of adsorption (GJ/mol)"
    elif x == "helium_void_fraction":
        return "Helium void fraction (frac.)"
    elif x == "number_density":
        return "Number density (atoms/$\AA^{3}$)"
    elif x == "regenerability":
        return "Regenerability (%)"
    elif x == "sorbent_selection_parameter":
        return "Sorbent selection parameter (dim.)"
    elif x == "vdw_interaction":
        return "van der Waals heat of adsorption (GJ/mol)"
    elif x == "volumetric_surface_area":
        return "Volumetric surface area (m$^2$/cm$^{3}$)"
    else:
        print("UNKNOWN PROPERTY TYPE")
        
def get_tick_labels(run_id, x, bins):
    min_, max_ = get_min_max(x)
    step = (max_ - min_) / bins
    return np.arange(min_, max_ + step, step)

def fmt(x, pos):
    a, b = "{:.1e}".format(x).split('e')
    b = int(b)
    return r"${} \times 10^{{{}}}$".format(a, b)

def annotate_plot(ax, x, y, run_id, n=None):
    # x/y-axis labels
    ax.set_xlabel(get_axis_label(x))
    ax.set_ylabel(get_axis_label(y))
    # x/y-axis ticks
    x_ticks = get_tick_labels(run_id, x, bins)
    ax.set_xticks(x_ticks[::2])
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    y_ticks = get_tick_labels(run_id, y, bins)
    ax.set_yticks(y_ticks[::2])
    # x/y limits
    x_min, x_max = get_min_max(x)
    x_step = (x_max - x_min) / bins
    #x_max = x_max +  x_step
    ax.set_xlim(x_min, x_max)
    y_min, y_max = get_min_max(y)
    y_step = (y_max - y_min) / bins
    #y_max = y_max + y_step
    ax.set_ylim(y_min, y_max)
    if n != None:
        ax.set_title(n)

def add_plot(fig, coord, run_id, bins, x, y, c, n=None, colorbar="on"):
    ax = fig.add_subplot(coord)
    annotate_plot(ax, x, y, run_id, n)
    x_ticks = get_tick_labels(run_id, x, bins)
    y_ticks = get_tick_labels(run_id, y, bins)
    dx, dy = [x_ticks[1] - x_ticks[0], y_ticks[1] - y_ticks[0]]
    #query data
    sql = text("""
    select {0}_bin, {1}_bin, avg({2}), count(*)
    from processed_data
    where run_id=:run_id
    group by {0}_bin, {1}_bin""".format(x, y, c))
    
    # gather min, max averaged bin-values
    rows = post_processed_engine.connect().execute(sql, run_id=run_id)
    min_color, max_color = 10 ** 6, 0
    min_count, max_count = 10 ** 6, 0
    for row in rows:
        if row[2] < min_color:
            min_color = row[2]
        if row[2] > max_color:
            max_color = row[2]
        if row[3] < min_count:
            min_count = row[3]
        if row[3] > max_count:
            max_count = row[3]
    
    # add patches for each bin
    colors, patches = [], []
    rows = post_processed_engine.connect().execute(sql, run_id=run_id)
    for row in rows:
        color = row[2]
        colors.append(color)
        facecolor = (color - min_color) / (max_color - min_color)
        if row[3] < threshold:
            alpha = row[3] / threshold
        else:
            alpha = 1.
        if x not in ["vdw_interaction", "coulombic_interaction", "heat_of_adsorption"]:
            x_pos = row[0] * dx
        else:
            x_pos = row[0] * -dx
        if y not in ["vdw_interaction", "coulombic_interaction", "heat_of_adsorption"]:
            y_pos = row[1] * dy
        else:
            y_pos = row[1] * -dy
        square = Rectangle((x_pos, y_pos), dx, dy)
        patches.append(square)
        #ax.add_patch(Rectangle((x_pos, y_pos), dx, dy, facecolor=cm.viridis(facecolor), alpha=alpha))
        ax.add_patch(Rectangle((x_pos, y_pos), dx, dy, facecolor=cm.seismic(facecolor), alpha=alpha))
    #p = PatchCollection(patches, cmap=cm.viridis)
    p = PatchCollection(patches, cmap=cm.seismic)
    p.set_array(np.array(colors))
    
    ax.grid()
    
    # add colorbar
    if colorbar != "off":
        width = 25
        cb = fig.colorbar(p, ax=ax, aspect=width)
        cb.set_label(get_axis_label(c))
        cb
    return ax

def add_plot_no_color(fig, coord, run_id, bins, x, y, n=None):
    ax = fig.add_subplot(coord)
    annotate_plot(ax, x, y, run_id, n)
    x_ticks = get_tick_labels(run_id, x, bins)
    y_ticks = get_tick_labels(run_id, y, bins)
    dx, dy = [x_ticks[1] - x_ticks[0], y_ticks[1] - y_ticks[0]]
    #query data
    sql = text("""
    select {0}_bin, {1}_bin, count(*)
    from processed_data
    where run_id=:run_id
    group by {0}_bin, {1}_bin""".format(x, y))
    
    # gather min, max averaged bin-values
    rows = post_processed_engine.connect().execute(sql, run_id=run_id)
    min_count, max_count = 10 ** 6, 0
    for row in rows:
        if row[2] < min_count:
            min_count = row[2]
        if row[2] > max_count:
            max_count = row[2]
    
    # add patches for each bin
    colors, patches = [], []
    rows = post_processed_engine.connect().execute(sql, run_id=run_id)
    for row in rows:
        if row[2] < threshold:
            alpha = row[2] / threshold
        else:
            alpha = 1.
        if x not in ["vdw_interaction", "coulombic_interaction", "heat_of_adsorption"]:
            x_pos = row[0] * dx
        else:
            x_pos = row[0] * -dx
        if y not in ["vdw_interaction", "coulombic_interaction", "heat_of_adsorption"]:
            y_pos = row[1] * dy
        else:
            y_pos = row[1] * -dy
        square = Rectangle((x_pos, y_pos), dx, dy)
        patches.append(square)
        ax.add_patch(Rectangle((x_pos, y_pos), dx, dy, facecolor='grey', alpha=alpha))
    #p = PatchCollection(patches, cmap=cm.viridis)
    p = PatchCollection(patches, cmap=cm.seismic)
    p.set_array(np.array(colors))
    ax.grid()
    return ax

# figures

r = "regenerability"
s = "co2_selectivity"
ss = "sorbent_selection_parameter"
wc = "co2_working_capacity"

vf = "helium_void_fraction"
sa = "volumetric_surface_area"
nd = "number_density"
sig = "average_vdw_radius"
eps = "average_well_depth"

vdw = "vdw_interaction"
cou = "coulombic_interaction"

def figure_two():
    fig = plt.figure(figsize=(7 , 6.222), facecolor='white')
    add_plot(fig, '221', run_id, bins, sig, s, sa, "A")
    add_plot(fig, '222', run_id, bins, sig, wc, sa, "B")
    add_plot(fig, '223', run_id, bins, sig, r, sa, "C")
    add_plot(fig, '224', run_id, bins, sig, ss, sa, "D")
    plt.tight_layout()
    plt.savefig("man003/figures/fig2.png", dpi=600)

def figure_three():
    fig = plt.figure(figsize=(7 , 6.222), facecolor='white')
    add_plot(fig, '221', run_id, bins, eps, s, sa, "A")
    add_plot(fig, '222', run_id, bins, eps, wc, sa, "B")
    add_plot(fig, '223', run_id, bins, eps, r, sa, "C")
    add_plot(fig, '224', run_id, bins, eps, ss, sa, "D")
    plt.tight_layout()
    plt.savefig("man003/figures/fig3.png", dpi=600)

def figure_four():
    fig = plt.figure(figsize=(7, 2.25), facecolor='white')
    add_plot(fig, '131', run_id, bins, vf, s, r, "A")
    add_plot(fig, '132', run_id, bins, sa, s, r, "B")
    add_plot(fig, '133', run_id, bins, nd, s, r, "C")
    plt.tight_layout()
    plt.savefig("man003/figures/fig4.png", dpi=600)

def figure_five():
    fig = plt.figure(figsize=(7, 2.25), facecolor='white')
    add_plot(fig, '131', run_id, bins, vf, wc, r, "A")
    add_plot(fig, '132', run_id, bins, sa, wc, r, "B")
    add_plot(fig, '133', run_id, bins, nd, wc, r, "C")
    plt.tight_layout()
    plt.savefig("man003/figures/fig5.png", dpi=600)

def figure_six():
    fig = plt.figure(figsize=(7, 3.111), facecolor='white')
    add_plot(fig, '121', run_id, bins, vdw, eps, nd, "A")
    add_plot(fig, '122', run_id, bins, sig, sa, nd, "B")
    plt.tight_layout()
    plt.savefig("man003/figures/fig6.png", dpi=600)

def figure_seven():
    fig = plt.figure(figsize=(7, 6.222), facecolor='white')
    add_plot(fig, '221', run_id, bins, vdw, s, vf, "A")
    add_plot(fig, '222', run_id, bins, vdw, wc, vf, "B")
    add_plot(fig, '223', run_id, bins, vdw, r, vf, "C")
    add_plot(fig, '224', run_id, bins, vdw, ss, vf, "D")
    plt.tight_layout()
    plt.savefig("man003/figures/fig7.png", dpi=600)

def figure_eight():
    fig = plt.figure(figsize=(7, 6.222), facecolor='white')
    add_plot(fig, '221', run_id, bins, cou, s, vf, "A")
    add_plot(fig, '222', run_id, bins, cou, wc, vf, "B")
    add_plot(fig, '223', run_id, bins, cou, r, vf, "C")
    add_plot(fig, '224', run_id, bins, cou, ss, vf, "D")
    plt.tight_layout()
    plt.savefig("man003/figures/fig8.png", dpi=600)

from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from os import listdir
import matplotlib.image as mpimg

def figure_nine():
    render_files = listdir("man003/renders/")

    for color_by in ["epsilon", "charge"]:
    
        fig = plt.figure(figsize=(7, 5), facecolor='white')
        ax = add_plot_no_color(fig, '231', run_id, bins, s, wc)
        ax.text(65, 5, "A", bbox=dict(facecolor='black', alpha=0.5), color='white')
        ax.text(325, 15, "B", bbox=dict(facecolor='black', alpha=0.5), color='white')
        ax.text(325, 55, "C", bbox=dict(facecolor='black', alpha=0.5), color='white')
        ax.text(325, 95, "D", bbox=dict(facecolor='black', alpha=0.5), color='white')
        ax.text(585, 85, "E", bbox=dict(facecolor='black', alpha=0.5), color='white')
    
        
        index = 2
        for r in ["A", "B", "C", "D", "E"]:
            for render_file in render_files:
                if r in render_file and color_by in render_file and "SI" not in render_file:
                    selected_file = render_file
            
            ax = fig.add_subplot("23{}".format(index))
            img = mpimg.imread("man003/renders/{}".format(selected_file))
            ax.imshow(img)
            ax.set_title(r)
            ax.xaxis.set_ticks_position('none')
            ax.yaxis.set_ticks_position('none')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            index += 1

        if color_by == "epsilon":
            min_, max_ = 1.258, 513.2
            label = "$\epsilon$ value (K)"
        else:
            min_, max_ = -1, 1
            label = "Partial charge"
        norm = Normalize(vmin=min_, vmax=max_)
#        fig.subplots_adjust(left=0.05, right=0.9, hspace=0.3)
        fig.subplots_adjust(left=0.1, right=0.85, hspace=0.45)
#        cb_ax = fig.add_axes([0.92, 0.25/2, 0.02, 0.75])
        cb_ax = fig.add_axes([0.87, 0.25/2, 0.02, 0.75])
        #cbar = ColorbarBase(cb_ax, cmap=cm.viridis, norm=norm, orientation='vertical')
        cbar = ColorbarBase(cb_ax, cmap=cm.seismic, norm=norm, orientation='vertical')
        cbar.set_label(label)

        plt.savefig("man003/figures/fig9_{}.png".format(color_by), dpi=600)

import matplotlib.gridspec as gridspec

def si_renders():
    render_files = listdir("man003/renders/")

    fig = plt.figure(figsize=(3.5, 3.5), facecolor='white')
    ax = add_plot_no_color(fig, '111', run_id, bins, s, wc)
    ax.text(65, 5, "A", bbox=dict(facecolor='black', alpha=0.5), color='white')
    ax.text(325, 15, "B", bbox=dict(facecolor='black', alpha=0.5), color='white')
    ax.text(325, 55, "C", bbox=dict(facecolor='black', alpha=0.5), color='white')
    ax.text(325, 95, "D", bbox=dict(facecolor='black', alpha=0.5), color='white')
    ax.text(585, 85, "E", bbox=dict(facecolor='black', alpha=0.5), color='white')
    plt.tight_layout()
    plt.savefig("man003/figures/fig_s2.png", dpi=600)

    gs = gridspec.GridSpec(1, 3, hspace=0, wspace=0.1, left=0.05, right=0.85)
    cb = gridspec.GridSpec(1, 1, left=0.875, right=0.9)

    fig_num = 3
    for q in ["A", "B", "C", "D", "E"]:
        for color_by in ["epsilon", "charge"]:
            fig = plt.figure(figsize=(7, 3), facecolor="white")
            index = 0
            for render_file in render_files:
                if q in render_file and color_by in render_file and "SI" in render_file:
                    print(render_file)
                    #ax = fig.add_subplot("13{}".format(index))
                    print(index)
                    ax = fig.add_subplot(gs[0, index])
                    img = mpimg.imread("man003/renders/{}".format(render_file))
                    ax.imshow(img)
                    if index == 0:
                        title = "A"
                    elif index == 1:
                        title = "B"
                    elif index == 2:
                        title = "C"
                    ax.set_title(title)
                    ax.xaxis.set_ticks_position('none')
                    ax.yaxis.set_ticks_position('none')
                    ax.set_xticklabels([])
                    ax.set_yticklabels([])
                    if color_by == "epsilon":
                        min_, max_ = 1.258, 513.2
                        label = "$\epsilon$ value (K)"
                    else:
                        min_, max_ = -1, 1
                        label = "Partial charge"
                    norm = Normalize(vmin=min_, vmax=max_)
                    #fig.subplots_adjust(left=0.1, right=0.85, hspace=0.45)
                    #cb_ax = fig.add_axes([0.87, 0.25/2, 0.02, 0.75])
                    cb_ax = fig.add_subplot(cb[0, 0])
                    cbar = ColorbarBase(cb_ax, cmap=cm.seismic, norm=norm, orientation='vertical')
                    #cbar = ColorbarBase(cb_ax, cmap=cm.viridis, norm=norm, orientation='vertical')
                    cbar.set_label(label)
                    index += 1

            plt.savefig("man003/figures/fig_s{}_{}_{}.png".format(fig_num, q, color_by))
            fig_num += 1
#            ax = fig.add_subplot("23{}".format(index))
#            img = mpimg.imread("man003/renders/{}".format(selected_file))
#            ax.imshow(img)
#            ax.set_title(r)
#            ax.xaxis.set_ticks_position('none')
#            ax.yaxis.set_ticks_position('none')
#            ax.set_xticklabels([])
#            ax.set_yticklabels([])
#            index += 1
#
#        if color_by == "epsilon":
#            min_, max_ = 1.258, 513.2
#            label = "$\epsilon$ value (K)"
#        else:
#            min_, max_ = -1, 1
#            label = "Partial charge"
#        norm = Normalize(vmin=min_, vmax=max_)
##        fig.subplots_adjust(left=0.05, right=0.9, hspace=0.3)
#        fig.subplots_adjust(left=0.1, right=0.85, hspace=0.45)
##        cb_ax = fig.add_axes([0.92, 0.25/2, 0.02, 0.75])
#        cb_ax = fig.add_axes([0.87, 0.25/2, 0.02, 0.75])
#        cbar = ColorbarBase(cb_ax, cmap=cm.viridis, norm=norm, orientation='vertical')
#        cbar.set_label(label)
#
#        plt.savefig("man003/figures/fig9_{}.png".format(color_by), dpi=600)


def figure_ten():
    fig = plt.figure(figsize=(7, 2.25), facecolor='white')
    add_plot(fig, '131', run_id, bins, s, wc, ss, "A")
    add_plot(fig, '132', run_id, bins, s, r, ss, "B")
    add_plot(fig, '133', run_id, bins, wc, r, ss, "C")
    plt.tight_layout()
    plt.savefig("man003/figures/fig10.png", dpi=600)

import os

import matplotlib.gridspec as gridspec

def no_labels(ax, image=True):
    if image == True:
        ax.xaxis.set_ticks_position("none")
        ax.yaxis.set_ticks_position("none")
        ax.axis("off")
    ax.set_xticklabels([])
    ax.set_yticklabels([])

def add_render(fig, position, filepath):
    ax = fig.add_subplot(position)
    no_labels(ax)
    img = mpimg.imread(filepath)
    ax.imshow(img)

def table_of_contents():
    fig = plt.figure(figsize=(1.75, 1.75), facecolor="white")
    gs = gridspec.GridSpec(2, 2, hspace=0, wspace=0, left=0, bottom=0, top=1, right=1)
    ax = add_plot(fig, gs[0,1], run_id, bins, s, wc, r, colorbar="off")
    no_labels(ax, image=False)
    ax = add_plot(fig, gs[1,0], run_id, bins, sa, wc, r, colorbar="off")
    no_labels(ax, image=False)
    add_render(fig, gs[0,0], "man003/renders/C_250982367_epsilon.png")
    add_render(fig, gs[1,1], "man003/renders/D_683739552_epsilon.png")
    plt.savefig("man003/figures/ToC.png", dpi=600)

def plot_all_figures():
    figure_two()
    figure_three()
    figure_four()
    figure_five()
    figure_six()
    figure_seven()
    figure_eight()
    figure_nine()
    figure_ten()
    table_of_contents()
    print("Completed plotting.")
