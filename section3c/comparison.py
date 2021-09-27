import pandas as pd
import numpy as np
import scipy.stats as stats
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Define files
file_best = "rsMultipleInput_func_06000_00000_04000_merged.csv"
file_worst = "rsMultipleInput_func_01000_08000_01000_merged.csv"

# Proccess data

vol_full = []
vol_merged = []
vol_merged_rel_error = []

data_best = pd.read_csv(file_best, header = 0, index_col = None, dtype = np.float64)
data_worst = pd.read_csv(file_worst, header = 0, index_col = None, dtype = np.float64)

vol_full_local = []
vol_merged_local = []

data_best.loc[:,"pressure"] = (1. / 133300.) * (data_best.loc[:,"eqResistance"] * data_best.loc[:,"flow"]) / ((data_best.loc[:,"radius"]) ** 4)
data_best.loc[:,"pressure"] = data_best.loc[:,"pressure"] - data_best.loc[:,"pressure"].max()
vol_full_local.append(np.array(data_best.loc[:,"length"] * data_best.loc[:,"radius"] * data_best.loc[:,"radius"] * np.pi).sum())

data_worst.loc[:,"pressure"] = (1. / 133300.) * (data_worst.loc[:,"eqResistance"] * data_worst.loc[:,"flow"]) / ((data_worst.loc[:,"radius"]) ** 4)
data_worst.loc[:,"pressure"] = data_worst.loc[:,"pressure"] - data_worst.loc[:,"pressure"].max()
vol_merged_local.append(np.array(data_worst.loc[:,"length"] * data_worst.loc[:,"radius"] * data_worst.loc[:,"radius"] * np.pi).sum())

# Volumes
vol_full_local_np = np.array(vol_full_local)
vol_merged_local_np = np.array(vol_merged_local)
vol_full.append(vol_full_local_np)
vol_merged.append(vol_merged_local_np)
vol_merged_rel = np.divide(np.absolute((vol_merged_local_np - np.mean(vol_full_local_np))), (np.mean(vol_full_local_np)))
vol_merged_rel_error.append(vol_merged_rel)


data_best_flat = data_best.loc[:,["level", "flow", "pressure", "radius", "length"]]
data_best_np = data_best_flat.to_numpy()
    
data_worst_flat = data_worst.loc[:,["level", "flow", "pressure", "radius", "length"]]
data_worst_np = data_worst_flat.to_numpy()

idxLevel = 0
idxFlow = 1
idxPressure = 2
idxRadius = 3
idxLength = 4

# Define plot functions

my_dpi = 300
my_dark_blue_rgb = [0, 0, 154./255]
my_dark_green_rgb = [0, 1, 180./255]

# linewitdh = 1.5
# #fontsize = "xx-large"
# fontsize = 20
# fontweight = "bold"
# fig_size = 10

linewitdh = 0.5
#fontsize = "xx-large"
fontsize = 12
fontweight = "bold"
cm = 1./2.54
fig_width=7.5*cm*3
fig_height=5.14*cm*3

def comparison_volumes_error_boxplot(filename, y, labels, ylog):
    fig, axs = plt.subplots()
    if (ylog):
        axs.set_yscale("log")
#    axs.set_title("Volume", fontsize = fontsize)
    axs.set_ylabel("Erro relativo de volume", fontsize = fontsize)
    axs.set_xlabel("$N_\mathrm{base}$", fontsize = fontsize)
    axs.boxplot(y, whis=(0, 100), labels=labelsAr)
    # axs.tick_params(axis="y", which="major", length=8, width=2, labelsize=fontsize)
    # axs.tick_params(axis="x", which="major", length=8, width=2, labelsize=fontsize,labelrotation=45)
    axs.tick_params(axis="y", which="major", labelsize=fontsize)
    axs.tick_params(axis="x", which="major", labelsize=fontsize,labelrotation=45)
    for axis in ['top','bottom','left','right']:
        axs.spines[axis].set_linewidth(linewitdh)
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_height)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

def comparative_profile(filename, y, xlabel, ylabel, ylog):
    lvl_size = max(len(y[0]), len(y[1]))
    fig, axs = plt.subplots()
    if(ylog):
        axs.set_yscale("log")
    dcco_style = dict(color="xkcd:orange", linestyle="-")
    pdcco_style = dict(color="xkcd:blue", linestyle="-")
    dcco_median_style = dict(color="xkcd:orange", linestyle="-")
    pdcco_median_style = dict(color="xkcd:blue", linestyle="-")
    dict_1 = axs.boxplot(y[0], whis=(0, 100), positions=range(2, (2 * lvl_size) + 1, 2),boxprops=dcco_style, whiskerprops=dcco_style, capprops=dcco_style, medianprops=dcco_median_style)
    # axs.set_title(title, fontsize=fontsize)
    axs.set_ylabel(ylabel, fontsize=fontsize)
    axs.set_xlabel(xlabel, fontsize=fontsize)
    dict_2 = axs.boxplot(y[1], whis=(0, 100), positions=range(1, (2 * lvl_size) + 1, 2), boxprops=pdcco_style, whiskerprops=pdcco_style, capprops=pdcco_style, medianprops=pdcco_median_style)
    axs.legend((dict_2["whiskers"][0], dict_1["whiskers"][0]), ("Best", "Worst"), fontsize=fontsize, loc='upper right')
    axs.xaxis.set_major_locator(ticker.FixedLocator(np.arange(1.5, (2 * lvl_size) + 1, 2)))
    axs.xaxis.set_major_formatter(ticker.FixedFormatter(range(0, lvl_size + 1, 1)))
    axs.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(2.5, (2 * lvl_size) + 1, 2)))
    axs.tick_params(axis="y", which="major", labelsize=fontsize)
    axs.tick_params(axis="x", which="major", labelsize=fontsize)
    axs.tick_params(axis="x", which="minor", grid_linestyle=":", length=0, width=0)
    axs.grid(True, axis="x", which="minor")
    fig.set_dpi(300)
    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)


# Plots

best_max = np.amax(data_best_np[:,0])
worst_max = np.amax(data_worst_np[:,0])
real_max = max(best_max, worst_max)
steps = np.arange(0, real_max + 1, 1, dtype=int)
size = np.size(steps)
bx_labels = steps

worst = []
best = []
for j in range(size):
    best_part = data_best_np[(data_best_np[:,0] == steps[j])]
    best.append(best_part)
    worst_part = data_worst_np[(data_worst_np[:,0] == steps[j])]
    worst.append(worst_part)

filename_base = "rsMultipleInput_"

radius = [[data[:,3] for data in worst], [data[:,3] for data in best]]
length = [[data[:,4] for data in worst], [data[:,4] for data in best]]
pressure = [[data[:,2] for data in worst], [data[:,2] for data in best]]
flow = [[data[:,1] for data in worst], [data[:,1] for data in best]]

filename = filename_base + "radius_comparison.pdf"
comparative_profile(filename, radius, "Nível", "Raio [cm]", True)

filename = filename_base + "pressure_comparison.pdf"
comparative_profile(filename, pressure, "Nível", "Pressão [mmHg]", False)

# filename = filename_base + "radius_log.pdf"
# profile(filename, radius, "", "Level", "Radius [cm]", bx_labels, True)

# filename = filename_base + "pressure.pdf"
# profile(filename, pressure, "", "Level", "Pressure [mmHg]", bx_labels, False)

# plot_type = "_s"
# filename = prefix + labels[i] + plot_type + "_rad_x_len" + ".png"
# comparison_scatter_1x2(filename, [data_worst_np[:, 3], data_best_np[:, 3]], [data_worst_np[:, 4], data_best_np[:, 4]], labelsArLatex[i],  "Raio [cm]", "Comprimento [cm]", True, True)

