import pandas as pd
import numpy as np
import scipy.stats as stats
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Define files

p = Path('.')

list_d025_full = list(p.glob("CSV/*d025*full*.csv"))
list_d025_merged = list(p.glob("CSV/*d025*merged*.csv"))
files_d025 = [list_d025_full, list_d025_merged]

list_d050_full = list(p.glob("CSV/*d050*full*.csv"))
list_d050_merged = list(p.glob("CSV/*d050*merged*.csv"))
files_d050 = [list_d050_full, list_d050_merged]

list_d075_full = list(p.glob("CSV/*d075*full*.csv"))
list_d075_merged = list(p.glob("CSV/*d075*merged*.csv"))
files_d075 = [list_d075_full, list_d075_merged]

list_d100_full = list(p.glob("CSV/*d100*full*.csv"))
list_d100_merged = list(p.glob("CSV/*d100*merged*.csv"))
files_d100 = [list_d100_full, list_d100_merged]

list_d150_full = list(p.glob("CSV/*d150*full*.csv"))
list_d150_merged = list(p.glob("CSV/*d150*merged*.csv"))
files_d150 = [list_d150_full, list_d150_merged]

list_d200_full = list(p.glob("CSV/*d200*full*.csv"))
list_d200_merged = list(p.glob("CSV/*d200*merged*.csv"))
files_d200 = [list_d200_full, list_d200_merged]

files = [files_d025, files_d050, files_d075, files_d100, files_d150, files_d200]
labels = ["d025", "d050", "d075", "d100", "d150", "d200"]
labelsAr = ["0.25", "0.50", "0.75", "1.00", "1.50", "2.00"]

# Proccess data

full = []
merged = []
vol_full = []
vol_merged = []
vol_merged_rel_error = []

for i in range(len(files)):
    
    data_full = [pd.read_csv(listed, header = 0, index_col = None, dtype = np.float64) for listed in files[i][0]]
    data_merged = [pd.read_csv(listed, header = 0, index_col = None, dtype = np.float64) for listed in files[i][1]]
 
    vol_full_local = []
    vol_merged_local = []

    for table in data_full:
        table.loc[:,"pressure"] = (1. / 133300.) * (table.loc[:,"eqResistance"] * table.loc[:,"flow"]) / ((table.loc[:,"radius"]) ** 4)
        table.loc[:,"pressure"] = table.loc[:,"pressure"] - table.loc[:,"pressure"].max()
        vol_full_local.append(np.array(table.loc[:,"length"] * table.loc[:,"radius"] * table.loc[:,"radius"] * np.pi).sum())
        
    for table in data_merged:
        table.loc[:,"pressure"] = (1. / 133300.) * (table.loc[:,"eqResistance"] * table.loc[:,"flow"]) / ((table.loc[:,"radius"]) ** 4)
        table.loc[:,"pressure"] = table.loc[:,"pressure"] - table.loc[:,"pressure"].max()
        vol_merged_local.append(np.array(table.loc[:,"length"] * table.loc[:,"radius"] * table.loc[:,"radius"] * np.pi).sum())

    # Volumes
    vol_full_local_np = np.array(vol_full_local)
    vol_merged_local_np = np.array(vol_merged_local)
    vol_full.append(vol_full_local_np)
    vol_merged.append(vol_merged_local_np)
    vol_merged_rel = np.divide(np.absolute((vol_merged_local_np - np.mean(vol_full_local_np))), (np.mean(vol_full_local_np)))
    vol_merged_rel_error.append(vol_merged_rel)
    
    
    data_full_flat = (pd.concat(data_full)).loc[:,["level", "flow", "pressure", "radius", "length"]]
    data_full_np = data_full_flat.to_numpy()
     
    data_merged_flat = (pd.concat(data_merged)).loc[:,["level", "flow", "pressure", "radius", "length"]]
    data_merged_np = data_merged_flat.to_numpy()
 
    full.append(data_full_np)
    
    merged.append(data_merged_np)   

# Define plot functions

my_dpi = 300
my_dark_blue_rgb = [0, 0, 154./255]
my_dark_green_rgb = [0, 1, 180./255]

linewitdh = 1.5
#fontsize = "xx-large"
fontsize = 20
fontweight = "bold"
fig_size = 10

def profile(filename, y, title, xlabel, ylabel, labels, ylog):
    fig, axs = plt.subplots(2, 1, sharey='all', sharex='all', squeeze=False)
    if (ylog):
        axs[0,0].set_yscale("log")
    axs[0,0].boxplot(y[0], whis=(0, 100), labels=labels)
    axs[0,0].set_title(title + " - DCCO", fontsize = fontsize)
    axs[0,0].set_ylabel(ylabel, fontsize = fontsize)
    axs[0,0].set_xlabel(xlabel, fontsize = fontsize)
    axs[0,0].tick_params(axis="both", which="major", length=8, width=2, labelsize=fontsize)
    axs[1,0].set_title("PDCCO", fontsize = fontsize)
    axs[1,0].boxplot(y[1], whis=(0, 100), labels=labels)
    axs[1,0].set_ylabel(ylabel, fontsize = fontsize)
    axs[1,0].set_xlabel(xlabel, fontsize = fontsize)
    axs[1,0].tick_params(axis="both", which="major", length=8, width=2, labelsize=fontsize)
    for axis in ['top','bottom','left','right']:
        axs[0, 0].spines[axis].set_linewidth(linewitdh)
        axs[1, 0].spines[axis].set_linewidth(linewitdh)
    fig.set_dpi(my_dpi)
    fig.set_figheight(1 * fig_size)
    fig.set_figwidth(1.5 * fig_size)
    fig.tight_layout()    
    fig.savefig(filename)
    plt.close(fig)

def comparison_scatter_1x2(filename, x, y, title, xlabel, ylabel, xlog, ylog):
    fig, axs = plt.subplots(1, 2, sharey='all', sharex='row', squeeze=False)
    if (xlog):
        axs[0, 0].set_xscale("log")
    if (ylog):
        axs[0, 0].set_yscale("log")
    axs[0, 0].scatter(x[0], y[0], c=[my_dark_green_rgb], alpha=0.1, linewidths=0)
    axs[0, 0].set_xlabel(xlabel, fontsize = 40)
    axs[0, 0].set_ylabel(ylabel, fontsize = 40)
#    axs[0, 0].set_title("DCCO", fontsize = 40)
    axs[0,0].tick_params(axis="both", which="major", length=8, width=2, labelsize=40)
    axs[0,1].scatter(x[1], y[1], c=[my_dark_blue_rgb], alpha=0.1, linewidths=0)
    axs[0,1].set_xlabel(xlabel, fontsize = 40)
#    axs[0,1].set_ylabel(ylabel, fontsize = 40)
#    axs[0,1].set_title("PDCCO", fontsize = 40)
    axs[0,1].tick_params(axis="both", which="major", length=8, width=2, labelsize=40)
    for axis in ['top','bottom','left','right']:
        axs[0, 0].spines[axis].set_linewidth(linewitdh)
        axs[0, 1].spines[axis].set_linewidth(linewitdh)
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_size)
    fig.set_figwidth(2 * fig_size)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

def comparison_volumes_error_boxplot(filename, y, labels, ylog):
    fig, axs = plt.subplots()
    if (ylog):
        axs.set_yscale("log")
#    axs.set_title("Volume", fontsize = fontsize)
    axs.set_ylabel("Relative Volume Error", fontsize = fontsize)
    axs.set_xlabel("$A_r$", fontsize = fontsize)
    axs.boxplot(y, whis=(0, 100), labels=labelsAr)
    axs.tick_params(axis="both", which="major", length=8, width=2, labelsize=fontsize)
    for axis in ['top','bottom','left','right']:
        axs.spines[axis].set_linewidth(linewitdh)
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_size)
    fig.set_figwidth(fig_size)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

def comparative_profile(filename, y, title, xlabel, ylabel, ylog):
    lvl_size = max(len(y[0]), len(y[1]))
    fig, axs = plt.subplots()
    if(ylog):
        axs.set_yscale("log")
    dcco_style = dict(color="xkcd:battleship grey", linestyle="--")
    pdcco_style = dict(color="xkcd:black", linestyle="-")
    dcco_median_style = dict(color="xkcd:orange", linestyle="--")
    pdcco_median_style = dict(color="xkcd:red", linestyle="-")
    dict_1 = axs.boxplot(y[0], whis=(0, 100), positions=range(2, (2 * lvl_size) + 1, 2),boxprops=dcco_style, whiskerprops=dcco_style, capprops=dcco_style, medianprops=dcco_median_style)
#    axs.set_title("Pressure profile", fontsize=fontsize)
    axs.set_ylabel(ylabel, fontsize=fontsize)
    dict_2 = axs.boxplot(y[1], whis=(0, 100), positions=range(1, (2 * lvl_size) + 1, 2), boxprops=pdcco_style, whiskerprops=pdcco_style, capprops=pdcco_style, medianprops=pdcco_median_style)
    axs.legend((dict_2["whiskers"][0], dict_1["whiskers"][0]), ("PDCCO", "DCCO"), fontsize=fontsize, loc='upper right')
    axs.xaxis.set_major_locator(ticker.FixedLocator(np.arange(1.5, (2 * lvl_size) + 1, 2)))
    axs.xaxis.set_major_formatter(ticker.FixedFormatter(range(0, lvl_size + 1, 1)))
    axs.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(2.5, (2 * lvl_size) + 1, 2)))
    axs.tick_params(axis="y", which="major", labelsize=fontsize)
    axs.tick_params(axis="x", which="major", labelsize=fontsize)
    axs.tick_params(axis="x", which="minor", grid_linestyle=":", length=0, width=0)
    axs.grid(True, axis="x", which="minor")
    fig.set_dpi(300)
    fig.set_figheight(fig_size)
    fig.set_figwidth(1.5*fig_size)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

# Plots

prefix = "IMG_updated/part_ar_"
filename = prefix + "comparison_volumes_err_lin.pdf"
comparison_volumes_error_boxplot(filename, vol_merged_rel_error, labels, False)

for i in range(len(files)):
    merged_max = np.amax(merged[i][:,0])
    full_max = np.amax(full[i][:,0])
    real_max = max(merged_max, full_max)
    steps = np.arange(0, real_max + 1, 1, dtype=np.int)
    size = np.size(steps)
    bx_labels = steps
    
    y_full = []
    y_merged = []
    for j in range(size):
        y_merged_part = merged[i][(merged[i][:,0] == steps[j])]
        y_merged.append(y_merged_part)
        y_full_part = full[i][(full[i][:,0] == steps[j])]
        y_full.append(y_full_part)
    
    filename_base = "IMG_updated/part_ar_" + labels[i] + "_"
    
    radius = [[data[:,3] for data in y_full], [data[:,3] for data in y_merged]]
    length = [[data[:,4] for data in y_full], [data[:,4] for data in y_merged]]
    pressure = [[data[:,2] for data in y_full], [data[:,2] for data in y_merged]]
    flow = [[data[:,1] for data in y_full], [data[:,1] for data in y_merged]]

    filename = filename_base + "radius_comparison.pdf"
    comparative_profile(filename, radius, "Radius profile", "Level", "Radius [cm]", True)

    filename = filename_base + "pressure_comparison.pdf"
    comparative_profile(filename, pressure, "Pressure profile", "Level", "Pressure [mmHg]", False)


    # filename = filename_base + "radius_log.pdf"
    # profile(filename, radius, "Radius profile", "Level", "Radius (cm)", bx_labels, True)

    # filename = filename_base + "pressure.pdf"
    # profile(filename, pressure, "Pressure profile", "Level", "Pressure (mmHg)", bx_labels, False)

    plot_type = "_s"
    filename = prefix + labels[i] + plot_type + "_rad_x_len" + ".png"
    comparison_scatter_1x2(filename, [full[i][:, 3], merged[i][:, 3]], [full[i][:, 4], merged[i][:, 4]], "",  "Radius [cm]", "Length [cm]", True, True)
