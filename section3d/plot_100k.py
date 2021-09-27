import pandas as pd
import numpy as np
import scipy.stats as stats
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Define files

file_p16 = "100k_merged_p16.csv"
file_p25 = "100k_merged_p25.csv"
file_p36 = "100k_merged_p36.csv"

files = [file_p16,file_p25,file_p36]
labels = ["p16","p25","p36"]

# Proccess data
merged = []
vol_merged = []

for i in range(len(files)):
    data_merged = pd.read_csv(files[i], header = 0, index_col = None, dtype = np.float64) 
    vol_merged_local = []
        
    data_merged.loc[:,"pressure"] = (1. / 133300.) * ((data_merged.loc[:,"eqResistance"] * data_merged.loc[:,"flow"]) / ((data_merged.loc[:,"radius"]) ** 4))
    data_merged.loc[:,"pressure"] = data_merged.loc[:,"pressure"] - data_merged.loc[:,"pressure"].max()
    vol_merged.append(np.array(data_merged.loc[:,"length"] * data_merged.loc[:,"radius"] * data_merged.loc[:,"radius"] * np.pi).sum())
 
    data_merged_np = (data_merged).loc[:,["level", "flow", "pressure", "radius", "length"]].to_numpy()

    merged.append(data_merged_np)   

# print(merged)
# Define plot functions

my_dpi = 300
my_dark_blue_rgb = [0, 0, 154./255]
my_dark_green_rgb = [0, 1, 180./255]

linewitdh = 1.5
#fontsize = "xx-large"
fontsize = 20
fontweight = "bold"
fig_size = 10

def comparison_scatter_1x2(filename, x, y, title, xlabel, ylabel, xlog, ylog):
    fig, axs = plt.subplots(1, 2, sharey='all', sharex='row', squeeze=False)
    if (xlog):
        axs[0, 0].set_xscale("log")
    if (ylog):
        axs[0, 0].set_yscale("log")
    axs[0, 0].scatter(x[0], y[0], c=[my_dark_green_rgb], alpha=0.1, linewidths=0)
    axs[0, 0].set_xlabel(xlabel, fontsize = fontsize)
    axs[0, 0].set_ylabel(ylabel, fontsize = fontsize)
    axs[0, 0].set_title("DCCO", fontsize = fontsize)
    axs[0,0].tick_params(axis="both", which="major", length=8, width=2, labelsize=fontsize)
    axs[0,1].scatter(x[1], y[1], c=[my_dark_blue_rgb], alpha=0.1, linewidths=0)
    axs[0,1].set_xlabel(xlabel, fontsize = fontsize)
#    axs[0,1].set_ylabel(ylabel, fontsize = fontsize)
    axs[0,1].set_title("PDCCO", fontsize = fontsize)
    axs[0,1].tick_params(axis="both", which="major", length=8, width=2, labelsize=fontsize)
    for axis in ['top','bottom','left','right']:
        axs[0, 0].spines[axis].set_linewidth(linewitdh)
        axs[0, 1].spines[axis].set_linewidth(linewitdh)
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_size)
    fig.set_figwidth(2 * fig_size)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

def comparison_volumes(filename, y):
    fig, axs = plt.subplots()
    # if (ylog):
    #     axs.set_yscale("log")
    # axs.set_title("Volume", fontsize = fontsize)
    axs.set_ylabel("Volume [mm$^{3}$]", fontsize = fontsize)
    axs.set_xlabel("$N_\mathrm{part}$", fontsize = fontsize)
    rect0=axs.bar(1,y[0],label="$N_\mathrm{part}=16$")
    # axs.bar_label(rect0)
    rect1=axs.bar(2,y[1],label="$N_\mathrm{part}=25$")
    # axs.bar_label(rect1)
    rect2=axs.bar(3,y[2],label="$N_\mathrm{part}=36$")
    # axs.bar_label(rect2)
    axs.xaxis.set_major_locator(ticker.FixedLocator([1,2,3]))
    axs.xaxis.set_major_formatter(ticker.FixedFormatter(["16","25","36"]))
    # axs.legend()
    axs.tick_params(axis="y", which="major", length=8, width=2, labelsize=fontsize)
    axs.tick_params(axis="x", which="major", length=8, width=2, labelsize=fontsize)
    for axis in ['top','bottom','left','right']:
        axs.spines[axis].set_linewidth(linewitdh)
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_size)
    fig.set_figwidth(fig_size)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

def comparative_profile(filename, y, xlabel, ylabel, ylog, xbins):
    lvl_size = max(len(y[0]),max(len(y[1]),len(y[2])))
    fig, axs = plt.subplots()
    if(ylog):
        axs.set_yscale("log")
    style_00 = dict(color="xkcd:black", linestyle="-")
    style_01 = dict(color="xkcd:greyish blue", linestyle="-")   
    style_02 = dict(color="xkcd:battleship grey", linestyle="-")
    median_style_00 = dict(color="xkcd:red", linestyle="-")
    median_style_01 = dict(color="xkcd:orange", linestyle="-")
    median_style_02 = dict(color="xkcd:yellowish orange", linestyle="-")
    dict_0 = axs.boxplot(y[0], whis=(0, 100), positions=range(1, (3 * lvl_size) + 1, 3),boxprops=style_00, whiskerprops=style_00, capprops=style_00, medianprops=median_style_00)
    dict_1 = axs.boxplot(y[1], whis=(0, 100), positions=range(2, (3 * lvl_size) + 1, 3), boxprops=style_01, whiskerprops=style_01, capprops=style_01, medianprops=median_style_01)
    dict_2 = axs.boxplot(y[2], whis=(0, 100), positions=range(3, (3 * lvl_size) + 1, 3), boxprops=style_02, whiskerprops=style_02, capprops=style_02, medianprops=median_style_02)
    axs.set_ylabel(ylabel, fontsize=fontsize)
    axs.set_xlabel(xlabel, fontsize=fontsize)
    axs.legend((dict_0["whiskers"][0], dict_1["whiskers"][0], dict_2["whiskers"][0]), ("$N_{\mathrm{part}}=$16", "$N_{\mathrm{part}}=$25", "$N_{\mathrm{part}}=$36"), fontsize=fontsize, loc='upper right')
    axs.xaxis.set_major_locator(ticker.FixedLocator(np.arange(2, (3 * lvl_size) + 1, 3)))
    axs.xaxis.set_major_formatter(ticker.FixedFormatter(xbins))
    # my_minor_ticks=np.array([])
    axs.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(3.5, (3 * lvl_size) + 1, 3)))
    axs.tick_params(axis="y", which="major", labelsize=fontsize)
    axs.tick_params(axis="x", which="major", labelsize=fontsize-8)
    axs.tick_params(axis="x", which="minor", grid_linestyle=":", length=0, width=0)
    axs.grid(True, axis="x", which="minor")
    fig.set_dpi(300)
    fig.set_figheight(fig_size)
    fig.set_figwidth(2*fig_size)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)


# Plots
# print(repr(merged[0]))
# print(repr(merged[1]))
# print(repr(merged[2]))
# print(merged[0][:,0])
# print(merged[1][:,0])
# print(merged[2][:,0])
real_max = max(np.amax(merged[0][:,0]),max(np.amax(merged[1][:,0]),np.amax(merged[2][:,0])))
print(real_max)
bin_size = 5
steps = np.arange(bin_size, real_max + bin_size, bin_size, dtype=np.int)
size = np.size(steps)
print(steps)
bx_labels = []
for i in range(size):
    bx_labels.append("{:d}-{:d}".format(steps[i]-bin_size,steps[i]-1))

   
#Binarizing data    
y_binarized_all = []
for i in range(3):
    y_merged = []
    for j in range(size):
        y_merged_part = merged[i][(merged[i][:,0] >= (steps[j]-bin_size)) & (merged[i][:,0] < steps[j])] 
        y_merged.append(y_merged_part)
    y_binarized_all.append(y_merged)
        
filename_base = "100k_network_"
    
radius = [[data[:,3] for data in y_binarized_all[0]], [data[:,3] for data in y_binarized_all[1]], [data[:,3] for data in y_binarized_all[2]]]
length = [[data[:,4] for data in y_binarized_all[0]], [data[:,4] for data in y_binarized_all[1]], [data[:,4] for data in y_binarized_all[2]]]
pressure = [[data[:,2] for data in y_binarized_all[0]], [data[:,2] for data in y_binarized_all[1]], [data[:,2] for data in y_binarized_all[2]]]
flow = [[data[:,1] for data in y_binarized_all[0]], [data[:,1] for data in y_binarized_all[1]], [data[:,1] for data in y_binarized_all[2]]]

filename = filename_base + "radius_comparison.pdf"
comparative_profile(filename, radius, "Level", "Radius [mm]", True, bx_labels)

filename = filename_base + "pressure_comparison.pdf"
comparative_profile(filename, pressure, "Level", "Pressure [mmHg]", False, bx_labels)

# filename = filename_base + "volumes.pdf"
# comparison_volumes(filename, vol_merged)

# plot_type = "_s"
# filename = prefix + labels[i] + plot_type + "_rad_x_len" + ".png"
# comparison_scatter_1x2(filename, [full[i][:, 3], merged[i][:, 3]], [full[i][:, 4], merged[i][:, 4]], "",  "Radius [cm]", "Length [cm]", True, True)
