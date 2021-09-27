import pandas as pd
import numpy as np
import scipy.stats as stats
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Define files

files = [["kidney_merged_gam30_strahler.csv"]]

# Proccess data
merged = []
vol_merged = []

data_merged = [pd.read_csv(listed, header = 0, index_col = None, dtype = np.float64) for listed in files[0]]
vol_merged_local = []
for table in data_merged:
    table.loc[:,"pressure"] = (1. / 133300.) * (table.loc[:,"eqResistance"] * table.loc[:,"flow"]) / ((table.loc[:,"radius"]) ** 4)
    table.loc[:,"pressure"] = table.loc[:,"pressure"] - table.loc[:,"pressure"].max()
    vol_merged_local.append(np.array(table.loc[:,"length"] * table.loc[:,"radius"] * table.loc[:,"radius"] * np.pi).sum())
    # Volumes
    vol_merged_local_np = np.array(vol_merged_local)
    vol_merged.append(vol_merged_local_np)
     
    data_merged_flat = (pd.concat(data_merged)).loc[:,["level","strahler", "subtreeVolume", "strahlerConnectivity", "flow", "pressure", "radius", "length"]]
    data_merged_np = data_merged_flat.to_numpy()

    merged.append(data_merged_np)

# Define plot functions

idxLevel = 0
idxStrahler=1
idxSubtreeVolume=2
idxStrahlerConnectivity=3
idxFlow=4
idxPressure=5
idxRadius=6
idxLength=7

my_dpi = 300
my_dark_blue_rgb = [0, 0, 154./255]
my_dark_green_rgb = [0, 1, 180./255]

linewitdh = 0.5
#fontsize = "xx-large"
fontsize = 12
fontweight = "bold"
cm = 1./2.54
fig_width=7.5*cm*3
fig_height=5.14*cm*3

def profile(filename, y, xlabel, ylabel, ylog):
    size = len(y)
    fig, axs = plt.subplots()
    if (ylog):
        axs.set_yscale("log")
    axs.boxplot(y, whis=(0, 100))
    # axs.set_title(title + " - DCCO", fontsize = fontsize)
    axs.set_ylabel(ylabel, fontsize = fontsize)
    axs.set_xlabel(xlabel, fontsize = fontsize)
    axs.xaxis.set_major_locator(ticker.FixedLocator(np.arange(1, (size) + 1, 5)))
    # axs.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(2, (size) + 1, 2)))
    axs.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    axs.xaxis.set_major_formatter(ticker.FixedFormatter(np.arange(0, (size) + 1, 5)))
    # axs.xaxis.set_minor_formatter(ticker.FixedFormatter(np.arange(1, (size) + 1, 2)))
    axs.tick_params(axis="x", which="major", labelsize=fontsize)
    axs.tick_params(axis="x", which="minor", labelsize=fontsize-5)
    # axs.tick_params(axis="y", which="major", length=8, width=2, labelsize=fontsize)
    axs.tick_params(axis="y", which="major", labelsize=fontsize)
    # axs.tick_params(axis="x", which="major", length=8, width=2, labelsize=12)
    
    # axs.grid(True, axis="x", which="minor")
    for axis in ['top','bottom','left','right']:
        axs.spines[axis].set_linewidth(linewitdh)
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

def scatter(filename, x, y, title, xlabel, ylabel, xlog, ylog):
    fig, axs = plt.subplots()
    if (xlog):
        axs.set_xscale("log")
    if (ylog):
        axs.set_yscale("log")
    axs.scatter(x, y, c="xkcd:blue", alpha=0.1, linewidths=0)
    axs.set_xlabel(xlabel, fontsize = fontsize)
    axs.set_ylabel(ylabel, fontsize = fontsize)
    # axs.set_title("PDCCO", fontsize = fontsize)
    # axs.tick_params(axis="both", which="major", length=8, width=2, labelsize=fontsize)
    axs.tick_params(axis="both", which="major", labelsize=fontsize)
    for axis in ['top','bottom','left','right']:
        axs.spines[axis].set_linewidth(linewitdh)
    # fig.suptitle(title)
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

def simplePlot(filename, x, y, xlabel, ylabel):
    fig, axs = plt.subplots()
    axs.plot(x, y, 'o-k')
    axs.set_xlabel(xlabel, fontsize = fontsize)
    axs.set_ylabel(ylabel, fontsize = fontsize)
    axs.tick_params(axis="both", which="major", labelsize=fontsize)
    for axis in ['top','bottom','left','right']:
        axs.spines[axis].set_linewidth(linewitdh)
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

def heatmap(data, size):
    fig, axs = plt.subplots()
    axs.matshow(data)
    # axs.xaxis.set_major_locator(ticker.FixedLocator(np.arange(-0.5, (size)+1,1)))
    # axs.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))
    # axs.yaxis.set_minor_locator(ticker.AutoMinorLocator(10))
    # axs.tick_params(axis="both", which="both", labelsize=fontsize)
    axs.grid(True, which="both", axis="both", color="w")
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

# Plots
prefix = "kidney_100k"

for i in range(len(files)):
    merged_max = np.amax(merged[i][:,idxStrahler])
    steps = np.arange(0, merged_max+1, 1, dtype=int)
    size = np.size(steps)
    bx_labels = steps
   
    y_merged = []
    for j in range(size):
        y_merged_part = merged[i][(merged[i][:,idxStrahler] == steps[j])]
        y_merged.append(y_merged_part)
   
    filename_base = "kidney_100k_"
    

    radius = [data[:,idxRadius] for data in y_merged]
    radius_mean = [np.mean(rad) for rad in radius]
    radius_std = [np.std(rad) for rad in radius]
    strahler_size = [np.size(rad) for rad in radius]
    length = [data[:,idxLength] for data in y_merged]
    length_mean = [np.mean(leng) for leng in length]
    length_std = [np.std(leng) for leng in length]
    
    filename = filename_base + "summary_strahler.csv"
    strahler_summary = pd.DataFrame(np.concatenate((np.reshape(steps, (size,1)), np.reshape(radius_mean, (size,1)),np.reshape(radius_std, (size,1)),
    np.reshape(length_mean, (size,1)),np.reshape(length_std, (size,1)), np.reshape(strahler_size, (size,1))),axis=1), index=None,
     columns=["Strahler order", "Mean radius [cm]", "Std radius [cm]", "Mean length [cm]", "Std length [cm]", "Size"])
    strahler_summary.to_csv(filename, sep=",", index=False)
    print(strahler_summary)
    strahler_summary.to_csv()

    crossArea = [data[:,idxRadius]*data[:,idxRadius]*np.pi for data in y_merged]
    # print(crossArea)
    totalCrossArea = [np.sum(area) for area in crossArea]
    # print(totalCrossArea)
    pressure = [data[:,idxPressure] for data in y_merged]
    flow = [data[:,idxFlow] for data in y_merged]
    strahlerConnectivity = [data[:,idxStrahlerConnectivity] for data in y_merged]
    # print(strahlerConnectivity)
    connectivity_matrix = np.zeros((size,size), dtype=int)
    for j in range(len(strahlerConnectivity)):
        for k in range(np.size(strahlerConnectivity[j])):
            if (strahlerConnectivity[j][k]>=0):
                connectivity_matrix[j,int(strahlerConnectivity[j][k])] = connectivity_matrix[j,int(strahlerConnectivity[j][k])]+1;
    connectivity_matrix_df = pd.DataFrame(connectivity_matrix)
    # filename = filename_base + "_connectivity_strahler.csv"
    # filename = filename_base + "_connectivity_heatmap_strahler.pdf"
    # print(connectivity_matrix)
    # heatmap(connectivity_matrix, size)
    # connectivity_matrix_df.to_csv(filename)
    # print(connectivity_matrix)
    subtreeVolume = [data[:,idxSubtreeVolume] for data in y_merged]

    filename = filename_base + "cross_area_strahler.pdf"
    simplePlot(filename, steps, totalCrossArea, "Strahler order", "Total cross-sectional area [$\mathrm{cm}^2$]")

    # print(idxRadius)    # print(merged[i].size)
    # print(merged[i][:,idxRadius])
    # print(merged[i][:,idxSubtreeVolume])
    filename = filename_base + "scatter_strahler"
    scatter(filename, merged[i][:,idxRadius], merged[i][:,idxSubtreeVolume], "", "Radius [cm]", "Subtree Volume [mL]", False, True)
    
    filename = filename_base + "radius_strahler.pdf"
    profile(filename, radius,  "Strahler order", "Radius [cm]", True)

    filename = filename_base + "pressure_strahler.pdf"
    profile(filename, pressure, "Strahler order", "Pressure [mmHg]", False)

    filename = filename_base + "length_strahler.pdf"
    profile(filename, length, "Strahler order", "length [cm]", False)

    filename = filename_base + "flow_strahler.pdf"
    profile(filename, flow, "Strahler order", "flow [mL/s]", True)

    filename = filename_base + "subtreeVolume_strahler.pdf"
    profile(filename, subtreeVolume, "Strahler order", "Subtree volume [mL]", True)

    filename = filename_base + "strahlerConnectivity_strahler.pdf"
    profile(filename, strahlerConnectivity, "Strahler order", "Connectivity", False)