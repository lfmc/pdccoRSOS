import pandas as pd
import numpy as np
# import scipy.stats as stats
# from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Define files

files = [["kidney_merged_gam30.csv"]]

# Proccess data
merged = []
vol_merged = []

data_merged = [pd.read_csv(listed, header=0, index_col=None, dtype=np.float64) 
               for listed in files[0]]
vol_merged_local = []
for table in data_merged:
    table.loc[:,"pressure"] = (1. / 133300.) * (table.loc[:,"eqResistance"] * table.loc[:,"flow"]) / ((table.loc[:,"radius"]) ** 4)
    table.loc[:,"pressure"] = table.loc[:,"pressure"] - table.loc[:,"pressure"].max()
    vol_merged_local.append(np.array(table.loc[:,"length"] * table.loc[:,"radius"] * table.loc[:,"radius"] * np.pi).sum())
    # Volumes
    vol_merged_local_np = np.array(vol_merged_local)
    vol_merged.append(vol_merged_local_np)
     
    data_merged_flat = (pd.concat(data_merged)).loc[:, ["level", "flow",
                                                        "pressure", "radius",
                                                        "length"]]
    data_merged_np = data_merged_flat.to_numpy()

    merged.append(data_merged_np)

# Define plot functions

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
    axs.set_ylabel(ylabel, fontsize=fontsize)
    axs.set_xlabel(xlabel, fontsize=fontsize)
    axs.xaxis.set_major_locator(ticker.FixedLocator(np.arange(1, (size) + 1, 2)))
    axs.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(2, (size) + 1, 2)))
    axs.xaxis.set_major_formatter(ticker.FixedFormatter(np.arange(0, (size) + 1, 2)))
    # axs.xaxis.set_minor_formatter(ticker.FixedFormatter(np.arange(1, (size) + 1, 2)))
    axs.tick_params(axis="x", which="major", labelsize=fontsize)
    axs.tick_params(axis="x", which="minor", labelsize=fontsize-5)
    # axs.tick_params(axis="y", which="major", length=8, width=2, labelsize=fontsize)
    axs.tick_params(axis="y", which="major", labelsize=fontsize)
    # axs.tick_params(axis="x", which="major", length=8, width=2, labelsize=12)
    
    # axs.grid(True, axis="x", which="minor")
    for axis in ['top', 'bottom', 'left', 'right']:
        axs.spines[axis].set_linewidth(linewitdh)
    fig.set_dpi(my_dpi)
    fig.set_figheight(fig_height)
    fig.set_figwidth(fig_width)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)

# Plots
prefix = "kidney_100k"

for i in range(len(files)):
    merged_max = np.amax(merged[i][:,0])
    real_max = merged_max
    steps = np.arange(0, real_max + 1, 1, dtype=int)
    size = np.size(steps)
    bx_labels = steps
   
    y_merged = []
    for j in range(size):
        y_merged_part = merged[i][(merged[i][:,0] == steps[j])]
        y_merged.append(y_merged_part)
   
    filename_base = "kidney_100k_"
    
    radius = [data[:,3] for data in y_merged]
    length = [data[:,4] for data in y_merged]
    pressure = [data[:,2] for data in y_merged]
    flow = [data[:,1] for data in y_merged]

    filename = filename_base + "radius.pdf"
    profile(filename, radius,  "Nível", "Raio [cm]", True)

    filename = filename_base + "pressure.pdf"
    profile(filename, pressure, "Nível", "Pressão [mmHg]", False)