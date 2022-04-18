import pandas as pd
import numpy as np
# import scipy.stats as stats
from pathlib import Path
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker

# Define files
p = Path('.')
files = list(p.glob("rsMultipleInput*.csv"))
std = []
max_dev = []
ratio = []
flows = []

for file in files:
    df = pd.read_csv(file, header=0, index_col=None)
    flow_pt = df.query("level == 4")["flow"] / 0.72
    flows.append(flow_pt)
    dev = flow_pt.std(ddof=0)
    std.append(dev)
    abs_dev = np.abs(flow_pt - 0.25)
    max_dev.append(np.amax(abs_dev))
    ratio_cur = flow_pt.min()/flow_pt.max()
    ratio.append(ratio_cur)

std_np = np.array(std)
sort_idx = np.argsort(std_np)
for j in range(10):
    i = sort_idx[sort_idx.size-j-1]
    print(flows[i])
    print("ratio = ", ratio[i])
    print("max = ", max_dev[i])
    print("std = ", std_np[i])
    print("file = ", files[i])
