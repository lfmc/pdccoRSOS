import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as curve_fit

# def lin(x,a,b):
#     return (a*x) + b

# def sqr(x,a,b,c):
#     return ((a*x*x)+(b*x)+c)

# def cub(x,a,b,c,d):
#     return ((a*x*x*x)+(b*x*x)+(c*x)+d)

# def my_exp(x,a,b):
#     return (a*np.exp(x/1000))+b

# def lin(x,a):
#     return (a*x)

# def sqr(x,a,b):
#     return ((a*x*x)+(b*x))

def cub(x,a,b,c):
    return ((a*x*x*x)+(b*x*x)+(c*x))

# def my_exp(x,a):
#     return (a*np.exp(x/1000))

data = pd.read_csv("time_processed.csv", sep=",", header=None, names=["Horas", "N"], index_col=False)
# print(data)
data = data.to_numpy(dtype=np.float128);

# p1, _ = curve_fit(lin,data[:,1],data[:,0])
# p2, _ = curve_fit(sqr,data[:,1],data[:,0])
p3, _ = curve_fit(cub,data[:,1],data[:,0])
# p4, _ = curve_fit(my_exp,data[:,1],data[:,0])
x = np.linspace(0,100000,num=100,endpoint=True)
# y1 = lin(x, *p1)
# y2 = sqr(x, *p2)
y3 = cub(x, *p3)
# y4 = my_exp(x, *p4)

# print("lin:",p1,"\n")
# print("sqr:",p2,"\n")
print("cub:",p3,"\n")
# print("exp:",p4,"\n")
# print("lin(100000)={:f}\n".format(lin(100000,*p1)))
# print("sqr(100000)={:f}\n".format(sqr(100000,*p2)))
print("cub(95000)={:f}\n".format(cub(100000,*p3)))
# print("exp(100000)={:f}\n".format(my_exp(100000,*p4)))

linewitdh = 1.5
fontsize = 20
fig, axs = plt.subplots()
for axis in ['top','bottom','left','right']:
        axs.spines[axis].set_linewidth(linewitdh)
axs.tick_params(axis="both", which="major", length=8, width=2, labelsize=fontsize)
axs.set_xlim([0,40000])
axs.set_ylim([0,100])
axs.plot(data[:,1],data[:,0], label="DCCO", marker="o", linestyle="None")
axs.set_xlabel("$N_{\mathrm{seq}}$", fontsize=fontsize)
axs.set_ylabel("Tempo de execução [h]", fontsize=fontsize)
# axs.plot(x, y1, label="lin")
# axs.plot(x, y2, label="sqr")
axs.plot(x, y3, label="Ajuste cúbico")
# axs.plot(x, y4, label="exp")
axs.legend(fontsize=15)
fig.set_dpi(300)
fig.set_figheight(10)
fig.set_figwidth(10)
fig.tight_layout()
fig.savefig("times_fit.pdf")
plt.close(fig)