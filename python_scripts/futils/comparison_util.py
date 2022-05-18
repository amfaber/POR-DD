import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
def Linreg(x, y, ax = None, coefpos = (0.8, 0.85), corpos = (0.2, 0.1), coef = False):
    if ax is None:
        ax = plt.gca()
    ax.plot(x, y, ".")
    linreg = scipy.stats.linregress(x, y)
    spearmanr = scipy.stats.mstats.spearmanr(x, y)[0]
    xlims = np.array(ax.get_xlim())
    ys = linreg[0]*xlims+linreg[1]
    ax.plot(xlims, ys)

    if coef:
        txt = "$y=m x + b$\n"
        txt += f"$m = {linreg[0]:.3g}$\n$b = {linreg[1]:.3g}$"
        ax.text(*coefpos, txt, transform = ax.transAxes, va = "center", ha = "center")
    
    txt = f"Pearson's r = {linreg[2]:.3g}"#\nspearmans r = {spearmanr:.3g}"
    ax.text(*corpos, txt, transform = ax.transAxes, ha = "center")