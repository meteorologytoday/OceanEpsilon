import matplotlib as mplt
#mplt.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm

from netCDF4 import Dataset

import sys, argparse
import numpy as np


models = [
    "NCAR_CESM2",
    "E3SM",
]

log_post = {}
post = {}
max_i = {}
max_eps = {}
max_post = {}

beg_y = 1
end_y = 2

init = True

for model in models:



    for y in range(beg_y,end_y+1):
        print("Doing year: %d" % (y,))
        for m in range(1, 13):
            with Dataset("result/%s/epsilon_log_posterior_%03d-%02d.nc" % (model, y, m), "r") as f:
                if init == True:
                    init = False
                    eps   = f.variables["eps"][:] * 86400.0
                    d_ep = eps[1] - eps[0]
                    _log_post = np.zeros((12, len(eps)), dtype=float)

                _log_post[m-1, :] += f.variables["log_post"][:]


    _post        = np.zeros(_log_post.shape, dtype=float)
    _max_i       = np.zeros(_post.shape[1], dtype=int)
    _max_eps = np.zeros(_post.shape[1], dtype=float)
    _max_post    = np.zeros(_post.shape[1], dtype=float)

    for m in range(12):
        _log_post[m, :]  -= np.amax(_log_post[m, :])
        _post[m, :]  = np.exp(_log_post[m, :])
        _post[m, :] /= np.sum(_post[m, :]) * d_ep
        _max_i[m] = np.argmax(_post[m, :])
        _max_eps[m] = eps[_max_i[m]]
        _max_post[m] = _post[m, _max_i[m]]


    
    post[model]  = _post
    max_i[model] = _max_i
    max_eps[model] = _max_eps
    max_post[model] = _max_post
    
#print(eps)


for model in models:

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))

    detail = ""

    for m in range(12):
        ax.plot(eps, post[model][m, :])
        ax.scatter(max_eps[model][m], max_post[model][m], s=20, marker='o', color='red', zorder=999)
        ax.text(max_eps[model][m], max_post[model][m] + 1, '%d' % (m+1,), ha="center")
        detail = '%s\n$\epsilon^*[%d] = %.2f \, \mathrm{day}^{-1}$' % (detail, m+1, max_eps[model][m],)


    #plt.text(0.9, 0.1, np.amax(max_post[model]), detail, ha="left", va="top", transform=fig.transFigure)
    plt.text(0.8, 0.9, detail, ha="left", va="top", transform=fig.transFigure)

    ax.set_title('[%s] Posterior of $ \epsilon $ of year %d to %d' % ( model, beg_y, end_y ) )
    ax.set_xlabel('$\epsilon$ [ $ \mathrm{day}^{-1} $ ]')
    ax.set_ylabel('PDF of $\epsilon$ [ $ \mathrm{day} $ ]')
    #ax.set_ylim([-0, 1])
    ax.grid(True)

    plt.subplots_adjust(right=.75)

    fig.savefig("graph/%s_monthly_%03d-%03d.png" % (model, beg_y, end_y) , dpi=200)

    plt.show()
    plt.close(fig)



