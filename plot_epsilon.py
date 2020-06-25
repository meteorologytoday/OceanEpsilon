import matplotlib as mplt
mplt.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm

from netCDF4 import Dataset

import sys, argparse
import numpy as np

beg_y = 1
end_y = 100

init = True
for y in range(beg_y,end_y+1):
    print("Doing year: %d" % (y,))
    for m in range(1, 13):
        f = Dataset("output/epsilon_log_posterior_%03d-%02d.nc" % (y, m), "r")

        if init == True:
            init = False
            epsilon   = f.variables["epsilon"][:] * 86400.0
            log_post  = np.zeros((12, len(epsilon)), dtype=float)

        log_post[m-1, :] += f.variables["log_post"][:]
        f.close()

print(epsilon)

d_ep = epsilon[1] - epsilon[0]
post = np.zeros(log_post.shape, dtype=float)
max_i = np.zeros(post.shape[1], dtype=int)
max_epsilon = np.zeros(post.shape[1], dtype=float)
max_post = np.zeros(post.shape[1], dtype=float)
for m in range(12):
    log_post[m, :]  -= np.amax(log_post[m, :])
    post[m, :]  = np.exp(log_post[m, :])
    post[m, :] /= np.sum(post[m, :]) * d_ep
    max_i[m] = np.argmax(post[m, :])
    max_epsilon[m] = epsilon[max_i[m]]
    max_post[m] = post[m, max_i[m]]


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
  


#ax.plot([max_epsilon, max_epsilon], [0, 1], dashes=(7,3), color='black')

detail = ""

for m in range(12):
    ax.plot(epsilon, post[m, :])
    ax.scatter(max_epsilon[m], max_post[m], s=20, marker='o', color='red', zorder=999)
    ax.text(max_epsilon[m], max_post[m] + 1, '%d' % (m+1,), ha="center")
    detail = '%s\n$\epsilon^*[%d] = %.2f \, \mathrm{day}^{-1}$' % (detail, m+1, max_epsilon[m],)


ax.text(0.5, np.amax(max_post), detail, ha="left", va="top", transform=ax.transData)

ax.set_title('Posterior of $ \epsilon $ of year %d - %d' % ( beg_y, end_y) )
ax.set_xlabel('$\epsilon$ [ $ \mathrm{day}^{-1} $ ]')
ax.set_ylabel('PDF of $\epsilon$ [ $ \mathrm{day} $ ]')
#ax.set_ylim([-0, 1])
ax.grid(True)


fig.savefig("result_monthly_%03d-%03d.png" % (beg_y, end_y) , dpi=200)
