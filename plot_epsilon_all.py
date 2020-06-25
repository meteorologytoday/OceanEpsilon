import matplotlib as mplt
mplt.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm

from netCDF4 import Dataset

import sys, argparse
import numpy as np


beg_y = 11
end_y = 100

init = True
for y in range(beg_y,end_y+1):
    print("Doing year: %d" % (y,))
    for m in range(1, 13):
        f = Dataset("output/epsilon_log_posterior_%03d-%02d.nc" % (y, m), "r")

        if init == True:
            init = False
            epsilon   = f.variables["epsilon"][:] * 86400.0
            log_post  = np.zeros(len(epsilon), dtype=float)

        log_post += f.variables["log_post"][:]
        f.close()

print(epsilon)

d_ep = epsilon[1] - epsilon[0]
log_post  -= np.amax(log_post)
posterior  = np.exp(log_post)
posterior /= np.sum(posterior) * d_ep


max_i = np.argmax(posterior)
max_epsilon = epsilon[max_i]
max_posterior = posterior[max_i]


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
  

ax.plot(epsilon, posterior)
ax.plot([max_epsilon, max_epsilon], [0, 1], dashes=(7,3), color='black')
ax.scatter(max_epsilon, max_posterior, s=30, marker='o', color='red', zorder=999)
ax.text(max_epsilon + 0.1, max_posterior, '$\epsilon^* = %.2f \, \mathrm{day}^{-1}$' % (max_epsilon,))

ax.set_title('Posterior of $ \epsilon $')
ax.set_xlabel('$\epsilon$ [ $ \mathrm{day}^{-1} $ ]')
ax.set_ylabel('PDF of $\epsilon$ [ $ \mathrm{day} $ ]')
#ax.set_ylim([-0, 1])
ax.grid(True)


fig.savefig("result.png", dpi=200)
