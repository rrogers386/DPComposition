from optimal_comp import *
import matplotlib.pyplot as plt
import numpy as np

##############################################################################################################
##############################################################################################################
# This code is used to generate families of curves in a single plot to show the ratio of the optimal DP comp
# bounds from Murtagh and Vadhan '16  to the bounds for bounded range algos derived in the following  paper:
# https://arxiv.org/abs/1905.04273 (where these plots can be found).
##############################################################################################################
##############################################################################################################

delta_end = 10 ** (-6)
k_range = range(1, 101)


# Function to generate plots comparing the ratio of optimal DP composition with our composition for BR mechanisms #
def plotter(eps_per, delta, k_array):
    opt_comp_vals = [opt_compMV(k, eps_per, delta) for k in k_array]
    adv_comp_vals = [advanced_comp_br(k, eps_per, delta) for k in k_array]
    adv_comp_eps = np.array([adv_comp_vals[j] for j in range(len(k_array))])
    opt_comp_eps = np.array([opt_comp_vals[j] for j in range(len(k_array))])
    return opt_comp_eps / adv_comp_eps


##############################################################################################################
# We will generate a family of curves in a single plot, each with its own color, for various eps_per values. #
##############################################################################################################

# Switch parameter values here
# eps_per1 = 0.005
# eps_per2 = 0.01
# eps_per3 = 0.025
# eps_per4 = 0.05
eps_per1 = 0.1
eps_per2 = 0.25
eps_per3 = 0.5
eps_per4 = 1.0
y1 = plotter(eps_per1, delta_end, k_range)
y2 = plotter(eps_per2, delta_end, k_range)
y3 = plotter(eps_per3, delta_end, k_range)
y4 = plotter(eps_per4, delta_end, k_range)

####### Plot the curves ################
plt.interactive(False)
plt.plot(k_range, y1, 'r', label="$\epsilon =$  %1.3f" % eps_per1)
plt.plot(k_range, y2, 'g', label="$\epsilon =$ %1.3f" % eps_per2)
plt.plot(k_range, y3, 'b', label="$\epsilon =$ %1.3f" % eps_per3)
plt.plot(k_range, y4, 'y', label="$\epsilon =$ %1.3f" % eps_per4)

plt.xlabel("k")
plt.ylabel("Ratio of Opt Comp to Range Bounded Comp")
plt.title("Ratio of Privacy Loss for $\delta< 10^{-6}$")
plt.legend(loc='best')
# plt.savefig('SAVE_FILE_HERE.pdf')
plt.show()
