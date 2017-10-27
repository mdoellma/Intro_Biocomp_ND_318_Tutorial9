import pandas as pd
import numpy as np
import scipy.stats
from scipy.optimize import minimize

ponzr = pd.read_csv("ponzr1.csv")

def nllike (p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    expected = B0 + B1 * obs.x
    nll = -1 * scipy.stats.norm(expected, sigma).logpdf(obs.ponzr1Counts).sum()
    return nll
    
def ttest (data, group1, group2):
    initialGuess = np.array([1,1,1])
    temp_df = data[(data.mutation==group1) | (data.mutation==group2)]
    temp_df["x"] = 0
    temp_df["x"][temp_df.mutation == group2] = 1

    fit1 = minimize(nllike, initialGuess, method = "Nelder-Mead", options={'disp': True}, args = temp_df[temp_df.mutation == group1])
    fit2 = minimize(nllike, initialGuess, method = "Nelder-Mead", options={'disp': True}, args = temp_df[temp_df.mutation == group2])
    
    D = fit1.x - fit2.x
    print scipy.stats.chi2.cdf(x = (2 * D), df = 1)

ttest(ponzr, "WT", "M124K")

