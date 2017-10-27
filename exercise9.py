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
    
def nllike_null (p, obs):
    B0 = p[0]
    sigma = p[1]
    expected = B0
    nll = -1 * scipy.stats.norm(expected, sigma).logpdf(obs.ponzr1Counts).sum()
    return nll

def ttest (data, group1, group2):
    initialGuess = np.array([1,1,1])
    temp_df = data[(data.mutation==group1) | (data.mutation==group2)]
    temp_df["x"] = 0
    temp_df["x"][temp_df.mutation == group2] = 1
    
    # y = B0 + B1*x + E
    model = minimize(nllike, [1, 1, 1], method = "Nelder-Mead", options={'disp': True}, args = temp_df)
    
    # y = B0 + E
    null_model = minimize(nllike_null, [1, 1], method = "Nelder-Mead", options={'disp': True}, args = temp_df)
    
    
    D = model.fun - null_model.fun
    print scipy.stats.chi2.cdf(x = (2 * D), df = 1)

ttest(ponzr, "WT", "I213N")


