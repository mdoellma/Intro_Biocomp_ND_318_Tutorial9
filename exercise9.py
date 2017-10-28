import pandas as pd
import numpy as np
import scipy.stats
from scipy.optimize import minimize
import plotnine

### Begin challenge 1 ------------------------------------------------------------------------------------
ponzr1 = pd.read_csv("ponzr1.csv")

#function for returning negative log likelihood for t-test model
def nllike (p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    expected = B0 + B1 * obs.x
    nll = -1 * scipy.stats.norm(expected, sigma).logpdf(obs.ponzr1Counts).sum()
    return nll

#function for returning negative log likelihood for null model
def nllike_null (p, obs):
    B0 = p[0]
    sigma = p[1]
    expected = B0
    nll = -1 * scipy.stats.norm(expected, sigma).logpdf(obs.ponzr1Counts).sum()
    return nll

#function for returning p-value for t-test
def ttest (data, group1, group2):
    #define a temporary slice of the data
    temp_df = data[(data.mutation==group1) | (data.mutation==group2)]
    
    #Make new column 'x' set as 0 or 1 based on group
    temp_df["x"] = 0
    temp_df["x"][temp_df.mutation == group2] = 1
    
    # y = B0 + B1*x + E
    model = minimize(nllike, [1, 1, 1], method = "Nelder-Mead", options={'disp': True}, args = temp_df)
    
    # y = B0 + E
    null_model = minimize(nllike_null, [1, 1], method = "Nelder-Mead", options={'disp': True}, args = temp_df)
    
    #Get differences in fit
    D = (null_model.fun - model.fun) * 2
    
    #Use chi3.sf() for returning p-value
    p = scipy.stats.chi2.sf(D,1)
    
    #Print results
    print "-----------------------------"
    print (group1 + " vs. " + group2)
    print ("p-value = " + str(p))
    if p <= 0.05:
        print "Significance"
    else:
        print "No significance"
    print "-----------------------------"

#Perform t-tests
ttest(ponzr1, group1="WT", group2="M124K")
ttest(ponzr1, group1="WT", group2="V456D")
ttest(ponzr1, group1="WT", group2="I213N")

#Just plotting to visualize the data
plotnine.ggplot(ponzr1, plotnine.aes(x="mutation", y="ponzr1Counts")) + plotnine.geom_point()

### End challenge 1 -------------------------------------------------------------------------------
