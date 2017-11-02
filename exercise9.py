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
    print ("-----------------------------")
    print (group1 + " vs. " + group2)
    print ("p-value = " + str(p))
    if p <= 0.05:
        print ("Significance")
    else:
        print ("No significance")
    print ("-----------------------------")

#Perform t-tests
ttest(ponzr1, group1="WT", group2="M124K")
ttest(ponzr1, group1="WT", group2="V456D")
ttest(ponzr1, group1="WT", group2="I213N")

#Just plotting to visualize the data
plotnine.ggplot(ponzr1, plotnine.aes(x="mutation", y="ponzr1Counts")) + plotnine.geom_point()

### End challenge 1 -------------------------------------------------------------------------------
### Begin Challenge 2 -----------------------------------------------------------------------------

#read in files
MG=pd.read_csv('MmarinumGrowth.csv')

#function for returning negative log likelihood for t-test model
def nllike_K (p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    expected = B0 * (obs.S/(obs.S + B1))
    nll = -1 * scipy.stats.norm(expected, sigma).logpdf(obs.u).sum()
    return nll

#initial guess variable, supposedly minimizes the negative log likelihood but I don't know how this works
initialGuess=np.array([1,1,1])
fit=minimize(nllike_K,initialGuess,method="Nelder-Mead",options={'disp': True},args=MG)
print(fit)
### End of Challenge 2 -----------------------------------------------------------------------------
###  Begin Challenge 3 -----------------------------------------------------------------------------
#read in file
LD=pd.read_csv('leafDecomp.csv')

#an nll function that can take any number of Betas
def nllike_v (p, obs):
    B = []
    x = obs.Ms
    
    #fill all B0 to Bn-1 because the last number is reserved for sigma
    for i in range(len(p) - 1):
        B.append(p[i])
    sigma = p[len(p) - 1]
    
    #iterally build the expected value:
    ### expected += B0 * x^0 = B0
    ### expected += B1 * x^1 = B0 + B1 * x
    ### expected += B2 * x^2 = B0 + B1 * x + B2 * x^2
    expected = 0
    for i in range(len(B)):
        expected = expected + B[i] * (x ** i)
        
    #calculate nll
    nll = -1 * scipy.stats.norm(expected, sigma).logpdf(obs.decomp).sum()
    return nll
    
#get the fit values
fit_1 = minimize(nllike_v, [1,1], method="Nelder-Mead", options={'disp': True}, args=LD)
fit_2 = minimize(nllike_v, [1,1,1], method="Nelder-Mead", options={'disp': True}, args=LD)
fit_3 = minimize(nllike_v, [1,1,1,1], method="Nelder-Mead", options={'disp': True}, args=LD)

#do the likelihood comparisons
scipy.stats.chi2.sf((fit_1.fun - fit_2.fun) * 2,1)
scipy.stats.chi2.sf((fit_2.fun - fit_3.fun) * 2,1)
scipy.stats.chi2.sf((fit_1.fun - fit_3.fun) * 2,2)

#constant rate maximum likelihood d=a
def nllike_LDC (p, obs):
    B0 = p[0]
    sigma = p[1]
    expected = B0
    nll = -1 * scipy.stats.norm(expected, sigma).logpdf(obs.decomp).sum()
    return nll
initialGuess=[1,1]
fitLDC=minimize(nllike_LDC,initialGuess,method="Nelder-Mead",options={'disp': True},args=LD)

#linear rate maximum likelihood d=a+bMs
def nllike_LDL (p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    expected = B0 + B1 * obs.Ms
    nll = -1 * scipy.stats.norm(expected, sigma).logpdf(obs.decomp).sum()
    return nll
initialGuess=[1,1,1]
fitLDL=minimize(nllike_LDL,initialGuess,method="Nelder-Mead",options={'disp': True},args=LD)

#quadratic rate maximum likelihood d=a+bMS+cMs^2
def nllike_LDQ (p, obs):
    B0 = p[0]
    B1 = p[1]
    B2 = p[2]
    sigma = p[3]
    
    expected = B0 + (B1 * obs.Ms) + B2 * (obs.Ms**(2))
    nll = -1 * scipy.stats.norm(expected, sigma).logpdf(obs.decomp).sum()
    return nll
initialGuess=[1,1,1,1]
fitLDQ=minimize(nllike_LDQ,initialGuess,method="Nelder-Mead",options={'disp': True},args=LD)

#likelihood comparisons
scipy.stats.chi2.sf((fitLDC.fun - fitLDL.fun) * 2,1)
scipy.stats.chi2.sf((fitLDL.fun - fitLDQ.fun) * 2,1)
scipy.stats.chi2.sf((fitLDC.fun - fitLDQ.fun) * 2,2)


