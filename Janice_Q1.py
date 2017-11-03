#Exercise 9 Q1 
#Authors: Janice Love and Melissa Stephens 
#Date: October 27, 2017 

#T-test and maximum liklihood ratio scripts 

import pandas as pd 
import numpy as np 
import scipy
from scipy.optimize import minimize 
from scipy.stats import norm 
from scipy.stats import chi2


file=pandas.read_csv("ponzr1.csv",header=0,sep=",")

#--from tutorial-------------------------------------------------------------------------
def nllike(arguments):
    unpack arguments, assign variables
    calc expected value (model equation)
    calc negative log likelihood 
    return negative log likelihood 

def nllike(p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    
    expected = B0 + B1*obs.x #make sure data import has a column called x 
    
    nll=-1*norm(expected,sigma).logpdf(obs.y)sum()
    return nll 
    
pval = 1-chi2.cdf(x=2*(nllalt-nullnull), df = diffInNumParam)    
#----------------------------------------------------------------------------

WT_M124=file.loc[file.mutation.isin(['WT', 'M124K']),:]
WT_V456D=file.loc[file.mutation.isin(['WT', 'V456D']),:]
WT_I213N=file.loc[file.mutation.isin(['WT', 'I213N']),:]

#change columns to x and y, and treatment to 0 (WT) or 1 (mutation) for x in equation down below 
WT_M124=pandas.DataFrame({'y':mut1.ponzr1Counts, 'x':0})
WT_M124.loc[mut1.mutation=='M124K', 'x']=1

WT_V456D=pandas.DataFrame({'y':mut2.ponzr1Counts, 'x':0})
WT_V456D.loc[mut2.mutation=='V456D', 'x']=1

WT_I213N=pandas.DataFrame({'y':mut3.ponzr1Counts, 'x':0})
WT_I213N.loc[mut3.mutation=='I213N', 'x']=1

#for WT 
def nullnull(p, obs):
    B0 = p[0]
    sigma = p[1]
    
    expected = B0 
    
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll 
    
initialGuess=np.array([1,1,1])

#------------------------------------------------------------------------------------------------------



#For WT_M124K -----------------------------------------------------------------------------------------
def nllike(p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    
    expected = B0 + B1*obs.x #make sure data import has a column called x 
    
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll 
    
initialGuess=np.array([1,1,1])
fit=minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp': True},args=WT_M124K)    
print (fit.x)

WT_M124K = fit.x 

nullfit=minimize(nullnull,initialGuess,method="Nelder-Mead",options={'disp': True},args=WT_M124)    
print (fit.x)

WT_null = fit.x

#--------------------------------------------------------------------------------------------------------
D = (2*(nullfit.fun - WT_M124K.fun))
chi_M124K = 1-scipy.stats.chi2.cdf(x=D,df=1)
print ("Value for M124K mutation:",chi_M124K)

#For WT_I213N ----------------------------------------------------------------------------------------
def nllikeI213N(p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    
    expected = B0 + B1*obs.x #make sure data import has a column called x 
    
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll 
    
initialGuess=np.array([1,1,1])
fit=minimize(nllikeI213N,initialGuess,method="Nelder-Mead",options={'disp': True},args=WT_I213N)
I213N_fit = fit.x

WT_null = minimize(nllike,initialGuess,method="Nelder-Mead", options = {'disp': True}, args = WT_I213N)
print (fit.x)
nullfit = fit.x 
#------------------------------------------------------------------------------------------------------
D = (2*(nullfit.fun - I213N_fit.fun))
chi_I213N = 1-scipy.stats.chi2.cdf(x=D,df=1)
print ("Answer for I213N mutation:", chi_I213N)

#For WT_V456D ----------------------------------------------------------------------------------------
def nllike3(p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    
    expected = B0 + B1*obs.x #make sure data import has a column called x 
    
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll 
    
initialGuess=np.array([1,1,1])
fit=minimize(nllike3,initialGuess,method="Nelder-Mead",options={'disp': True},args=WT_V456D)    
print (fit.x)
V456D_fit = fit.x

WT_null = minimize(nllike,initialGuess,method="Nelder-Mead", options = {'disp': True}, args = WT_V456D)
print (fit.x)
nullfit = fit.x 
#------------------------------------------------------------------------------------------------------
D = (2*(nullfit.fun - V456D_fit.fun))
chi_V456N = 1-scipy.stats.chi2.cdf(x=D,df=1)
print ("Answer for V456D mutation:", chi_V456N)


