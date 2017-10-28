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


infile = open ("ponzr1.csv", 'r')

#---------------------------------------------------------------------------
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

df = pd.read_csv("ponzr1.csv", delimiter= ',')

WT = []
M124K = []
I213N = [] 
V456D = []

#place gene counts into respective list: WT or specific mutation 
for i in range(0,len(df),1):
    if df.mutation[i]=="WT":
        WT.append(df.ponzr1Counts[i])
    elif df.mutation[i] == "M124K":
        M124K.append(df.ponzr1Counts[i])
    elif df.mutation[i] == "I213N":
        I213N.append(df.ponzr1Counts[i])
    elif df.mutation[i] == "V456D":
        V456D.append(df.ponzr1Counts[i])

#create df with conditions to be compared 
WT_df = pd.DataFrame({'WT':WT})
WT_df['M124K']= M124K

WT_I213N = pd.DataFrame({'WT':WT})
WT_I213N['I213N']= I213N

WT_V456D = pd.DataFrame({'WT':WT})
WT_V456D['V456D']= V456D

#For WT_I213N ----------------------------------------------------------------------------------------
def nllike(p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    
    expected = B0 + B1*WT_I213N.I213N #make sure data import has a column called x 
    
    nll=-1*norm(expected,sigma).logpdf(WT_I213N.WT).sum()
    return nll 
    
initialGuess=np.array([1,1,1])
fit=minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp': True},args=df)    
print (fit.x)
I213N_fit = fit.x
#------------------------------------------------------------------------------------------------------

#For WT_V456D ----------------------------------------------------------------------------------------
def nllike(p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    
    expected = B0 + B1*WT_V456D.V456D #make sure data import has a column called x 
    
    nll=-1*norm(expected,sigma).logpdf(WT_V456D.WT).sum()
    return nll 
    
initialGuess=np.array([1,1,1])
fit=minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp': True},args=df)    
print (fit.x)
V456D_fit = fit.x
#------------------------------------------------------------------------------------------------------


#for WT 
def nullnull(p, obs):
    B0 = p[0]
    sigma = p[1]
    
    expected = B0 + WT_df.WT #make sure data import has a column called x 
    
    nll=-1*norm(expected,sigma).logpdf(WT_df.WT).sum()
    return nll 
    
initialGuess=np.array([1,1,1])
fit=minimize(nullnull,initialGuess,method="Nelder-Mead",options={'disp': True},args=df)    
print (fit.x)

WT_null = fit.x
#------------------------------------------------------------------------------------------------------

#For WT_M124K -----------------------------------------------------------------------------------------
def nllike(p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    
    expected = B0 + B1*WT_df.M124K #make sure data import has a column called x 
    
    nll=-1*norm(expected,sigma).logpdf(WT_df.WT).sum()
    return nll 
    
initialGuess=np.array([1,1,1])
fit=minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp': True},args=df)    
print (fit.x)

WT_M124K = fit.x 
#--------------------------------------------------------------------------------------------------------

chi_M124K = 1-scipy.stats.chi2.cdf(x=D,df=1)

diffInNumParam = 1

M124K_pval = 1-chi2.cdf(x=2*(WT_M124K-WT_null), df = diffInNumParam)
V456D_pval = 1-chi2.cdf(x=2*(V456D_fit-WT_null), df = diffInNumParam)

