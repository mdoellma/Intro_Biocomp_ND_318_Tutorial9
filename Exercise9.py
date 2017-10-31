#import os and packages
import os
import fileinput
import string
import sys
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.stats import chi2
from scipy.optimize import minimize
from scipy.stats import norm
from plotnine import *

#set working directory
os.chdir('/Users/omneelay/Desktop/ExerciseNine/Intro_Biocomp_ND_318_Tutorial9/')

###PART 1

#Custom likelihood functions null and alt
def nllikealt(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    
    expected=B0+B1*obs.x
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll
    
def nllikenull(p,obs):
    B0=p[0]
    sigma=p[1]
    
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll
    
initialGuessnull=np.array([1,1])
initialGuessalt=np.array([1,1,1])

data=pd.read_csv("ponzr1.csv")

##Creating data frames for WT and each mutation where the columns are x and y. x is 0 for wild type and 1 for mutations. The y values are the ponzr1 counts.
#WT data frame
WT=data.loc[0:9]
WT=WT.ponzr1Counts
WT=WT.values.tolist()
a=[0,0,0,0,0,0,0,0,0,0]
WT=pd.DataFrame({'x': a, 'y': WT})

#M124K data frame
M124K=data.loc[10:19]
M124K=M124K.ponzr1Counts
M124K=M124K.values.tolist()
a=[1,1,1,1,1,1,1,1,1,1]
M124Kdata=pd.DataFrame({'x': a, 'y': M124K})
M124Kdata=M124Kdata.append(WT)

#V456D data frame
V456D=data.loc[20:29]
V456D=V456D.ponzr1Counts
V456D=V456D.values.tolist()
a=[1,1,1,1,1,1,1,1,1,1]
V456Ddata=pd.DataFrame({'x': a, 'y': V456D})
V456Ddata=V456Ddata.append(WT)

#I213N data frame
I213N=data.loc[30:39]
I213N=I213N.ponzr1Counts
I213N=I213N.values.tolist()
a=[1,1,1,1,1,1,1,1,1,1]
I213Ndata=pd.DataFrame({'x': a, 'y': I213N})
I213Ndata=I213Ndata.append(WT)

##minimizing the nll function for 3 mutations null and alt. And calculating D:
#M124K NLL minimize
fit=minimize(nllikenull,initialGuessnull,method="Nelder-Mead",options={'disp': True},args=M124Kdata)
M124K_null=fit.fun
fit=minimize(nllikealt,initialGuessalt,method="Nelder-Mead",options={'disp': True},args=M124Kdata)
M124K_alt=fit.fun
M124K_D=-2*(M124K_alt-M124K_null)

#V456D NLL minimize
fit=minimize(nllikenull,initialGuessnull,method="Nelder-Mead",options={'disp': True},args=V456Ddata)
V456D_null=fit.fun
fit=minimize(nllikealt,initialGuessalt,method="Nelder-Mead",options={'disp': True},args=V456Ddata)
V456D_alt=fit.fun
V456D_D=-2*(V456D_alt-V456D_null)

#I213N NLL minimize
fit=minimize(nllikenull,initialGuessnull,method="Nelder-Mead",options={'disp': True},args=I213Ndata)
I213N_null=fit.fun
fit=minimize(nllikealt,initialGuessalt,method="Nelder-Mead",options={'disp': True},args=I213Ndata)
I213N_alt=fit.fun
I213N_D=-2*(I213N_alt-I213N_null)

##Test for statistical significance for each mutation
print("Part 1:")
print("p-value M124K: ",1-scipy.stats.chi2.cdf(x=M124K_D,df=1)," NOT SIGNIFICANT")
print("p-value V456D: ",1-scipy.stats.chi2.cdf(x=V456D_D,df=1)," SIGNIFICANT")
print("p-value I213N: ",1-scipy.stats.chi2.cdf(x=I213N_D,df=1)," NOT SIGNIFICANT")


### PART 2

#Define Function
def nllikebacteria(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    
    expected=B0*((obs.S)/(obs.S+B1))
    nll=-1*norm(expected,sigma).logpdf(obs.u).sum()
    return nll
#read in data
data2=pd.read_csv("MmarinumGrowth.csv")
initialGuess=np.array([1,1,1])
#solve function
fit=minimize(nllikebacteria,initialGuess,method="Nelder-Mead",options={'disp': True},args=data2)
print("Part 2:")
print("Max Growth Rate=",fit.x[0])
print("Half Saturation Constant=",fit.x[1])



### PART 3

#Define Functions
def nllikequadratic(p,obs):
    B0=p[0]
    B1=p[1]
    B2=p[2]
    sigma=p[3]
    
    expected=B0+B1*obs.Ms+B2*obs.Ms*obs.Ms
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

def nllikelinear(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    
    expected=B0+B1*obs.Ms
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll
    
def nllikeconstant(p,obs):
    B0=p[0]
    sigma=p[1]
    
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll
    
initialGuessquadratic=np.array([1,1,1,1])
initialGuesslinear=np.array([1,1,1])
initialGuessconstant=np.array([1,1])

#read data
soildata=pd.read_csv("leafDecomp.csv")

#solve functions
fitquadratic=minimize(nllikequadratic,initialGuessquadratic,method="Nelder-Mead",options={'disp': True},args=soildata)
fitlinear=minimize(nllikelinear,initialGuesslinear,method="Nelder-Mead",options={'disp': True},args=soildata)
fitconstant=minimize(nllikeconstant,initialGuessconstant,method="Nelder-Mead",options={'disp': True},args=soildata)
1-scipy.stats.chi2.cdf(x=-2*(fitquadratic.fun-fitconstant.fun),df=2)
1-scipy.stats.chi2.cdf(x=-2*(fitlinear.fun-fitconstant.fun),df=1)
1-scipy.stats.chi2.cdf(x=-2*(fitquadratic.fun-fitlinear.fun),df=1)

print("Part 3:")
print("Constant Model: ","B0=",fitconstant.x[0], ", ", "sigma=",fitconstant.x[1])
print("Linear Model: ","B0=",fitlinear.x[0], ", ", "B1=",fitlinear.x[1], ", ", "sigma=", fitlinear.x[2])
print("Quadratic Model: ","B0=",fitquadratic.x[0], ", ", "B1=",fitquadratic.x[1], ", ", "B2=", fitquadratic.x[2], ", ", "sigma=", fitquadratic.x[3])
print("Quadratic is best. Linear better than Constant. p-values for comparisons are basically 0.0")



