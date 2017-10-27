#Exercise 9 Q1 
#Authors: Janice Love and Melissa Stephens 
#Date: October 27, 2017 

#T-test and maximum liklihood ratio scripts 

import pandas as pd 
import numpy as np 
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
M124K - []
I213N = [] 

for i in df:
    if mutation[i]=="WT":
        WT.append()
        
    
    