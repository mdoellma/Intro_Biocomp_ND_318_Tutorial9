# Exercise 9 
# Author: Grant Keller and Kathleen Nicholson
# Import Packages

import numpy
import pandas
from scipy.optimize import minimize
from scipy.stats import norm
from plotnine import *


# Q1
## This code will perform a likelihood ratio test to determine which mutation 
##      of three in the ponzr1 gene causes a significant reduction in mRNA levels
#           setup
#               packages: look at code from wednesday's class to find this.
#               read data file
#           function (maximum likelihood)
#               def nllike(arguments):
#                   "unpack" arguments, assign variables
#                   calculate expected value (model equation). This defines model.
#                   calculate negative log likelihood (result=nll calc)
#                   return negative log likelihood
#               def nllike(p,obs): 
#                   B0=p[0]
#                   B1=p[1]
#                   sigma=p[2]
#                   expected=B0+B1*obs.x
#                   nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
#                   return nll
#           evaluation of function
#                   fit nll to data 
#                   Likelihood ratio test, graph.



#Q2




#Q3