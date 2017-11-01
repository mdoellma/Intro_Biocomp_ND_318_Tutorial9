# Exercise 9 
# Author: Grant Keller and Kathleen Nicholson

import numpy
import pandas
from scipy.optimize import minimize
from scipy.stats import norm, chi2
from plotnine import *

def nllike(p, obs, model, x, y):
    return -1 * norm(eval(model), p[2]).logpdf(obs[y]).sum()

def test_null(data, control, treat):
    guess = numpy.array([1, 1, 1])
    control_x = [0 for v in data.mutation if v == control]
    control_y = [data.iloc[i][1] for i in range(len(data)) if data.iloc[i][0] == control]
    treat_x = [1 for v in data.mutation if v == treat]
    treat_y = [data.iloc[i][1] for i in range(len(data)) if data.iloc[i][0] == treat]
    
    d = {'x': control_x + treat_x,
         'y': control_y + treat_y}
    df = pandas.DataFrame(d)

    model0 = 'p[0]' # null hypothesis - same model describes control & experiment well
    args0 = (df, model0, 'x', 'y')
    
    model1 = 'p[0] + p[1] * obs[x]' # model treating two conditions (exp. vs. control) differently
    args1 = (df, model1, 'x', 'y')
    
    nullFit = minimize(nllike, guess, method="Nelder-Mead", args=args0)
    treatFit = minimize(nllike, guess, method="Nelder-Mead", args=args1)
    diff = abs(nllike(treatFit.x, df, model1, 'x', 'y') - nllike(nullFit.x, df, model0, 'x', 'y'))
    return 1 - chi2.cdf(x=2*diff, df=1)

#Q1#######################################################
## This code will perform a likelihood ratio test to determine which mutation 
##      of three in the ponzr1 gene causes a significant reduction in mRNA levels
data = pandas.read_csv('ponzr1.csv')
control = 'WT'
muts = []
#creates lists of data for future comparison
for item in data.mutation.unique():
    if item != control:
        muts.append(item)

#runs tests for p-values and determines if mutations significantly affected expression
for condition in muts:
    p_val = test_null(data, control, condition)
    if p_val <= 0.05:
        msg = "significantly affected"
    else:
        msg = "did not significantly affect"
    print "The {0} mutation {1} ponzr1 expression (p-value: {2:.2})".format(condition, msg, p_val)

###########################################################

#Q2#######################################################
#This code will determine the maximum growth rate and half-maximal growth 
## concentration of cell growth data. 

#reads csv file
growth_data = pandas.read_csv('MmarinumGrowth.csv')
#creates the model we are interested in
model = 'p[0] * obs[x] / (obs[x] + p[1])'
guess = numpy.array([1,1,1])
fit = minimize(nllike, guess, method="Nelder-Mead", options={'disp': True}, args=(growth_data, model, 'S', 'u'))

print "The maximum growth rate is: {:.6}.".format(fit.x[0])
print "The half maximal growth concentration is: {}.".format(int(round(fit.x[1])))

###########################################################
#This code will create three models for one set of data and determine the p-values
# when comparing them. 

#Q3
#reads csv file
decomposition = pandas.read_csv('leafDecomp.csv')
#This estimates the negative log likelihood using a constant rate model.
def nllike(p,obs): 
        B0=p[0]
        sigma=p[1]
        expected=B0
        nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
        return nll
Guess=numpy.array([1,1])
fit=minimize(nllike,Guess,method="Nelder-Mead",options={'disp': True},args=decomposition)
print(fit.x)
#This estimates the negative log likelihood using a linear model.
def nllike(p,obs): 
        B0=p[0]
        B1=p[1]
        sigma=p[2]
        expected=B0+B1*obs.Ms
        nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
        return nll
Guess=numpy.array([1,1,1])
fit2=minimize(nllike,Guess,method="Nelder-Mead",options={'disp': True},args=decomposition)
print(fit2.x)
#This estimates the negative log likelihood using a quadratic model.
def nllike(p,obs): 
        B0=p[0]
        B1=p[1]
        B2=p[2]
        sigma=p[3]
        expected=B0+B1*obs.Ms+B2*obs.Ms*obs.Ms
        nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
        return nll
Guess=numpy.array([200,10,-.02,1])
fit3=minimize(nllike,Guess,method="Nelder-Mead",options={'disp': True},args=decomposition)
print(fit3.x)

#This saves the nll functions as variables to be called during the comparison of models below. 
nllalt=fit.fun
nllnul=fit2.fun
nllnul2=fit3.fun

#This compares the constant rate model and linear model and returns a p-value.
pval= 1 - chi2.cdf(x=2*(nllalt-nllnul), df=1)
print "The p-value when comparing the constant rate model and linear model is {0}".format(pval)
#This compares the linear model and quadratic model and returns a p-value.
pval=1 - chi2.cdf(x=2*(nllnul-nllnul2), df=1)
print "The p-value when comparing the linear model and quadratic model is {0}".format(pval)
#This compares constant rate model and quadratic model and returns a p-value.
pval=1-chi2.cdf(x=2*(nllalt-nllnul2), df=2)
print "The p-value when comparing the constant rate model and quadratic model is {0}".format(pval)
"""

"""
def main():
    pass

# if __name__ == '__main__':
    