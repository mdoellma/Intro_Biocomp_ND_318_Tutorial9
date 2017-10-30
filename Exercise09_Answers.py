# Exercise 9 
# Author: Grant Keller and Kathleen Nicholson
# Import Packages

import numpy
import pandas
from scipy.optimize import minimize
from scipy.stats import norm
from plotnine import *

def nllike(p, obs, x, y):
    expected = p[0] + p[1] * obs[x]
    return -1 * norm(expected, p[2]).logpdf(obs[y]).sum()

# Q1
## This code will perform a likelihood ratio test to determine which mutation 
##      of three in the ponzr1 gene causes a significant reduction in mRNA levels
"""
1. In contrast to the study we read this week (Kelly et al. 2014), 
many laboratory experiments are designed to test the speciﬁc eﬀect 
of a treatment with all other variables controlled. 

These experiments will have replicated treatment and control sets 
of experimental units.
a t-test is used to evaluate whether there is a statistically 
signiﬁcant diﬀerence in the mean of measurements taken from the
treatment and control experimental units. 

In this approach, the likelihood of a null model with a single
parameter describing the mean behavior of all experimental units
(y=B0+error) is compared to the likelihood of a model with an
additional parameter that describes the diﬀerence between treatment
and control groups (y=B0+B1*treat+error). 

If the diﬀerence in likelihood between the two models is large
enough, we can say that our treatment had a statistically signiﬁcant
eﬀect. 

We test for statistical signiﬁcance by asking whether two times
the diﬀerence in likelihoods (D) is large relative to a chi-squared
distribution with one degree of freedom. 

This test can be accomplished with 1-scipy.stats.chi2.cdf(x=D,df=1) in Python.

Use a likelihood ratio test to determine which of three mutations
signiﬁcantly reduced the expression of ponzr1, a gene involved in
the formation of the glomerulus in a developing zebraﬁsh kidney.

ANSWERS

M124K: p-value ~ 0.72 (no effect of treatment)
V456D: p-value ~ 5.6e-6 (effect of treatment)
I213N: p-value ~ 0.88 (no effect of treatment)
"""
# this isn't correct, but it's a (working) start
data = pandas.read_csv('ponzr1.csv')
data['x1'] = [0 for v in data.mutation]
data['x2'] = [['WT', 'M124K', 'V456D', 'I213N'].index(v) for v in data.mutation]
data['y'] = data.ponzr1Counts
initGuess=numpy.array([1,1,1])
fit0 = minimize(nllike, initGuess, method="Nelder-Mead", options={'disp': True}, args=(data,'x1','y'))
fit1 = minimize(nllike, initGuess, method="Nelder-Mead", options={'disp': True}, args=(data,'x2','y'))
D1 = nllike(fit0.x, data)
D2 = nllike(fit1.x, data)

diff = abs(D2 - D1)
1-scipy.stats.chi2.cdf(x=2*diff, df=1)


#           evaluation of function
#                   fit nll to data 
#                   Likelihood ratio test, graph.

#Q2

def nlllike_growth(p, obs, x, y):
    expected = p[0] * obs[x] / (obs[x] + p[1])
    return -1 * norm(expected, p[2]).logpdf(obs[y]).sum()

growth_data = pandas.read_csv('MmarinumGrowth.csv')
fit = minimize(nllike2, initGuess, method="Nelder-Mead", options={'disp': True}, args=(growth_data,'S','u'))

print "The maximum growth rate is: {:.6}.".format(fit.x[0])
print "The half maximal growth concentration is: {}.".format(int(round(fit.x[1])))
print "Std. dev. (sigma) is: {:.1}.".format(fit.x[2])

#Q3

"""
3. As you’ll see next week, when recreating biological processes in
simulation models we often make strong simplifying assumptions.

The degree to which we simplify the representation of a biological
process in a simulation model is often guided by observational or
experimental data. One way to quantify the simplicity, or conversely
complexity, of a simulation model is the number of parameters used.

Imagine we want to model decomposition of leaves (d) in the soils of
a forest as a function of water availability (measured as soil
moisture, Ms).

For this aspect of our simulation model we want to keep the model as
simple (i.e. fewest parameters) as possible, but also want to capture
the “shape” of the relationship between decomposition and soil
moisture.

Using the observations in “leafDecomp.csv” determine whether we
should use a constant rate for all soil moistures (d = a), a linear
response of decomposition to soil moisture (d = a + bMs), or a
hump-shaped response of decomposition to soil moisture
(d = a+bMs +cM2 s). Hint: You can use the likelihood ratio test to
compare these models too. This is because the likelihood ratio test
can be used to compare any set of models that are subsets of each
other. You know one model is a subset of another if you can set a
parameter equal to zero or some non-zero constant and have the same
equation for both models. Also note that the degrees of freedom used
in the chi-squared distribution is determined by the diﬀerence in
the number of parameters between the two models.

In the t-test example above, that will always be 1, but in this case
it will be 1 or 2 depending on which models you are comparing.

ANSWERS:
Constant fit: B0 ~ 589.7, sigma ~ 164
Linear fit: B0 ~ 318, B1 ~ 6.3, sigma ~ 54
Quadratic fit: B0 ~ 180, B1 ~ 15.7, B2 ~ -0.11, sigma ~ 10.7
The quadratic model is by for the best and the linear model is 
a lot better than the constant or null model. The p-values for 
the likelihood ratio tests should all be essentially zero (~1e-20 
or less).
"""
