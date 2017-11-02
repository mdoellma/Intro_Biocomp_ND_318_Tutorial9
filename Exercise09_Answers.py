"""
Exercise 9
Authors: Grant Keller and Kathleen Nicholson
"""

from numpy import array
import pandas
from scipy.optimize import minimize
from scipy.stats import norm, chi2

METHOD = 'Nelder-Mead'

def nllike(p, obs, x='x', y='y', model='p[0]'):
    """
    Returns negative log likelihood of model (default constant rate)
        with parameters p to describe the 2D data obs with ind. and dep.
        variables x and y.
    """
    return -1 * norm(eval(model), p[-1]).logpdf(obs[y]).sum()

def test_null(data, control, treat):
    """
    Returns the chi2 p value of the null hypothesis that the segment of
        data represented by control is not significantly different
        than the segment represented by treat.
    """
    control_x = [0 for v in data.mutation if v == control]
    control_y = [data.iloc[i][1] for i in range(len(data)) if data.iloc[i][0] == control]
    treat_x = [1 for v in data.mutation if v == treat]
    treat_y = [data.iloc[i][1] for i in range(len(data)) if data.iloc[i][0] == treat]

    df = pandas.DataFrame({'x': control_x + treat_x,
                           'y': control_y + treat_y})

    model1 = 'p[0] + p[1] * obs[x]' # model treating two conditions (exp. vs. control) differently
    args1 = (df, 'x', 'y', model1)

    nullfit = minimize(nllike, array([1, 1]), method=METHOD, args=(df))
    treatfit = minimize(nllike, array([1, 1, 1]), method=METHOD, args=args1)
    diff = abs(nllike(treatfit.x, df, model=model1) - nllike(nullfit.x, df))
    return 1 - chi2.cdf(x=2*diff, df=1)

def q1solution():
    """
    This code will perform a likelihood ratio test to determine which mutation
        of three in the ponzr1 gene causes a significant reduction in mRNA levels
    """
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
        str1 = "The {0} mutation {1} ponzr1 expression".format(condition, msg)
        print str1 + " (p-value: {:.2})".format(p_val)

def q2solution():
    """
    This code will determine the maximum growth rate and half-maximal growth
        concentration of cell growth data.
    """

    #reads csv file
    growth_data = pandas.read_csv('MmarinumGrowth.csv')
    #creates the model we are interested in
    model = 'p[0] * obs[x] / (obs[x] + p[1])'
    fit = minimize(nllike, array([1, 1, 1]), method=METHOD, args=(growth_data, 'S', 'u', model))

    print "The maximum growth rate is: {:.6}.".format(fit.x[0])
    print "The half maximal growth concentration is: {}.".format(int(round(fit.x[1])))

def q3solution():
    """
    This code will create three models for one set of data and determine the p-values
        when comparing them.
    """
    #reads csv file
    decomposition = pandas.read_csv('leafDecomp.csv')
    #This estimates the negative log likelihood.
    #   Assumes constant rate, linear or quadratic model.
    models = ['p[0]+p[1]*obs[x]',
              'p[0]+p[1]*obs[x]+p[2]*obs[x]**2']

    args1 = (decomposition, 'Ms', 'decomp')
    fit1 = minimize(nllike, array([1, 1]), method=METHOD, args=args1)
    print fit1.x
    #This estimates the negative log likelihood using a linear model.

    args2 = (decomposition, 'Ms', 'decomp', models[0])
    fit2 = minimize(nllike, array([1, 1, 1]), method=METHOD, args=args2)
    print fit2.x
    #This estimates the negative log likelihood using a quadratic model.

    args3 = (decomposition, 'Ms', 'decomp', models[1])
    fit3 = minimize(nllike, array([200, 10, -.02, 1]), method=METHOD, args=args3)
    print fit3.x

    #This saves the nll functions as variables to be called during the comparison of models below.
    nllalt = fit1.fun
    nllnul = fit2.fun
    nllnul2 = fit3.fun

    #This compares the constant rate model and linear model and returns a p-value.
    pval = 1 - chi2.cdf(x=2*(nllalt-nllnul), df=1)
    print "The p-value when comparing the constant rate model and linear model is %s." %pval
    #This compares the linear model and quadratic model and returns a p-value.
    pval = 1 - chi2.cdf(x=2*(nllnul-nllnul2), df=1)
    print "The p-value when comparing the linear model and quadratic model is %s." %pval
    #This compares constant rate model and quadratic model and returns a p-value.
    pval = 1-chi2.cdf(x=2*(nllalt-nllnul2), df=2)
    print "The p-value when comparing the constant rate model and quadratic model is %s." %pval
    print "The quadratic model best fits the shape of this data."
if __name__ == '__main__':
    print "Question 1:"
    q1solution()
    print "Question 2:"
    q2solution()
    print "Question 3:"
    q3solution()
