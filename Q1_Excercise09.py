#Exercise 9
#Authors: Janice Love and Melissa Stephens

#Load data and import packages
import numpy
import pandas
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.stats import chi2
from plotnine import *

#need to read-in/populate dataframe with existing data
MyDataFrame=pandas.read_csv("ponzr1.csv",sep=",",header=0)

#make 3 data frames with the control and treatment we are interested in. Initialize x=0
#designate treatment group as x=1, and x=0 for control group

Var1=MyDataFrame.loc[MyDataFrame.mutation.isin(["WT","M124K"]),:]
Var2=pandas.DataFrame({'y':Var1.ponzr1Counts,'x':0})
Var2.loc[Var1.mutation=="M124K","x"]=1

Var3=MyDataFrame.loc[MyDataFrame.mutation.isin(["WT","V456D"]),:]
Var4=pandas.DataFrame({'y':Var3.ponzr1Counts,'x':0})
Var4.loc[Var3.mutation=="V456D","x"]=1

Var5=MyDataFrame.loc[MyDataFrame.mutation.isin(["WT","I213N"]),:]
Var6=pandas.DataFrame({'y':Var5.ponzr1Counts,'x':0})
Var6.loc[Var5.mutation=="I213N","x"]=1

#p is your guesses. One for each parameter (B0,B1,alpha). Use mean of control group for B0 and mean of control-mean of B1 group for B1. 1 for alpha
#obs refers to my data frame
#need to define 2 custom functions. One for our T-test and one for the Null model
def nllike(p,obs):
	B0=p[0]
	B1=p[1]
	sigma=p[2]
	expected=B0+B1*obs.x
	nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
	return nll

def nllikeNull(p,obs):
        B0=p[0]
        expected=B0
        sigma=p[1]
	nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
        return nll

#need to use the minimize function to check wheither the parameters(p) we guessed/passed are true
#need a fit function for each comparison (ie WT vs M, WT vs V, WT vs I) changing only the function name for each iteration and changing the dataframe for each comparison
#makes guesses as to the best "p" one for each parameter
guessTest=[2395,56,1]
guessNull=[2395,1]
fit1=minimize(nllike,guessTest,method="Nelder-Mead",options={'disp':True},args=Var2)
fit2=minimize(nllikeNull,guessNull,method="Nelder-Mead",options={'disp':True},args=Var2)
fit3=minimize(nllike,guessTest,method="Nelder-Mead",options={'disp':True},args=Var4)
fit4=minimize(nllikeNull,guessNull,method="Nelder-Mead",options={'disp':True},args=Var4)
fit5=minimize(nllike,guessTest,method="Nelder-Mead",options={'disp':True},args=Var6)
fit6=minimize(nllikeNull,guessNull,method="Nelder-Mead",options={'disp':True},args=Var6)

#Optional- print both for the nll alt model and the nll null model. This atrribues 'x' to a list of the most likely parameter values
##print(fit1.x)#parameter estimates for t-test for Var2
##print(fit2.x)#parameter estimates for null model for Var2 etc ..

#need to output the nll for each of the null model and the t-test model 
print(fit1.fun)
print(fit2.fun)
print(fit3.fun)
print(fit4.fun)
print(fit5.fun)
print(fit6.fun)
	
#likelihood ratio test. Need one for each of the 3 models for Alt Test vs Null
##D=2*(nllAlt-nllNuull)
##pval=1-chi2.cdf(x=2*(),df=1)

#pvalWTvM
print(1-chi2.cdf(x=2*(fit2.fun-fit1.fun),df=1))
#pvalWTvV
print(1-chi2.cdf(x=2*(fit4.fun-fit3.fun),df=1))
#pvalWTvI
print(1-chi2.cdf(x=2*(fit6.fun-fit5.fun),df=1))


