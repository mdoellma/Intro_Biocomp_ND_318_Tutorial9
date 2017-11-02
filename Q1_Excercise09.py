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
def nllike(guess=[2395,56,1],Var2):
	B0=p[0]
	B1=p[1]
	sigma=p[2]
	expected=B0+B1*Var2.x
	nll=-1*norm(expected,sigma).logpdf(Var2.y).sum()
	return nll

def nllikeNull(guess=2395],obs):
        B0=p[0]
        expected=B0
        nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
        return nll

#need to use the minimize function to check wheither the parameters(p) we guessed/passed are true
#need a fit function for each comparison (ie WT vs M, WT vs V, WT vs I) changing only the function name for each iteration and changing the dataframe for each comparison

fit=minimize(function,initialVals,method="Nelder-Mead",options={'disp':True},args=observedData)
fit=minimize(function,initialVals,method="Nelder-Mead",options={'disp':True},args=observedData)
fit=minimize(function,initialVals,method="Nelder-Mead",options={'disp':True},args=observedData)
fit=minimize(function,initialVals,method="Nelder-Mead",options={'disp':True},args=observedData)
fit=minimize(function,initialVals,method="Nelder-Mead",options={'disp':True},args=observedData)
fit=minimize(function,initialVals,method="Nelder-Mead",options={'disp':True},args=observedData)

	
#likelihood ratio test
pval=1-chi2.cdf(x=2*(),df=1)
print pval
