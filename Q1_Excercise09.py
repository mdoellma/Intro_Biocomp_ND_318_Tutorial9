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

def nllike(p,obs):
	B0=p[0]
	B1=p[1]
	sigma=p[2]
	expected=B0+B1*obs.x
	nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
	return nll
	
#likelihood ratio test
pval=1-chi2.cdf(x=2*(),df=1)
print pval
