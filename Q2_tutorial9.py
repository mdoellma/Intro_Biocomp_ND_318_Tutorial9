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
MyDataFrame=pandas.read_csv("MmarinumGrowth.csv",sep=",",header=0)

#the model
##Monod Equation-u=umax*S/(S=Ks)

#write a custom function for the growth model
#How many parameters does the model have? 3. umax, Ks, and sigma (need to include sigma because it is a normal distribution, and sigma will always be included)
#assign starting guesses
#fit the model
#plot results
#a custom function below. Change parameters to fit the model we are given, using obs.x and obs.y to call on the column names of the csv input file.
def nllike(p,obs):
	umax=p[0]
	Ks=p[1]
	sigma=p[2]
	expected=umax*obs.S/(obs.S+Ks)
	nll=-1*norm(expected,sigma).logpdf(obs.u).sum()
	return nll
#starting guess
guess=[1,1,1]

#fit the model
fitModel=minimize(nllike,guess,method="Nelder-Mead",options={'disp':True},args=MyDataFrame)

#print
print(fitModel.x)
print(fitModel.fun)
#plot results
ggplot(MyDataFrame,aes(x='u',y='S'))+geom_point()+theme_classic()


