#set working directory
os.chdir('/Users/brittnibertolet/Desktop/bcTutorials/Intro_Biocomp_ND_318_Tutorial9/')

#add packages
import numpy
import pandas
import re
import scipy
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.stats import chi2
from plotnine import *

#1
#add data 'ponzr1.csv'
mrna=pandas.read_csv("ponzr1.csv",header=0)
#check out data
mrna.head()
#visualize data
ggplot(mrna, aes(x="mutation", y="ponzr1Counts"))+geom_point()+theme_classic()

#write likelihood function for null model
def null(p,obs):
    B0=p[0]
    sigma=p[1]
    
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll
    
#write likelihood function for model where a mutation effects expression
def mut(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    
    expected=B0+B1*obs.x
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll
################# First t-test between control and M124K mutation
#data for control vs first mutation (M124K)
data1=mrna.loc[(mrna.mutation == "WT") | (mrna.mutation == "M124K")]
data1.columns=['x', 'y']
data1['x'] = data1['x'].map({'WT': 0, 'M124K': 1})

#estimate parameters with null model
initialGuess=numpy.array([2000,1])
null_fit=minimize(null, initialGuess, method="Nelder-Mead", options={'disp': True}, args=data1)
print(null_fit.x) #print parameters
print(null_fit.fun) #print negative log likelihood

#estimate parameters with mutation model 
initialGuess=numpy.array([2000,1000, 1])
mut_fit=minimize(mut, initialGuess, method="Nelder-Mead", options={'disp': True}, args=data1)
print(mut_fit.x) #print parameters
print(mut_fit.fun) #print negative log likelihood

#calculate the different in negative log likelihood
D=2*(null_fit.fun-mut_fit.fun)
#test for statistical significance
1-scipy.stats.chi2.cdf(x=D,df=1)

#2 
#add data 'MmarinumGrowth.csv'
mar=pandas.read_csv("MmarinumGrowth.csv",header=0)
#check out data
mar.head()
#visualize data
ggplot(mar,aes(x='S',y='u'))+geom_point()+theme_classic()

#3
#add data 'leafDecomp.csv'
leaf=pandas.read_csv("leafDecomp.csv",header=0)
#check out data
leaf.head()
#visualize data
ggplot(leaf,aes(x='Ms',y='decomp'))+geom_point()+theme_classic()


