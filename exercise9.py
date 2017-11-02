#set working directory
os.chdir('/Users/Chloe/Desktop/data-shell/Intro_Biocomp_ND_318_Tutorial9/')

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

#calculate the difference in negative log likelihood
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
mar.columns=['x', 'y']

### Custom likelihood function 
def nllike(p,obs):
    umax=p[0]
    Ks=p[1]
    sigma=p[2]
    
    expected=umax*(obs.x/(obs.x-Ks))
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll

### estimate parameters by minimizing the negative log likelihood
mar.columns=['x', 'y']
initialGuess=numpy.array([1,1,1])
fit=minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp': True},args=mar)

print(fit.x)

#3
#add data 'leafDecomp.csv'
leaf=pandas.read_csv("leafDecomp.csv",header=0)
#check out data
leaf.head()
#visualize data
ggplot(leaf,aes(x='Ms',y='decomp'))+geom_point()+theme_classic()
#make new dataframe with x and y as the headers
leaves=leaf
leaves.columns=['x', 'y']
leaves.head()
#define custom liklihood function for constant decomp
def constant(p,obs):
    B0=p[0]
    sigma=p[1]
    
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll
#set intial guesses
constantguess=numpy.array([600,1])
#estimate parameters
constant_fit=minimize(constant,constantguess,method="Nelder-Mead",options={'disp':True},args=leaves)
print(constant_fit.x)
print(constant_fit.fun)#nll
#define custom liklihood function for linear decomp
def linear(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    
    expected=B0+B1*obs.x
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll
#set intial guesses
linearguess=numpy.array([10,6,1])
#estimate parameters
linear_fit=minimize(linear,linearguess,method="Nelder-Mead",options={'disp':True},args=leaves)
print(linear_fit.x)
print(linear_fit.fun)#nll
#define custom liklihood function for hump-shaped decomp
def hump(p,obs):
    B0=p[0]
    B1=p[1]
    B2=p[2]
    sigma=p[3]
    
    expected=B0+B1*obs.x+B2*(obs.x)*(obs.x)
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll
#set initial guesses
humpguess=numpy.array([200,10,-.2,1])
#estimate parameters
hump_fit=minimize(hump,humpguess,method="Nelder-Mead",options={'disp':True},args=leaves)
print(hump_fit.x)
print(hump_fit.fun)#nll

#calculate the difference in negative log likelihood constant vs linear
first_D=2*(constant_fit.fun-linear_fit.fun)
#test for statistical significance
1-scipy.stats.chi2.cdf(x=first_D,df=1)

#calculate the difference in negative log likelihood linear vs hump
second_D=2*(linear_fit.fun-hump_fit.fun)
#test for statistical significance
1-scipy.stats.chi2.cdf(x=second_D,df=1)

#calculate the difference in negative log likelihood constant vs hump
third_D=2*(constant_fit.fun-hump_fit.fun)
#test for statistical significance
1-scipy.stats.chi2.cdf(x=third_D,df=2)

#plot
ggplot(leaves,aes(x='x',y='y'))+geom_point()+theme_classic()+geom_line(aes(y=constant_fit.x[0]))+geom_line(aes(x='x',y=linear_fit.x[0]+linear_fit.x[1]*leaves.x))+geom_line(aes(x='x',y=hump_fit.x[0]+hump_fit.x[1]*leaves.x+hump_fit.x[2]*((leaves.x)*(leaves.x))))


