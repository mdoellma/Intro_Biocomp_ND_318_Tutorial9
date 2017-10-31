#Exercise 09

#Q1 Zebrafish

#load file
import pandas
file=pandas.read_csv("ponzr1.csv",header=0,sep=",")

#subset each different mutation values vs. WT values
mut1=file.loc[file.mutation.isin(['WT', 'M124K']),:]
mut2=file.loc[file.mutation.isin(['WT', 'V456D']),:]
mut3=file.loc[file.mutation.isin(['WT', 'I213N']),:]

#mutation 1 group into a new dataframe and changes the x column to 0's and 1's
mut1a=pandas.DataFrame({'y':mut1.ponzr1Counts, 'x':0})
mut1a.loc[mut1.mutation=='M124K', 'x']=1

#mutation 2 group into a new dataframe and changes the x column to 0's and 1's
mut2a=pandas.DataFrame({'y':mut2.ponzr1Counts, 'x':0})
mut2a.loc[mut2.mutation=='V456D', 'x']=1

#mutation 3 group into a new dataframe and changes the x column to 0's and 1's
mut3a=pandas.DataFrame({'y':mut3.ponzr1Counts, 'x':0})
mut3a.loc[mut3.mutation=='I213N', 'x']=1

#import packages
import numpy
import pandas
from scipy.optimize import minimize
from scipy.stats import norm
from plotnine import *

#plot the values for WT (0) & mutation (1)
ggplot(mut1a,aes(x='x',y='y'))+geom_point()+theme_classic()

#Null hypothesis likelihood equation
def nllike(p,obs):
    B0=p[0]
    sigma=p[1]
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll

#Alternative hypothesis likelihood equation
def nllike2(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    expected=B0+B1*obs.x
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll

#estimating parameters by minimizing the nll 
initialVals1=numpy.array([1,1,1])

#for mutation M124K
fitNull=minimize(nllike,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut1a)
fitAlt=minimize(nllike2,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut1a)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitNull.fun-fitAlt.fun))
mut1aanswer=(1-chi2.cdf(x=D,df=1))
print('M124K p value')
print(mut1aanswer)


#for mutation V456D
fitNull=minimize(nllike,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut2a)
fitAlt=minimize(nllike2,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut2a)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitNull.fun-fitAlt.fun))
mut2aanswer=1-chi2.cdf(x=D,df=1)
print('V456D p value')
print(mut2aanswer)

#for mutation I213N

fitNull=minimize(nllike,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut3a)
fitAlt=minimize(nllike2,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut3a)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitNull.fun-fitAlt.fun))
mut3aanswer=1-chi2.cdf(x=D,df=1)
print('I213N p value')
print(mut3aanswer)



#Q2

#load file
import pandas
file2=pandas.read_csv("MmarinumGrowth.csv",header=0,sep=",")

#plot data
ggplot(file2,aes(x='S',y='u'))+geom_point()+theme_classic()

def nllike(p,obs):
    umax=p[0]
    Ks=p[1]
    sigma=p[2]
    expected=umax*obs.S/(obs.S+Ks)
    nll=-1*norm(expected,sigma).logpdf(obs.u).sum()
    return nll

guess2=numpy.array([1, 1, 1])

fitNull2=minimize(nllike,guess2, method="Nelder-Mead",options={'disp': True}, args=file2)

print(fitNull2.x)



#Q3

#load file
import pandas
file3=pandas.read_csv("leafDecomp.csv",header=0,sep=",")

#plot data
ggplot(file3,aes(x='Ms',y='decomp'))+geom_point()+theme_classic()

#constant fit: d=a
def nllike(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

guess3=numpy.array([1, 1, 1])

fitNull3constant=minimize(nllike,guess3, method="Nelder-Mead",options={'disp': True}, args=file3)

print('constant fit')
print(fitNull3constant.x)

#linear fit: d = a + bMs
def nllike(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    expected=B0+B1*obs.Ms 
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

guess3=numpy.array([1, 1, 1])

fitNull3linear=minimize(nllike,guess3, method="Nelder-Mead",options={'disp': True}, args=file3)

print('linear fit')
print(fitNull3linear)

from scipy.stats import chi2
D=(2*(fitNull3constant.fun-fitNull3linear.fun))
linearVSconstant=(1-chi2.cdf(x=D,df=1))
print('linear vs constant p value')
print(linearVSconstant)

#quadratic fit: d = a +bMs + cMs^2
def nllike(p,obs):
    a=p[0]
    b=p[1]
    c=p[2]
    sigma=p[3]
    expected=a+b*obs.Ms+c*((obs.Ms)*(obs.Ms))
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

guess3=numpy.array([180, 15, -0.1, 10])

fitNull3quad=minimize(nllike,guess3, method="Nelder-Mead",options={'disp': True}, args=file3)

print('quadratic fit')
print(fitNull3quad)

from scipy.stats import chi2
D=(2*(fitNull3constant.fun-fitNull3quad.fun))
quadVSconstant=(1-chi2.cdf(x=D,df=1))
print('quadratic vs constant p value')
print(quadVSconstant)