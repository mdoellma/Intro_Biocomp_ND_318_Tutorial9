#question 3
#Melissa Stephens and Janice Love
#Date: 10/31/2017

import pandas
file3=pandas.read_csv("leafDecomp.csv",header=0,sep=",")

ggplot(file3,aes(x='Ms',y='decomp'))+geom_point()+theme_classic()
#constant fit 
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

### linear fit
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

### quad fit --- doesn't work
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
