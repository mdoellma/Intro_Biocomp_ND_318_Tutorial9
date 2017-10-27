#set working directory

#add packages
import numpy
import pandas
import re
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
ggplot(mrna,aes(x="mutation",y="ponzr1Counts"))+geom_point()+theme_classic()


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


