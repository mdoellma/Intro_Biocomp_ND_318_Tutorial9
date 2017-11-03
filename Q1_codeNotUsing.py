#Code that didn't work 

WT = []
M124K = []
I213N = [] 
V456D = []

#place gene counts into respective list: WT or specific mutation 
for i in range(0,len(df),1):
    if df.mutation[i]=="WT":
        WT.append(df.ponzr1Counts[i])
    elif df.mutation[i] == "M124K":
        M124K.append(df.ponzr1Counts[i])
    elif df.mutation[i] == "I213N":
        I213N.append(df.ponzr1Counts[i])
    elif df.mutation[i] == "V456D":
        V456D.append(df.ponzr1Counts[i])


#create df with conditions to be compared 
WT_df = pd.DataFrame({'WT':WT})
WT_df['M124K']= M124K

WT_mut1 = pd.DataFrame(['WT','M124k']),:]

WT_I213N = pd.DataFrame({'WT':WT})
WT_I213N['I213N']= I213N

WT_V456D = pd.DataFrame({'WT':WT})
WT_V456D['V456D']= V456D
