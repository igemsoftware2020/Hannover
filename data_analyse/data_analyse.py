# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 20:00:51 2020

@author: jonas
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


'''Get the data from a file'''
data = pd.read_json('C:/Users/jonas/Python Projekte/data/10start 75steps.json')#reades the data
print(data['BACTERIA'])
#print(len(data.index))

'''Clean the data'''
newdata=pd.DataFrame(data['BACTERIA'][0:-19])#get the interesting data(labeld with bacteria)
li=[]#temp. storrage for the data
for i in range(len(data.index)-19):
    li.append(pd.DataFrame(newdata['BACTERIA'][i]))#fill list with dataframe  
index=data['BACTERIA'].index[0:len(data.index)-19]#get the correct index for each bacteria
clean_data=pd.concat(li,keys=index)#add all dataframes in the list to one dataframe
print(clean_data)


'''Analyse velocity'''
velocity=pd.Series(clean_data['velocity']).apply(np.array)#convert lists to np.arrays
velocity=velocity.apply(np.linalg.norm)#get the magnitude 
velocity=velocity.unstack(level=0)#get a dataframe of velocitys 
velocity=velocity.transform(lambda x: sorted(x,key=pd.isnull,reverse=True))#move NaNs to beginning 
print(velocity)
plt.plot(velocity.iloc[:,0:].mean(axis=1),label='mean')#plot means for each iteration
plt.legend()
plt.show()


'''Analyse length'''
print(clean_data['length'].unstack(level=0).transform(lambda x: sorted(x,key=pd.isnull,reverse=True)))#dataframe of length
print(clean_data['length'].unstack(level=0).transform(lambda x: sorted(x,key=pd.isnull,reverse=True)).mean(axis=1))#means for each row 

'''Bacteria growth'''
living=clean_data['living'].unstack(level=0).transform(lambda x: sorted(x,key=pd.isnull,reverse=True))
grow=[]#number of bacteria in iteration steps
for i in range(len(living.index)):
    a=living.iloc[i,:].count()#count True values in i'th row 
    print(a)
    grow.append(a)
plt.plot(range(len(living.index)),grow)#plot number of bacteria as a function of iteration steps  
plt.show()