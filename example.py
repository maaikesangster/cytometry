#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:57:19 2021

@author: msangste
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 13:01:43 2021

@author: msangste

An example of the cytometry functions
"""

#import methods
import cytoflow as flow
from cytometry import set_figure, set_names, fractions, number_without_dust
from os import listdir
import pandas as pd

#improve figures 
set_figure(axes=4, xtick=4, ytick=4)

#default some names of often used channels
FSC_H, FSC_A, SSC, Orange_G, Yellow_G, Red_G = set_names(area=True)    

#create empty list of tubes
tubes = []

#load first data set
#specify main folder
main_folder = '220114/2022-01-14_at_04-49-28pm'
#list all flow cytometry files in main folder
files = listdir(f'{main_folder}')
files.sort()
if 'output' in files:
    files.remove('output')
    
#give a list of sample names corresponding to the tubes
samples = ['1', '2', '3', '4', '10', 'scarlet', 'cerulean', 'M9']

##import cytometry files in tubes
for file, sample, in zip(files, samples, ) :
    tube = flow.Tube(file = f'{main_folder}/{file}', 
                              conditions={
                                          'sample': sample,
                                
                                         
                                          })
    

    tubes.append(tube)

#create a cytoflow experiment from tubes 
import_op = flow.ImportOp(conditions= { 'sample':'str', }, 
                          tubes=tubes)
ex = import_op.apply()

#set threshold of forward scatter 
thresh_op = flow.ThresholdOp(name = 'Threshold',
                            channel = FSC_H,
                            threshold = 40)

ex2 = thresh_op.apply(ex)

ex2 = ex2.subset('Threshold', True)


#apply kmeans to select scarlet
km_op = flow.KMeansOp(name = 'KMeans',
                    channels = [Orange_G, FSC_H],
                       scale = {Orange_G : 'log', FSC_H : 'log'},
                      num_clusters = 2,
                     
                     
                      )

km_op.estimate(ex2,
              )

ex3 = km_op.apply(ex2)

#calculate scarlet fractions based on KMeans. KMeans_2 are scarlet cells, KMeans_1 are not
#need to figure out which is which before applying this function
df_scarlet = fractions(ex3, 'KMeans', 'KMeans_2', 'KMeans_1')


