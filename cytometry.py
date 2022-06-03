#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:07:03 2021

@author: msangste
"""
import pandas as pd
from os import listdir
import cytoflow as flow
from cytoflow.utility.cytoflow_errors import CytoflowError


def set_figure(dpi=200, axes=5, xtick=5, ytick=5):
    """
    improve figure as created by cytoflow
    """
    import matplotlib
    matplotlib.rc('figure', dpi = dpi)
    matplotlib.rc('axes', labelsize=axes)
    matplotlib.rc('xtick', labelsize=xtick)
    matplotlib.rc('ytick', labelsize=ytick)
    
def set_names(area = False):
   """
   since FC strings are very long and complicated,
   set some common channel names as variables
   """
   FSC_H = 'Forward Scatter (FSC-HLin)'
   if area:
       FSC_A = 'Forward Scatter Area (FSC-ALin)'
   SSC = 'Side Scatter (SSC-HLin)'
   Orange_G = 'Orange-G Fluorescence (ORG-G-HLin)'
   Yellow_G = 'Yellow-G Fluorescence (YEL-G-HLin)'
   Red_G = 'Red-G Fluorescence (RED-G-HLin)'
   
   if area:
       return FSC_H, FSC_A, SSC, Orange_G, Yellow_G, Red_G
   else:
  
       return FSC_H, SSC, Orange_G, Yellow_G, Red_G


def number(sample, column, value):
    """
    count the number of events that have the given value for the given column
    e.g. calculate: how many events with 'Threshold' True
    set value to None to take all the events
    """
    if value != None:
        try:
            subset = sample.subset(column, value)
        except CytoflowError:
            print(f'no events with this value: {value}')
    else:
        print('no value given')
        subset = sample
    try:
        number = len(subset.data)
    except UnboundLocalError:
        number = 0
    return number

def number_without_dust(sample, sample_dust, column, value):
    number_sample = number(sample, column, value)
    number_dust = number(sample_dust, column, value)
    return number_sample - number_dust

def fraction(sample, sample_dust, column, value, nvalue):
    """
    calculate the fraction of events with the given value for the given column
    e.g. the fraction of events with 'Threshold' True
    """
    true = number_without_dust(sample, sample_dust, column, value)
    false = number_without_dust(sample, sample_dust, column, nvalue)
    try:
        fraction = true/(true+false)
    except ZeroDivisionError:
        fraction = 0
    return fraction

def fractions(experiment, column, value, nvalue, name_dust='M9', expected=False):
    """
    calculate the fraction of events with given value for each sample and
    append to a dataframe 
    """
    
    #determine dust sample
    sample_dust = experiment.subset('sample', name_dust)
    #create empty dataframe
    df = pd.DataFrame()
    for i, condition in enumerate(experiment.conditions['sample']):
        sample = experiment.subset('sample', condition)
        fraction_sample = fraction(sample, sample_dust, column, value, nvalue)
        appendix = {}
        
        #get a list of conditions, without fraction criteria
        conditions = list(experiment.conditions.keys())
        conditions.remove(column)
        
        for condition in conditions:
            appendix[condition] = sample[condition].iloc[0]
        
        appendix['scarlet fraction'] = fraction_sample
        df = df.append(appendix, ignore_index=True)

    return df


def numbers(experiment, column, value, nvalue, name_dust='M9'):
    """
    calculate the fraction of events with given value for each sample and
    append to a dataframe 
    """
    
    #determine dust sample
    sample_dust = experiment.subset('sample', name_dust)
    #create empty dataframe
    df = pd.DataFrame()
    for i, condition in enumerate(experiment.conditions['sample']):
        sample = experiment.subset('sample', condition)
        fraction_sample = number_without_dust(sample, sample_dust, column, value)
        appendix = {}
        
        #get a list of conditions, without fraction criteria
        conditions = list(experiment.conditions.keys())
        conditions.remove(column)
        
        for condition in conditions:
            appendix[condition] = sample[condition].iloc[0]
        
        appendix['scarlet number'] = fraction_sample
        df = df.append(appendix, ignore_index=True)
        
    return df

def create_experiment(folder, conditions):
    """
    folder: location of the files
    conditions: dictionary of condition name and list of names/values
        e.g. {'sample':['scarlet', 'WT', 'half scarlet'], 'volume scarlet':[1, 0, 0.5]}
    list should be in the same order as files after sorting 
    """
    tubes = []

    files = listdir(folder)
    files.sort()
    if 'output' in files:
        files.remove('output')
    
    
    for i, file in enumerate(files) :
        conditions2 = {}
        for key in conditions.keys():
            conditions2[key] = conditions[key][i]
        tube = flow.Tube(file = f'{folder}/{file}', 
                                  conditions=conditions2)
        
    
        tubes.append(tube)
    
    
    types = {}
    for key in conditions.keys():
        types[key] = type(conditions[key][0]).__name__
    
    import_op = flow.ImportOp(conditions=types, 
                              tubes=tubes)
    ex = import_op.apply()
    return ex

def expected_fraction(experiment, column, value, nvalue, scarlet, cerulean, dust, df_scarlet):
    """
    add a column of expected scarlet fraction to a df with already a column 'volume_scarlet'
    """
    sample_scarlet = experiment.subset('sample', scarlet)
    sample_WT = experiment.subset('sample', cerulean)
    sample_dust = experiment.subset('sample', dust)
    scarlet = number_without_dust(sample_scarlet, sample_dust, column, value)
    WT = number_without_dust(sample_WT, sample_dust, column, nvalue)
    
    df_scarlet['expected'] = (df_scarlet['volume_scarlet'] * 
     scarlet / (df_scarlet['volume_scarlet'] * scarlet + 
                (1-df_scarlet['volume_scarlet'])*WT))
    return df_scarlet


def load_bioreactor_cytometry(main_folder, df):
    tubes = []
    files = listdir(main_folder)
    if 'output' in files:
        files.remove('output')
    for file in files:
        try:
            bioreactor = int(df.loc[df['file']==file]['bioreactor'])
            time = float(df.loc[df['file']==file]['time'])
            tube = flow.Tube(file = f'{main_folder}/{file}',
                         conditions={'sample':bioreactor, 'time':time})
            tubes.append(tube)
        except TypeError:
            print(f'{file} not known')
        
       
        
    
    import_op = flow.ImportOp(tubes=tubes,conditions={'sample':'int','time':'float'})
    ex = import_op.apply()
    
    return ex
        
        
 