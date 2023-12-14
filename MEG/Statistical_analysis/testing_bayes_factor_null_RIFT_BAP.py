# -*- coding: utf-8 -*-
"""
This script will calculate the Bayes Factor
for the null behavioural and RIFT restuls
from the perceptual_load project.

the input is the Pearson rs from MATLAB
the output is bayes factor.

@author: GhafariT
"""

import pingouin as pg
import numpy as np
import pandas as pd
import os.path as op

# locate the LV and dependent variables (BAP or HLM(RIFT))
result_fpath = 'Z:\Programming\MATLAB\Perceptual_Load\FieldTrip\Results\group_level\Lateralization_indices'
LV_fpath = op.join(result_fpath, 'LV_tbl.csv')
BAP_fpath = op.join(result_fpath, 'BAP_RT_IES.csv')
HLM_RIFT_fpath = op.join(result_fpath, 'HLM_RIFT.csv' )

LV_df = pd.read_csv(LV_fpath)
LV_df = LV_df.drop(columns='sub #')  # drop the subject number column
BAP_df = pd.read_csv(BAP_fpath)
HLM_RIFT_df = pd.read_csv(HLM_RIFT_fpath)


# Ensure that the DataFrames have the same number of rows
if len(LV_df) != len(BAP_df) or len(LV_df) != len(HLM_RIFT_df):
    raise ValueError("DataFrames must have the same number of rows.")

# Iterate through columns and calculate Bayes factor for correlation
HLM_RIFT_data = HLM_RIFT_df['HLM_RIFT'].values
BAP_data = BAP_df['BAP'].values

for lv_column in LV_df.columns:
    # Extract columns as arrays
    lv_data = LV_df[lv_column].values  

    # Perform Bayesian correlation test
    bf_result_HLM = pg.corr(lv_data, HLM_RIFT_data, method='pearson')
    bf_result_BAP = pg.corr(lv_data, BAP_data, method='pearson')

    
    # Print the Bayes factor for each correlation
    print(f"Bayes Factor for Correlation between {lv_column} and RIFT HLM: {bf_result_HLM['BF10'].values[0]}")
    print(f"Bayes Factor for Correlation between {lv_column} and BA: {bf_result_BAP['BF10'].values[0]}")
    print("=============================")