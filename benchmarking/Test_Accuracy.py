#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import binom_test 
import os

import seaborn as sns
sns.set(font_scale= 1.5)
sns.set_style('whitegrid') 

import sys
sys.path.append('/work-zfs/abattle4/heyuan/old_work_files/yuan/tools/python_lib/lib/python2.7/site-packages')

from matplotlib_venn import venn3, venn2

from statsmodels.stats.multitest import multipletests
sns.set(font_scale= 1.5)
sns.set_style('whitegrid')

from scipy.stats import ranksums


# In[52]:


discover_number = []
precision_matric = []
for sample in ['HG00108']:
    #print(sample)
    try:
        for minDP in range(2,11):
            bowtie_matrix = pd.read_csv('confusion_matrices/%s_bowtie_GATK_DPmin%d.txt' % (sample, minDP), 
                                        sep='\t', index_col = 0)
            discover_number.append(list(np.diag(bowtie_matrix)) + [minDP, sample])
            sensitivity_arr = np.diag(np.array(bowtie_matrix) / np.reshape(np.array(bowtie_matrix.sum(axis=1)), [1,3]))
            specificity_arr = np.diag(np.array(bowtie_matrix) / np.reshape(np.array(bowtie_matrix.sum(axis=0)), [1,3]))
            precision_matric.append(list(sensitivity_arr) + list(specificity_arr) + [minDP, sample])
    except:
        print(sample)
        continue

        
        
discover_df = pd.DataFrame(discover_number)
discover_df.columns = ['AA', 'AB', 'BB', 'minDP', 'sample']

precision_df = pd.DataFrame(precision_matric)
precision_df.columns = ['Sens_AA','Sens_AB','Sens_BB','Spec_AA', 'Spec_AB', 'Spec_BB', 'minDP', 'sample']


# In[51]:


plt.figure(figsize = (25,4))

plt.subplot(131)
for col in precision_df.columns[:3]:
    plt.plot(precision_df['minDP'], precision_df[col], 'o',label = col)
plt.legend()

#plt.figure()
plt.subplot(132)
for col in precision_df.columns[3:6]:
    plt.plot(precision_df['minDP'], precision_df[col], 'P',label = col)
plt.legend()


discover_df = discover_df.iloc[3:]
plt.subplot(133)
for col in discover_df.columns[:3]:
    plt.plot(discover_df['minDP'], discover_df[col], 'v',label = col)
plt.legend()


plt.show()
plt.close()


# In[ ]:





# In[66]:


discover_number = []
precision_matric = []
for sample in ['HG00108']:
    if 1:
        for minDP in range(3,11):
            bowtie_matrix = pd.read_csv('confusion_matrices/%s_bowtie_GATK_DPmin%d_SNPs_noAmbig.txt' % (sample, minDP), 
                                        sep='\t', index_col = 0)

            discover_number.append(list(np.diag(bowtie_matrix)) + [minDP, sample])
            sensitivity_arr = np.diag(np.array(bowtie_matrix) / np.reshape(np.array(bowtie_matrix.sum(axis=1)), [1,3]))
            specificity_arr = np.diag(np.array(bowtie_matrix) / np.reshape(np.array(bowtie_matrix.sum(axis=0)), [1,3]))
            precision_matric.append(list(sensitivity_arr) + list(specificity_arr) + [minDP, sample])
        
        
discover_df = pd.DataFrame(discover_number)
discover_df.columns = ['AA', 'AB', 'BB', 'minDP', 'sample']

precision_df = pd.DataFrame(precision_matric)
precision_df.columns = ['Sens_AA','Sens_AB','Sens_BB','Spec_AA', 'Spec_AB', 'Spec_BB', 'minDP', 'sample']

precision_df


# In[67]:


plt.figure(figsize = (25,4))

plt.subplot(131)
for col in precision_df.columns[:3]:
    plt.plot(precision_df['minDP'], precision_df[col], 'o',label = col)
plt.legend()

#plt.figure()
plt.subplot(132)
for col in precision_df.columns[3:6]:
    plt.plot(precision_df['minDP'], precision_df[col], 'P',label = col)
plt.legend()


discover_df = discover_df.iloc[3:]
plt.subplot(133)
for col in discover_df.columns[:3]:
    plt.plot(discover_df['minDP'], discover_df[col], 'v',label = col)
plt.legend()


plt.show()
plt.close()


# In[ ]:




