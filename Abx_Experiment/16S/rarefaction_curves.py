# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 10:29:57 2022

@author: etb
"""

# Code adapted from https://forum.qiime2.org/t/how-to-convert-qzv-figure-into-pdf/15101/8
# To plot rarefaction curve output from QIIME2

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Formatting the csv file for plotting
path = 'C:/Users/etb/Documents/McGill/Experiments/Abx/abx1_2_controls/revised_diversity/alpha/nsnc-otus-rare.csv'
otus = pd.read_csv(path,index_col=0, sep=',')

cols = [col for col in otus.columns if 'iter' not in col]
mean,data = otus[cols],pd.DataFrame(columns=cols)
depths = [col.split('_')[0] for col in otus.columns if 'depth' in col]
otus = otus.drop(cols,axis=1)
otus.columns = depths
for depth in depths:
    mean['ASV'] = otus[depth].mean(axis=1)
    mean['depth']= depth.split('-')[-1]
    data = pd.concat([data,mean])

# here provide colors for each item that will be plotted
pal={'Control':'blue','Initial':'green'} 
#pal={'Control':'blue','Initial':'green','Van':'red','Pmb':'black'}
fig,ax = plt.subplots(figsize=(2.5,2),dpi=600,tight_layout=True)
sns.set(style='ticks',rc={"lines.linewidth":0.7,"axes.linewidth":0.5})

# use your column name to plot here instead of 'treatment'
sns.lineplot(x='depth',y='ASV',data=data,palette=pal,hue='treatment',sort=False,err_style='bars',\
             dashes=True,style='treatment',ci=67)
ax.set_xlabel('Sequencing depth',fontsize=8)
ax.set_ylabel('Observed Features',fontsize=8)
ax.tick_params(axis='x', labelrotation=90)
ax.tick_params(axis='both',which='major',length=2,pad=0.5,labelsize=6)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[1:],labels=labels[1:],fontsize=5,frameon=False,numpoints=4,borderaxespad=0,handletextpad=0.2,loc=2,)
plt.savefig('abx1_2_controls_otus.png', bbox_inches='tight')
#plt.savefig('abx2_otus.png', bbox_inches='tight')
