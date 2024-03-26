# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 08:08:24 2022

@author: laszl

Follow the white rabbit
"""

from __future__ import division
import pandas as pd
from scipy.sparse import csr_matrix
import numpy as np
import numpy.matlib
from sklearn.manifold import MDS
from matplotlib import pyplot as plt

"""Jaccard similarity calculation"""
dfm = pd.read_csv('Adj_matrix_T.csv')
t = csr_matrix(dfm)
T = np.multiply(t.T,t)
s = sum(t,0)
ss = s.toarray()
S = ss
l = S.shape[1]
for i in range(1,l,1): S = np.append(S,ss,axis=0)   
AB = S + S.T
Jaccard = T/(AB-T)

"""MDS coordinates calculation"""
mds = MDS(random_state=0,dissimilarity='precomputed',max_iter=2500,eps=1e-10)
X = mds.fit_transform(1-Jaccard)
coord = pd.DataFrame(X,columns={'X','y'})

"""Graph data"""
events = pd.read_csv('Seqs_w_colors.csv')
conf = pd.read_csv('Sequence_attributes.csv')

events['X'] = coord['X']
events['y'] = coord['y']
events['poidx'] = 0
events['paidx'] = 0
events['Conf'] = conf['C_conf']

for i in range(0,len(events),1):
    events.at[i,'Sequence'] = events.at[i,'Sequence'][:-1]
    p = str(events.at[i,'Sequence']).rfind('_')
    events.at[i,'Parent'] = events.at[i,'Sequence'][:p] if p > 0 else ''  
    events.at[i,'poidx'] = i
    for k in range(0,len(events),1):
        if events.at[k,'Sequence'] == events.at[i,'Parent']:
            events.at[i,'paidx'] = k
            k = len(events)

events['startswith'] = events.startswith.astype('str')
events.to_csv('seq.csv',index=False)

par = events[events['Parent'] != '']
parents = pd.DataFrame(list(set(par['Parent'].astype(str).values.tolist())))
events['Parent'] = events.Parent.astype('str')
events['startswith'] = events.startswith.astype('str')
events = events.replace('nan','')

"""Filtering options"""
filt_ev = ['2','3']
#filt_ev = list(set(events['startswith'].astype(str).values.tolist()))

ann_ind = []
annotations = []
idxs = events.index.tolist()
for i in idxs:
    if events.at[i,'startswith'] in filt_ev:
        ann_ind.append(i)
        annotations.append(events.at[i,'Sequence'])

"""Plot graph"""
fig, ax = plt.subplots(figsize=(72,29))
ax.set(xlim=(-1, 1), ylim=(-1, 1), aspect='equal')
# ax.tick_params(axis='both', labelsize=15)
# ax.spines['bottom'].set_position(('data',-1))
# ax.spines['left'].set_position(('data',-1))
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
maxConf = max(events['Conf'])
minConf = maxConf*0.00001 #Set minimum confidence of transitions to represent with an edge
ymin = events['y'].min()
ymax = events['y'].max()
dy = ymax-ymin
plt.axis('off')
for i in ann_ind: 
    #annX = events.at[i,'X'] + events.at[i,'X']*0.0015
    annX = events.at[i,'X']-int(events.at[i,'startswith'])*0.018
    mod = i % 4
    anny = events.at[i,'y'] + (mod-2)*dy*0.01
    #anny = events.at[i,'y']-int(events.at[i,'startswith'])*0.02
    color = events.at[i,'color']
    plt.scatter(events.at[i,'X'], events.at[i,'y'], s=200, c=color, cmap="RdYlGn")
    plt.annotate(events.at[i,'Sequence'], (annX, anny),fontsize=30)
    if events.at[i,'Parent'] != '' and events.at[i,'Conf'] > minConf:
        width = 0.01*events.at[i,'Conf']
        plt.arrow(events.at[i,'X'],events.at[i,'y'],events.at[events.at[i,'paidx'],'X']-events.at[i,'X'],events.at[events.at[i,'paidx'],'y']-events.at[i,'y'],width=width,color=color,head_length=0.0,head_width=0.0)

#plt.title('NbVFS-TM MDS projection')
#plt.savefig('MDS_T' + str(filt_ev) + '.svg', bbox_inches='tight') if len(filt_ev) < len(set(events['startswith'].astype(str).values.tolist())) else plt.savefig('MDS_T.svg', bbox_inches='tight')
plt.savefig('Fig5.svg', bbox_inches='tight')