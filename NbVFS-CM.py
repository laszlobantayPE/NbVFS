# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 12:24:11 2022

@author: laszl

Follow the white rabbit
"""

from __future__ import division
import pandas as pd
import numpy as np
import numpy.matlib
from sklearn.manifold import MDS
from matplotlib import pyplot as plt

"""Transitive distance calculation"""
A = pd.read_csv('Adj_matrix_C.csv')
S = np.array(A)
S = np.fmax(S,S.T) 
N = len(A)
for j in range(0,N,1):
    for k in range(0,N,1): 
        d1 = np.asmatrix(S[:,k]).reshape(N,1)
        d2 = np.asmatrix(S[k,:])
        dum1 = np.matlib.repmat(d1,1,N)
        dum2 = np.matlib.repmat(d2,N,1)
        DUM = np.multiply(dum1,dum2)
        S = np.fmax(S,DUM)

"""MDS coordinates calculation"""
mds = MDS(random_state=0,dissimilarity='precomputed',max_iter=2500,eps=1e-10)
X = mds.fit_transform(1-S)
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
par = events[events['Parent'] != '']
parents = pd.DataFrame(list(set(par['Parent'].astype(str).values.tolist())))
events['Parent'] = events.Parent.astype('str')
events['startswith'] = events.startswith.astype('str')
events = events.replace('nan','')

"""Filtering options"""
filt_ev = ['2','3','4']
#filt_ev = list(set(events['startswith'].astype(str).values.tolist()))

ann_ind = []
annotations = []
idxs = events.index.tolist()
for i in idxs:
    if events.at[i,'startswith'] in filt_ev:
        ann_ind.append(i)
        annotations.append(events.at[i,'Sequence'])

"""Plot graph"""
#fig = plt.figure(figsize=(72,29),frameon=False)
fig, ax = plt.subplots(figsize=(72,29))
ax.set(xlim=(-1, 1), ylim=(-1, 1), aspect='equal')
#ax.tick_params(axis='both', labelsize=15)
#ax.spines['bottom'].set_position(('data',-1))
#ax.spines['left'].set_position(('data',-1))
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
maxConf = max(events['Conf'])
minConf = maxConf*0.00001 #Set minimum confidence of transitions to represent with an edge
ymin = events['y'].min()
ymax = events['y'].max()
dy = ymax-ymin
plt.axis('off')
for i in ann_ind: 
    annX = events.at[i,'X'] + events.at[i,'X']*0.015
    mod = i % 4
    anny = events.at[i,'y'] + (mod-2)*dy*0.01
    color = events.at[i,'color']
    plt.scatter(events.at[i,'X'], events.at[i,'y'], s=200, c=color, cmap="RdYlGn")
    plt.annotate(events.at[i,'Sequence'], (annX, anny),fontsize=30)
    if events.at[i,'Parent'] != '' and events.at[i,'Conf'] > minConf:
        width = 0.01*events.at[i,'Conf']
        plt.arrow(events.at[i,'X'],events.at[i,'y'],events.at[events.at[i,'paidx'],'X']-events.at[i,'X'],events.at[events.at[i,'paidx'],'y']-events.at[i,'y'],width=width,color=color,head_length=0.0,head_width=0.0)

#plt.title('NbVFS-CM MDS projection')
#plt.savefig('MDS_C' + str(filt_ev) + '.svg', bbox_inches='tight') if len(filt_ev) < len(set(events['startswith'].astype(str).values.tolist())) else plt.savefig('MDS_C.svg', bbox_inches='tight')
plt.savefig('Fig4.svg', bbox_inches='tight')