# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 09:59:38 2022

@author: laszl

Follow the white rabbit
"""

from __future__ import division
import pandas as pd
from pyvis.network import Network
import networkx as nx
import random

"""Weighted graph"""
df = pd.read_csv('Weighted_network_data.csv')
filt = 1 #Filter on first event, 0:no, 1:yes
filt_ev = 2 #ID of first event of filtered sequences
if filt == 1:
    filter_event = str(filt_ev)
    k = 0
    k_max = len(df)
    while k < k_max:
        p = str(df.at[k,'Children']).find('_')
        if str(df.at[k,'Children'][:p]) != filter_event:
            df = df.drop([k])
        k += 1
    df = df.reset_index(drop=True)

df['SeqLen'] = df.SeqLen.astype('int')
df['Confidence'] = round(df['Confidence'],3)

"""Parent nodes"""
pnodes = pd.DataFrame(data=df['Parent'])
pnodes['Support'] = df['P_sup']
pnodes['SeqLen'] = df['SeqLen']
pnodes['Confidence'] = df['P_conf']
pnodes.rename(columns={'Parent':'Node'},inplace=True)
pnodes['SeqLen'] = pnodes['SeqLen']-1
pnodes = pnodes.drop_duplicates(subset=['Node'])

"""Children nodes"""
cnodes = pd.DataFrame(data=df['Children'])
cnodes['Support'] = df['C_sup']
cnodes['SeqLen'] = df['SeqLen']
cnodes['Confidence'] = df['C_conf']
cnodes.rename(columns={'Children':'Node'},inplace=True)
cnodes = cnodes.drop_duplicates(subset=['Node'])

"""All nodes"""
nodes = pnodes.append(cnodes)
nodes['Confidence'] = round(nodes['Confidence'],3)
nodes = nodes.drop_duplicates(subset=['Node'])
nodes = nodes.reset_index(drop=True)

"""Generating the graph and allocating colors to sequence length"""
G = nx.Graph()
sl = set(df['SeqLen'].values)
sl.add(1)
sl = list(sl)
colors = list(["#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)]) for j in range(len(sl))])

"""Add nodes"""
j = 0
j_max = len(nodes)
while j < j_max:
    idx = nodes.at[j,'SeqLen']-1
    value = nodes.at[j,'Confidence']
    node = str(nodes.at[j,'Node'])
    title = node + ' Support:' + str(nodes.at[j,'Support']) + ' Confidence:' + str(nodes.at[j,'Confidence'])
    G.add_node(nodes.at[j,'Node'], title=title, color=colors[idx], value=value)
    j += 1

"""Add edges"""
i = 0
i_max = len(df)
while i < i_max:
    title = 'Confidence: ' + str(df.at[i,'Confidence']) + ' Added event: ' + str(df.at[i,'Added_event'])
    G.add_edge(df.at[i,'Parent'], df.at[i,'Children'], value=df.at[i,'Confidence']*1000, title=title)
    i += 1

"""Plot graph"""
nt = Network('900px', '1300px')
nt.from_nx(G)
nt.show_buttons(filter_=['layout'])
nt.show('NbVFS-WN_weighted_network.html')