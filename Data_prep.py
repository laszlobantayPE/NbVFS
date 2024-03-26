# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 08:59:19 2022

@author: laszl

Follow the white rabbit
"""

from __future__ import division
import os
import pandas as pd
import csv
import random

"""Basic data"""
basic = pd.read_csv('Transactional_database.csv',sep=',')
transids = basic['SessionID'].unique().tolist()
eventids = basic['service_detail_EN'].unique().tolist()
eventcodes = pd.DataFrame({'Event':eventids})
i = 0
while i < len(eventcodes): 
    eventcodes.at[i,'Code'] = str(i)
    i += 1
evc = dict(zip(eventcodes.Event, eventcodes.Code))
eventcodes.to_csv('Event_codes.csv',index=False)

transactions = pd.DataFrame()
for i in transids:
    trans = basic[basic['SessionID'] == i]
    events = pd.DataFrame([trans['service_detail_EN'].tolist()])
    transactions = pd.concat([transactions,events])

#transactions.index = transids
transactions = transactions.reset_index(drop=True)
transactions.replace(evc,inplace=True)
transactions.to_csv('Transactions.csv',index=False)

"""Sequence database to SPMF format"""
command = "java -jar spmf.jar run Convert_a_sequence_database_to_SPMF_format Transactions.csv Transactions_spmf.csv CSV_INTEGER 100000"  #The command needs to be a string
os.system(command) #The command can also be passed as a string, instead of a variable

"""Frequent sequence mining"""
command = 'java -jar spmf.jar run CM-SPAM Transactions_spmf.csv Freqseqs.csv 0.01 1 10 "" 1 true'   #The command needs to be a string
os.system(command)

"""Information extraction and sequence attribute generation"""
with open('Freqseqs.csv', 'r') as csvfile:
     reader = csv.reader(csvfile, delimiter=' ')
     fc = [ row[0] for row in reader ]

c = 0
c_max = len(fc)
while c < c_max:
    data = pd.read_csv('Freqseqs.csv',sep=" ", skiprows=range(0,c), header=None, nrows=1)
    data = data.replace(-1,'_')
    data = data.astype({0:'str'})
    n = 1
    n_max = len(data.columns)
    while n < n_max:
        if data.at[0,n] == '#SID:':
            k = 1
            while k < n-2:
                data = data.astype({k:'str'})
                data.at[0,0] = str(data.at[0,0]) + str(data.at[0,k])
                data = data.drop([k], axis=1)
                k += 1
            data = data.drop([n-2], axis=1)
            data = data.iloc[:,:2]
            n = n_max
        else:
            n += 1
    data.columns = ['Sequence','Support']
    if c == 0:
        df = data
    else:
        df = df.append(data,ignore_index=True)
    c += 1

events = pd.DataFrame(sorted(list(set(fc))),columns=['Events'])
events.to_csv('Events.csv',index=False)

i = 0
while i < len(df):
    df.at[i,'SeqLen'] = str(df.at[i,'Sequence']).count('_')
    i += 1

df = df.reset_index(drop=True)
df.to_csv('Sequences.csv',index=False)

"""Generating network data"""
value = df.index.values.T
df['SeqID'] = value
df['Parent'] = ''
df['Confidence'] = 0
df['Confidence'] = df.Confidence.astype('float')
df['P_sup'] = 0
df['P_conf'] = 0
df['P_conf'] = df.P_conf.astype('float')
df['C_conf'] = 0
df['C_conf'] = df.C_conf.astype('float')
df.rename(columns={'Sequence':'Children','Support':'C_sup'},inplace=True)
df['Added_event'] = ''

i = 0
i_max = len(df)
while i < i_max:
    df.at[i,'Children'] = df.at[i,'Children'][:-1]
    p = str(df.at[i,'Children']).rfind('_')
    df.at[i,'Parent'] = '' if df.at[i,'SeqLen'] == 1 else df.at[i,'Children'][:p]   
    idx = df.index[df['Children'] == df.at[i,'Parent']]
    df.at[i,'Confidence'] = 1 if df.at[i,'SeqLen'] == 1 else df['C_sup'][i]/df['C_sup'][idx[0]]
    df.at[i,'P_sup'] = df['C_sup'][idx[0]] if df.at[i,'SeqLen'] > 1 else df['C_sup'][i]
    df.at[i,'P_conf'] = 1 if df.at[i,'SeqLen'] == 1 else df['C_conf'][idx[0]]
    df.at[i,'C_conf'] = 1 if df.at[i,'SeqLen'] == 1 else df.at[i,'P_conf']*df.at[i,'Confidence']
    p = str(df.at[i,'Children']).rfind('_')
    df.at[i,'Added_event'] = df.at[i,'Children'][p+1:]
    i += 1
df.to_csv('Sequence_attributes.csv',index=False)

df = df[df['Parent'] != '']
cols = list(df)
cols[0],cols[3] = cols[3], cols[0]
df = df[cols]
cols[1],cols[4] = cols[4], cols[1]
df = df[cols]
cols[2],cols[6] = cols[6], cols[2]
df = df[cols]
df = df.reset_index(drop=True)
df.to_csv('Weighted_network_data.csv',index=False)

"""Collecting supporting transactions"""
with open('Freqseqs.csv', 'r') as csvfile:
     reader = csv.reader(csvfile, delimiter=' ')
     fc = [ row[0] for row in reader ]

c = 0
c_max = len(fc)
traces = []
while c < c_max:
    data = pd.read_csv('Freqseqs.csv',sep=" ", skiprows=range(0,c), header=None, nrows=1)
    n = 1
    n_max = len(data.columns)
    while n < n_max:
        if data.at[0,n] == '#SID:':
            tt = data.iloc[:,n+1:].values.flatten().tolist()
            traces.extend(tt)
            n = n_max
        else:
            n += 1
    c += 1
sid = pd.DataFrame(list(set(traces)))

"""Allocation of frequent sequences to supporting transactions"""
sid.rename(columns={0:'Trace'},inplace=True)
sid['Trace'] = sid.Trace.astype('int')
sd = len(sid)
seq = pd.read_csv('Sequences.csv')
sq = len(seq)

i = 0
while i < sq: 
    column_name = 'Sequence' + str(i+1)
    sid[column_name] = ''
    i += 1

with open('freqseqs.csv', 'r') as csvfile:
     reader = csv.reader(csvfile, delimiter=' ')
     fc = [ row[0] for row in reader ]
c = 0
c_max = len(fc)
while c < c_max:
    print('c='+str(c)+'/'+str(c_max))
    data = pd.read_csv('Freqseqs.csv',sep=" ", skiprows=range(0,c), header=None, nrows=1)
    data = data.replace(-1,'_')
    data = data.astype({0:'str'})
    n = 0
    n_max = len(data.columns)
    while n < n_max:
        if data.at[0,n] == '#SID:':
            tr = data.iloc[:,n+1:]
            tr = tr.transpose()
            tr = tr[0].tolist()
            k = 1
            while k < n-2:
                data = data.astype({k:'str'})
                data.at[0,0] = str(data.at[0,0]) + str(data.at[0,k])
                data = data.drop([k], axis=1)
                k += 1
            data = data.drop([n-2], axis=1)
            data = data.iloc[:,:1]
            sequence = data.at[0,0] #sequence identification
            n = n_max
        else:
            n += 1        
    j = 0
    while j < sd:
        s = sid.at[j,'Trace']
        if s in tr: #if trace is in tr
            m = 1
            while m < len(sid.columns):
                col = 'Sequence' + str(m)
                if sid.at[j,col] == '': #locating the first empty spot
                    sid.at[j,col] = sequence
                    m = len(sid.columns)
                else:
                    m += 1
            j += 1
        else:
            j += 1
    c += 1
sid = sid.drop(columns='Trace',axis=1)
sid2 = sid.dropna(axis=1,how='all')
sid2.to_csv('Traces.csv',index=False)

seq0 = pd.read_csv('Sequence_attributes.csv')
seq = seq0['Children']

"""Creation of the adjacency matrix based on confidence of transitions"""
A = pd.DataFrame(data=0,index=seq,columns=seq,dtype=float)
for i in range(0,len(seq),1):
    sq1 =seq[i]
    for k in range(0,len(seq),1):
        sq2 = seq[k]
        p = sq2.rfind('_')
        if sq1 == sq2:
            A.at[sq1,sq2] = 0
        elif p > 0 and sq2[:p] == sq1:
            A.at[sq1,sq2] = seq0.at[k,'C_conf']
A.to_csv('Adj_matrix_C.csv',index=False)

"""Creation of the adjacency matrix based on supporting transactions"""
matrix = pd.read_csv('Traces.csv')
i = 0
i_max = len(matrix.columns)
while i < i_max:
    colnameold = 'Sequence' + str(i+1)
    colnamenew = seq[i] + '_'
    matrix.rename(columns={colnameold:colnamenew},inplace=True)
    matrix[colnamenew] = 0
    i += 1

cseqs = matrix.columns.values.tolist()
sid = pd.read_csv('Traces.csv')
i = 0
i_max = len(sid)
j_max = len(matrix.columns)
while i < i_max:
    seqs = sid.iloc[i,:]
    seqs = seqs.dropna()
    seqs = seqs.values.tolist()
    j = 0
    while j < j_max:
        checkseq = cseqs[j]
        if checkseq in seqs:
            matrix.at[i,checkseq] = 1
        j += 1
    i += 1           
matrix.to_csv('Adj_matrix_T.csv',index=False) 

"""Visualization of data"""
seqs = pd.read_csv('Sequences.csv')
for i in range(0,len(seqs),1):
    p = str(seqs.at[i,'Sequence']).find('_')
    seqs.at[i,'startswith'] = seqs.at[i,'Sequence'][:p]
seqs['color'] = seqs['startswith']
sqs = list(seqs['startswith'].unique())
colors = list(["#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)]) for j in range(len(sqs))])
seq_colors = {sqs[i]: colors[i] for i in range(len(sqs))}
seqs['color'].replace(seq_colors, inplace=True)
seqs.to_csv('Seqs_w_colors.csv',index=False)
