import numpy as np
import pandas as pd
import networkx as nx
import pygraphviz as pgv
import itertools
import subprocess
import time
import ast

def MIPUP_Acessment(G_Dir,new_mutation_in_clone,n_mutation,columns,tree):
    f = open(columns).readlines()
    f.pop(0);muta_group = {}
    for line in f:
        k = line.split(';')[0] 
        v = line.split(';')[1:]
        new_list = []
        for x in v:
            new_list.append(int(x.split(':')[0]))
        muta_group[k] = new_list
    
    G_outp = pgv.AGraph(tree)
    G_outp_nx = nx.DiGraph();muta_fraction = 0
    for edge in G_outp.edges():
        label = edge.attr['label'].encode('utf8')
        if len(label) != 0:
            G_outp_nx.add_edge(*edge,label = label[0])
            muta_fraction = muta_fraction + len(muta_group[label[0]])
        else:
            G_outp_nx.add_edge(*edge,label = '')
    
    muta_relation_inp = np.zeros((n_mutation,n_mutation))     
    node_com = itertools.combinations(G_Dir.nodes(),2)
    node_com = list(node_com)
    for pair in node_com:
        node1 = pair[0];node2 = pair[1]
        if nx.has_path(G_Dir, node1, node2):
            for i in new_mutation_in_clone[node1]:
                for j in new_mutation_in_clone[node2]:
                    muta_relation_inp[i,j] = 1
        if nx.has_path(G_Dir, node2, node1): 
            for i in new_mutation_in_clone[node1]:
                for j in new_mutation_in_clone[node2]:
                    muta_relation_inp[j,i] = 1                            
        if G_Dir.predecessors(node1) == G_Dir.predecessors(node2):
            for i in new_mutation_in_clone[node1]:
                for j in new_mutation_in_clone[node2]:
                    muta_relation_inp[i,j] = 2
    
    muta_relation_outp = np.zeros((n_mutation,n_mutation))
    edge_com = itertools.combinations(G_outp_nx.edges(),2)
    edge_com = list(edge_com)
    for pair in edge_com:
        edge1 = pair[0];label1 = G_outp_nx.get_edge_data(*edge1)['label']
        edge2 = pair[1];label2 = G_outp_nx.get_edge_data(*edge2)['label']
        if len(label1)!=0 and len(label2)!=0:
            if nx.has_path(G_outp_nx, edge1[1], edge2[0]):
                for i in muta_group[label1]:
                    for j in muta_group[label2]:
                        muta_relation_outp[i,j] = 1
            if nx.has_path(G_outp_nx, edge2[1], edge1[0]):
                for i in muta_group[label1]:
                    for j in muta_group[label2]:
                        muta_relation_outp[j,i] = 1
            if edge1[0] == edge2[0]:
                for i in muta_group[label1]:
                    for j in muta_group[label2]:
                        muta_relation_outp[i,j] = 2
            if nx.has_path(G_outp_nx, edge1[0], edge2[0]):
                p = nx.shortest_path(G_outp_nx,edge1[0],edge2[0])
                lab_len = 0
                for i in range(len(p)-1):
                    lab = G_outp_nx.get_edge_data(p[i],p[i+1])['label']
                    lab_len = lab_len + len(lab)
                if lab_len == 0:
                    for i in muta_group[label1]:
                        for j in muta_group[label2]:
                            muta_relation_outp[i,j] = 2
            if nx.has_path(G_outp_nx, edge2[0], edge1[0]):
                p = nx.shortest_path(G_outp_nx,edge2[0],edge1[0])
                lab_len = 0
                for i in range(len(p)-1):
                    lab = G_outp_nx.get_edge_data(p[i],p[i+1])['label']
                    lab_len = lab_len + len(lab)
                if lab_len == 0:
                    for i in muta_group[label1]:
                        for j in muta_group[label2]:
                            muta_relation_outp[i,j] = 2
           
    MIPUP_acc = pd.DataFrame(np.zeros((1,6)),index=None,columns=['%Corr(AD+Sib)','%Corr AD','%AD reverse','%Corr Sib','%AD->Sib','%Sib->AD'])    
    muta_com = itertools.combinations(range(n_mutation),2)
    muta_com = list(muta_com)
    for pair in muta_com:
        muta1 = pair[0];muta2 = pair[1]
        if muta_relation_inp[muta1,muta2] == 1:
            if muta_relation_outp[muta1,muta2] == 1:MIPUP_acc['%Corr AD']+=1;MIPUP_acc['%Corr(AD+Sib)']+=1
            if muta_relation_outp[muta2,muta1] == 1:MIPUP_acc['%AD reverse']+=1
            if muta_relation_outp[muta1,muta2] == 2 or muta_relation_outp[muta2,muta1] == 2:MIPUP_acc['%AD->Sib']+=1
        if muta_relation_inp[muta2,muta1] == 1:
            if muta_relation_outp[muta2,muta1] == 1:MIPUP_acc['%Corr AD']+=1;MIPUP_acc['%Corr(AD+Sib)']+=1
            if muta_relation_outp[muta1,muta2] == 1:MIPUP_acc['%AD reverse']+=1
            if muta_relation_outp[muta1,muta2] == 2 or muta_relation_outp[muta2,muta1] == 2:MIPUP_acc['%AD->Sib']+=1
        if muta_relation_inp[muta1,muta2] == 2 or muta_relation_inp[muta2,muta1] == 2:
            if muta_relation_outp[muta1,muta2] == 1 or muta_relation_outp[muta2,muta1] == 1:MIPUP_acc['%Sib->AD']+=1
            if muta_relation_outp[muta1,muta2] == 2 or muta_relation_outp[muta2,muta1] == 2:MIPUP_acc['%Corr Sib']+=1;MIPUP_acc['%Corr(AD+Sib)']+=1
   
    unique, counts = np.unique(muta_relation_inp, return_counts=True)
    relation_counts = dict(zip(unique, counts))
    MIPUP_acc['%SSNVs'] = muta_fraction

    if 2 in relation_counts:
        MIPUP_acc['%Corr(AD+Sib)'] = MIPUP_acc['%Corr(AD+Sib)']/(relation_counts[1.0]+relation_counts[2.0])
        MIPUP_acc['%Corr Sib'] = MIPUP_acc['%Corr Sib']/relation_counts[2.0]
        MIPUP_acc['%Sib->AD'] = MIPUP_acc['%Sib->AD']/relation_counts[2.0]
    else:
        MIPUP_acc['%Corr(AD+Sib)'] = MIPUP_acc['%Corr(AD+Sib)']/relation_counts[1.0]
        MIPUP_acc['%Corr Sib'] = 0.0
        MIPUP_acc['%Sib->AD'] = 0.0
    MIPUP_acc['%Corr AD'] = MIPUP_acc['%Corr AD']/relation_counts[1.0]
    MIPUP_acc['%AD reverse'] = MIPUP_acc['%AD reverse']/relation_counts[1.0]
    MIPUP_acc['%AD->Sib'] = MIPUP_acc['%AD->Sib']/relation_counts[1.0]
    
    MIPUP_acc['%SSNVs'] = muta_fraction/(n_mutation+0.0)
    
    Acc2 = pd.DataFrame(np.zeros((1,4)),index=None,columns=['Precision','Recall','True Negative Rate','Accuracy'])  
    muta_com = itertools.combinations(range(n_mutation),2)
    muta_com = list(muta_com)
    tp = 0;fn = 0;fp = 0;tn = 0
    for pair in muta_com:
        muta1 = pair[0];muta2 = pair[1]
        if muta_relation_inp[muta1,muta2] == 1:
            if muta_relation_outp[muta1,muta2] == 1:tp = tp +1
            else:fn = fn + 1
        elif muta_relation_outp[muta1,muta2] == 1:fp = fp + 1
        else:tn = tn + 1
        if muta_relation_inp[muta2,muta1] == 1:
            if muta_relation_outp[muta2,muta1] == 1:tp = tp +1
            else:fn = fn + 1
        elif muta_relation_outp[muta2,muta1] == 1:fp = fp + 1
        else:tn = tn + 1
                              
    if tp+fp == 0:
        Acc2['Precision'] = None
    else:
        Acc2['Precision'] = tp/(tp+fp+0.0)
    if tp+fn == 0:
        Acc2['Recall'] = None
    else:
        Acc2['Recall'] = tp/(tp+fn+0.0)
    if tn+fp == 0:
        Acc2['True Negative Rate'] = None
    else:
        Acc2['True Negative Rate'] = tn/(tn+fp+0.0)
    if tp+fn+fp+tn == 0:
        Acc2['Accuracy'] = None
    else:
        Acc2['Accuracy'] = (tp+tn)/(tp+fn+fp+tn+0.0)

    result = pd.concat([MIPUP_acc, Acc2], axis=1)
    return result

######################################################################################
def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", type=str)
    parser.add_argument("textfile", type=str)
    parser.add_argument("csvfile", type=str)
    parser.add_argument("dotfile", type=str)
    
    args = parser.parse_args()
    
    graph = args.graph
    textfile = args.textfile
    csvfile = args.csvfile
    dotfile = args.dotfile
    
    return graph, textfile, csvfile, dotfile

########################################################################################
#parse the command line arguments
graph, textfile, csvfile, dotfile = parse_arguments()

#read in ground truth tree topology and mutations in its clones
G_Dir = nx.read_gpickle(graph)
f = open(textfile,'r')             
lines = f.readlines()                        
new_mutation_in_clone = ast.literal_eval(lines[2].rstrip('\n').lstrip('\r'))
f.close()
n_mutation = 0
for c in new_mutation_in_clone:
    n_mutation = n_mutation + len(new_mutation_in_clone[c])

#acessment by comparing the ground truth tree with the predicted tree
result = MIPUP_Acessment(G_Dir,new_mutation_in_clone,n_mutation,columns=csvfile,tree=dotfile)           

#print out results            
print result.to_string(index=False)
