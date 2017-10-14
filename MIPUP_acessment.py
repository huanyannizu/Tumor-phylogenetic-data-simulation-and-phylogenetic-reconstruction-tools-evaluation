import numpy as np
import pandas as pd
import networkx as nx
import pygraphviz as pgv
import itertools
import argparse
import ast

def MIPUP_Acessment(G_Dir,new_mutation_in_clone,columns,tree):
    n_mutation = 0
    for k in new_mutation_in_clone.keys():
        n_mutation = n_mutation + len(new_mutation_in_clone[k])
    
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
            
    MIPUP_acc = pd.DataFrame(np.zeros((1,1)),index=None,columns=['%Corr AD'])    
    muta_com = itertools.combinations(range(n_mutation),2)
    muta_com = list(muta_com)
    for pair in muta_com:
        muta1 = pair[0];muta2 = pair[1]
        if muta_relation_inp[muta1,muta2] == 1:
            if muta_relation_outp[muta1,muta2] == 1:MIPUP_acc['%Corr AD']+=1
        if muta_relation_inp[muta2,muta1] == 1:
            if muta_relation_outp[muta2,muta1] == 1:MIPUP_acc['%Corr AD']+=1
            
    unique, counts = np.unique(muta_relation_inp, return_counts=True)
    relation_counts = dict(zip(unique, counts))
    
    MIPUP_acc['%Corr AD'] = MIPUP_acc['%Corr AD']/relation_counts[1.0]
    
    return MIPUP_acc

def parse_arguments():   
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
##################################################################################################################################################################################
#parse the command line arguments
graph, textfile, csvfile, dotfile = parse_arguments()

#read in ground truth tree topology and mutations in its clones
G_Dir = nx.read_gpickle(graph)
f = open(textfile,'r')             
lines = f.readlines()                        
new_mutation_in_clone = ast.literal_eval(lines[2].rstrip('\n').lstrip('\r'))
f.close()

#acessment by comparing the ground truth tree with the predicted tree
result_ip = MIPUP_Acessment(G_Dir,new_mutation_in_clone,columns=csvfile,tree=dotfile)
result_ipd = MIPUP_Acessment(G_Dir,new_mutation_in_clone,columns=csvfile,tree=dotfile)            

#print out results            
print 'ip result:'
print result_ip.to_string(index=False)
print 'ipd result:'
print result_ipd.to_string(index=False)















