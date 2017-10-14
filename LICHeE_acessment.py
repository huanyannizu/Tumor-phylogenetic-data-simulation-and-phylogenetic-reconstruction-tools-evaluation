import numpy as np
import pandas as pd
import networkx as nx
import itertools
import argparse
import ast

def LICHeE_Acessment(G_Dir,new_mutation_in_clone,tree):
    n_mutation = 0
    for k in new_mutation_in_clone.keys():
        n_mutation = n_mutation + len(new_mutation_in_clone[k])
    
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

    muta_relation_inp = np.r_[np.zeros((1,n_mutation)),muta_relation_inp]
    muta_relation_inp = np.c_[np.zeros((n_mutation+1,1)),muta_relation_inp]

    f = open(tree).readlines()
    f.pop(0);muta_group = {};muta_group[0] = []
    for line in f:
        if line == '\n':break
        line = line.rstrip('\n')
        line = line.split('\t')
        v = []
        for vv in line[3:]:
            v.append(int(vv[3:]))
        muta_group[int(line[0])] = v
    
    muta_fraction = 0
    for k in muta_group.keys():
        muta_fraction = muta_fraction + len(muta_group[k])
    
    G_outp_nx = nx.DiGraph()    
    for line in f:
        if '->' in line: #and line[0] != '0':
            nodes = line.rstrip('\n').split(' ')
            G_outp_nx.add_edge(int(nodes[0]),int(nodes[2]))
        
    muta_relation_outp = np.zeros((n_mutation+1,n_mutation+1))#only read from row,1 is AD, 2 is Sib,0 is no relation
    node_com = itertools.combinations(G_outp_nx.nodes(),2)
    node_com = list(node_com)
    for pair in node_com:
        node1 = pair[0];node2 = pair[1]
        if nx.has_path(G_outp_nx, node1, node2):
            for i in muta_group[node1]:
                for j in muta_group[node2]:
                    muta_relation_outp[i,j] = 1
        if nx.has_path(G_outp_nx, node2, node1):
            for i in muta_group[node1]:
                for j in muta_group[node2]:
                    muta_relation_outp[j,i] = 1
            
    Acc = pd.DataFrame(np.zeros((1,1)),index=None,columns=['%Corr AD'])    
    muta_com = itertools.combinations(range(1,n_mutation+1),2)
    muta_com = list(muta_com)
    for pair in muta_com:
        muta1 = pair[0];muta2 = pair[1]
        if muta_relation_inp[muta1,muta2] == 1:
            if muta_relation_outp[muta1,muta2] == 1:Acc['%Corr AD']+=1
            
        if muta_relation_inp[muta2,muta1] == 1:
            if muta_relation_outp[muta2,muta1] == 1:Acc['%Corr AD']+=1
            
        
    unique, counts = np.unique(muta_relation_inp, return_counts=True)
    relation_counts = dict(zip(unique, counts))
    
    Acc['%Corr AD'] = Acc['%Corr AD']/relation_counts[1.0]
    
    return Acc

def parse_arguments():   
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", type=str)
    parser.add_argument("textfile", type=str)
    parser.add_argument("treefile", type=str)
    
    args = parser.parse_args()
    
    graph = args.graph
    textfile = args.textfile
    treefile = args.treefile
    
    return graph, textfile, treefile
##################################################################################################################################################################################
#parse the command line arguments
graph, textfile, treefile = parse_arguments()

#read in ground truth tree topology and mutations in its clones
G_Dir = nx.read_gpickle(graph)
f = open(textfile,'r')             
lines = f.readlines()                        
new_mutation_in_clone = ast.literal_eval(lines[2].rstrip('\n').lstrip('\r'))
f.close()

#acessment by comparing the ground truth tree with the predicted tree
result = LICHeE_Acessment(G_Dir,new_mutation_in_clone,tree=treefile)           

#print out results            
print result.to_string(index=False)

