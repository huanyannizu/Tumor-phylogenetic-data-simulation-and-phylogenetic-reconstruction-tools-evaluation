import numpy as np
import pandas as pd
import networkx as nx
import itertools
import ast
import argparse

def Treeomics_Acessment(G_Dir,new_mutation_in_clone,columns):
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
    
    unique, counts = np.unique(muta_relation_inp, return_counts=True)
    relation_counts = dict(zip(unique, counts))
        
    f = open(columns).readlines()
    muta_group = {}
    
    G_outp_nx = nx.DiGraph()
    for i in range(len(f)):
        if 'end{itemize}' in f[i]:break
        if 'to' in f[i]:
            k1 = f[i].split(' to ')[-2].split('{')[-1]
            k2 = f[i].split(' to ')[-1].split('}')[0]
            if k1 != 'Germline Data':G_outp_nx.add_edge(k1,k2)
            v = []
            for i in f[i+1].split(' '):
                if '>' in i:
                    v.append(int(i.split('__')[1]))
            muta_group[k2] = v
    
    muta_fraction = 0
    for k in muta_group.keys():
        muta_fraction = muta_fraction + len(muta_group[k])            
  
    muta_relation_outp = np.zeros((n_mutation,n_mutation))
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
    muta_com = itertools.combinations(range(n_mutation),2)
    muta_com = list(muta_com)
    for pair in muta_com:
        muta1 = pair[0];muta2 = pair[1]
        if muta_relation_inp[muta1,muta2] == 1:
            if muta_relation_outp[muta1,muta2] == 1:Acc['%Corr AD']+=1
            
        if muta_relation_inp[muta2,muta1] == 1:
            if muta_relation_outp[muta2,muta1] == 1:Acc['%Corr AD']+=1
    
    Acc['%Corr AD'] = Acc['%Corr AD']/relation_counts[1.0]
    
    return Acc

def parse_arguments():   
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", type=str)
    parser.add_argument("textfile", type=str)
    parser.add_argument("outfile", type=str)
    
    args = parser.parse_args()
    
    graph = args.graph
    textfile = args.textfile
    outfile = args.outfile
    
    return graph, textfile, outfile
##################################################################################################################################################################################
#parse the command line arguments
graph, textfile, outfile = parse_arguments()

#read in ground truth tree topology and mutations in its clones
G_Dir = nx.read_gpickle(graph)
f = open(textfile,'r')             
lines = f.readlines()                        
new_mutation_in_clone = ast.literal_eval(lines[2].rstrip('\n').lstrip('\r'))
f.close()

#acessment by comparing the ground truth tree with the predicted tree
result = Treeomics_Acessment(G_Dir,new_mutation_in_clone,columns=outfile)      

#print out results            
print result.to_string(index=False)
