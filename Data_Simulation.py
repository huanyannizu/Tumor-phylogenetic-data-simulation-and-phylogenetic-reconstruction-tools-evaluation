import numpy as np
import pandas as pd
import networkx as nx
import random
import os
import copy

def parse_arguments():    
    import argparse
    import ast
    parser = argparse.ArgumentParser()
    parser.add_argument("repSaveData", type=str)
    parser.add_argument("n_clone", type=int)
    parser.add_argument("n_mutation", type=int)
    parser.add_argument("cell_size_min", type=int)
    parser.add_argument("cell_size_max", type=int)
    parser.add_argument("n_sample", type=int)
    parser.add_argument("read_cov", type=int)
    parser.add_argument("d",type=int)
    
    args = parser.parse_args()
    
    repSaveData = args.repSaveData
    n_clone = args.n_clone
    n_mutation = args.n_mutation
    cell_size_min = args.cell_size_min
    cell_size_max = args.cell_size_max
    n_sample = args.n_sample
    read_cov = args.read_cov
    d = args.d
    
    return repSaveData,n_clone,n_mutation,cell_size_min,cell_size_max,n_sample,read_cov,d

###Functions of creating VAF values################################################################################################################################################################
def prufer_to_tree(a):
    tree = [];T = range(0, len(a)+2);deg = [1]*len(T)
    for i in a: deg[i] += 1
    for i in a:
        for j in T:
            if deg[j] == 1:
                tree.append((i,j));deg[i] -= 1;deg[j] -= 1
                break

    last = [x for x in T if deg[x] == 1];tree.append((last[0],last[1]))
    return tree

def clone_size(n, total):
    dividers = sorted(random.sample(xrange(1, total), n - 1))
    return [a - b for a, b in zip(dividers + [total], [0] + dividers)]

def assign_mutation(mutation,cz):
    Clone = {};i = -1 #clones are marked from 0 to 9.
    for s in cz:
        i = i + 1
        clone = np.random.choice(mutation,s,replace=False)
        Clone[i] = list(clone)
        mutation = set(mutation) - set(clone)
        mutation = np.array(list(mutation))
    
    return Clone

def get_samples(n_sample,n_clone):	
    samples = {} 
    for i in range(n_sample):
        c = np.random.randint(2,5) 
        selected_clones = random.sample(range(n_clone), c)
        samples[i] = selected_clones
    
    return samples

def tree(n_clone,n_mutation,cell_size_min,cell_size_max): 
    mutation = range(n_mutation);mutation = np.array(mutation)        
    ps = np.random.randint(0,n_clone,size = n_clone-2) #creating a prufer sequence
    edges = prufer_to_tree(ps)
    root = np.random.randint(0,n_clone)
    
    cell_size = np.random.randint(cell_size_min,cell_size_max,size = n_clone)
    cell_size = dict(enumerate(cell_size))

    return root,edges,cell_size

def VAF(n_sample,n_clone,cell_size,n_mutation,read_cov,root,edges,d):
    samples = get_samples(n_sample,n_clone)
    mutation = range(n_mutation);mutation = np.array(mutation) 
    cz = clone_size(n_clone, n_mutation)

    #variable new_mutation_in_clone saved new mutations occured in each clone
    #variable Clone saved mutations in each clone after mutation inheritance
    Clone = assign_mutation(mutation,cz) 
    new_mutation_in_clone = copy.deepcopy(Clone)
    
    G=nx.Graph();G.add_nodes_from(range(n_clone));G.add_edges_from(edges)
    p=nx.shortest_path(G,source=root) #paths from root to each node
    edges_with_label = {}
       
    for node in G.nodes(): #mutation inheritance
        pn = p[node] 
        for k in pn[:-1]:
            Clone[node] = Clone[node] + Clone[k]
        Clone[node] = list(set(Clone[node]))        
        for i in range(1,len(pn)):
            edges_with_label[(pn[i-1],pn[i])] = len(new_mutation_in_clone[pn[i]])
    
    G_Dir=nx.DiGraph();edges_with_label_dir = {}
    for node in G.nodes():
        pn = p[node]
        for j in range(len(pn) - 1):
            G_Dir.add_edge(pn[j], pn[j+1]) 
            edges_with_label_dir[(pn[j], pn[j+1])] = len(new_mutation_in_clone[pn[j+1]])
    
    nodes = range(n_clone);nodes.remove(root)    
    D = np.random.choice(nodes, d, replace=False);mutas = []
    for node in D:
        muta=np.random.choice(list(set(Clone[node])-set(new_mutation_in_clone[node])))
        mutas.append(muta)
        #print 'muta:',muta
        #print muta
        Clone[node].remove(muta)
        descendants = nx.descendants(G_Dir, node)
        for de in descendants:
            if muta in Clone[de]:
                Clone[de].remove(muta)
    
    U = np.zeros((n_sample,n_clone))
    for s in samples.keys():
        c = samples[s]
        for x in c:
            U[s,x] = cell_size[x]
        U[s] = U[s]/sum(U[s])
            
    B = np.zeros((n_clone,n_clone))        
    for c in Clone.keys():
        for i in new_mutation_in_clone.keys():
            if set(new_mutation_in_clone[i]).issubset(set(Clone[c])):
                B[c,i] = 1
    
    F = np.dot(U, B)*0.5
    
    F_unpack = np.zeros((n_sample,n_mutation))
    for i in range(n_clone):
        m = new_mutation_in_clone[i]
        for j in m:
            for s in range(n_sample):
                F_unpack[s,j] = F[s,i]
        
    n = np.random.poisson(read_cov, n_mutation*n_sample).reshape((n_mutation, n_sample))
    x = np.zeros((n_mutation,n_sample))
    F_unpack_noise = np.zeros((n_mutation,n_sample))
    for i in mutation:
        for p in samples.keys():
            x[i,p] = np.random.binomial(n[i,p], F_unpack[p,i])
            F_unpack_noise[i,p] = x[i,p] / (n[i,p]+0.0)
              
    n_x = np.subtract(n, x)
    	
    read_count = np.zeros((n_mutation,1))
    for j in range(n_sample):
        read_count = np.c_[read_count,n_x[:,j]]
        read_count = np.c_[read_count,x[:,j]]
    read_count = read_count[:,1:]
    
    return U,B,F,F_unpack,n,x,F_unpack_noise,samples,G_Dir,new_mutation_in_clone,Clone


###Functions of creating input for each tool#####################################################################################################################################
def Data_MIPUP(F,n_mutation,n_sample):
    Data_MIPUP = F
    row_name = ['#chr'] + range(n_mutation)
    row_name = np.array(row_name).reshape(n_mutation+1,1)
    norm = np.array(['normal'] + [0]*n_mutation).reshape(n_mutation+1,1)
    
    col_name = []
    for i in range(n_sample):
        col_name.append("sample%d" % i)
    col_name = np.array([col_name])
    
    Data_MIPUP = np.append(col_name,Data_MIPUP,axis=0)
    Data_MIPUP = np.append(norm,Data_MIPUP,axis=1)   
    Data_MIPUP = np.append(row_name,Data_MIPUP,axis=1)
    Data_MIPUP = pd.DataFrame(Data_MIPUP) 
    return Data_MIPUP 

def Data_LICHeE(F,n_mutation,n_sample):
    Data_LICHeE = F
    row_name = ['#chr'] + range(1,n_mutation+1)
    row_name = np.array(row_name).reshape(n_mutation+1,1)
    norm = np.array(['normal'] + [0]*n_mutation).reshape(n_mutation+1,1)
    description = np.array(['description'] + [1]*n_mutation).reshape(n_mutation+1,1)
    position = np.array(['position'] + [1]*n_mutation).reshape(n_mutation+1,1)
    col_name = []
    for i in range(n_sample):
        col_name.append("sample%d" % i)
    col_name = np.array([col_name])
    
    Data_LICHeE = np.append(col_name,Data_LICHeE,axis=0)
    Data_LICHeE = np.append(norm,Data_LICHeE,axis=1)
    Data_LICHeE = np.append(description,Data_LICHeE,axis=1)
    Data_LICHeE = np.append(position,Data_LICHeE,axis=1)    
    Data_LICHeE = np.append(row_name,Data_LICHeE,axis=1)
    Data_LICHeE = pd.DataFrame(Data_LICHeE)  
    return Data_LICHeE

def Data_AncesTree(n,x,n_mutation,n_sample):
    n_x = np.subtract(n, x)
    	
    read_count = np.zeros((n_mutation,1))
    for j in range(n_sample):
        read_count = np.c_[read_count,n_x[:,j]]
        read_count = np.c_[read_count,x[:,j]]
    read_count = read_count[:,1:]
    
    Data_AncesTree = read_count
    Data_AncesTree = Data_AncesTree.astype(int)    
    row_name = ['gene_id'] + [i for i in range(n_mutation)]
    row_name = np.array(row_name).reshape(n_mutation+1,1)
    col_name = []
    for i in range(n_sample):
        col_name.append("sample%d" % i)
        col_name.append("sample%d" % i)
    col_name = np.array([col_name])
    
    Data_AncesTree = np.append(col_name,Data_AncesTree,axis=0)
    Data_AncesTree = np.append(row_name,Data_AncesTree,axis=1)
    Data_AncesTree = pd.DataFrame(Data_AncesTree)  
    return Data_AncesTree

def Data_CITUP(F,n_mutation,n_sample,error_rate):
    Data_CITUP = F * 2  
    return Data_CITUP

def Data_Treeomics(n,x,n_mutation,n_sample):
    cov = n.astype(int)
    mu = x.astype(int)
    
    row_name = ['Change']
    for i in range(n_mutation):
        m = np.random.choice(['G>A','G>C','G>T','A>G','A>C','A>T','C>G','C>A','C>T','T>G','T>A','T>C'])
        row_name = row_name + [m]
    
    row_name = np.array(row_name).reshape(n_mutation+1,1)
    col_name = []
    for i in range(n_sample):
        col_name.append("sample%d" % i)
    col_name = np.array([col_name])
    
    Gene = np.array(['Gene'] + [1]*n_mutation).reshape(n_mutation+1,1)
    Chromosome = np.array(['Chromosome'] + [1]*n_mutation).reshape(n_mutation+1,1)
    Position = np.array(['Position'] + [i for i in range(n_mutation)]).reshape(n_mutation+1,1)
       
    cov = np.append(col_name,cov,axis=0)
    cov = np.append(Gene,cov,axis=1)
    cov = np.append(row_name,cov,axis=1)
    cov = np.append(Position,cov,axis=1)
    cov = np.append(Chromosome,cov,axis=1)   
    cov = pd.DataFrame(cov)  
     
    mu = np.append(col_name,mu,axis=0)
    mu = np.append(Gene,mu,axis=1)
    mu = np.append(row_name,mu,axis=1)
    mu = np.append(Position,mu,axis=1)
    mu = np.append(Chromosome,mu,axis=1)   
    mu = pd.DataFrame(mu)  
    return cov,mu


######################################################################################################################################################################
repSaveData,n_clone,n_mutation,cell_size_min,cell_size_max,n_sample,read_cov,d = parse_arguments()
mutation = range(n_mutation); mutation = np.array(mutation)

root,edges,cell_size = tree(n_clone,n_mutation,cell_size_min,cell_size_max)
U,B,F,F_unpack,n,x,F_unpack_noise,samples,G_Dir,new_mutation_in_clone,Clone = VAF(n_sample,n_clone,cell_size,n_mutation,read_cov,root,edges,d)

data_mipup=Data_MIPUP(F_unpack_noise,n_mutation,n_sample)
data_lichee=Data_LICHeE(F_unpack_noise,n_mutation,n_sample)
data_ancestree=Data_AncesTree(n,x,n_mutation,n_sample) 
if read_cov == 100:
    error_rate = 0.08
    data_citup=Data_CITUP(F_unpack_noise,n_mutation,n_sample,error_rate)
if read_cov == 1000:
    error_rate = 0.05
    data_citup=Data_CITUP(F_unpack_noise,n_mutation,n_sample,error_rate)
if read_cov == 10000:
    error_rate = 0.03
    data_citup=Data_CITUP(F_unpack_noise,n_mutation,n_sample,error_rate)    
cov,mu=Data_Treeomics(n,x,n_mutation,n_sample)
                       
data_mipup.to_csv(repSaveData + '/Data_MIPUP.txt',sep = '\t',header = False,index = False)
data_lichee.to_csv(repSaveData + '/Data_LICHeE.txt',sep = '\t',header = False,index = False)
data_ancestree.to_csv(repSaveData + '/Data_AncesTree.txt',sep = '\t',header = False,index = False)
header='Num_mutations: %d\nNum_samples: %d\nError_rate: %f\n'%(n_mutation,n_sample,error_rate)
np.savetxt(repSaveData + '/Data_CITUP.txt', data_citup, delimiter=' ', newline='\n', header=header,comments='')
cov.to_csv(repSaveData + '/Data_Treeomics_n.txt',sep = '\t',header = False,index = False)
mu.to_csv(repSaveData + '/Data_Treeomics_x.txt',sep = '\t',header = False,index = False)

            
                        