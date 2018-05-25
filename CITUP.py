import numpy as np
import pandas as pd
import networkx as nx
import pygraphviz as pgv
import itertools
import subprocess
import time
import shutil
import os
import ast

def Data_CITUP(F,n_mutation,n_sample,error_rate,name):
    Data_CITUP = F*2 
    header='Num_mutations: %d\nNum_samples: %d\nError_rate: %f\n'%(n_mutation,n_sample,error_rate)
    os.mkdir('/u/85/lix10/unix/Documents/%s'%name)
    np.savetxt('/u/85/lix10/unix/Documents/%s/Frequencies.txt'%name, Data_CITUP, delimiter=' ', newline='\n', header=header,comments='')
    return Data_CITUP
########################################################################################    
def CITUP_Acessment(G_Dir,new_mutation_in_clone,n_mutation,dot,muta): 
    muta_relation_inp = np.zeros((n_mutation,n_mutation))#remember the table must be read from rows for correct relation     
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
    
    unique, counts = np.unique(muta_relation_inp, return_counts=True)
    relation_counts = dict(zip(unique, counts))
    
    G_outp = pgv.AGraph(dot)
    G_outp_nx = nx.DiGraph()

    for edge in G_outp.edges():
        #if edge[0][0] != 's' and edge[0][0] != 'd' and edge[1][0] != 's' and edge[1][0] != 'd':
        G_outp_nx.add_edge(int(edge[0]),int(edge[1]))
        
    f = open(muta).readlines()
    for line in f:
        if len(line) == 2+n_mutation+(n_mutation-1)*3+2:break
    
    line = line[2:-2].split('  ')
    line = map(int, line)
        
    muta_group = dict([(key, []) for key in G_outp_nx.nodes()])  
    for i in range(len(line)):
        muta_group[line[i]].append(i)
            
    muta_fraction = 0
    for k in muta_group.keys():
        muta_fraction = muta_fraction + len(muta_group[k])
        
    muta_relation_outp = np.zeros((n_mutation,n_mutation))#only read from row,1 is AD, 2 is Sib,0 is no relation
    node_com = itertools.combinations(G_outp_nx.nodes()[1:],2)
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
        if G_outp_nx.predecessors(node1) == G_outp_nx.predecessors(node2):
            for i in muta_group[node1]:
                for j in muta_group[node2]:
                    muta_relation_outp[i,j] = 2

    if 1 not in muta_relation_outp:return 'no AD'   
        
    Acc = pd.DataFrame(np.zeros((1,7)),columns=['%Corr(AD+Sib)','%Corr AD','%AD reverse','%Corr Sib','%AD->Sib','%Sib->AD','%SSNVs'])    
    muta_com = itertools.combinations(range(n_mutation),2)
    muta_com = list(muta_com)
    for pair in muta_com:
        muta1 = pair[0];muta2 = pair[1]
        if muta_relation_inp[muta1,muta2] == 1:
            if muta_relation_outp[muta1,muta2] == 1:Acc['%Corr AD']+=1;Acc['%Corr(AD+Sib)']+=1
            if muta_relation_outp[muta2,muta1] == 1:Acc['%AD reverse']+=1
            if muta_relation_outp[muta1,muta2] == 2 or muta_relation_outp[muta2,muta1] == 2:Acc['%AD->Sib']+=1
        if muta_relation_inp[muta2,muta1] == 1:
            if muta_relation_outp[muta2,muta1] == 1:Acc['%Corr AD']+=1;Acc['%Corr(AD+Sib)']+=1
            if muta_relation_outp[muta1,muta2] == 1:Acc['%AD reverse']+=1
            if muta_relation_outp[muta1,muta2] == 2 or muta_relation_outp[muta2,muta1] == 2:Acc['%AD->Sib']+=1
        if muta_relation_inp[muta1,muta2] == 2 or muta_relation_inp[muta2,muta1] == 2:
            if muta_relation_outp[muta1,muta2] == 1 or muta_relation_outp[muta2,muta1] == 1:Acc['%Sib->AD']+=1
            if muta_relation_outp[muta1,muta2] == 2 or muta_relation_outp[muta2,muta1] == 2:Acc['%Corr Sib']+=1;Acc['%Corr(AD+Sib)']+=1
            
    Acc['%Corr AD'] = Acc['%Corr AD']/relation_counts[1.0]
    Acc['%AD reverse'] = Acc['%AD reverse']/relation_counts[1.0]
    Acc['%AD->Sib'] = Acc['%AD->Sib']/relation_counts[1.0]    
    Acc['%SSNVs'] = muta_fraction/(n_mutation+0.0)
    if 2 in relation_counts:
        Acc['%Corr(AD+Sib)'] = Acc['%Corr(AD+Sib)']/(relation_counts[1.0]+relation_counts[2.0])
        Acc['%Corr Sib'] = Acc['%Corr Sib']/relation_counts[2.0]
        Acc['%Sib->AD'] = Acc['%Sib->AD']/relation_counts[2.0]
    else:
        Acc['%Corr(AD+Sib)'] = 0
        Acc['%Corr Sib'] = 0
        Acc['%Sib->AD'] = 0   
        
    Acc2 = pd.DataFrame(np.zeros((1,4)),columns=['Precision','Recall','True Negative Rate','Accuracy'])  
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
#    print tp,fn,fp,tn     
    result = pd.concat([Acc, Acc2], axis=1)
    return result


##Results for CITUP###########################################################################################  
'''       
n_clone = 10;n_mutation = 100;n_trees = 3;mutation = range(n_mutation);mutation = np.array(mutation)
Results_best_CITUP = pd.DataFrame();Results_avg_CITUP = pd.DataFrame();Results_ignored = pd.DataFrame()
start_time = time.time()
for n_sample in [5]:
    for read_cov in [100]: 
        CITUP_avg_ACC = pd.DataFrame();CITUP_best_ACC = pd.DataFrame();ignored = []       
        for i in range(1,n_trees+1):
            print('loop%d'%i)
            F_unpack_noise = np.loadtxt('/u/85/lix10/unix/Documents/%d_%d/%d_%d_%d_F_unpack_noise.txt'%(n_sample,read_cov,n_sample,read_cov,i))
            G_Dir = nx.read_gpickle('/u/85/lix10/unix/Documents/%d_%d/%d_%d_%d_G_Dir.gpickle'%(n_sample,read_cov,n_sample,read_cov,i))
            f = open('/u/85/lix10/unix/Documents/%d_%d/%d_%d_%d.txt'%(n_sample,read_cov,n_sample,read_cov,i),'r')             
            lines = f.readlines()                        
            new_mutation_in_clone = ast.literal_eval(lines[2].rstrip('\n').lstrip('\r'))
            #Clone = ast.literal_eval(lines[3].rstrip('\n').lstrip('\r'))
            #samples = ast.literal_eval(lines[4].rstrip('\n').lstrip('\r'))
            f.close()
            
            if read_cov == 100:
                error_rate = 0.08
                Data_CITUP(F_unpack_noise,n_mutation,n_sample,error_rate)
            if read_cov == 1000:
                error_rate = 0.05
                Data_CITUP(F_unpack_noise,n_mutation,n_sample,error_rate)
            if read_cov == 10000:
                error_rate = 0.03
                Data_CITUP(F_unpack_noise,n_mutation,n_sample,error_rate)
            
            in_start_time = time.time()
            subprocess.call('bash runCITUP.sh /u/85/lix10/unix/Documents/for-citup 10',shell=True,cwd='/u/85/lix10/unix/Documents/citup/bin',timeout=1000)             
            
            #outp=subprocess.check_output('bash runCITUP.sh /u/85/lix10/unix/Documents/for-citup 10',shell=True,cwd='/u/85/lix10/unix/Documents/citup/bin',timeout=1000)             
            #subprocess.call('perl ./visualizeResults.pl -r /u/85/lix10/unix/Documents/for-citup -g /u/85/lix10/unix/Documents/citup/GammaAdjMatrices -o /u/85/lix10/unix/Documents/for-citup/results_summary',shell=True,cwd='/u/85/lix10/unix/Documents/citup/bin')
            print('software run time: %s seconds'%(time.time() - in_start_time))
            
            
            shutil.rmtree('/u/85/lix10/unix/Documents/for-citup')
            
            best=[]
            for i in range(len(outp)):
                if outp[i] == 'R':
                    best.append(outp[i+8:i+11])
            
            no_AD = 0;onedata = pd.DataFrame()
            onedata_best = pd.DataFrame(np.zeros((1,1)),columns=['%Corr AD'])
            
            for t in best:
                dot='/u/85/lix10/unix/Documents/for-citup/results_summary/dots/node%d_tree%d.0.dot'%(int(t[0]),int(t[2]))
                muta='/u/85/lix10/unix/Documents/for-citup/Results_%d_%d.txt'%(int(t[0]),int(t[2]))
                CITUP_acc = CITUP_Acessment(G_Dir,new_mutation_in_clone,n_mutation,dot=dot,muta=muta)
                if type(CITUP_acc) != str:
                    onedata = pd.concat([onedata, CITUP_acc])
                    if CITUP_acc['%Corr AD'].item() >= onedata_best['%Corr AD'].item():onedata_best = CITUP_acc
                else:no_AD = no_AD+1 
            
            shutil.rmtree('/u/85/lix10/unix/Documents/for-citup')
            
            no_AD = no_AD/(len(best)+0.0)
            if no_AD == 1:continue
            else:onedata_avg = onedata.groupby(level=0).mean()    
                        
            CITUP_best_ACC = pd.concat([CITUP_best_ACC, onedata_best]) 
            CITUP_avg_ACC = pd.concat([CITUP_avg_ACC, onedata_avg])
            ignored.append(no_AD)
            
        CITUP_best_ACC = CITUP_best_ACC.groupby(level=0).mean()
        CITUP_avg_ACC = CITUP_avg_ACC.groupby(level=0).mean()
        ignored = sum(ignored)/(len(ignored)+0.0)
       
        #Results_best_CITUP = pd.concat([Results_best_CITUP, CITUP_best_ACC])
        #Results_avg_CITUP = pd.concat([Results_avg_CITUP, CITUP_avg_ACC])
        #Results_ignored = pd.concat([Results_ignored, ignored])
        
print 'total run time: %s seconds'%(time.time() - start_time)
print CITUP_best_ACC
print CITUP_avg_ACC
print ignored
#print Results_best_CITUP
#print Results_avg_CITUP
#print Results_ignored

'''