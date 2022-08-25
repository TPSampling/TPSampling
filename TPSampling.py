#-*- coding:Utf-8 -*-
__author__ = "DIOP Lamine BSF"

import os
import time
import sys
from random import randrange, random
import  collections
import statistics 
import math
from asizeof.asizeof import asizeof


#######################################################################################################

class TrieNodeBis():
    def __init__(self):
        self.children = {}
        self.weight= {}


def combin(n, k):
    if k>n: return 0 #LDIOPBSF
    if k==0 and n==0: return 1 #LDIOPBSF
    if k > n//2:
        k = n-k
    x = 1
    y = 1
    i = n-k+1
    while i <= n:
        x = (x*i)//y
        y += 1
        i += 1
    return x



def trouver(tab,i,j,x):
    m=int((i+j)/2)
    if m==0 or (tab[m-1]<x and x<=tab[m]):
        return m
    if tab[m]<x:
        return trouver(tab,m+1,j,x)
    return trouver(tab,i,m,x)


def utility(name,alpha, norme):
    if name=="Freq":
        return 1
    elif name=="Area":
        return norme
    elif name=="Decay":
        return pow(alpha,norme)
    elif name=="Gauss":
        return Gauss(m,M,norme)
    else:
        print("Error of utility definition!")
        sys.exit(1)
        
def Gauss(m,M,l):
    mu=(M+m)/2 # loi normale centré en mu
    sigmaCarre=mu-m#(M-m)/2
    
    x=1/(math.sqrt(sigmaCarre*2*math.pi))
    y=(l-mu)**2/(2*sigmaCarre)
    return x*math.exp(-y)


#######################################################################################################

def constructionCnk(maxL,M):
    tabCnk=[]
    for n in range(maxL+1):
        tabCnk.append([combin(n, k) for k in range(M+1)])
    return tabCnk


def sortAccordingtoItems(list_a, list_b):
    return sorted(list_a, key=lambda x: list_b.index(x))

def lexicOrderItems(dataset):
    s=[]
    maxL=0
    with open(dataset, 'r') as base:
            line=base.readline()
            while line:
                itemset = line.replace(" \n","").replace("\n","").split(" ")
                n=len(itemset)
                if maxL<n:
                    maxL=n
                for item in itemset:
                    if item not in s:
                        s.append(item)
                line=base.readline()
    base.close()
    del base
    s.sort()
    return s,maxL


def freqOrderItems(dataset):
    dictionnary = {}
    maxL=0
    with open(dataset, 'r') as base:
            line=base.readline()
            while line:
                itemset = line.replace(" \n","").replace("\n","").split(" ")
                n=len(itemset)
                if maxL<n:
                    maxL=n
                for item in itemset:
                    try:
                        dictionnary[item]+=1
                    except KeyError:
                        dictionnary[item]=1
                line=base.readline()
    base.close()
    del base
    s=collections.OrderedDict(sorted(dictionnary.items(), key=lambda t: t[1]))
    s=list(s.keys())
    s.reverse()
    del dictionnary
    return s,maxL

def addTrans(root, trans, m, M, tabCnk, nbNode):
    M=min(M, len(trans))
    node = root
    l = len(trans)
    for i in range(m, min(M,l+1)+1):
        try:
            node.weight[i][0]+=tabCnk[l][i]
        except KeyError:
            node.weight[i]=[tabCnk[l][i]]
    l-=1
    for item in trans:
        if not item in node.children:
            node.children[item] = TrieNodeBis()
            nbNode+=1
        node = node.children[item]
        for i in range(m, min(M,l+1)+1):
            try:
                node.weight[i][0]+=tabCnk[l][max(0,i-1)]
                node.weight[i][1]+=tabCnk[l][i]
            except KeyError:
                node.weight[i]=[tabCnk[l][max(0,i-1)], tabCnk[l][i]]
        l-=1
    return nbNode
    

def TPSpace(dataset, m, Ml, ordre):
    nbNode=0
    litterals, maxL = [], 0
    if ordre==0:
        litterals, maxL = lexicOrderItems(dataset) # freqOrderItems  lexicOrderItems
    else:
        litterals, maxL = freqOrderItems(dataset) 
    #print(items)
    M=min(Ml, maxL)
    tabCnk = constructionCnk(maxL,M)
    root = TrieNodeBis()
    with open(dataset, 'r') as base:
            line=base.readline()
            while line:
                trans=line.replace(" \n","").replace("\n","").split(" ")
                nbNode=addTrans(root, sortAccordingtoItems(trans, litterals), m, M,tabCnk, nbNode)
                line=base.readline()
    base.close()
    del litterals
    del tabCnk
    return root, nbNode

#######################################################################################################
    
def TPSampling(root, x, normPattern):
    node = root
    pattern=[]
    while not normPattern==0:
        child=""
        for child in node.children.keys():
            try:
                w=sum(node.children[child].weight[normPattern])
                if w < x:
                    x -= w
                else:
                    break
            except KeyError:
                pass
        node = node.children[child]
        if node.weight[normPattern][1] < x:
            pattern.append(child)
            x -= node.weight[normPattern][1]
            normPattern -= 1
    return pattern


def sampleBasedTimer(sample, runningTime, weightedNorm, tabNorme):
    beginTime =time.process_time()
    while time.process_time() < beginTime + runningTime:
        alea=weightedNorm[-1]*random()
        j = trouver(weightedNorm, 0, len(weightedNorm), alea)
        l=tabNorme[j]
        x= randrange(1,root.weight[l][0]+1)
        pattern = TPSampling(root, x, l)
        if str(pattern) in sample:
            sample[str(pattern)]+=1
        else:
            sample[str(pattern)]=1
    return sample
            
#######################################################################################################
            
            
if __name__ == "__main__":
    print("====== Thanks for using TPSampling ======")
    ################################################################
    database=["adult", "connect", "mushroom", "chess", "pumsb", "uscensus", "susy"]#
    measures=[ "Area", "Decay", "Freq", "Gauss"]
    ################################################################
    N = 1000
    m, M = 1, 6
    alpha = 0.1
    measure = measures[2]
    dataset = database[0]
    ordre = 1 #0=lexicalOrder, 1=freqOrder
    runningTime = 0
    
    print("Dataset :",dataset)
    print("Interestingness measure :",measure)
    print("Norm constraint :","["+str(m)+".."+str(M)+"]")
    if measure=="Decay":
        print("Exponential decay : alpha="+str(alpha))
    dataset = "Datasets/"+dataset+".num"
    if ordre==1:
        print("Total order relation :", "Frequency order")
    elif ordre==0:
        print("Total order relation :", "Lexicographic order")
        
    beginTime =time.process_time()
    root, nbNode = TPSpace(dataset, m, M, ordre)
    weightedNorm = []
    w=0
    for norme in root.weight:
        w+=root.weight[norme][0]*utility(measure, alpha, norme)
        weightedNorm.append(w)
    del w
    tabNorme=list(root.weight.keys())
    endTime =time.process_time()-beginTime
    print("Preprocessing time :",endTime,"s")
    print("Number of nodes in the trie :",nbNode)
    print("Trie size in main memory :", round(asizeof(root)/(1024*1024),3),"MB") 
        
    sample ={} 
    if runningTime==0:
        print("Sample size :",N)
        beginTime =time.process_time()
        for k in range(N): 
            alea=weightedNorm[-1]*random()
            j = trouver(weightedNorm, 0, len(weightedNorm), alea)
            l=tabNorme[j]
            x= randrange(1,root.weight[l][0]+1)
            pattern = TPSampling(root, x, l)
            if str(pattern) in sample:
                sample[str(pattern)]+=1
            else:
                sample[str(pattern)]=1
        endTime = time.process_time() - beginTime
    else:
            beginTime =time.process_time()
            sample = sampleBasedTimer(sample, runningTime, weightedNorm, tabNorme)
            endTime = time.process_time() - beginTime
    print("Sampling time :",endTime,"s")
    print("Distinct sampled patterns :",len(sample))
    #print(sample)
    del sample
    del root
    print("================== END ==================")
    
    
