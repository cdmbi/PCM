# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 12:18:21 2014

@author: Fujitsu
"""

def PCA(X, Expect_ext):
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    Xtrans = pca.fit_transform(X)

    labels = []
    for Var in Xtrans:
        if Var[0] >= 0 and Var[1] >= 0:
            labels.append(0)
        elif Var[0] < 0 and Var[1] >= 0:
            labels.append(1)
        elif Var[0] < 0  and Var[1] < 0:
            labels.append(2)
        elif Var[0] >= 0 and Var[1] < 0:
            labels.append(3)
    return labels
 

def KMean(X, Expect_ext):
    from sklearn.cluster import KMeans       
    k_means = KMeans(n_clusters=Expect_ext)
    k_means.fit(X)
    return k_means.labels_
    
def HierachicalClustering(X, Expect_ext):
    from sklearn.cluster import Ward       
    HC = Ward(n_clusters=Expect_ext)
    HC.fit(X)
    return HC.labels_
    
def CluterAnalysis(labels, Criteria, Expect_ext):
    import numpy as np
    import random, collections
    Ind_ext = []
    counter = collections.Counter(labels)
    c_high_f = counter.most_common()
    print c_high_f
    
    while len(Ind_ext) <= Expect_ext:
        for i in c_high_f:
            Segment = [ind for ind,val in enumerate(labels) if val == i[0]]
            R_segment = random.sample(Segment, int(np.round(Criteria*len(Segment))))
            Ind_ext.extend(R_segment)
            
        for j in Ind_ext:
            labels = [val for ind,val in enumerate(labels) if ind != j]
        
    if len(Ind_ext) > Expect_ext:
        Ind_ext = random.sample(Ind_ext, Expect_ext)
    return Ind_ext
    
def CV_determination(Y, Method):
    from sklearn.cross_validation import KFold, LeaveOneOut    
    if Method == 'loo':
        kf = LeaveOneOut(len(Y))
    else:
        ind_k = [ind for ind,val in enumerate(list(Method)) if val == '-']
        k = int(Method[:ind_k[0]])
        kf = KFold(len(Y), k, indices=True, shuffle=True, random_state=1)
    return kf