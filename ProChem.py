# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 09:55:28 2013

@author: wiwat
"""

def GetProteinData():
    import os, csv
    fileName = os.getcwd()+'/'+'Protein_Descriptors.csv'
    
    with open(fileName,'rb') as csvfile:
        dialect = csv.Sniffer().has_header(csvfile.read())
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect) 
        h = next(reader)
        Data = []
        for row in reader:
            Data.append(row)
    Abbre, Des=[],[]
    for i in range(len(Data)):
        Protein = Data[i]
        Des.append(Protein[3:])
        Abbre.append(Protein[1])
    return Des, Abbre, h
    
def GetFileProteinSequence(f, Des, Abbre):
    import csv
    import numpy as np
    
    with open(f,'rb') as csvfile:
        dialect = csv.Sniffer().has_header(csvfile.read())
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect) 
        data = []
        for row in reader:
            data.append(row)
            
    listProtein = []
    for i in range(len(data)):
        Protein = data[i]
        ProteinSeq = Protein[0]
        
        Pseq = []
        for j in range(len(ProteinSeq)):
            index_Protein = [ind for ind,val in enumerate(Abbre) if ProteinSeq[j]==Abbre[ind]]
            Seq = Des[index_Protein[0]]
            
            for k in range(len(Seq)):
                Pseq.append(Seq[k])
        listProtein.append(Pseq)
    return np.array((listProtein), dtype=np.float)
    
def GetProteinSeq(Protein_var, Des, Abbre):
    import numpy as np
    Protein_sequence = []
    for i in range(len(Protein_var)):
        ProSeq = []
        for ii in Protein_var[i]:
            if ii == '-':
                Seq = ['0','0','0','0','0']
            else:
                index_Protein = [ind for ind,val in enumerate(Abbre) if ii==Abbre[ind]]
                Seq = Des[index_Protein[0]]
            ProSeq.append(Seq)
        Protein_sequence.append(ProSeq)
        
    P = np.array(Protein_sequence, dtype=np.float)
    return np.reshape(P, (P.shape[0], P.shape[1]*P.shape[2]))
    
def GetLigand(f):
    import csv 
    import numpy as  np
    
    with open(f,'rb') as csvfile:
        dialect = csv.Sniffer().has_header(csvfile.read())
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect) 
        data = []
        for row in reader:
            data.append(row)
            
    D = np.array(data[1:], dtype=np.float)
    return D[:,:-1]
    
def Combine_Protein_Protein(protein, header_protein):
    import numpy as np
    nr, nc = np.shape(protein)
    Cross_protein = np.zeros((nr,1))
    H_pro = []
    for i in range(nc):
        initial_h = header_protein[0]
        initial_protein = protein[:,0]
        for ii in range(len(header_protein)-1):
            H_pro.append(initial_h+'*'+header_protein[ii+1])
            cross = np.multiply(initial_protein, protein[:,ii+1])
            Cross_protein = np.append(Cross_protein,np.reshape(cross,(nr,1)),axis=1)
        protein = np.delete(protein,0,axis=1)
        header_protein.pop(0)
    return np.delete(Cross_protein,0,axis=1), H_pro
    

def Combine_Ligand_Ligand(ligand, header_ligand):
    import numpy as np
    nr, nc = np.shape(ligand)
    Cross_ligand = np.zeros((nr,1))
    H_li = []
    for i in range(nc):
        initial_h = header_ligand[0]
        initial_ligand =ligand[:,0]
        for ii in range(len(header_ligand)-1):
            H_li.append(initial_h+'*'+header_ligand[ii+1])
            cross = np.multiply(initial_ligand, ligand[:,ii+1])
            Cross_ligand = np.append(Cross_ligand,np.reshape(cross,(nr,1)),axis=1)
        ligand = np.delete(ligand,0,axis=1)
        header_ligand.pop(0)
    return np.delete(Cross_ligand,0,axis=1), H_li
    
def Combine_Ligand_Protein(ligand, header_ligand, protein, header_protein):
    import numpy as np
    R, Cl = np.shape(ligand)
    R, Cp = np.shape(protein)
    Cross_lp = np.zeros((R,1))
    H_lp = []
    for j in range(Cl):
        for jj in range(Cp):
            H_lp.append(header_ligand[j]+'*'+header_protein[jj])
            cross = np.multiply(ligand[:,j], protein[:,jj])
            Cross_lp = np.append(Cross_lp, np.reshape(cross,(ligand.shape[0],1)),axis=1)
    return np.delete(Cross_lp,0,axis=1), H_lp
            
def Normalized(var):    #mean centering and scaling
    from sklearn import preprocessing
    return preprocessing.scale(var)
    
    
def Inputs_2(First, Header_First, Second, Header_Second):
    import numpy as np
    r, c1 = np.shape(First)
    r, c2 = np.shape(Second)
    
    norm_1 = np.sqrt(c1)*np.ones([r,c1])
    First = np.divide(First,norm_1)
    norm_2 = np.sqrt(c2)*np.ones([r,c2])
    Second = np.divide(Second,norm_2)

    print 'Sum variance of b1: %0.1f, b2: %0.1f' %(sum(np.var(First, axis=0)), 
                                                   sum(np.var(Second, axis=0)))
            
    Descriptors = np.zeros((r, c1+c2))
    Descriptors[:,0:c1] = First
    Descriptors[:,c1:c1+c2] = Second

    HEADER = []
    for i in range(c1):
        HEADER.append(Header_First[i])
    for i in range(c2):
        HEADER.append(Header_Second[i])
           
    return Descriptors, HEADER
    
def Inputs_3(First, Header_First, Second, Header_Second, Third, Header_Third):
    import numpy as np
    r, c1 = np.shape(First)
    r, c2 = np.shape(Second)
    r, c3 = np.shape(Third)
    
    norm_1 = np.sqrt(c1)*np.ones([r,c1])
    First = np.divide(First,norm_1)
    norm_2 = np.sqrt(c2)*np.ones([r,c2])
    Second = np.divide(Second,norm_2)
    norm_3 = np.sqrt(c3)*np.ones([r,c3])
    Third = np.divide(Third,norm_3)
    
    print 'Sum variance of b1: %0.1f, b2: %0.1f, b3: %0.1f' %(sum(np.var(First, axis=0)),
                                                              sum(np.var(Second, axis=0)), 
                                                              sum(np.var(Third, axis=0)))
        
    Descriptors = np.zeros((r, c1+c2+c3))
    Descriptors[:,0:c1] = First
    Descriptors[:,c1:c1+c2] = Second
    Descriptors[:,c1+c2:c1+c2+c3] = Third

    HEADER = []
    for i in range(c1):
        HEADER.append(Header_First[i])
    for i in range(c2):
        HEADER.append(Header_Second[i])
    for i in range(c3):
        HEADER.append(Header_Third[i])
            
    return Descriptors, HEADER
        
def Inputs_4(First, Header_First, Second, Header_Second,Third, Header_Third, Fouth, Header_Fouth):
    import numpy as np
    r, c1 = np.shape(First)
    r, c2 = np.shape(Second)
    r, c3 = np.shape(Third)
    r, c4 = np.shape(Fouth)
    
    norm_1 = np.sqrt(c1)*np.ones([r,c1])
    First = np.divide(First,norm_1)
    norm_2 = np.sqrt(c2)*np.ones([r,c2])
    Second = np.divide(Second,norm_2)
    norm_3 = np.sqrt(c3)*np.ones([r,c3])
    Third = np.divide(Third,norm_3)
    norm_4 = np.sqrt(c4)*np.ones([r,c4])
    Fouth = np.divide(Fouth,norm_4)

    print 'Sum variance of b1: %0.1f, b2: %0.1f, b3: %0.1f, b4: %0.1f' %(sum(np.var(First, axis=0)),
                                                                            sum(np.var(Second, axis=0)), 
                                                                            sum(np.var(Third, axis=0)),
                                                                            sum(np.var(Fouth, axis=0)))
    Descriptors = np.zeros((r, c1+c2+c3+c4))
    Descriptors[:,0:c1] = First
    Descriptors[:,c1:c1+c2] = Second
    Descriptors[:,c1+c2:c1+c2+c3] = Third
    Descriptors[:,c1+c2+c3:c1+c2+c3+c4]= Fouth
        
    HEADER=[]
    for i in range(c1):
        HEADER.append(Header_First[i])
    for i in range(c2):
        HEADER.append(Header_Second[i])
    for i in range(c3):
        HEADER.append(Header_Third[i])
    for i in range(c4):
        HEADER.append(Header_Fouth[i])
            
    return Descriptors, HEADER
        
def Inputs_5(First, Header_First, Second, Header_Second,Third, Header_Third, Fouth, Header_Fouth, Fifth, Header_Fifth):
        
    import numpy as np
    r, c1 = np.shape(First)
    r, c2 = np.shape(Second)
    r, c3 = np.shape(Third)
    r, c4 = np.shape(Fouth)
    r, c5 = np.shape(Fifth)
    norm_1 = np.sqrt(c1)*np.ones([r,c1])
    norm_2 = np.sqrt(c2)*np.ones([r,c2])
    norm_3 = np.sqrt(c3)*np.ones([r,c3])
    norm_4 = np.sqrt(c4)*np.ones([r,c4])
    norm_5 = np.sqrt(c5)*np.ones([r,c5])
    First = np.divide(First,norm_1)
    Second = np.divide(Second,norm_2)
    Third = np.divide(Third,norm_3)
    Fouth = np.divide(Fouth,norm_4)
    Fifth = np.divide(Fifth,norm_5)

    print 'Sum variance of b1: %0.1f, b2: %0.1f, b3: %0.1f, b4: %0.1f, b5: %0.1f' %(sum(np.var(First, axis=0)),
                                                                                    sum(np.var(Second, axis=0)), 
                                                                                    sum(np.var(Third, axis=0)),
                                                                                    sum(np.var(Fouth, axis=0)),
                                                                                    sum(np.var(Fifth, axis=0)))
        
    Descriptors = np.zeros((r, c1+c2+c3+c4+c5))
    Descriptors[:,0:c1] = First
    Descriptors[:,c1:c1+c2] = Second
    Descriptors[:,c1+c2:c1+c2+c3] = Third
    Descriptors[:,c1+c2+c3:c1+c2+c3+c4]= Fouth
    Descriptors[:,c1+c2+c3+c4:c1+c2+c3+c4+c5] = Fifth
        
    HEADER=[]
    for i in range(c1):
        HEADER.append(Header_First[i])
    for i in range(c2):
        HEADER.append(Header_Second[i])
    for i in range(c3):
        HEADER.append(Header_Third[i])
    for i in range(c4):
        HEADER.append(Header_Fouth[i])
    for i in range(c5):
        HEADER.append(Header_Fifth[i])
        
    return Descriptors, HEADER