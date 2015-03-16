# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 10:08:00 2014

@author: Fujitsu
"""
def UserDefined(Rawfile,Indicator,Ligand_index,Protein_index, Model_index,SpiltCriteria,
               CV_Method,FeatureSelectionMode, Iteration, NumPermute, SpiltMethod):
    import os
    user = {}
    user['Root'] = os.getcwd()
    user['Rawfile'] = Rawfile
    user['Indicator'] = Indicator
    user['Ligand_index'] = Ligand_index[1:-1].split(',')
    user['Protein_index'] = Protein_index[1:-1].split(',')
    user['Model_index'] = Model_index[1:-1].split(',')
    user['Spiltcriteria'] = SpiltCriteria
    user['CV_Method'] = CV_Method
    user['SelectionMode'] = FeatureSelectionMode
    user['Iteration'] = Iteration
    user['NumPermute'] = NumPermute
    user['SpiltMethod'] = SpiltMethod
    
    from time import gmtime, strftime
    user['Date Started'] = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    return user
    
def AnalysisInputfile(user):
    Root = user['Root']
    Rawfile = user['Rawfile']  
    Proteingroup = user['Protein_index']  
    Ligandgroup = user['Ligand_index'] 
    
    import csv
    import numpy as np
    fileName = Root+'/'+Rawfile+'.csv'
    with open(fileName,'rb') as csvfile:
        dialect = csv.Sniffer().has_header(csvfile.read())
        csvfile.seek(0)
        reader = csv.reader(csvfile, dialect)
        h = next(reader)
        data = []
        for row in reader:
            data.append(row)
        data_array = np.array(data)

    Yname = h.pop(-1)
    Y_array = np.append(np.reshape(np.array(Yname),(1,1)),data_array[:,-1])
    
    if len(np.unique(data_array[:,-1])) > 3:  #regression
        user['Datatype'] = 'Regression'
    elif len(np.unique(data_array[:,-1])) == 3:
        user['Datatype'] = 'Classification 3 classes'
    elif len(np.unique(data_array[:,-1])) == 2:
        user['Datatype'] = 'Classification 2 classes'
        
    psmiles = [ind for ind, val in enumerate(h) if val[0:5] == 'Smile' or val[0:5] == 'smile']
    psequence = [ind for ind, val in enumerate(h) if val[0:8] == 'Sequence' or val[0:8] == 'sequence']
    
    if len(psmiles) == 0 and len(psequence) == 0:
        print 'All Descriptors were prepared by user'
        ise = [ind for ind,val in enumerate(h) if val == '']
        
        if len(ise) != 0: 
            hi = np.reshape(np.array(h), (1,len(h)))
            hf = np.reshape(np.array(Yname), (1,1))
            h = np.append(hi,hf, axis=1)
        
            data_array = np.append(h,data_array,axis=0)
            Array_ligand = data_array[:,:ise[0]]
            Array_Pro = data_array[:,ise[0]+1:-1]
        else:
            if user['Ligand_index'] == []:
                Array_ligand = []
                Array_Pro = data_array[:,:-1]
            elif user['Protein_index'] == []:
                Array_ligand = data_array[:,:-1]
                Array_Pro = []
       
    elif len(psmiles) == 1 and len(psequence) == 0:
        print 'Ligand descriptors will be generated'
        import Descriptors_Extraction as DE
        data = data_array[:,psmiles[0]]
        Array_ligand = DE.Ligand_gen(data, Ligandgroup)
        px = [ind for ind,val in enumerate(h) if ind!=psmiles[0]]
        hx = np.array(h)[px]
        Array_Pro = np.append(np.reshape(hx,(1,len(hx))),data_array[:,px],axis=0)
        
    elif len(psmiles) == 0 and len(psequence) == 1:
        print 'Protein descriptors will be generated'
        import Descriptors_Extraction as DE
        data = data_array[:,psequence[0]]
        Array_Pro = DE.Protein_gen(data,Proteingroup)
        px = [ind for ind,val in enumerate(h) if ind!=psequence[0]]
        hx = np.array(h)[px]
        Array_ligand = np.append(np.reshape(hx,(1,len(hx))),data_array[:,px],axis=0)
        
    elif len(psmiles) == 1 and len(psequence) == 1:
        print 'Ligand & Protein descriptors will be generated'
        import Descriptors_Extraction as DE
        data1 = data_array[:,psmiles[0]]
        data2 = data_array[:,psequence[0]]
        Array_ligand = DE.Ligand_gen(data1,Ligandgroup)
        Array_Pro = DE.Ligand_gen(data2,Proteingroup)
        
    elif len(psmiles) == 2 and len(psequence) == 0:
        print 'Two different Ligand descriptors will be generated'
        import Descriptors_Extraction as DE
        data1 = data_array[:,psmiles[0]]
        data2 = data_array[:,psmiles[1]]
        Array_ligand = DE.Ligand_gen(data1,Ligandgroup)
        Array_Pro = DE.Ligand_gen(data2,Proteingroup)
        
    elif len(psmiles) == 0 and len(psequence) == 2:
        print 'Two different Protein descriptors will be generated'
        import Descriptors_Extraction as DE
        data1 = data_array[:,psequence[0]]
        data2 = data_array[:,psequence[1]]
        Array_ligand = DE.Protein_gen(data1,Ligandgroup)
        Array_Pro = DE.Protein_gen(data2,Proteingroup)
        
################## Comnbine All array for saving ##############
    emp = np.array([None for i in range(Array_Pro.shape[0])])
    emp = np.reshape(emp, (emp.shape[0],1))
    
    Array = np.append(Array_ligand, emp, axis=1)
    Array = np.append(Array, Array_Pro, axis=1)
    Array = np.append(Array, np.reshape(Y_array,(len(Y_array),1)), axis=1)

    path = user['Root']
    raw = user['Rawfile']
    Indica = user['Indicator']
    
    import os
    try:
        os.makedirs(path+'/'+Indica)
    except OSError:
        pass
    
    with open(path+'/'+Indica+'/'+raw+'_complete'+'.csv', 'wb') as csvfile:
        spam = csv.writer(csvfile,delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL )
        for k in range(len(Array)):
            spam.writerow(Array[k])
            
    return Array_ligand, Array_Pro, Y_array, user
   
    
def Ligand_gen(data, Ligandgroup):
    import os
    import numpy as np
    os.chdir('C:\Users\Fujitsu\Anaconda\envs\Rdkit')
    from pydpi.pydrug import PyDrug
    drug=PyDrug()

    HL_list, D_list = [], []
        
    for i in range(len(data)):
        drug.ReadMolFromSmile(data[i])
        keys, values = [],[]
    
        for j in Ligandgroup:
            if j == '0':    #all descriptors   615
                res = drug.GetAllDescriptor()
            elif j == '1':    # constitution   30
                res = drug.GetConstitution()
            elif j == '2':    # topology       25
                res = drug.GetTopology()
            elif j == '3':    #connectivity    44
                res = drug.GetConnectivity()
            elif j == '4':    #E-state         237
                res = drug.GetEstate()
            elif j == '5':    #kappa            7
                res = drug.GetKappa()
            elif j == '6':    #Burden           64
                res = drug.GetBurden()
            elif j == '7':    #information      21
                res = drug.GetBasak()
            elif j == '8':    #Moreau-Boto      32
                res = drug.GetMoreauBroto()
            elif j == '9':    #Moran            32
                res = drug.GetMoran()
            elif j == '10':   #Geary            32
                res = drug.GetGeary()
            elif j == '11':   #charge           25
                res = drug.GetCharge()
            elif j == '12':   #property          6
                res = drug.GetMolProperty()
            elif j == '13':   #MOE-type          60
                res = drug.GetMOE()
            
            keys.extend(res.viewkeys())
            values.extend(res.viewvalues())
        
        if i == 0:
            HL_list = keys
            D_list.append(values)
        else:
            D_list.append(values)

    D_ligand = np.zeros((len(data),len(HL_list)), dtype=float)
    for k in range(len(data)):
        D_ligand[k,:] = D_list[k]   
    
    
    #Variance threshold       std > 0.01  
    import Descriptors_Selection as DesSe
    ind_var = DesSe.VarinceThreshold(D_ligand)
    D_ligand = D_ligand[:,ind_var]
    HL_list = np.array(HL_list)[ind_var]
        
#    #Intra pearson's correlation           p-value > 0.05
#    ind_corr = DesSe.Correlation(D_ligand, Y.astype(np.float))
#    D_ligand = D_ligand[:,ind_corr]
#    HL_list = np.array(HL_list)[ind_corr]
        
    H_ligand = np.reshape(HL_list,(1,len(HL_list)))
    Array_ligand = np.append(H_ligand, D_ligand, axis=0) 
    return Array_ligand

def Protein_gen(data, Proteingroup):
    import numpy as np
    from pydpi.pypro import PyPro
    protein = PyPro() 
    
    HP_list, D_list = [], []     
    for ii in range(len(data)):
        p = data[ii]
        protein.ReadProteinSequence(p)
        keys, values = [],[]
        for jj in Proteingroup:
            if jj == '0':    #All descriptors          2049
                res = protein.GetALL()
            elif jj == '1':    #amino acid composition   20
                res = protein.GetAAComp()
            elif jj == '2':    #dipeptide composition    400
                res = protein.GetDPComp()
            elif jj == '3':    #Tripeptide composition   8000
                res = protein.GetTPComp()
            elif jj == '4':    #Moreau-Broto autocorrelation  240
                res = protein.GetMoreauBrotoAuto()   
            elif jj == '5':    #Moran autocorrelation       240
                res = protein.GetMoranAuto()
            elif jj == '6':    #Geary autocorrelation       240
                res = protein.GetGearyAuto()
            elif jj == '7':    #composition,transition,distribution  21+21+105
                res = protein.GetCTD()
            elif jj == '8':    #conjoint triad features     343
                res = protein.GetTriad()
            elif jj == '9':    #sequence order coupling number  60
                res = protein.GetSOCN(30)
            elif jj == '10':   #quasi-sequence order descriptors   100
                res = protein.GetQSO()
            elif jj == '11':    #pseudo amino acid composition   50
                res = protein.GetPAAC(30)
                    
            keys.extend(res.viewkeys())
            values.extend(res.viewvalues())  
        if ii == 0:
            HP_list = keys
            D_list.append(values)
        else:
            D_list.append(values)
            
    D_Pro = np.zeros((len(D_list),len(HP_list)), dtype=float)
    for k in range(len(D_list)):
        D_Pro[k,:] = D_list[k]
        
    #Variance threshold       std > 0.01  
    import Descriptors_Selection as DesSe
    ind_var = DesSe.VarinceThreshold(D_Pro)
    D_Pro = D_Pro[:,ind_var]
    HP_list = np.array(HP_list)[ind_var]

    H_Pro = np.reshape(HP_list,(1,len(HP_list)))
    Array_Pro = np.append(H_Pro, D_Pro, axis=0) 
    
    return Array_Pro

#def LigandProteinExtraction(user):
#    Rawfile = user['Rawfile'] 
#    IndicatorName = user['Indicator'] 
#    Proteingroup = user['Protein_index']  
#    Ligandgroup = user['Ligand_index'] 
#    
#    import os, csv, sys
#    import numpy as np
#    
#    path = os.path.dirname(os.path.abspath(sys.argv[0]))
#    fileName = path+'/'+ Rawfile +'.csv'
#    
#    with open(fileName,'rb') as csvfile:
#        dialect = csv.Sniffer().has_header(csvfile.read())
#        csvfile.seek(0)
#        reader = csv.reader(csvfile, dialect)
#        h = next(reader)
#        data = []
#        for row in reader:
#            data.append(row)
#        data_array = np.array(data)
#        
#    
#    ind = [ind for ind, val in enumerate(h) if val == 'SMILES' or val == 'Sequence']
#        
#    if len(ind) == 0:   
#        print 'Descriptors were prepared by user'
#        
#        Y = data_array[:,-1]
#        data = np.delete(data_array,0,axis=1)
#        data = np.delete(data,-1,axis=1)
#    
#        h.pop(0)
#        Yname = h.pop(-1)
#    
#        idx_pt = [inx for inx,val in enumerate(data[0,:]) if val.isdigit() == True]
#        Protein = np.transpose(np.transpose(data)[idx_pt])
#        Ligand = np.delete(data, idx_pt, axis=1) 
#        
#        h_pt = np.array(h)[idx_pt]
#        h_li = np.delete(np.array(h), idx_pt)
#        
#        Array_Pro = np.append(np.reshape(h_pt,(1,len(h_pt))),Protein,axis=0)
#        Array_ligand = np.append(np.reshape(h_li,(1,len(h_li))),Ligand,axis=0)
#        Y_array = np.append(np.reshape(np.array(Yname),(1,1)),Y)
#
#    else:
#        import sys
#        import Descriptors_Selection as DesSe
#        print '### No Descriptors. Automatic preparation is being processed ###'
#        
#        if len(Ligandgroup) == 0 or len(Proteingroup) == 0:
#            sys.exit('Index groups of Descriptor must be determined')
#        
#        Smile = data_array[:,1]
#        Sequence = data_array[:,3]
#        Y = data_array[:,-1]
#        
#        ############ Protein desciptors ################
#        import Descriptors_Extraction as DesEx
#        D_Pro, HP_list = DesEx.ProteinDescriptors(Sequence,Proteingroup)
#        
#        #Variance threshold       std > 0.01  
#        ind_var = DesSe.VarinceThreshold(D_Pro)
#        D_Pro = D_Pro[:,ind_var]
#        HP_list = np.array(HP_list)[ind_var]
#        
#        #Intra pearson's correlation           p-value > 0.05
#        ind_corr = DesSe.Correlation(D_Pro, Y.astype(np.float))
#        D_Pro = D_Pro[:,ind_corr]
#        HP_list = np.array(HP_list)[ind_corr]
#        
#        H_Pro = np.reshape(HP_list,(1,len(HP_list)))
#        Array_Pro = np.append(H_Pro, D_Pro, axis=0) 
#            
#        ############## Ligand descriptors ################
#        from pydpi.pydrug import PyDrug
#        drug=PyDrug()
#
#        HL_list, D_list = [], []
#        
#        for i in range(len(Smile)):
#            drug.ReadMolFromSmile(Smile[i])
#            keys, values = [],[]
#    
#            for j in Ligandgroup:
#                if j == '0':    #all descriptors   615
#                    res = drug.GetAllDescriptor()
#                elif j == '1':    # constitution   30
#                    res = drug.GetConstitution()
#                elif j == '2':    # topology       25
#                    res = drug.GetTopology()
#                elif j == '3':    #connectivity    44
#                    res = drug.GetConnectivity()
#                elif j == '4':    #E-state         237
#                    res = drug.GetEstate()
#                elif j == '5':    #kappa            7
#                    res = drug.GetKappa()
#                elif j == '6':    #Burden           64
#                    res = drug.GetBurden()
#                elif j == '7':    #information      21
#                    res = drug.GetBasak()
#                elif j == '8':    #Moreau-Boto      32
#                    res = drug.GetMoreauBroto()
#                elif j == '9':    #Moran            32
#                    res = drug.GetMoran()
#                elif j == '10':   #Geary            32
#                    res = drug.GetGeary()
#                elif j == '11':   #charge           25
#                    res = drug.GetCharge()
#                elif j == '12':   #property          6
#                    res = drug.GetMolProperty()
#                elif j == '13':   #MOE-type          60
#                    res = drug.GetMOE()
#            
#                keys.extend(res.viewkeys())
#                values.extend(res.viewvalues())
#        
#            if i == 0:
#                HL_list = keys
#                D_list.append(values)
#            else:
#                D_list.append(values)
#
#        D_ligand = np.zeros((len(Smile),len(HL_list)), dtype=float)
#        for k in range(len(Smile)):
#            D_ligand[k,:] = D_list[k]
#            
#        #Variance threshold       std > 0.01  
#        ind_var = DesSe.VarinceThreshold(D_ligand)
#        D_ligand = D_ligand[:,ind_var]
#        HL_list = np.array(HL_list)[ind_var]
#        
#        #Intra pearson's correlation           p-value > 0.05
#        ind_corr = DesSe.Correlation(D_ligand, Y.astype(np.float))
#        D_ligand = D_ligand[:,ind_corr]
#        HL_list = np.array(HL_list)[ind_corr]
#        
#        H_ligand = np.reshape(HL_list,(1,len(HL_list)))
#        Array_ligand = np.append(H_ligand, D_ligand, axis=0) 
#        
#         
#    ################## Comnbine All array for saving ##############
#    
#        Array = np.append(Array_Pro, Array_ligand, axis=1)
#        
#        Y_array = np.append(np.array(h[-1]), Y)
#        Y_array = np.reshape(Y_array, (len(Y_array),1))
#        Array = np.append(Array, Y_array, axis=1)
#
#        try:
#            os.makedirs(path+'/'+IndicatorName)
#        except OSError:
#            pass
#        
#        with open(path+'/'+IndicatorName+'/'+IndicatorName+'.csv', 'wb') as csvfile:
#            spam = csv.writer(csvfile,delimiter=',',quotechar='|', 
#                          quoting=csv.QUOTE_MINIMAL )
#            for k in range(len(Array)):
#                spam.writerow(Array[k])
#            
#    
#    return Array_ligand, Array_Pro, Y_array
#    
#def ProteinExtraction(user):
#    RawfileName = user['Rawfile']
#    IndicatorName = user['Indicator']
#    Proteingroup = user['Protein_index']
#    
#    import os, csv, sys
#    import numpy as np
#    import Descriptors_Selection as DesSe
#    
#    path = os.path.dirname(os.path.abspath(sys.argv[0]))
#    fileName = path+'/'+ RawfileName +'.csv'
#    
#    with open(fileName,'rb') as csvfile:
#        dialect = csv.Sniffer().has_header(csvfile.read())
#        csvfile.seek(0)
#        reader = csv.reader(csvfile, dialect)
#        h = next(reader)
#        data = []
#        for row in reader:
#            data.append(row)
#        data_array = np.array(data)
#    Y = data_array[:,-1]
#        
#    ind_P1 = [ind for ind,val in enumerate(h) if val == 'Protein1_Sequence']
#    ind_P2 = [ind for ind,val in enumerate(h) if val == 'Protein2_Sequence']
#    
#    data_P1 = data_array[:,ind_P1]
#    data_P2 = data_array[:,ind_P2]
#    
#    
#    S1,S2 = [],[]    
#    for i in range(len(data_P1)):
#        A, AA = data_P1[i][0], data_P2[i][0]
#        if A and AA == []:
#            pass
#        S1.append(A)
#        S2.append(AA)
#
#    import Descriptors_Extraction as DesEx
##    ########## Performing Protein Block 1 #################    
#    PS1, H1 = DesEx.ProteinDescriptors(S1,Proteingroup)
#    #Variance threshold       std > 0.01  
#    ind_var = DesSe.VarinceThreshold(PS1)
#    PS1 = PS1[:,ind_var]
#    H1 = np.array(H1)[ind_var]
#            
##        #Intra pearson's correlation           p-value > 0.05
###        ind_corr = DesSe.Correlation(D_Pro, Y.astype(np.float))
###        D_Pro = D_Pro[:,ind_corr]
###        HP_list = np.array(HP_list)[ind_corr]
##        
#    H_Pro = np.reshape(H1,(1,len(H1)))
#    Array_Pro = np.append(H_Pro, PS1, axis=0) 
#    Array_ligand = Array_Pro 
#               
#    ########## Performing Protein Block 2 #################    
#    PS2, H2 = DesEx.ProteinDescriptors(S2,Proteingroup)
#    #Variance threshold       std > 0.01  
#    ind_var = DesSe.VarinceThreshold(PS2)
#    PS2 = PS2[:,ind_var]
#    H2 = np.array(H2)[ind_var]
#            
##        #Intra pearson's correlation           p-value > 0.05
###        ind_corr = DesSe.Correlation(D_Pro, Y.astype(np.float))
###        D_Pro = D_Pro[:,ind_corr]
###        HP_list = np.array(HP_list)[ind_corr]
##        
#    H_Pro = np.reshape(H2,(1,len(H2)))
#    Array_Pro = np.append(H_Pro, PS2, axis=0) 
#    
#    ####################### Feature selection using VIP ################
#    import Descriptors_Selection as DS
#    
#    X1, Y1, H1 = DS.VIP_origin(PS1, Y, H_Pro)
#    H1 = np.reshape(H1,(1,len(H1)))
#    Array_ligand = np.append(H1, X1, axis=0)
#    
#    X2, Y2, H2 = DS.VIP_origin(PS2, Y, H_Pro)
#    H2 = np.reshape(H2,(1,len(H2)))
#    Array_Pro = np.append(H2, X2, axis=0)
#    
#    Y_array = np.append(np.array(h[-1]), Y2)
#    Y_array = np.reshape(Y_array, (len(Y_array),1))
#    
#
#        ################## Comnbine All array for saving ##############
#    
##    Array = np.append(Array_Pro, Array_ligand, axis=1)
##    Y_array = np.append(np.array(h[-1]), Y)
##    Y_array = np.reshape(Y_array, (len(Y_array),1))
##    Array = np.append(Array, Y_array, axis=1)
############################################################
##    try:            
##        os.makedirs(path+'/'+IndicatorName)        
##    except OSError:
##        pass
##        
##    with open(path+'/'+IndicatorName+'/'+IndicatorName+'.csv', 'wb') as csvfile:
##        spam = csv.writer(csvfile,delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL )
##        for k in range(len(Array)):
##            spam.writerow(Array[k])
#            
#    return Array_ligand, Array_Pro, np.squeeze(Y_array)
#            