#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.append('../pkg_mod')

import db_connector as db
import chemistry_psql as cp
import importlib
importlib.reload(cp)
import psycopg2
#import pandas as pd
from rdkit import Chem
#from rdkit.Chem import SaltRemover as rdsr


############ Open Database Connection and Pull Structures ############


db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
cur = db_chem.cursor()
sql = 'select * from structure.structures where chiral is NULL limit 1000;'
#print(sql)
cur.execute(sql)
molsql = cur.fetchall()
#print(molsql)
mollist = []
molrev = []

############ Create List of Multifragment Structures to Review ############

print('===== Structures to Review ===== ')
print()
for i in molsql:
    j = i[1].split('.')
    #print(j)
    if len(j) == 1:
        mollist.append(i)
    else:
        print(j)
        molrev.append(i)
        
print()
print(len(molrev), ' Structures need review - Full line entry is following')
print()
print(molrev)
        

############ Generate Chirality Information in Database Structures Table ############


print()
print('===== Structures that need chirality treatment =====')
print()
print(len(mollist), ' Structures to assign chirality')
print()


for i in mollist:
    m = Chem.MolFromSmiles(i[1])
    mc = Chem.FindMolChiralCenters(m,force=True,includeUnassigned=True)
    mct = Chem.FindMolChiralCenters(m,force=True,includeUnassigned=False)
    mcl = ''.join([str(e) for e in mc])
    #print(i[0],i[1],mc)
    ### No chiral centers ###
    if len(mc) == 0:
        #print()
        sql = 'UPDATE structure.structures SET chiral = 0, chiral_assigned = 0 WHERE structure_id = '+str(i[0])+';'
        print(sql)
        cur.execute(sql)
        db_chem.commit()
    else:
    ### No defined chiral centers ###
        mct = Chem.FindMolChiralCenters(m,force=True,includeUnassigned=False)
        mcl = ''.join([str(e) for e in mc])
        print()
        print(i[0],len(mc), i[1],mc)
        print(i[0],len(mct),i[1],mct)
        mcl = mcl.replace(chr(39),'')
        mcl = mcl.replace(' ','')
        print(mcl)
        if len(mct) == 0:
            sql = "UPDATE structure.structures SET chiral = "+str(len(mc))+", chiral_assigned = 0, chiral_list = '"+mcl+"' WHERE structure_id = "+str(i[0])+";"
            print(sql)
            cur.execute(sql)
            db_chem.commit()
        else:
    ### Defined chiral centers ###
    ### Flatten structures ### 
            mf = m
            Chem.RemoveStereochemistry(mf)
            mfc = Chem.FindMolChiralCenters(mf,force=True,includeUnassigned=True)
            mfct = Chem.FindMolChiralCenters(mf,force=True,includeUnassigned=False)
            mfs = Chem.MolToSmiles(mf)
            mfcl = ''.join([str(e) for e in mfc])
            print(len(mfc), mfs, mfc)
            print(len(mfct), mfs, mfct)
            mfcl = mfcl.replace(chr(39),'')
            mfcl = mfcl.replace(' ','')
            print(mfcl)
     ### Does flattened structures exist ? ###       
            #sql = 'select * from structure.structures where molecule::structure.mol @= structure.mol_from_smiles(\''+mfs+'\'::cstring);'
            sql = 'select * from structure.structures where substruct(molecule, \''+mfs+'\'::mol) and mol_inchi(molecule)::text = mol_inchi(\''+mfs+'\'::mol)::text;'
            #print(sql)
            cur.execute(sql)
            molflat = cur.fetchall()
            print(molflat)
            if len(molflat) == 1:
     ### One flattened structure exists - update flattened structure record and original record ###
                #print(molflat)
                sid = molflat[0][0]
                sql = "UPDATE structure.structures SET chiral = "+str(len(mc))+", chiral_assigned = 0, chiral_list = '"+mfcl+"' WHERE structure_id = "+str(sid)+";"
                print(sql)
                cur.execute(sql)
                sql = "UPDATE structure.structures SET chiral = "+str(len(mc))+", chiral_assigned = "+str(len(mct))+", chiral_list = '"+mcl+"' WHERE structure_id = "+str(i[0])+";"
                print(sql)
                cur.execute(sql)
                db_chem.commit()
            elif len(molflat) == 0:
    ### No flattened structure - insert flattened structure and update original record ###
                smi = (mfs,)
                print(smi)
                cp.SmilesInsert(smi)
                sql = 'select max(structure_id) from structure.structures'
                cur.execute(sql)
                smax = cur.fetchone()
                sm = smax[0]
                print(smax, sm)
                sql = "UPDATE structure.structures SET chiral = "+str(len(mc))+", chiral_assigned = 0, chiral_list = '"+mfcl+"' WHERE structure_id = "+str(sm)+";"
                print(sql)
                cur.execute(sql)
                sql = "UPDATE structure.structures SET chiral = "+str(len(mc))+", chiral_assigned = "+str(len(mct))+", chiral_list = '"+mcl+"' WHERE structure_id = "+str(i[0])+";"
                print(sql)
                cur.execute(sql)
                db_chem.commit()
            else:
     ### More than one flattened structure --- Alert ###
                print('More than one match!!!!!')
        print()





