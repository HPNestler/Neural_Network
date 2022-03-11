import psycopg2

# expects a tuple or list of tuples

def SmilesInsert(smi_in):
    
    #db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
    db_chem = psycopg2.connect(host = "192.168.86.31", dbname="Chemistry", user="postgres", password="postgres")
    cur = db_chem.cursor()
    sql = 'INSERT INTO public.structures (smiles) VALUES(%s);'
    #print(sql)
    #smi_in = smi_t
    #print(smi_in)

    for i in smi_in:   
        smi = (i,)
        #print(smi)
        try:
            cur.execute(sql, smi)
            db_chem.commit()
            print('Inserted:', i)
        except Exception as error:
            print("Structure Exist: ", i) #, error)
            cur.close()
            db_chem.close()
            #db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
            db_chem = psycopg2.connect(host = "192.168.86.31", dbname="Chemistry", user="postgres", password="postgres")
            cur = db_chem.cursor() 
            sql = 'INSERT INTO public.structures(smiles) VALUES(%s);'
            
# insert new Morgan Fingerprints - expect fpn as dict, smi as string
            
def FPInsertMorgan(fpn, smi):
           
    #db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
    db_chem = psycopg2.connect(host = "192.168.86.31", dbname="Chemistry", user="postgres", password="postgres")
    cur = db_chem.cursor()
    sql = 'INSERT INTO public.structures (smiles) VALUES(%s);'
    #print(sql)
    
# expects a tuple or list of tuples

def StructPropInsert(sp):
    
    #db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
    db_chem = psycopg2.connect(host = "192.168.86.31", dbname="Chemistry", user="postgres", password="postgres")
    cur = db_chem.cursor()
    sql = 'INSERT INTO public.structure_properties (property_type_id, structure_id, property_value) VALUES(%s, %s, %s);'
    #print(sql)
    #smi_in = smi_t
    #print(smi_in)

    print(sp)
    try:
        cur.execute(sql, (sp[0], sp[1], sp[2]))
        db_chem.commit()
        print('Inserted:', sp)
    except Exception as error:
        print("Property Exist: ", sp) #, error)
        cur.close()
        db_chem.close()
        #db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
        db_chem = psycopg2.connect(host = "192.168.86.31", dbname="Chemistry", user="postgres", password="postgres")
        cur = db_chem.cursor() 
        sql = 'INSERT INTO public.structure_properties (property_type_id, structure_id, property_value) VALUES(%s, %s, %s);'

# expects a tuple or list of tuples

def StructCalcPropInsert(sp):
    
    #db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
    db_chem = psycopg2.connect(host = "192.168.86.31", dbname="Chemistry", user="postgres", password="postgres")
    cur = db_chem.cursor()
    sql = 'INSERT INTO public.structure_prop_calc (property_type_id, structure_id, property_value) VALUES(%s, %s, %s);'
    #print(sql)
    #smi_in = smi_t
    #print(smi_in)

    print(sp)
    try:
        cur.execute(sql, (sp[0], sp[1], sp[2]))
        db_chem.commit()
        print('Inserted:', sp)
    except Exception as error:
        print("Property Exist: ", sp) #, error)
        cur.close()
        db_chem.close()
        #db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
        db_chem = psycopg2.connect(host = "192.168.86.31", dbname="Chemistry", user="postgres", password="postgres")
        cur = db_chem.cursor() 
        sql = 'INSERT INTO public.structure_prop_calc (property_type_id, structure_id, property_value) VALUES(%s, %s, %s);'

# find structure id for SMIlES

def SmilesID(smi):

    db_chem = psycopg2.connect(host = "192.168.86.31", dbname="Chemistry", user="postgres", password="postgres")
    cur = db_chem.cursor()
    sql = 'select structure_id, smiles from public.structures where molecule =\' '+smi+' \'::mol;'
    #print(sql)
    cur.execute(sql)
    mollist = cur.fetchall()
    return mollist[0][0]



 