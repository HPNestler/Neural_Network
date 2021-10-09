import psycopg2

# expects a tuple or list of tuples

def SmilesInsert(smi_in):
    
    db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
    cur = db_chem.cursor()
    sql = 'INSERT INTO public."Structures"("SMILES") VALUES(%s);'
    #print(sql)
    #smi_in = smi_t
    #print(smi_in)

    for i in smi_in:   
        smi = (i,)
        # print(smi)
        try:
            cur.execute(sql, smi)
            db_chem.commit()
            print('Inserted:', i)
        except:
            print("Structure Exist: ", i)
            cur.close()
            db_chem.close()
            db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
            cur = db_chem.cursor()
            sql = 'INSERT INTO public."Structures"("SMILES") VALUES(%s);'
            
# insert new Morgan Fingerprints - expect fpn as dict, smi as string
            
def FPInsertMorgan(fpn, smi):
           
    db_chem = psycopg2.connect(host = "localhost", dbname="Chemistry", user="postgres", password="postgres")
    cur = db_chem.cursor()
    sql = 'INSERT INTO public."Structures"("SMILES") VALUES(%s);'
    #print(sql)