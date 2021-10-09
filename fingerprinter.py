# Generate Fingerprint SMILES

from rdkit import Chem
from rdkit.Chem import AllChem

def FingerprintToSmiles(m, s):
    fp_sm = []
    bi = {}
    fp = AllChem.GetMorganFingerprint(m,s, bitInfo=bi)
    # print('FPSM : ', bi)
    for f in bi:
        # print('K:', f,' V:', bi[f])
        a = bi[f][0][0]
        r = bi[f][0][1]
        # print(f, a, r)
        if r > 0:
            env = Chem.FindAtomEnvironmentOfRadiusN(m,r,a)
            amap={}
            submol=Chem.PathToSubmol(m,env,atomMap=amap)
            sm =Chem.MolToSmiles(submol)
        else:
            am = m.GetAtomWithIdx(a)
            sm = am.GetSymbol()
            if am.GetIsAromatic():
                sm = sm.lower()
        fp_sm.append((f,sm))
        # print(f,' - ',sm)
        # print(f,' - ',len(v),' - ',v,' - ',a,' - ',r,' - ',sm)
    return fp_sm