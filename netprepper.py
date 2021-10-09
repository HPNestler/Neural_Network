import pandas as pd
from scipy import stats as scistat

def shannon(x, n):
    pv = 0
    ps = n - len(x)
    xf = x  + [pv] * ps
    pd_series = pd.Series(xf)
    # print(pd_series)
    counts = pd_series.value_counts()
    # print(counts)
    entropy = scistat.entropy(counts)
    return xf, entropy
    
print('Net_Prepper_Imported')

