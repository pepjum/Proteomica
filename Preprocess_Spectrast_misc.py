import re
import pandas as pd
import multiprocessing as mp


def bining(panda,Dalton):
    new = 0
    old = len(panda)
    pattern = r'^[0-9]+'
    cPattern = re.compile(pattern)
    panda.loc[:, 'M'] = map(lambda x: cPattern.findall(x)[0], panda.loc[:,'M'])
    panda = pd.DataFrame( {'M' : pd.to_numeric(panda.loc[:,'M']), "Z" : pd.to_numeric(panda.loc[:,"Z"])})
    while old != new:
        tmp = pd.DataFrame( {'M' : [], 'Z' : []})
        SKIP_LINE = False
        iteration = range(0,len(panda)-2)
        for i in iteration:
            if SKIP_LINE:
                SKIP_LINE = False
                next
            panda.sort_values(by=['M'])
            if pd.to_numeric(panda.loc[:,'M'][i+1]) - pd.to_numeric(panda.loc[:,'M'][i]) <= Dalton:
                maxM = max(pd.to_numeric(panda.loc[:,'M'][i+1]), pd.to_numeric(panda.loc[:,'M'][i]))
                maxZ = max(pd.to_numeric(panda.loc[:,'Z'][i+1]), pd.to_numeric(panda.loc[:,'Z'][i]))
                tmp = tmp.append({'M': maxM, 'Z': maxZ}, ignore_index=True)
                SKIP_LINE = True
            else:
                tmp = tmp.append({ 'M': pd.to_numeric(panda.loc[:,'M'][i]), 'Z': pd.to_numeric(panda.loc[:,'Z'][i])}, ignore_index=True)
                SKIP_LINE = False
        tmp = tmp.append({ 'M': pd.to_numeric(panda.loc[:,'M'][i+1]), 'Z': pd.to_numeric(panda.loc[:,'Z'][i+1])}, ignore_index=True)
        new = len(tmp)
        old = len(panda)
        panda = tmp
    if pd.to_numeric(panda.loc[:,'M'][len(panda)-1]) - pd.to_numeric(panda.loc[:,'M'][len(panda)-2]) <= Dalton:
        maxM = max(pd.to_numeric(panda.loc[:,'M'][len(panda)-1]), pd.to_numeric(panda.loc[:,'M'][len(panda)-2]))
        maxZ = max(pd.to_numeric(panda.loc[:,'Z'][len(panda)-1]), pd.to_numeric(panda.loc[:,'Z'][len(panda)-2]))
        panda2 = panda
        panda2 = panda2.append({'M': maxM, 'Z': maxZ}, ignore_index=True)
    else:
        panda2 = panda
    panda2 = panda2.drop_duplicates()
    return panda2
