#/usr/bin/python

import argrparse
import re
import pandas as pd
import multiprocessing as mp
import Preprocess_Spectrast_misc as misc

parser = argparse.ArgumentParser(description='Options for creating database')
parser.add_argument('database', type=str, help='database to create')
parser.add_argument('spread_out_value', type=int, help='value to spread out peaks')
parser.add_argument('Daltons', type=int, help='Value of bining peaks in Daltons')
parser.add_argument('n_peaks',type=int, help='Minimun number of peaks for each spectra in the new database')
parser.add_argument('intensity', type=int, help="Minimum intensity of spectra")

args = parser.parse_args()

def readSpectraSTDB(database):
    with open(database, 'rb') as origin_file:
        Names = []
        Lib_id = []
        Status = []
        Full_name = []
        Comments = []
        NumPeaks = []
        PrecursorMZ = []
        processed_file =  []
        MW = []
        for line in origin_file:
            if not line.startswith('###'):
                processed_file.append(line)
            if line.startswith('Name'):
                Names.append(line)
            if line.startswith('LibID'):
                Lib_id.append(line)
            if line.startswith('PrecursorMZ'):
                PrecursorMZ.append(line)
            if line.startswith('Status'):
                Status.append(line)
            if line.startswith('FullName'):
                Full_name.append(line)
            if line.startswith('Comment'):
                Comments.append(line)
            if line.startswith('NumPeaks'):
                NumPeaks.append(line)
            if line.startswith('MW'):
                MW.append(line)
        # pattern = r'NumPeaks: [0-9]+\n'
        # cPattern = re.compile(pattern)
    processed_file = '###'.join(processed_file)
    processed_file = processed_file.replace('NumPeaks:', 'NumPeaks:@')
    b = processed_file.split('\n###Name:')
    # c = b.split('@')
    c = map(lambda x: x.split('@')[-1], b)
    pattern2 = r' [0-9]+\n###'
    cPattern = re.compile(pattern2)
    d = map(lambda x: cPattern.split(x)[-1],c)
    pattern3 = r'(?<=[ ])([0-9A-Z\/]+)'
    cPattern = re.compile(pattern3)
    names_clean = map(lambda x: cPattern.findall(x), Names)
    lib_id_clean = map(lambda x: cPattern.findall(x), Lib_id)
    pattern4 = r'(?<=[ ])([0-9.\/]+)'
    cPattern = re.compile(pattern4)
    Precursor_clean = map(lambda x: cPattern.findall(x), PrecursorMZ)
    MW_clean = map(lambda x: cPattern.findall(x), MW)
    pattern5 = r'(?<=[ ])([A-Za-z\/]+)'
    cPattern = re.compile(pattern5)
    Status_clean = map(lambda x: cPattern.findall(x), Status)
    pattern6 = r'(?<=[ ])([A-Za-z.[0-9\/]+)'
    cPattern = re.compile(pattern6)
    Fullname_clean = map(lambda x: cPattern.findall(x), Full_name)
    pattern7 = r'(?<=[ ])([0-9\/]+)'
    cPattern = re.compile(pattern7)
    Numpeaks_clean = map(lambda x: cPattern.findall(x), NumPeaks)
    df_panda = pd.DataFrame( {'Name': names_clean, 'LibID' : lib_id_clean, 'MW' : MW_clean, 'PrecursorMZ': Precursor_clean, 'Status' :  Status_clean, 'FullName' : Fullname_clean, 'Comment' : Comments, 'NumPeaks':  Numpeaks_clean, 'Peaks': d})
    return df_panda

db_specta = readSpectraSTDB(database)

def filteringSpectraByNumberOfPeaks(db_spectra , n_of_peaks):
    db_spectra_filt = db_spectra.loc[db_spectra['NumPeaks'] > n_of_peaks]
    return db_spectra_filt

db_spectra_f = filteringSpectraByNumberOfPeaks(db_specta, n_peaks)

def filteringSpectraByIntensity(db_spectra , intensity):
    db_spectra_filt_2 = db_spectra.loc[db_spectra['PrecursorMZ'] > intensity]
    return db_spectra_filt_2

db_spectra_filt_2 = filteringSpectraByIntensity(db_spectra_f, intensity)

list_all_peaks = db_spectra_filt_2.loc[:,"Peaks"]

pandas=[]
for i in list_all_peaks:
    tmp = i
    tmp_split = tmp.split('\n###')
    if tmp_split[-1] == '':
        tmp_split = tmp_split[:-1]
    M = map(lambda x: x.split('\t')[0], tmp_split)
    Z = map(lambda x: x.split('\t')[1], tmp_split)
    df_panda = pd.DataFrame( {'M': M , 'Z' : Z})
    pandas.append(df_panda)

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
        print '====================='
        #iteration = range(0,len(panda)-2)
        # for i, row in panda.iterrows():
        for i in xrange(len(panda)-2):
            if SKIP_LINE:
                SKIP_LINE = False
                next
            panda = panda.sort_values(by=['M'])
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

pandas_binned = []
for panda in pandas:
    binned = bining(panda, Dalton)
    pandas_binned.append(binned)

#spread out
pandas_spread = []
for i in pandas_binned:
    tmp = i
    tmp['percentage'] = [0] * len(tmp)
    tmp['peak_less'] = [0] * len(tmp)
    tmp['peaks_final'] = [0] * len(tmp)  #ojo aqui
    tmp['peaks_final'][1]<-
    for row in tmp.iterrows():
        tmp['percentage'][row] = (pd.to_numeric(tmp.loc[:'Z'][row]) * spread_value)/2
        tmp['peak_less'][row] = (pd.to_numeric(tmp.loc[:'Z'][row]) - (pd.to_numeric(tmp.loc[:'percentage'][row]))
