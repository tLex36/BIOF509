#import plotly.graph_objects as go
import pandas as pn
import numpy as np
from sklearn.cluster import KMeans
import json

data = pn.ExcelFile(r'HA NLGN3 developmental mass spec.xlsx')
p6xl = pn.read_excel(data, sheet_name='P6<0.05',usecols="B,C,D,E,F")
p12xl = pn.read_excel(data, sheet_name='P12<0.05',usecols="B,C,D,E,F")
p18xl = pn.read_excel(data, sheet_name='P18<0.05',usecols="B,C,D,E,F")
p21xl = pn.read_excel(data, sheet_name='P21<0.05',usecols="B,C,D,E,F")
p35xl = pn.read_excel(data, sheet_name='P35<0.05',usecols="B,C,D,E,F")


tmpArray = [p6xl.size, p12xl.size, p18xl.size, p21xl.size, p35xl.size]
max = 0
length = len(tmpArray)

for i in range(1, length):
    if(tmpArray[i] > max):
        max = tmpArray[i]

print(p6xl.columns.values)


# this creates the array of every unique 'Gene Symbol' from every spreadsheet
#Change indexArray into a set
indexArray = []
for i in range(1, max):

    if((i < len(p6xl.index)) and (not (p6xl.at[i, 'Gene Symbol'] in indexArray))):
        indexArray.append(p6xl.at[i, 'Gene Symbol'])

    if((i < len(p12xl.index)) and (not (p12xl.at[i, 'Gene Symbol'] in indexArray))):
        indexArray.append(p12xl.at[i, 'Gene Symbol'])

    if((i < len(p18xl.index)) and (not (p18xl.at[i, 'Gene Symbol'] in indexArray))):
        indexArray.append(p18xl.at[i, 'Gene Symbol'])

    if((i < len(p21xl.index)) and (not (p21xl.at[i, 'Gene Symbol'] in indexArray))):
        indexArray.append(p21xl.at[i, 'Gene Symbol'])

    if((i < len(p35xl.index)) and (not (p35xl.at[i, 'Gene Symbol'] in indexArray))):
        indexArray.append(p35xl.at[i, 'Gene Symbol'])

print(len(indexArray))


tmp_hash = {}
NaN = np.nan
#change length to loop to the end of the gene list
for i in range(0, len(indexArray)):
    tmp = p6xl['Gene Symbol'] == indexArray[i]
    p6tmp = p6xl[tmp]
    #check here if p6tmp is empty and then fill it w/ na values as placeholder
    if p6tmp.empty:
        df1 = pn.DataFrame([[NaN] * len(p6tmp.columns)], columns=p6tmp.columns)
        p6tmp = df1.append(p6tmp, ignore_index=True)

    tmp = p12xl['Gene Symbol'] == indexArray[i]
    p12tmp = p12xl[tmp]
    if p12tmp.empty:
        df1 = pn.DataFrame([[NaN] * len(p12tmp.columns)], columns=p12tmp.columns)
        p12tmp = df1.append(p12tmp, ignore_index=True)

    tmp = p18xl['Gene Symbol'] == indexArray[i]
    p18tmp = p18xl[tmp]
    if p18tmp.empty:
        df1 = pn.DataFrame([[NaN] * len(p18tmp.columns)], columns=p18tmp.columns)
        p18tmp = df1.append(p18tmp, ignore_index=True)

    tmp = p21xl['Gene Symbol'] == indexArray[i]
    p21tmp = p21xl[tmp]
    if p21tmp.empty:
        df1 = pn.DataFrame([[NaN] * len(p21tmp.columns)], columns=p21tmp.columns)
        p21tmp = df1.append(p21tmp, ignore_index=True)

    tmp = p35xl['Gene Symbol'] == indexArray[i]
    p35tmp = p35xl[tmp]
    if p35tmp.empty:
        df1 = pn.DataFrame([[NaN] * len(p35tmp.columns)], columns=p35tmp.columns)
        p35tmp = df1.append(p35tmp, ignore_index=True)

    #print(p6tmp.values[0])
    #need to make this append to my hash

    tmp_hash.update({indexArray[i]: {
        'Description': p6tmp['Description'].values[0],
        'P6': 
            {
            'AR': p6tmp['Abundance Ratio: (P6) / (Ctrl)'].values[0],
            'P': p6tmp['Abundance Ratio P-Value: (P6) / (Ctrl)'].values[0]
            },
        'P12': 
            {
            'AR': p12tmp['Abundance Ratio: (P12) / (Ctrl)'].values[0],
            'P': p12tmp['Abundance Ratio P-Value: (P12) / (Ctrl)'].values[0]
            },
        'P18': 
            {
            'AR': p18tmp['Abundance Ratio: (P18) / (Ctrl)'].values[0],
            'P': p18tmp['Abundance Ratio P-Value: (P18) / (Ctrl)'].values[0]
            },
        'P21': 
            {
            'AR': p21tmp['Abundance Ratio: (P21) / (Ctrl)'].values[0],
            'P': p21tmp['Abundance Ratio P-Value: (P21) / (Ctrl)'].values[0]
            },
        'P35': 
            {
            'AR': p35tmp['Abundance Ratio: (P35) / (Ctrl)'].values[0],
            'P': p35tmp['Abundance Ratio P-Value: (P35) / (Ctrl)'].values[0]
            },
        #here I'm adding 4 different selected intervals for the algorithm to inspect the slope
        #these slopes need to be normalized.
        'Factors': 
            {
            'P6toP35Factor': p35tmp['Abundance Ratio: (P35) / (Ctrl)'].values[0]/p6tmp['Abundance Ratio: (P6) / (Ctrl)'].values[0],
            'P6toP18Factor': p18tmp['Abundance Ratio: (P18) / (Ctrl)'].values[0]/p6tmp['Abundance Ratio: (P6) / (Ctrl)'].values[0],
            'P6toP12Factor': p12tmp['Abundance Ratio: (P12) / (Ctrl)'].values[0]/p6tmp['Abundance Ratio: (P6) / (Ctrl)'].values[0],
           },
            
        }
    })

#writes the JSON
fh = open('NLGN3DevHash.json','w')
json.dump(tmp_hash,fh)
fh.close()
