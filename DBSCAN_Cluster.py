#import required libraries
import plotly.express as pl
import json
import numpy
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
import math




#create the dataframe that stores the gene name with the coordinate
def importJson(fileName):
    fh = open(fileName,"r")
    hash = json.load(fh)
    fh.close()
    return hash

def createDataFrame(hash):
    print("createDataFrame() called")
    geneCoor_df = pd.DataFrame()
    #columns = ['gene','x','y','z','cluster','cellular component','biological process']
    tmpdata = []
    for k,v in hash.items() :
        P6toP12 = v["Factors"]["P6toP12Factor"]
        P6toP18 = v["Factors"]["P6toP18Factor"]
        P6toP35 = v["Factors"]["P6toP35Factor"]
        
        if( not ( math.isnan(P6toP12) or  math.isnan(P6toP18) 
        or math.isnan(P6toP35))):
            if (P6toP12<10 and P6toP18<10 and P6toP35<10): #this line is to remove the outliers
                tmpdata.append([k,P6toP12,P6toP18,P6toP35])
                
                #tmpdf = pd.DataFrame({"Gene Name":k,
               # "x":P6toP12,"y":P6toP18,"z":P6toP35},index = i)
                #print (tmpdf)
    #this line returns the final dataFrame
    geneCoor_df = pd.DataFrame(tmpdata,columns = ["Gene Name","P6 to P12","P6 to P18","P6 to P35"])
    print (geneCoor_df)
    return geneCoor_df

def addClustersToDataFrame(origDF, clusterArray):
    print("addClusterstoDataFrame() called")
    origDF["Cluster"] = clusterArray
    return origDF


def clusterData(geneCoor_df):
    print("clusterData() called")
    arrayToBeClustered = []
    print(geneCoor_df)
    for row in geneCoor_df.itertuples():
        arrayToBeClustered.append([row[2],row[3],row[4]])


    print(arrayToBeClustered)
    clustering = DBSCAN(eps=0.1, min_samples=4).fit(arrayToBeClustered)
    clusterArray = clustering.labels_
    return clusterArray

def removeOutliers(gene_df):
    gene_df = gene_df[gene_df.Cluster != -1]
    return  gene_df

def generateGraph(gene_DF):
    print("generateGraph() called")
    #coorArray = convertDFCoorIntoArrays(gene_DF)
    fig = pl.scatter_3d(gene_DF, x= "P6 to P12",
                y="P6 to P18",
                z="P6 to P35",
                color="Cluster", 
                hover_data = ['Gene Name'])
    fig.show()

def convertDFCoorIntoArrays(gene_DF):
    print("convertDFCoorIntoArrays() called")
    coorArray = [[],[],[],[]]
    #I want the coordinate array to have 4 indivual arrays inside it. one for the x, one for the y coor, one for z, and one for the cluster array
    for row in gene_DF.itertuples():
        coorArray[0].append(row.x)
        coorArray[1].append(row.y)
        coorArray[2].append(row.z)
        coorArray[3].append(row.Cluster)
    return coorArray

hash = importJson("NLGN3DevHash.json")
gene_DF = createDataFrame(hash)
clusters = clusterData(gene_DF)
gene_DF = addClustersToDataFrame(gene_DF,clusters)
gene_DF = removeOutliers(gene_DF)
print(gene_DF)

generateGraph(gene_DF)

 