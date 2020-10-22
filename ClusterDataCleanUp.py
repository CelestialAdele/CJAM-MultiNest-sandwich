import pandas as pd

def UpdateBaumgardtDB():
    listingurl = "https://people.smp.uq.edu.au/HolgerBaumgardt/globular/veldis.html"
    dataFrames = pd.read_html(listingurl)
    
    dataFrame = dataFrames[0]
    
    indexOfPreviousName = 0
    for i in range(dataFrame.shape[0]):
        rowVelocity = dataFrame.at[i, "σ"]
        rowCluster = dataFrame.at[i, "Cluster"]
        if str(rowVelocity) == "nan":
            dataFrame = dataFrame.drop([i])
        elif str(rowCluster) == "nan":
            dataFrame.at[i, "Cluster"] = dataFrame.at[indexOfPreviousName, "Cluster"]
        else:
            indexOfPreviousName = i
    
    dataFrame.reset_index(inplace=True)
    
    dataFrame = dataFrame.drop(columns={"index"})
    
    dataFrame.to_pickle("Updated_Baumgardt_GC_DB.pkl")
    dataFrame.to_csv("Updated_Baumgardt_GC_DB.csv")
    
