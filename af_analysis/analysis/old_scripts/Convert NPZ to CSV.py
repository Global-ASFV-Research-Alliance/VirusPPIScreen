import os as os
import glob
import numpy as np
import pandas as pd
###  Set Directory ###
Directory = "C:\\Users\\Edward.Spinard\\Desktop\\JAcob\\clustalo-pMSAs_filt95ID_50cov2\\outputs\\*.npz"
######################################################################
OutPutDirectory = Directory.split("*")[0]
Files = glob.glob(Directory)
######################################################################
for npz in Files:
    OutPutFile = OutPutDirectory + npz.split('\\')[-1].split(".")[0] + '.csv'
    #ConvertArray to DataFrame
    NParray = np.load(npz)
    for key in NParray.keys():
        Dataframe = pd.DataFrame(NParray[key])
        Dataframe.to_csv(OutPutFile)
######################################################################
print('finished')
