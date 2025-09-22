import pandas as pd

source = pd.read_csv("bladedat.txt", sep=r"\s+",header=None)
mylist = []
for row in range(0, 18):
    fname = f"FFA_W3-{round(source.iloc[row,3] / 100, 4)}.csv"
    mylist.append(fname)

source[4] = mylist



# use join to join column to dataframe


dummy= 0;