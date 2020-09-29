import pandas as pd
import numpy as np
import os
import glob

data_dir = '/pylon5/mc3bggp/vhayot/SAMPLE-ICM/wrapper/1597346843'

def load_fort22(data_dir):
    f22_paths = glob.glob(os.path.join(data_dir, "**", "fort.22"))
    run_names = [os.path.basename(os.path.dirname(f)) for f in f22_paths]
    df = pd.concat([pd.read_csv(f22, delim_whitespace=True, header=2) for f22 in f22_paths], keys=run_names) \
           .reset_index(level=1, drop=True).rename_axis('run').reset_index()
    df['crevliq'] = df['run'].apply(lambda x: int(x.split('_')[1].strip('C')))
    df['cliffvmax'] = df['run'].apply(lambda x: int(x.split('_')[2].strip('V')))
    df = df[df['time'] != 'MARK']
    df = df.astype({'time': int, 'cliffvmax': int, 'crevliq': int})
    return df[['crevliq', 'cliffvmax', 'time', 'esl(m)']]

def nearest_neighbor_gradient(df):

    def get_gradient(row, df):
        gradients = []
        print(row)

        nn_crevliq_min = df[(df["crevliq"]==row["crevliq"]) & (df["cliffvmax"] < row["cliffvmax"])].tail(1)['esl(m)'].to_list()
        nn_crevliq_max = df[(df["crevliq"]==row["crevliq"]) & (df["cliffvmax"] < row["cliffvmax"])].head(1)['esl(m)'].to_list()

        nn_cliffvmax_min = df[(df["cliffvmax"]==row["cliffvmax"]) & (df["crevliq"] < row["crevliq"])].tail(1)['esl(m)'].to_list()
        nn_cliffvmax_max = df[(df["cliffvmax"]==row["cliffvmax"]) & (df["crevliq"] < row["crevliq"])].head(1)['esl(m)'].to_list()

        #print(nn_crevliq_min, nn_crevliq_max, nn_cliffvmax_min, nn_cliffvmax_max)

        #gradients.append((nn_crevliq_max-nn_crevliq_min)/2)
        #gradients.append((nn_cliffvmax_max-nn_cliffvmax_min)/2)

        print(gradients)

        return gradients 

    df.head().apply(lambda row: get_gradient(row, df), axis=1)
    #len_c = len(df['crevliq'].unique())
    #len_v = len(df['cliffvmax'].unique())
    #arr = np.zeros((len_c, len_v))

    #df['crevliq'] = df['crevliq'].apply(lambda x: np.where(df['crevliq'].unique()==x)[0][0])
    #df['cliffvmax'] = df['cliffvmax'].apply(lambda x: np.where(df['cliffvmax'].unique()==x)[0][0])
    
    #df = df.pivot(index='crevliq', columns='cliffvmax', values='esl(m)')

    #for index, row in df.iterrows():
    #    #print(row['crevliq'], row['cliffvmax'])
    #    arr[int(row['crevliq']), int(row['cliffvmax'])] = row['esl(m)']
    #return arr

def high_gradient(arr, cond, crevliq_r, cliffvmax_r):
    arr = np.gradient(arr)
     

df = load_fort22(data_dir)
df = df[df['time'] == 500]
df = df.sort_values(by=['crevliq', 'cliffvmax'])
nearest_neighbor_gradient(df)
#arr = numpy_setup(df)
#arr = np.abs(np.gradient(arr))
#print(arr)

