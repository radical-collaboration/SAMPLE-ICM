import pandas as pd
import numpy as np
import os
import glob
import click

data_dir = '/pylon5/mc3bggp/vhayot/SAMPLE-ICM/wrapper/1597346843'
ps_csv = 'parameter_space.csv' 
n_samples = 5
output = 'adaptive.csv'

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

        is_edge = lambda x: True if len(x) == 0 else False

        nn_crevliq_min = df[(df["crevliq"]==row["crevliq"]) & (df["cliffvmax"] < row["cliffvmax"])].tail(1)
        nn_crevliq_max = df[(df["crevliq"]==row["crevliq"]) & (df["cliffvmax"] > row["cliffvmax"])].head(1)
        nn_crevliq_min_m = nn_crevliq_min['esl(m)'].to_list()
        nn_crevliq_max_m = nn_crevliq_max['esl(m)'].to_list()

        nn_cliffvmax_min = df[(df["cliffvmax"]==row["cliffvmax"]) & (df["crevliq"] < row["crevliq"])].tail(1)
        nn_cliffvmax_max = df[(df["cliffvmax"]==row["cliffvmax"]) & (df["crevliq"] > row["crevliq"])].head(1)
        nn_cliffvmax_min_m = nn_cliffvmax_min['esl(m)'].to_list()
        nn_cliffvmax_max_m = nn_cliffvmax_max['esl(m)'].to_list()

        if is_edge(nn_crevliq_max_m) and  is_edge(nn_crevliq_min_m):
            row['grad_crevliq'] = row['esl(m)']
        elif is_edge(nn_crevliq_min_m):
            row['grad_crevliq'] = nn_crevliq_max_m[0] - row['esl(m)']
        elif is_edge(nn_crevliq_max_m):
            row['grad_crevliq'] = row['esl(m)'] - nn_crevliq_min_m[0]
        else:
            row['grad_crevliq'] = (nn_crevliq_max_m[0]-nn_crevliq_min_m[0])/2

        row['nn_crevliq'] = (nn_crevliq_min[['crevliq', 'cliffvmax']].values.tolist(), nn_crevliq_max[['crevliq', 'cliffvmax']].values.tolist())

        if is_edge(nn_cliffvmax_max_m) and  is_edge(nn_cliffvmax_min_m):
            row['grad_cliffvmax'] = row['esl(m)']
        elif is_edge(nn_cliffvmax_min_m):
            row['grad_cliffvmax'] = nn_cliffvmax_max_m[0] - row['esl(m)']
        elif is_edge(nn_cliffvmax_max_m):
            row['grad_cliffvmax'] = row['esl(m)'] - nn_cliffvmax_min_m[0]
        else:
            row['grad_cliffvmax'] = (nn_cliffvmax_max_m[0]-nn_cliffvmax_min_m[0])/2
        row['nn_cliffvmax'] = (nn_cliffvmax_min[['crevliq', 'cliffvmax']].values.tolist(), nn_cliffvmax_max[['crevliq', 'cliffvmax']].values.tolist())

        return row

    df = df.apply(lambda row: get_gradient(row, df), axis=1)

    return df

def select_candidates(ps_csv, df, output, n_samples, threshold=0.2):

    df = df[(df['grad_cliffvmax'] > threshold) | (df['grad_crevliq'] > threshold)]

    ps = pd.read_csv(ps_csv)

    def high_gradient(ps_row, df, threshold, param):
        if param == 'crevliq':
            df = df[(df['grad_crevliq'] > threshold)]
        else:
            df = df[(df['grad_cliffvmax'] > threshold)]
        
        fix_param_range = lambda x: (x[0][0] if len(x[0]) > 0 else None, x[1][0] if len(x[1]) > 0 else None ) 

        for index, row in df.iterrows():
            nn_name = 'nn_{}'.format(param)
            param_range = fix_param_range(row[nn_name])
            if param_range[0] is None and param_range[1] is None:
                return True
                #return None
            elif param_range[0] is None and row['crevliq'] == param_range[1][0] and row['cliffvmax'] == param_range[1][1]:
                return True
            elif param_range[1] is None and row['crevliq'] == param_range[0][0] and row['cliffvmax'] == param_range[0][1]:
                return True 
            elif (param_range[0] is not None and param_range[1] is not None and 
              row['crevliq'] >= param_range[0][0] and
              row['cliffvmax'] >= param_range[0][1] and
              row['crevliq'] <= param_range[1][0] and
              row['cliffvmax'] <= param_range[1][1]):
                return True

                #if len(crevliq_range[0]) == 0 and ps_row['crevliq'] row['crevliq'] > crevliq_range[0]:
        return False

    ps = ps[~(ps['crevliq'].isin(df['crevliq']) & ps['cliffvmax'].isin(df['cliffvmax']))]
    ps_crevliq_filter = ps.apply(lambda row: high_gradient(row, df, threshold, 'crevliq'), axis=1)
    ps_cliffvmax_filter = ps.apply(lambda row: high_gradient(row, df, threshold, 'cliffvmax'), axis=1)


    ps_crevliq = ps[ps_crevliq_filter]
    ps_cliffvmax = ps[ps_cliffvmax_filter]
    
    ps = pd.concat([ps_crevliq, ps_cliffvmax])

    ps.to_csv(ps_csv)

    selection = ps.sample(n_samples)
    selection.to_csv(output)


@click.command()
@click.argument('ps_csv', type=str, required=True)
@click.argument('n_samples', type=int, required=True)
@click.argument('output', type=str, default="selection.csv")
@click.argument('threshold', type=float, default=0.2)
def main(ps_csv, n_samples, output, threshold):
    df = load_fort22(data_dir)
    #df = df[df['time'] == 500]
    df = df.sort_values(by=['crevliq', 'cliffvmax'])
    df = nearest_neighbor_gradient(df)
    select_candidates(ps_csv, df, output, n_samples, threshold)

if __name__='__main__':
    main()
#arr = numpy_setup(df)
#arr = np.abs(np.gradient(arr))
#print(arr)

