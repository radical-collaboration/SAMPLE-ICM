#!/usr/bin/env python

import pandas as pd
import click
import pyDOE as doe
from scipy.stats import randint

@click.command()
@click.argument('ps_csv', type=str, required=True)
@click.argument('n_samples', type=int, required=True)
@click.argument('output', type=str, default="lhs.csv")
def lhsampling(ps_csv, n_samples, output):

    print("lhs sampling")

    param_space = pd.read_csv(ps_csv)
    if n_samples > len(param_space.index):
        n_samples = len(param_space.index)

    unit_lhs = pd.DataFrame(
        doe.lhs(
            n=len(param_space.columns),
            samples=n_samples,
            criterion="maximin",
        ),
        columns=param_space.columns,
    )
    lhs_idx = pd.DataFrame(index=unit_lhs.index)
    lhs = pd.DataFrame(index=unit_lhs.index)

    for i in unit_lhs.columns:
        lhs_idx[i] = (
            randint(0, len(param_space[i].unique()))
            .ppf(unit_lhs[i])
            .astype(int)
        )
        lhs[i] = lhs_idx[i].apply(lambda x: param_space[i].unique()[x])

    lhs[["crevliq", "cliffvmax"]].to_csv(output)

    # save updated param_space                                                   
    updated_ps = param_space[~param_space.isin(lhs)]                             
    updated_ps.to_csv(ps_csv)

if __name__ == '__main__':
    lhsampling()
