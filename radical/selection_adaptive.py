import pandas as pd
import click

@click.command()
@click.argument('param_space', type=str, required=True)
@click.argument('n_samples', type=int, required=True)
@click.argument('output', type=str, default="lhs.csv")
def lhsampling(param_space, n_samples, output):

    param_space = pd.read_csv(param_space)
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

    param_space.to_csv(output)

if __name__ == '__main__':
    lhsampling()
