#!/usr/bin/env python

import pandas as pd
import sys
import subprocess as sp


out_fldr = sys.argv[1]
row_index = int(sys.argv[2])
select_file = sys.argv[3]

model = "/pylon5/mc3bggp/vhayot/SAMPLE-ICM/wrapper/doruniceeg"

params = pd.read_csv(sys.argv[3])[["crevliq", "cliffvmax"]].iloc[row_index]

p = sp.Popen([model, out_fldr, *[str(i) for i in params.to_list()]], stderr=sp.PIPE, stdout=sp.PIPE)

out, err = p.communicate()

print(str(out))
print(str(err))
