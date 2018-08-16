#!/usr/bin/env python

import sys
import pandas as pd

df = []
for i in sys.argv[1:]:
    df.append(pd.read_csv(i).assign(ntrial = lambda x: i.split('.')[-1]))
pd.concat(df).to_csv('/dev/stdout', index=False)
