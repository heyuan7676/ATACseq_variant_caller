import numpy as np
import pandas as pd
import sys

filename = sys.argv[1]
dat = pd.read_csv('%s' % filename, sep='\t', comment='#', header=None)
dat.loc[np.array(dat[4]) == '<NON_REF>', 4] = dat.loc[np.array(dat[4]) == '<NON_REF>', 3]
dat.to_csv(filename, sep = '\t', index = False)
