
# load packages
import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import pandas as pd

pandas2ri.activate()

# import R packages
base = importr('base')
utils = importr('utils')
stats = importr('stats')
fastFMM = importr('fastFMM')

# load csv file
df = pd.read_csv("/Users/ecui/Documents/Research/Projects/photometry_Loewinger/pyfastFMM/example_data.csv")

# convert it to an R variable
robjects.globalenv['df'] = df

# fit an FUI model
m1 = fastFMM.fui(stats.as_formula('photometry ~ iri + (1 | id)'), data = base.as_symbol('df'))

# estimated coefficient
m1.rx2('betaHat')