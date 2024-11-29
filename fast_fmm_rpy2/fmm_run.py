from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.rinterface as rinterface  # type: ignore
from rpy2 import robjects as ro  # type: ignore
from rpy2.robjects import pandas2ri  # type: ignore
from rpy2.robjects.conversion import localconverter  # type: ignore
from rpy2.robjects.packages import importr  # type: ignore
from rpy2.rinterface_lib.sexp import NULLType  # type: ignore

from fast_fmm_rpy2.ingest import read_csv_in_pandas_pass_to_r

# import R packages
base = importr("base")
utils = importr("utils")
stats = importr("stats")
fastFMM = importr("fastFMM")

local_rules = ro.default_converter + pandas2ri.converter


@local_rules.rpy2py.register(rinterface.FloatSexpVector)
def rpy2py_floatvector(obj):
    x = np.array(obj)
    try:
        # if names is assigned, convert to pandas series
        return pd.Series(x, obj.names)
    except Exception:
        # if dimnames assigned, it's a named matrix,
        # convert to pandas dataframe
        try:
            rownames, colnames = obj.do_slot("dimnames")
            if not isinstance(rownames, NULLType) and not isinstance(
                colnames, NULLType
            ):
                x = pd.DataFrame(x, index=rownames, columns=colnames)
            else:
                x = pd.DataFrame(x, columns=colnames)

        finally:
            # plain vector/matrix
            return x


def run_with_pandas_dataframe(csv_filepath: Path, import_rules=local_rules):
    read_csv_in_pandas_pass_to_r(
        csv_filepath=csv_filepath, r_var_name="py_dat"
    )
    with localconverter(import_rules):
        mod = fastFMM.fui(
            stats.as_formula("photometry ~ cs + (1 | id)"),
            data=base.as_symbol("py_dat"),
        )
    return mod


def run_with_r_dataframe(csv_filepath: Path, import_rules=local_rules):
    ro.r(f'dat = read.csv("{str(csv_filepath.absolute())}")')
    ro.r("mod = fui(photometry ~ cs + (1 | id), data = dat, parallel = TRUE)")
    with localconverter(import_rules):
        r_mod = ro.r("mod")
    return r_mod


def fui(
    csv_filepath: Path,
    formula: str,
    parallel: bool = True,
    import_rules=local_rules,
):
    read_csv_in_pandas_pass_to_r(
        csv_filepath=csv_filepath, r_var_name="py_dat"
    )
    with localconverter(import_rules):
        mod = fastFMM.fui(
            stats.as_formula(formula),
            data=base.as_symbol("py_dat"),
            parallel=parallel,
        )
    return mod
