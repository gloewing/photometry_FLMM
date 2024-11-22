from fast_flmm_rpy2.flmm_run import (
    run_with_pandas_dataframe,
    run_with_r_dataframe,
    fui,
)
from pathlib import Path
import pytest
import numpy as np
from numpy import ndarray
from pandas import DataFrame
from rpy2 import robjects as ro
from rpy2.robjects.vectors import BoolVector
from rpy2.rinterface_lib.sexp import NULLType
import rpy2.rinterface as rinterface  # type: ignore
from rpy2 import robjects as ro  # type: ignore
from rpy2.robjects import pandas2ri  # type: ignore
from rpy2.robjects.conversion import localconverter  # type: ignore
from fast_flmm_rpy2.flmm_run import fui
import pandas as pd


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
            x = pd.DataFrame(x, index=rownames, columns=colnames)
        finally:
            # plain vector/matrix
            return x


def compare_models(mod, r_mod) -> None:
    for key in dict(mod).keys():
        if isinstance(mod[key], ndarray):
            # check for nans
            mod_data_is_nan = np.isnan(mod[key])
            mod_nan_cols = np.any(mod_data_is_nan, axis=0)
            if np.any(mod_nan_cols):
                mod_data = mod[key][:, ~mod_nan_cols]
                r_mod_data = r_mod[key][:, ~mod_nan_cols]
            else:
                mod_data = mod[key]
                r_mod_data = r_mod[key]
            mod_flat = mod_data.flatten()
            r_mod_flat = r_mod_data.flatten()
        elif isinstance(mod[key], DataFrame):
            mod_flat = mod[key].to_numpy().flatten()
            r_mod_flat = r_mod[key].to_numpy().flatten()
        elif isinstance(mod[key], BoolVector):
            mod_flat = np.array(mod[key]).flatten()
            r_mod_flat = np.array(mod[key]).flatten()
        elif isinstance(mod[key], NULLType):
            assert isinstance(mod[key], NULLType) == isinstance(
                r_mod[key], NULLType
            ), f"{key} is NULLType for mod but not r_mod!"
            continue
        else:
            raise ValueError(f"{key} is a {type(mod[key])} variable!")
        try:
            assert all(mod_flat == r_mod_flat), (
                "R dataframe vs Pandas DataFrame resulted"
                + f" in different models for {key}"
            )
        except ValueError:
            raise ValueError(f"{key} is a {type(mod[key])} variable!")


@pytest.mark.parametrize(
    "csv_filepath,formula,parallel,import_rules",
    [
        (
            Path(r"Tutorials/Photometry FLMM Guide Part I/binary.csv"),
            "photometry ~ cs + (1 | id)",
            True,
            local_rules,
        ),
        (
            Path(r"Tutorials/Photometry FLMM Guide Part I/binary.csv"),
            "photometry ~ cs + (cs | id)",
            True,
            local_rules,
        ),
        (
            Path(r"Tutorials/Photometry FLMM Guide Part III/corr_data.csv"),
            "photometry ~ cs + (1 | id)",
            True,
            local_rules,
        ),
        (
            Path(r"Tutorials/Photometry FLMM Guide Part III/corr_data.csv"),
            "photometry ~ cs + (cs | id)",
            True,
            local_rules,
        ),
    ],
)
def test_fui_compare(csv_filepath, formula, parallel, import_rules) -> None:
    bool_map: dict = {True: "TRUE", False: "FALSE"}
    ro.r(f'dat <- read.csv("{str(csv_filepath.absolute())}")')
    ro.r("library(fastFMM)")
    ro.r(f"mod <- fui({formula}, data = dat, parallel = {bool_map[parallel]})")
    with localconverter(import_rules):
        r_mod = ro.r("mod")
    mod = fui(csv_filepath, formula, parallel, import_rules)
    compare_models(mod, r_mod)
