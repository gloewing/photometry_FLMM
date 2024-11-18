from fast_flmm_rpy2.flmm_run import (
    run_with_pandas_dataframe,
    run_with_r_dataframe,
)
from pathlib import Path
import pytest
import numpy as np
from numpy import ndarray
from pandas import DataFrame
from rpy2.robjects.vectors import BoolVector
from rpy2.rinterface_lib.sexp import NULLType


@pytest.fixture
def binary_filepath() -> Path:
    return Path(r"Tutorials/Photometry FLMM Guide Part I/binary.csv")


def test_flmm_runs(binary_filepath: Path) -> None:
    mod = run_with_pandas_dataframe(binary_filepath)
    r_mod = run_with_r_dataframe(binary_filepath)
    for key in dict(mod).keys():
        print(key)
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
            mod_flat = np.array(mod[key])
            r_mod_flat = np.array(mod[key])
        elif isinstance(mod[key], NULLType):
            assert isinstance(mod[key], NULLType) == isinstance(
                r_mod[key], NULLType
            ), f"{key} is NULLType for mod but not r_mod!"
            continue
        else:
            raise ValueError(f"{key} is a {type(mod[key])} variable!")
        assert all(mod_flat == r_mod_flat), (
            "R dataframe vs Pandas DataFrame resulted"
            + f" in different models for {key}"
        )
