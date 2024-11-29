from pathlib import Path
from typing import Callable

import pandas as pd
import pytest
from rpy2 import robjects as ro

from fast_fmm_rpy2.ingest import (
    compare_df_dat,
    pandas_read_in_csv_roundtrip,
    r_read_in_csv_rpy2_convert,
    compare_df_dat_in_r,
)


@pytest.fixture
def binary_filepath() -> Path:
    return Path(r"Tutorials/Photometry FLMM Guide Part I/binary.csv")


@pytest.fixture
def example_filepath() -> Path:
    return Path(r"Tutorials/example_data.csv")


@pytest.fixture
def int_filepath() -> Path:
    return Path(r"Tutorials/example_data.csv")


@pytest.fixture
def anova_filepath() -> Path:
    return Path(r"Tutorials/Photometry FLMM Guide Part IV/anova_data.csv")


@pytest.fixture
def corr_filepath() -> Path:
    return Path(r"Tutorials/Photometry FLMM Guide Part III/corr_data.csv")


def test_binary_csv_exists(binary_filepath):
    assert binary_filepath.exists()


def test_example_csv_exists(example_filepath):
    assert example_filepath.exists()


def test_int_csv_exists(int_filepath):
    assert int_filepath.exists()


def test_anova_csv_exists(anova_filepath):
    assert anova_filepath.exists()


def test_corr_csv_exists(corr_filepath):
    assert corr_filepath.exists()


def test_roundtrip_compare_binary(binary_filepath: Path):
    df: pd.DataFrame = pandas_read_in_csv_roundtrip(binary_filepath)
    dat: ro.vectors.DataFrame = r_read_in_csv_rpy2_convert(binary_filepath)
    results_df: pd.DataFrame = compare_df_dat(df=df, dat=dat)
    assert results_df["float_match"].all()


def test_roundtrip_compare_example(example_filepath: Path):
    df: pd.DataFrame = pandas_read_in_csv_roundtrip(example_filepath)
    dat: ro.vectors.DataFrame = r_read_in_csv_rpy2_convert(example_filepath)
    results_df: pd.DataFrame = compare_df_dat(df=df, dat=dat)
    assert results_df["float_match"].all()


def test_roundtrip_compare_int(int_filepath: Path):
    df: pd.DataFrame = pandas_read_in_csv_roundtrip(int_filepath)
    dat: ro.vectors.DataFrame = r_read_in_csv_rpy2_convert(int_filepath)
    results_df: pd.DataFrame = compare_df_dat(df=df, dat=dat)
    assert results_df["float_match"].all()


def test_roundtrip_compare_anova(anova_filepath: Path):
    df: pd.DataFrame = pandas_read_in_csv_roundtrip(anova_filepath)
    dat: ro.vectors.DataFrame = r_read_in_csv_rpy2_convert(anova_filepath)
    results_df: pd.DataFrame = compare_df_dat(df=df, dat=dat)
    assert results_df["float_match"].all()


def test_roundtrip_compare_corr(corr_filepath: Path):
    df: pd.DataFrame = pandas_read_in_csv_roundtrip(corr_filepath)
    dat: ro.vectors.DataFrame = r_read_in_csv_rpy2_convert(corr_filepath)
    results_df: pd.DataFrame = compare_df_dat(df=df, dat=dat)
    assert results_df["float_match"].all()


def test_compare_df_dat_in_r_binary(binary_filepath: Path):
    result = compare_df_dat_in_r(binary_filepath)
    assert result


def test_compare_df_dat_in_r_example(example_filepath: Path):
    result = compare_df_dat_in_r(example_filepath)
    assert result


def test_compare_df_dat_in_r_int(int_filepath: Path):
    result = compare_df_dat_in_r(int_filepath)
    assert result


def test_compare_df_dat_in_r_anova(anova_filepath: Path):
    result = compare_df_dat_in_r(anova_filepath)
    assert result


def test_compare_df_dat_in_r_corr(corr_filepath: Path):
    result = compare_df_dat_in_r(corr_filepath)
    assert result
