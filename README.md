# photometry_FLMM

Code to reproduce analyses and figures from the manuscript: "A Statistical Framework for Analysis of Trial-Level Temporal Dynamics in Fiber Photometry Experiments". This repository contains tutorials in both R and Python Jupyter notebooks.

## Installation
In order to run this tutorial the following must be installed:
1. The R Project for Statistical Computing (R)
2. `fastFMM` R Package
3. Optionally Python dependencies if using the Jupyter notebook versions of the tutorials

### 1. Install R
For more information on installing R and system requirements see the official R [documentation](http://r-project.org/) and the [user guide](https://github.com/gloewing/photometry_FLMM/blob/main/Tutorials/Python%20rpy2%20installation/R%20and%20rpy2%20installation%20guide.ipynb) for installing R, `fastFMM` and `rpy2` for Python in the `Tutorials` folder.

### 2. Install R packages

Due to challenges with cross platform R environment management the R package dependencies are not managed with Renv. Download the $\texttt{R}$ Package `fastFMM` by running the following command within $\texttt{R}$ or $\texttt{RStudio}$:

```{R}
install.packages("fastFMM", dependencies = TRUE)
```

For more information see the [fastFMM R package repo](https://github.com/gloewing/fastFMM).

Additionally part IV of the Tutorial requires the following R packages to be installed:
```R
install.packages('lmerTest')
install.packages('emmeans')
```

### 3. Python dependencies
The tutorials are written in R and also Python. The [fast-fmm-rpy2 package](https://pypi.org/project/fast-fmm-rpy2/) was created to wrap the fastFMM R package. This can be installed directly from PyPI, for example `pip install fast-fmm-rpy2`. However for convenience all the Python dependencies required for these tutorials are managed using [uv](https://astral.sh/blog/uv). The following instructions will get you started with uv.

1. Install uv

uv can be installed from PyPI (ex. `pip install uv`) or as a standalone installer. See uv [docs](https://docs.astral.sh/uv/getting-started/installation/) for more info.

2. Install Photometry_FLMM python dependencies

Use uv to create a virtual environment with Python installed and all Python dependencies as defined in the pyproject.toml file.
```bash
uv sync
```

3. Add the uv managed venv as a Jupyter notebook kernel

```bash
uv run --with jupyter jupyter lab
```

See [Using uv with Jupyter](https://docs.astral.sh/uv/guides/integration/jupyter/) for more info.

> [!NOTE]
> As the name implies `fast-fmm-rpy2` uses the Python package `rpy2` to wrap the R package. Refer to `rpy2` [documentation](https://rpy2.github.io/doc/v3.0.x/html/overview.html#installation) for troubleshooting or any [issues loading shared C libraries](https://github.com/rpy2/rpy2?tab=readme-ov-file#issues-loading-shared-c-libraries).

### Known issues using Python Tutorials on Windows

1. Versions 3.6.0 and 3.6.1 of rpy2 have a known issue with Windows. See the [GitHub issue](https://github.com/rpy2/rpy2/issues/1188) for a full overview. The bug has been fixed and will ship in the next rpy2 release (3.6.2). In the mean time the rpy2 dependency has been rolled back to version <3.6.

2. If you encounter an issue with rpy2 not being able to find the shared C libraries [rpy2 suggests](https://github.com/rpy2/rpy2?tab=readme-ov-file#issues-loading-shared-c-libraries) setting the `LD_LIBRARY_PATH` environment variable:
```
export LD_LIBRARY_PATH="$(uv run python -m rpy2.situation LD_LIBRARY_PATH)":${LD_LIBRARY_PATH}
```
Note if you are not using uv to manage your environment remove `uv run` from the command above.

3. If you encounter additional issues it may help to install [RTools](https://cran.r-project.org/bin/windows/Rtools/).

4. Check current rpy2 [Windows issues on GitHub](https://github.com/rpy2/rpy2/issues?q=is%3Aissue%20state%3Aopen%20label%3AWindows). Windows is not officially supported by rpy2 however the community works hard to make it functional and bug fixes are made.

## Tutorials


###  Package Usage

For the fastFMM usage and a tutorial on package functions, please refer to [fastFMM's Vignette](https://rpubs.com/gloewinger/1110512).

### Photometry Analysis Guide in R
- $\textbf{Part 1}$: [Binary Variables](https://rpubs.com/gloewinger/1159094) 
- $\textbf{Part 2}$: [Testing changes within a trial between 2 periods (baseline vs. cue period)](https://rpubs.com/gloewinger/1159127)
- $\textbf{Part 3}$: [Associations with continuous variables](https://rpubs.com/gloewinger/1159129)
- $\textbf{Part 4}$: [Testing Factor Variables](https://rpubs.com/gloewinger/1159411)
- $\textbf{Part 5}$: [Testing how signal–covariate associations change across trials/sessions](https://rpubs.com/gloewinger/1159601)

See the `Tutorials` folder above for the datasets and Rmarkdown files used to generate the above guides.

### Photometry Analysis Guide in Python
- $\textbf{Part 1}$: [Binary Variables](https://github.com/gloewing/photometry_FLMM/blob/main/Tutorials/Photometry%20FLMM%20Guide%20Part%20I/fastFMM-photometry-binary.ipynb)
- $\textbf{Part 2}$: [Testing changes within a trial between 2 periods (baseline vs. cue period)](https://github.com/gloewing/photometry_FLMM/blob/main/Tutorials/Photometry%20FLMM%20Guide%20Part%20II/fastFMM-photometry-withinTrial.ipynb)
- $\textbf{Part 3}$: [Associations with continuous variables](https://github.com/gloewing/photometry_FLMM/blob/main/Tutorials/Photometry%20FLMM%20Guide%20Part%20III/fastFMM-photometry-Correlation.ipynb)
- $\textbf{Part 4}$: [Testing Factor Variables](https://github.com/gloewing/photometry_FLMM/blob/main/Tutorials/Photometry%20FLMM%20Guide%20Part%20IV/fastFMM-photometry-ANOVA.ipynb)
- $\textbf{Part 5}$: [Testing how signal–covariate associations change across trials/sessions](https://github.com/gloewing/photometry_FLMM/blob/main/Tutorials/Photometry%20FLMM%20Guide%20Part%20V/fastFMM-photometry-Interaction.ipynb)

## References
- Lawrimore, J., Loewinger, G., & Ptinis, A. (2025). fast-fmm-rpy2 (v1.0.2). Zenodo. https://doi.org/10.5281/zenodo.15653132
- Cui E, Loewinger G (2025). fastFMM: Fast Functional Mixed Models using Fast Univariate Inference. R package version 0.4.0, https://github.com/gloewing/fastfmm
- Cui et al. (2022) Implementation of the fast univariate inference approach
- Loewinger et al. (2024) A statistical framework for analysis of trial-level temporal dynamics in fiber photometry experiments.
