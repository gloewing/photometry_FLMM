# photometry_FLMM

Code to reproduce analyses and figures from the manuscript: "A Statistical Framework for Analysis of Trial-Level Temporal Dynamics in Fiber Photometry Experiments"

## `fastFMM` R Package

For more information see the `fastFMM` R package repo: https://github.com/gloewing/fastFMM

### Installation

Download the $\texttt{R}$ Package `fastFMM` by running the following command within $\texttt{R}$ or $\texttt{RStudio}$:

```{R}
install.packages("fastFMM", dependencies = TRUE)
```

###  Package Usage

For the usage and a tutorial on package functions, please refer to [fastFMM's Vignette](https://rpubs.com/gloewinger/1110512). 

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

### Calling `fastFMM` from Python
See the ipynb Python version of all user guides in the [Tutorials folder](https://github.com/gloewing/photometry_FLMM/tree/main/Tutorials), which provides examples of using `fastFMM` in Python through the Python packages `rpy2` and `fast_fmm_rpy2`. The tutorials assume the `fastFMM` R package (and all its dependenices), and the `rpy2` Python package have already been installed, but we provide a [user guide](https://github.com/gloewing/photometry_FLMM/blob/main/Tutorials/Python%20rpy2%20installation/R%20and%20rpy2%20installation%20guide.ipynb) for installing `rpy2`, `fastFMM` and `fast_fmm_rpy2` from Python in the `Tutorials` folder.

## Python pip package coming soon!
The pip package should be released shortly. In the meantime, the Python package above can be used. 



