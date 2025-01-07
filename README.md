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

### Photometry Analysis Guide
- $\textbf{Part 1}$: [Binary Variables](https://rpubs.com/gloewinger/1159094) 
- $\textbf{Part 2}$: [Testing changes within a trial between 2 periods (baseline vs. cue period)](https://rpubs.com/gloewinger/1159127)
- $\textbf{Part 3}$: [Associations with continuous variables](https://rpubs.com/gloewinger/1159129)
- $\textbf{Part 4}$: [Testing Factor Variables](https://rpubs.com/gloewinger/1159411)
- $\textbf{Part 5}$: [Testing how signal–covariate associations change across trials/sessions](https://rpubs.com/gloewinger/1159601)

See the `Tutorials` folder above for the datasets and Rmarkdown files used to generate the above guides.

### Calling `fastFMM` from Python

See the Python version of [Photometry FLMM Guide Part I](https://github.com/gloewing/photometry_FLMM/blob/main/Tutorials/Photometry%20FLMM%20Guide%20Part%20I/fastFMM-photometry-binary.ipynb) for an example of using `fastFMM` in Python through the Python packages `rpy2` and `fast_fmm_rpy2`. The tutorial assumes the `fastFMM` R package (and all its dependenices), and the `rpy2` Python package have already been installed, but we are working on more documentation for how to install and set up these packages for Python users.  Even if you intend to use the package purely within Python, it may be helpful to first install `fastFMM` from within RStudio to ensure all package dependenices are installed automatically. Finally, see 'python_fastFMM_vignette.py' in the [Github repo](https://github.com/gloewing/photometry_FLMM/tree/main/Tutorials) for a very brief example of using `fastFMM` on Python through the Python package `rpy2`. 



