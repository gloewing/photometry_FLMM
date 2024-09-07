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
$\textbf{Part 1}$: [Binary Variables](https://rpubs.com/gloewinger/1159094) \newline
Part 2: Testing changes within a trial between 2 periods (baseline vs. cue period) – https://rpubs.com/gloewinger/1159127
Part 3: Associations with continuous variables – https://rpubs.com/gloewinger/1159129
Part 4: Testing Factor Variables – https://rpubs.com/gloewinger/1159411
Part 5: Testing how signal–covariate associations change across time https://rpubs.com/gloewinger/1159601

### Calling `fastFMM` from Python

See 'python_fastFMM_vignette.py' in the Tutorials folder for a brief example of using `fastFMM` on Python through the Python package `rpy2`. We are working on more documentation. The tutorial assumes the `fastFMM` R package (and all its dependenices), and the `rpy2` Python package have already been installed.  Even if you intend to use the package purely within Python, it may be helpful to first install `fastFMM` from within RStudio to ensure all package dependenices are installed automatically.


