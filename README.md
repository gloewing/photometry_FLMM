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

### Calling `fastFMM` from Python

See 'python_fastFMM_vignette.py' for a brief example of how to `fastFMM` on Python. We are working on more documentation. This assumes the package (and all its dependenices) have already been installed. Even if you intend to use the package purely within Python, it may be helpful to first install `fastFMM` from within RStudio to ensure all package dependenices are installed automatically.


