import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec


def plot_fui(
    fuiobj,
    num_row=None,
    xlab="Functional Domain",
    title_names=None,
    ylim=None,
    align_x=None,
    x_rescale=1,
    y_val_lim=1.1,
    y_scal_orig=0.05,
    return_data=False,
):
    """
    Plot fixed effects from a functional univariate inference object.

    Parameters:
    -----------
    fuiobj : dict
        A dictionary containing the following keys:
        - betaHat: numpy array of shape (num_vars, num_points) containing coefficient estimates
        - betaHat_var: numpy array of shape (num_points, num_points, num_vars) containing variance estimates (optional)
        - argvals: numpy array of domain points
        - qn: numpy array of quantiles for joint confidence bands (if variance is included)
    num_row : int, optional
        Number of rows for subplot grid
    xlab : str, optional
        Label for x-axis
    title_names : list of str, optional
        Names for each coefficient plot
    ylim : tuple, optional
        Y-axis limits (min, max)
    align_x : float, optional
        Point to align x-axis to (useful for time domain)
    x_rescale : float, optional
        Scale factor for x-axis
    y_val_lim : float, optional
        Factor to extend y-axis limits
    y_scal_orig : float, optional
        Factor to adjust bottom y-axis limit
    return_data : bool, optional
        Whether to return the plotting data

    Returns:
    --------
    matplotlib.figure.Figure or tuple
        If return_data=False, returns the figure
        If return_data=True, returns (figure, list of dataframes)
    """

    num_var = fuiobj["betaHat"].shape[0]  # number of variables to plot
    if num_row is None:
        num_row = int(np.ceil(num_var / 2))
    num_col = int(np.ceil(num_var / num_row))

    align = 0 if align_x is None else align_x * x_rescale

    if title_names is None:
        try:
            title_names = fuiobj["var_names"]
        except KeyError:
            title_names = [f"Variable {i}" for i in range(num_var)]

    # Create figure and subplots
    fig = plt.figure(figsize=(5 * num_col, 4 * num_row))
    gs = GridSpec(num_row, num_col)

    res_list = []

    for r in range(num_var):
        row = r // num_col
        col = r % num_col
        ax = fig.add_subplot(gs[row, col])

        # Create plotting dataframe
        if "betaHat_var" not in fuiobj or fuiobj["betaHat_var"] is None:
            beta_hat_plt = pd.DataFrame(
                {"s": fuiobj["argvals"], "beta": fuiobj["betaHat"][r, :]}
            )

            # Plot estimate
            ax.plot(
                beta_hat_plt["s"] / x_rescale
                - align / x_rescale
                - 1 / x_rescale,
                beta_hat_plt["beta"],
                color="black",
                label="Estimate",
                linewidth=1,
            )

        else:
            var_diag = np.diag(fuiobj["betaHat_var"][:, :, r])
            beta_hat_plt = pd.DataFrame(
                {
                    "s": fuiobj["argvals"],
                    "beta": fuiobj["betaHat"][r, :],
                    "lower": fuiobj["betaHat"][r, :] - 2 * np.sqrt(var_diag),
                    "upper": fuiobj["betaHat"][r, :] + 2 * np.sqrt(var_diag),
                    "lower_joint": fuiobj["betaHat"][r, :]
                    - fuiobj["qn"][r] * np.sqrt(var_diag),
                    "upper_joint": fuiobj["betaHat"][r, :]
                    + fuiobj["qn"][r] * np.sqrt(var_diag),
                }
            )

            # Plot confidence bands
            x_vals = (
                beta_hat_plt["s"] / x_rescale
                - align / x_rescale
                - 1 / x_rescale
            )
            ax.fill_between(
                x_vals,
                beta_hat_plt["lower_joint"],
                beta_hat_plt["upper_joint"],
                color="gray",
                alpha=0.2,
            )
            ax.fill_between(
                x_vals,
                beta_hat_plt["lower"],
                beta_hat_plt["upper"],
                color="gray",
                alpha=0.4,
            )
            ax.plot(
                x_vals,
                beta_hat_plt["beta"],
                color="black",
                label="Estimate",
                linewidth=1,
            )

        # Add horizontal line at y=0
        ax.axhline(y=0, color="red", linestyle="--", alpha=0.75)

        # Set labels and title
        ax.set_xlabel(xlab)
        ax.set_ylabel(f"β{r}(s)")
        ax.set_title(title_names[r], fontweight="bold")

        # Set y limits
        if ylim is not None:
            ax.set_ylim(ylim)
        else:
            if "betaHat_var" not in fuiobj or fuiobj["betaHat_var"] is None:
                y_range = [
                    beta_hat_plt["beta"].min(),
                    beta_hat_plt["beta"].max(),
                ]
            else:
                y_range = [
                    beta_hat_plt["lower_joint"].min(),
                    beta_hat_plt["upper_joint"].max(),
                ]

            y_adjust = y_scal_orig * (y_range[1] - y_range[0])
            y_range[0] -= y_adjust
            y_range = [y * y_val_lim for y in y_range]
            ax.set_ylim(y_range)

        # Add vertical line at x=0 if aligned
        if align_x is not None:
            ax.axvline(
                x=0, color="black", linestyle="--", alpha=0.75, linewidth=0.5
            )

        res_list.append(beta_hat_plt)

    plt.tight_layout()

    if return_data:
        return fig, res_list
    return fig
