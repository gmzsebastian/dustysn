import os
import corner
from matplotlib import gridspec
from .utils import compute_rhat
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.family': 'serif'})


def plot_corner(object_name, samples_crop, results, n_components, output_dir='.'):
    """
    Plot and save a corner plot of the MCMC samples.

    Parameters
    ----------
    object_name : str
        Name of the object for the output filename
    samples_crop : numpy.ndarray
        The MCMC chain with shape (n_walkers, n_steps, n_dim)
    results : dict
        Dictionary containing results for labeling
    n_components : int
        Number of components in the model (1 or 2)
    output_dir : str
        Directory to save the plot
    """

    print('Plotting Corner...')

    # Define parameter labels based on number of components
    if n_components == 1:
        labels = [r"$\log(M/M_\odot)$", r"$T$ [K]"]
    elif n_components == 2:
        labels = [
            r"$\log(M_{\rm cold}/M_\odot)$",
            r"$T_{\rm cold}$ [K]",
            r"$\log(M_{\rm hot}/M_\odot)$",
            r"$T_{\rm hot}$ [K]"
        ]
    else:
        # For any other number of components, generate generic labels
        labels = []
        for i in range(n_components):
            comp_name = f"_{i+1}" if n_components > 1 else ""
            labels.extend([
                rf"$\log(M{comp_name}/M_\odot)$",
                rf"$T{comp_name}$ [K]"
            ])

    # Get median parameters as truths
    median_params = [results[key][0] for key in results.keys()][:-1]

    # Create the corner plot
    corner.corner(
        samples_crop,
        labels=labels,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 12},
        smooth=2,  # Apply smoothing as requested
        title_fmt='.3f',
        truths=median_params
    )

    # Save the figure
    output_filename = os.path.join(output_dir, f"{object_name}_{n_components}_corner.png")
    plt.savefig(output_filename, bbox_inches='tight')
    plt.clf()
    plt.close('all')


def plot_trace(param_chain, param_values, param_values_log, min_val, max_val,
               title_name, param, log, n_steps, burn_in, repeats, object_name, n_components,
               output_dir='.'):
    '''
    This function plots the trace of a parameter chain with autocorrelation time.

    Parameters
    ----------
    param_chain : np.ndarray
        The chain of the parameter with shape (n_walkers, n_steps).
    param_values : np.ndarray
        The median, upper and lower limits of the parameter.
    param_values_log : np.ndarray
        The median, upper and lower limits of the log of the parameter.
    min_val : float
        The minimum value of the parameter.
    max_val : float
        The maximum value of the parameter.
    title_name : str
        The name of the parameter.
    param : str
        The name of the parameter.
    log : bool
        Whether the parameter is in log scale.
    n_steps : int
        The number of steps in the chain.
    burn_in : float
        The fraction of steps to burn in.
    object_name : str
        The name of the object being fitted.
    n_components : int
        Number of components being fit
    output_dir : str
        The directory to save the plot.
    repeats: int
        Number of times the emcee was repeated
    '''
    # Average walker position
    print(f'Plotting {param} Trace...')
    Averageline = np.average(param_chain.T, axis=1)

    # The chain only after the last sigma-clipping
    n_walkers = param_chain.shape[0]

    # Plot Trace
    plt.subplots_adjust(wspace=0)
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
    ax0 = plt.subplot(gs[0])
    ax0.axhline(param_values[0], color='r', lw=2.0, linestyle='--', alpha=0.75)
    ax0.axhline(param_values[0] - param_values[2], color='r', lw=1.0, linestyle='--', alpha=0.50)
    ax0.axhline(param_values[0] + param_values[1], color='r', lw=1.0, linestyle='--', alpha=0.50)
    ax0.plot(Averageline, lw=1.0, color='b', alpha=0.75)
    ax0.plot(param_chain.T, '-', color='k', alpha=0.2, lw=0.5)

    # Add Gelman-Rubin statistic for multiple chains
    if param_chain.shape[0] > 1:
        try:
            # Split chains into two halves to compute Gelman-Rubin
            if repeats > 1:
                # Reshape to (n_walkers, repeats, n_steps)
                reshaped_chain = param_chain.reshape(n_walkers, repeats, n_steps)
                # Extract post burn-in samples for each repeat
                post_burn_chains = reshaped_chain[:, :, int(n_steps * burn_in):]

                rhat_values = []

                for walker_idx in range(n_walkers):
                    # Extract chains for this walker across repeats
                    walker_chains = post_burn_chains[walker_idx]  # shape: (repeats, post_burn_steps)

                    if walker_chains.shape[0] >= 2:  # Need at least 2 chains
                        try:
                            walker_rhat = compute_rhat(walker_chains)
                            if not np.isnan(walker_rhat):
                                rhat_values.append(walker_rhat)
                        except Exception:
                            pass

                if rhat_values:
                    mean_rhat = np.mean(rhat_values)
                    ax0.text(0.65, 0.05, r'$\hat{{R}}$ = {:.3f}'.format(mean_rhat),
                             transform=ax0.transAxes,
                             bbox=dict(facecolor='white', alpha=0.7))
                else:
                    ax0.text(0.65, 0.05, r'$\hat{{R}}$ calculation failed',
                             transform=ax0.transAxes,
                             bbox=dict(facecolor='white', alpha=0.7))
            else:
                # Just one run, so extract post burn-in samples directly
                post_burn_chains = param_chain[:, int(n_steps * burn_in):]
                # Add a dimension for consistency with multiple repeats case
                post_burn_chains = post_burn_chains.reshape(n_walkers, 1, -1)

                half = n_walkers // 2
                # Use the second half of the burn-in period for more stable chains
                chains_group1 = post_burn_chains[:half, 0, :]  # First half of walkers
                chains_group2 = post_burn_chains[half:, 0, :]  # Second half of walkers

                try:
                    # Compute R-hat between these two groups
                    combined_chains = np.vstack([
                        np.mean(chains_group1, axis=0, keepdims=True),
                        np.mean(chains_group2, axis=0, keepdims=True)
                    ])
                    rhat = compute_rhat(combined_chains)
                    ax0.text(0.65, 0.05, r'$\hat{{R}}$ = {:.3f}'.format(rhat),
                             transform=ax0.transAxes,
                             bbox=dict(facecolor='white', alpha=0.7))
                except Exception:
                    ax0.text(0.65, 0.05, r'$\hat{{R}}$ calculation failed',
                             transform=ax0.transAxes,
                             bbox=dict(facecolor='white', alpha=0.7))
        except Exception as e:
            print(f"Error calculating Gelman-Rubin statistic: {e}")
            pass

    plt.xlim(0, (repeats * n_steps) - 1)
    if log:
        plt.ylim(np.log10(min_val), np.log10(max_val))
    else:
        plt.ylim(min_val, max_val)

    title_string = r"$%s^{+%s}_{-%s}$" % (np.round(param_values[0], 5), np.round(param_values[1], 5),
                                          np.round(param_values[2], 5))
    if log:
        title_string += '  = log(' + r"$%s^{+%s}_{-%s}$" % (np.round(param_values_log[0], 5),
                                                            np.round(param_values_log[1], 5),
                                                            np.round(param_values_log[2], 5)) + ')'
    plt.title(title_string)
    plt.ylabel(title_name)
    plt.xlabel("Step")

    # Plot Histogram
    ax1 = plt.subplot(gs[1])
    plt.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False,
                    labeltop=False, labelright=False, labelbottom=False)
    if log:
        plt.ylim(np.log10(min_val), np.log10(max_val))
    else:
        plt.ylim(min_val, max_val)
    if repeats > 1:
        ax1.hist(np.ndarray.flatten(param_chain[:, -n_steps:]), bins='auto',
                 orientation="horizontal", color='k', range=(min_val, max_val))
    else:
        ax1.hist(np.ndarray.flatten(param_chain[:, -int(n_steps*(1-burn_in)):]), bins='auto',
                 orientation="horizontal", color='k', range=(min_val, max_val))
    ax1.axhline(Averageline[-1], color='b', lw=1.0, linestyle='-', alpha=0.75)
    ax1.axhline(param_values[0], color='r', lw=2.0, linestyle='--', alpha=0.75)
    ax1.axhline(param_values[0] - param_values[2], color='r', lw=1.0, linestyle='--', alpha=0.50)
    ax1.axhline(param_values[0] + param_values[1], color='r', lw=1.0, linestyle='--', alpha=0.50)

    # Create output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f"{object_name}_{param}_{n_components}_Trace.jpg")
    plt.savefig(output_path, bbox_inches='tight', dpi=200)
    plt.clf()
    plt.close('all')
