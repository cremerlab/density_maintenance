
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import bokeh.plotting
import bokeh.io
import bokeh.palettes
import bokeh.themes
import corner
import altair as alt
from bokeh.models import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import AnchoredText
import seaborn as sns


def lit_mapper():
    """
    Returns a dictionary mapping sources of the E. coli data with standard colors
    and glyphs. This ensures constant marking of data across plots.
    """
    rng = np.random.default_rng(666)
    colors, _ = get_colors()
    mapper = {
        'Bremer & Dennis 2008': {'m': 'X'},
        'Brunschede et al. 1977': {'m': 's'},
        'Dai et al. 2016': {'m': 'o'},
        'Forchhammer & Lindahl 1971': {'m': 'v'},
        'Li et al. 2014': {'m': 'd'},
        'Schmidt et al. 2016': {'m': '8'},
        'Scott et al. 2010': {'m': '^'},
        'Wu et al. 2021': {'m': '<'},
        'Bremer & Dennis, 1996': {'m': '>'},
        'Dalbow & Young, 1975': {'m': 'P'},
        'Young & Bremer, 1976': {'m': 'h'},
        'Skjold et al. 1973': {'m': '*'},
        'Dong et al. 1996': {'m': 'p'},
        'Dong et al. 1995': {'m': 'v'},
        'Bentley et al. 1990': {'m': 'X'},
        'Erickson et al. 2017': {'m': 'o'},
        'Oldewurtle et al. 2021': {'m': 's'},
        'Mori et al. 2017': {'m': '*'},
        'Sloan and Urban, 1976': {'m': 'h'},
        'Li et al. 2018': {'m': '>'},
        'Korem Kohanim et al. 2018': {'m': 'd'},
        'Panlilio et al. 2021': {'m': 'p'},
        'Basan et al. 2015': {'m': '8'},
        'You et al. 2013': {'m': 'h'},
        'Hernandez & Bremer 1993': {'m': 'X'},
        'Farewell & Neidhart 1998': {'m': 'h'},
        'Kepes & Beguin, 1966': {'m': 'o'},
        'Coffman et al. 1971':  {'m': 's'},
        'Morris & Hansen, 1973': {'m': '*'},
        'Schleif et al. 1973': {'m': 'v'},
        'Lacroute & Stent 1968': {'m': 'p'},
        'Dennis & Bremer 1974': {'m': 's'},
        'Albertson & Nyström 1994': {'m': '^'},
        'Gausing 1972': {'m': '>'},
        'Schleif 1967': {'m': '<'},
        'Hernandez & Bremer 1993': {'m': 'v'},
        'Pedersen 1984': {'m': 'X'},
        'Taheri-Araghi et al. 2015': {'m': 'o'},
        'Pierucci 1978': {'m': 'v'},
        'Grossman et al. 1982': {'m': 'X'},
        'Zaritsky & Woldringh 1978': {'m': '<'},
        'Rueba & Woldringh 1980': {'m': 's'},
        'Zaritsky et al. 1993': {'m': '>'},
        'Si et al. 2017': {'m': '^'},
        'Caglar et al. 2017': {'m': 'h'},
        'Li et al. 2014': {'m': 'p'},
        'Mori et al. 2021': {'m': 'P'},
        'Peebo et al. 2015': {'m': '*'},
        'Schmidt et al. 2016': {'m': 'o'},
        'Soufi et al. 2015': {'m': '8'},
        'Valgepea et al. 2013': {'m': 'd'},
        'Basan et al. 2015': {'m': '^'},
        'Dai et al. 2016': {'m': 'o'},
        'Churchward et al. 1982': {'m': '>'},
        'Neidhardt et al. 1992': {'m': 'v'},
        'Wright & Lockhart 1964': {'m': '<'},
        'Chohji et al. 1976': {'m': 'h'},
        'Dennis & Bremer 1987': {'m': 's'},
        'Dennis & Bremer 1974': {'m': '*'},
        'Arike et al. 2012': {'m': 'X'},
        'Balakrishnan et al. 2022': {'m': 'X'},
        'Kubitschek et al. 1983': {'m': 'h'},
        'Mir et al. 2011': {'m': 'p'},
        'Watson et al. 1976': {'m': '^'},
        'Woldringh et al. 1981': {'m': 'v'},
        'Poole 1977': {'m': '>'},
        'Zhu & Dai 2019': {'m': '.'}
    }
    # Set colors rooted in blue
    cmap = sns.color_palette(f"light:{colors['black']}",
                             n_colors=len(mapper)).as_hex()
    rng.shuffle(cmap)
    counter = 0
    for k, _ in mapper.items():
        mapper[k]['c'] = cmap[counter]
        counter += 1
    return mapper


def get_colors(all_palettes=False):
    """
    Generates a dictionary of standard colors and returns a sequential color
    palette.

    Parameters
    ----------
    all_palettes : bool
        If True, lists of `dark`, `primary`, and `light` palettes will be returned. If
        False, only the `primary` palette will be returned.
    """
    # Define the colors
    colors = {
        'dark_black': '#2b2b2a',
        'black': '#3d3d3d',
        'primary_black': '#4c4b4c',
        'light_black': '#8c8c8c',
        'pale_black': '#afafaf',
        'dark_blue': '#154577',
        'blue': '#005da2',
        'primary_blue': '#3373ba',
        'light_blue': '#5fa6db',
        'pale_blue': '#8ec1e8',
        'dark_green': '#356835',
        'green': '#488d48',
        'primary_green': '#5cb75b',
        'light_green': '#99d097',
        'pale_green': '#b8ddb6',
        'dark_red': '#79302e',
        'red': '#a3433f',
        'primary_red': '#d8534f',
        'light_red': '#e89290',
        'pale_red': '#eeb3b0',
        'dark_gold': '#84622c',
        'gold': '#b1843e',
        'primary_gold': '#f0ad4d',
        'light_gold': '#f7cd8e',
        'pale_gold': '#f8dab0',
        'dark_purple': '#43355d',
        'purple': '#5d4a7e',
        'primary_purple': '#8066ad',
        'light_purple': '#a897c5',
        'pale_purple': '#c2b6d6'
    }

    # Generate the sequential color palettes.
    keys = ['black', 'blue', 'green', 'red', 'purple', 'gold']
    dark_palette = [colors[f'dark_{k}'] for k in keys]
    primary_palette = [colors[f'primary_{k}'] for k in keys]
    light_palette = [colors[f'light_{k}'] for k in keys]

    # Determine what to return.
    if all_palettes:
        palette = [dark_palette, primary_palette, light_palette]
    else:
        palette = primary_palette

    return [colors, palette]


def matplotlib_style(return_colors=True, return_palette=True, **kwargs):
    """
    Assigns the plotting style for matplotlib generated figures.

    Parameters
    ----------
    return_colors : bool
        If True, a dictionary of the colors is returned. Default is True.
    return_palette: bool
        If True, a sequential color palette is returned. Default is True.
    """
    # Define the matplotlib styles.
    rc = {
        # Axes formatting
        "axes.facecolor": "#f0f3f7",
        "axes.edgecolor": "#ffffff",  # 5b5b5b",
        "axes.labelcolor": "#5b5b5b",
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "axes.axisbelow": True,
        "axes.linewidth": 0.15,
        "axes.grid": True,

        # Formatting of lines and points.
        "lines.linewidth": 0.5,
        "lines.dash_capstyle": "butt",
        "patch.linewidth": 0.25,
        "lines.markeredgecolor": '#f0f3f7',
        "lines.markeredgewidth": 0.5,

        # Grid formatting
        "grid.linestyle": '-',
        "grid.linewidth": 0.5,
        "grid.color": "#FFFFFF",

        # Title formatting
        "axes.titlesize": 8,
        "axes.titleweight": 700,
        "axes.titlepad": 3,
        "axes.titlelocation": "left",

        # Axes label formatting.
        "axes.labelpad": 0,
        "axes.labelweight": 700,
        "xaxis.labellocation": "center",
        "yaxis.labellocation": "center",
        "axes.labelsize": 8,
        "axes.xmargin": 0.03,
        "axes.ymargin": 0.03,

        # Legend formatting
        "legend.fontsize": 6,
        "legend.labelspacing": 0.25,
        "legend.title_fontsize": 6,
        "legend.frameon": True,
        "legend.edgecolor": "#5b5b5b",

        # Tick formatting
        "xtick.color": "#5b5b5b",
        "ytick.color": "#5b5b5b",
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "xtick.major.width": 0.25,
        "ytick.major.width": 0.25,
        "xtick.major.pad": 2,
        "ytick.major.pad": 2,
        "xtick.minor.size": 0,
        "ytick.minor.size": 0,

        # General Font styling
        "font.family": "sans-serif",
        "font.family": "Lato",
        "font.weight": 400,  # Weight of all fonts unless overriden.
        "font.style": "normal",
        "text.color": "#3d3d3d",  # "#5b5b5b",

        # Higher-order things
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "figure.dpi": 300,
        "errorbar.capsize": 1,
        "savefig.bbox": "tight",
        "mathtext.default": "regular",
    }
    matplotlib.style.use(rc)

    # Load the colors and palettes.
    colors, palette = get_colors(**kwargs)
    sns.set_palette(palette)

    # Determine what, if anything should be returned
    out = []
    if return_colors == True:
        out.append(colors)
    if return_palette == True:
        out.append(palette)

    if len(out) == 1:
        return out[0]
    else:
        return out


def cell_gallery(biometrics,
                 cells,
                 anatomy,
                 fname,
                 suptitle=None,
                 cols=10):
    """
    Creates and saves a figure of each cell segmentation mask, its labeled
    contours, and its dimensions, grouped by the cell ID and the image from
    which it was segmented.

    Parameters
    -----------
    biometrics : pandas DataFrame
        A DataFrame containing dimension measurements for each cell.
    cells : pandas DataFrame
        A Dataframe containing the cell intensity images.
    anatomy :  pandas DataFrame
        A DataFrame containing contour coordinates for each individual cell and
        their labeled regions. This DataFrame can be the default output from
        `size.image.assign_anatomy`
    fname : str
        The filename of the plot output.
    suptitle : str or None
        The title of the entire plot. If `None`, no title is added.
    cols : int
        The number of columns in the displayed plot. Default is 10

    Returns
    -------
    [fig, ax] : matplotlib figure canvas and axes
        The figure canvas and axes artists.
    """
    cor, _ = matplotlib_style()
    # Determine the total number of cells to display
    n_cells = biometrics.groupby(['cell_id', 'image']).ngroups
    n_rows = int(np.ceil(n_cells/cols))
    figsize = (8.5, n_rows)
    n_blank = (n_rows * cols) - n_cells
    fig, ax = plt.subplots(n_rows, cols, figsize=figsize)
    ax = ax.ravel()
    for a in ax:
        a.grid(False)
        a.set_xticks([])
        a.set_yticks([])
    if n_blank != 0:
        for i in range(n_blank):
            ax[-(i+1)].axis('off')
    _idx = 0
    for g, d in biometrics.groupby(['cell_id', 'image']):
        # Get the cell image

        image = cells[(cells['cell_id'] == g[0]) & (
            cells['image'] == g[1])]['cell_image'].values[0]

        # Get the contours
        anat = anatomy[(anatomy['cell_id'] == g[0]) &
                       (anatomy['image'] == g[1])]

        ax[_idx].imshow(image, cmap='Greys_r')
        caps = anat[(anat['component'] == 'top') |
                    (anat['component'] == 'bottom')]
        sides = anat[(anat['component'] == 'left') |
                     (anat['component'] == 'right')]
        ax[_idx].plot(caps['x_coords'], caps['y_coords'], '.', ms=2, color=cor['light_blue'],
                      markeredgewidth=0)

        ax[_idx].plot(sides['x_coords'], sides['y_coords'], '.', ms=2, color=cor['primary_green'],
                      markeredgewidth=0)

        # Assign bar coords for width
        median_left_x = sides[sides['component']
                              == 'left']['x_coords'].median()
        median_right_x = sides[sides['component']
                               == 'right']['x_coords'].median()
        median_left_y = sides[sides['component']
                              == 'left']['y_coords'].median()
        median_right_y = sides[sides['component']
                               == 'right']['y_coords'].median()

        # Assign bar coords for length
        median_top_x = caps[caps['component']
                            == 'top']['x_coords'].median()
        median_bottom_x = caps[caps['component']
                               == 'bottom']['x_coords'].median()
        median_top_y = caps[caps['component']
                            == 'top']['y_coords'].max()
        median_bottom_y = caps[caps['component']
                               == 'bottom']['y_coords'].min()

        # Plot measurement bars
        ax[_idx].plot([median_left_x, median_right_x], [
            median_left_y, median_right_y], '-', color='white')
        ax[_idx].plot([median_left_x, median_right_x], [
            median_left_y, median_right_y], '|', ms=3, color='white')
        ax[_idx].plot([median_top_x, median_bottom_x], [
            median_top_y, median_bottom_y], '-', color='white')
        ax[_idx].plot([median_top_x, median_bottom_x], [
            median_top_y, median_bottom_y], '_', ms=3, color='white')

        # Add text label
        ax[_idx].set_xlabel(
            f'$w$={d["width_median"].values[0]:0.2f} µm', fontsize=5)
        ax[_idx].set_ylabel(
            f'$\ell$={d["length"].values[0]:0.2f} µm', fontsize=5)
        _idx += 1

    plt.tight_layout()
    plt.subplots_adjust(wspace=0)
    if suptitle != None:
        fig.text(0, 1,  suptitle, fontsize=10)
    plt.savefig(fname, dpi=200)
    return [fig, ax]


def diagnostic_size_viz(samples,
                        data,
                        metadata,
                        dst,
                        hypervars=['width_mu',
                                   'vol_mu',
                                   'length_alpha',
                                   'length_beta',
                                   'SAV_mu'],
                        thin=10
                        ):
    cor, _ = matplotlib_style()
    out = f'{dst}/{metadata["strain"]}_{metadata["carbon"]}_{metadata["temp"]}C_{metadata["oe"]}OE_{metadata["ind"]}'
    # Generate the corner plot
    fig = plt.figure(figsize=(4, 4))
    _ = corner.corner(samples,
                      var_names=hypervars,
                      divergences=True,
                      smooth=1,
                      fig=fig,
                      divergences_kwargs={
                          'color': cor['light_blue'],
                          'ms': 4,
                          'markeredgewidth': 0})
    fig.text(
        0, 1, f'{metadata["strain"]}; {metadata["carbon"]} ; {metadata["temp"]} °C ; overexpression: {metadata["oe"]} @ {metadata["ind"]} ng/mL')
    plt.savefig(f'{out}_hyperparameter_joint_dists.pdf', bbox_inches='tight')

    # Generate the posterior predictive plots
    n_reps = data['idx'].max()
    fig, ax = plt.subplots(n_reps, 4, figsize=(6, n_reps * 1))

    # Format axes
    width_rep = samples.posterior.width_rep.to_dataframe().reset_index()
    length_rep = samples.posterior.length_rep.to_dataframe().reset_index()
    volume_rep = samples.posterior.volume_rep.to_dataframe().reset_index()
    sav_rep = samples.posterior.SAV_rep.to_dataframe().reset_index()
    titles = ['width', 'length', 'volume', 'SA/V']
    labs = ['µm', 'µm', 'µm$^3$', 'µm$^{-1}$']
    if n_reps == 1:
        for i in range(4):
            ax[i].set_title(titles[i], fontsize=6)
            ax[i].set_xlabel(labs[i], fontsize=6)
        ax[0].set_ylabel('ECDF')
    else:
        for i in range(4):
            ax[0, i].set_title(titles[i], fontsize=6)
            ax[-1, i].set_xlabel(labs[i], fontsize=6)
        for i in range(n_reps):
            ax[i, 0].set_ylabel('ECDF')

    # Generate ECDFs
    for i in range(n_reps):
        # Get width draws for each day
        loc = np.where(data['idx'].values == i + 1)[0]

        for j, ppc in enumerate([width_rep, length_rep, volume_rep, sav_rep]):
            k = 0
            if n_reps == 1:
                _ax = ax[j]
            else:
                _ax = ax[i, j]
            for g, _d in ppc.groupby(['chain', 'draw']):
                if k % thin == 0:
                    _x = np.sort(_d[_d[_d.keys()[-2]].isin(loc)]
                                 [_d.keys()[-1]].values)
                    _y = np.arange(len(_x)) / len(_x)
                _ax.plot(_x, _y, 'k-', lw=0.1, alpha=0.1)
                k += 1
        rep = data[data['idx'] == (i + 1)]
        width_x = np.sort(rep['width_median'].values)
        length_x = np.sort(rep['length'].values)
        vol_x = np.sort(rep['volume'].values)
        sav_x = np.sort(rep['surface_to_volume'].values)
        y = np.arange(len(sav_x)) / len(sav_x)
        for j, val in enumerate([width_x, length_x, vol_x, sav_x]):
            if n_reps == 1:
                _ax = ax[j]
            else:
                _ax = ax[i, j]
            _ax.plot(val, y, '-', lw=1.5, color=cor['primary_red'])
    fig.text(
        0, 1, f'{metadata["strain"]}; {metadata["carbon"]} ; {metadata["temp"]} °C ; overexpression: {metadata["oe"]} @ {metadata["ind"]} ng/mL')
    plt.tight_layout()

    # save figure
    plt.savefig(f'{out}_ppc.pdf', bbox_inches='tight')


def diagnostic_growth_viz(samples,
                          data,
                          metadata,
                          dst,
                          thin=10,
                          hypervars=['mu', 'tau', 'sigma']):
    cor, _ = matplotlib_style()
    out = f'{dst}/{metadata["strain"]}_{metadata["carbon"]}_{metadata["temp"]}C_{metadata["oe"]}OE_{metadata["ind"]}'
    # Generate the corner plot
    fig = plt.figure(figsize=(4, 4))
    _ = corner.corner(samples,
                      var_names=hypervars,
                      divergences=True,
                      smooth=1,
                      fig=fig,
                      divergences_kwargs={
                          'color': cor['light_blue'],
                          'ms': 4,
                          'markeredgewidth': 0})
    fig.text(
        0, 1, f'{metadata["strain"]}; {metadata["carbon"]} ; {metadata["temp"]} °C ; overexpression: {metadata["oe"]} @ {metadata["ind"]} ng/mL')
    plt.savefig(f'{out}_growth_hyperparameter_joint_dists.pdf',
                bbox_inches='tight')
    plt.close()
    # Generate the ppc plot
    n_reps = data['idx'].values.max()
    optical_density_rep = samples.posterior.optical_density_rep.to_dataframe().reset_index()

    n_cols = 2
    if n_reps > n_cols:
        n_rows = int(np.ceil(n_reps / n_cols))
    else:
        n_rows = 1
    # Instantiate canvas for ppc
    fig, ax = plt.subplots(n_rows, n_cols, figsize=(6, n_reps))
    axes = ax.ravel()
    extra_ax = (n_cols * n_rows) - n_reps
    for i in range(np.abs(extra_ax)):
        axes[-(i+1)].axis('off')
    for a in axes:
        a.set_xlabel('elapsed time [hr$^{-1}$]')
        a.set_ylabel('log$_{10}$ OD$_{600nm}$')

    for i in range(n_reps):
        # Plot the data
        _d = data[data['idx'] == i+1]
        time = _d['elapsed_time_hr'].values
        # Get the ppc
        loc = np.where(data['idx'].values == i+1)[0]
        _ppc = optical_density_rep[optical_density_rep['optical_density_rep_dim_0'].isin(
            loc)]
        j = 0
        for _g, __d in _ppc.groupby(['chain', 'draw']):
            if j % thin == 0:
                axes[i].plot(time, np.log10(__d['optical_density_rep'].values),
                             '-', lw=0.1, color=cor['primary_black'])
            j += 1
        axes[i].plot(_d['elapsed_time_hr'], np.log10(_d['od_600nm']),
                     '-o', ms=4, color=cor['primary_red'], lw=1)
    fig.text(
        0, 1, f'{metadata["strain"]}; {metadata["carbon"]} ; {metadata["temp"]} °C ; overexpression: {metadata["oe"]} @ {metadata["ind"]} ng/mL')
    plt.tight_layout()

    # save figure
    plt.savefig(f'{out}_growth_ppc.pdf', bbox_inches='tight')
    plt.close()


def compute_percentiles(df,
                        quantity,
                        groupby,
                        lower_bounds=[5, 10, 15, 20, 25, 30, 35, 40, 45, ],
                        upper_bounds=[95, 90, 85, 80, 75, 70, 65, 60, 55],
                        interval_labels=['90%', '80%', '70%', '60%',
                                         '50%', '40%', '30%', '20%', '10%']):

    # Allow flexibility in what quantities are being supplied
    if type(quantity) != str:
        if type(quantity) != list:
            raise TypeError("`quantity` must be a `str` or list of `str.`")
    else:
        quantity = [quantity]

    # Instantiate the dataframe and loop through every group
    perc_df = pd.DataFrame([])
    if type(groupby) != list:
        groupby = [groupby]

    for g, d in df.groupby(groupby):
        if type(g) != tuple:
            g = (g,)
        # Compute the percentiles for different quantities
        for q in quantity:
            lower = np.percentile(d[f'{q}'].values, lower_bounds)
            upper = np.percentile(d[f'{q}'].values, upper_bounds)
            _df = pd.DataFrame(np.array([lower, upper]).T, columns=[
                               'lower', 'upper'])
            _df['quantity'] = q
            _df['interval'] = interval_labels

            # Add the grouping informaton
            for i in range(len(groupby)):
                _df[groupby[i]] = g[i]
            perc_df = pd.concat([perc_df, _df], sort=False)
    return perc_df


def style_point(source, alpha=0.75):
    mapper = lit_mapper()
    cor, _ = matplotlib_style()
    style = {'marker': mapper[source]['m'],
             'color': mapper[source]['c'],
             'markeredgecolor': cor['primary_black'],
             'markeredgewidth': 0.25,
             'alpha': alpha,
             'label': source,
             'linestyle': 'none'}
    return style
