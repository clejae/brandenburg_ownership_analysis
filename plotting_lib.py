#### Plotting
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def plot_corr(df, size=10):
    """Function plots a graphical correlation matrix for each pair of columns in the dataframe.

    Input:
        df: pandas DataFrame
        size: vertical and horizontal size of the plot
    """
    import matplotlib.pyplot as plt

    corr = df.corr()
    fig, ax = plt.subplots(figsize=(size, size))
    ax.matshow(corr)
    plt.xticks(range(len(corr.columns)), corr.columns)
    plt.yticks(range(len(corr.columns)), corr.columns)
    ax.set_xticklabels(corr.columns, rotation=90)


def scatterplot_two_columns(df, col1, col2, out_pth, xminmax=None, yminmax=None,
                            x_label=None, y_label=None, title=None, hue=None,
                            x_log=False, y_log=False, normal_scatter=None):
    import matplotlib.pyplot as plt
    import seaborn as sns
    print(f'Plotting scatterplot of {col1} against {col2}.\n\tSaving at {out_pth}')

    fig, ax = plt.subplots(1, 1, figsize=cm2inch(10, 10))
    if hue:
        sns.scatterplot(x=col1, y=col2, data=df, hue=hue, s=4, linewidth=0.1, ax=ax)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size': 8})
    elif normal_scatter:
        sns.scatterplot(x=col1, y=col2, data=df, s=4, linewidth=0.1, ax=ax)
    else:
        sns.regplot(x=col1, y=col2, data=df, ax=ax, robust=True, n_boot=100,
                    scatter_kws=dict(s=4, linewidths=.7, edgecolors='none'),
                    line_kws={"color": "black", "linewidth": 1}, truncate=True)

    if xminmax:
        plt.xlim(xminmax[0], xminmax[1])
    if yminmax:
        plt.ylim(yminmax[0], yminmax[1])
        # Decorations
        # gridobj.set(xlim=(0.5, 7.5), ylim=(0, 50))

    if title:
        ax.set_title(title)
    else:
        ax.set_title(f"{col1} vs. {col2}")

    if x_log:
        ax.set(xscale="log")
    if y_log:
        ax.set(yscale="log")

    if x_label:
        plt.xlabel(x_label)
    if y_label:
        plt.ylabel(y_label)

    fig.tight_layout()
    plt.savefig(out_pth, dpi=600)
    plt.close()


def scatter_distribution_plot_two_columns(df, col1, col2, out_pth, xminmax=None, yminmax=None,
                            x_label=None, y_label=None, title=None, hue=None,
                            x_log=False, y_log=False, s=4):

    import matplotlib.pyplot as plt
    import seaborn as sns
    print(f'Plotting scatterplot of {col1} against {col2}.\n\tSaving at {out_pth}')

    fig, ax = plt.subplots(1, 1, figsize=cm2inch(10, 10))

    if hue:
        sns.jointplot(x=col1, y=col2, data=df, hue=hue, s=s, linewidth=0.1, ax=ax)
        # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size': 8})
    else:
        sns.jointplot(x=col1, y=col2, data=df, ax=ax, s=s, linewidths=0.1, edgecolor=None)

    if xminmax:
        plt.xlim(xminmax[0], xminmax[1])
    if yminmax:
        plt.ylim(yminmax[0], yminmax[1])

    if title:
        ax.set_title(title)
    else:
        ax.set_title(f"{col1} vs. {col2}")

    if x_log:
        ax.set(xscale="log")
    if y_log:
        ax.set(yscale="log")

    if x_label:
        plt.xlabel(x_label)
    if y_label:
        plt.ylabel(y_label)

    plt.tight_layout()
    plt.savefig(out_pth, dpi=600)
    plt.close()


def lorenz_curve(arr_lst, label_lst, color_lst, legend_title, out_pth, xlabel=None, ylabel=None, title=None,
                 scatter=False, annotation=None):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np

    print("Plot lorenz curve")
    plt.ioff()
    fig, ax = plt.subplots(figsize=cm2inch((20, 20)))

    if len(color_lst) < len(arr_lst):
        print('\tColor list too short, choosing other colors.')
        color_lst = ["#000000",
            "#191970", "#006400", "#bc8f8f", "#ff4500", "#00ff00", "#00ffff", "#0000ff", "#ffff54", "#ff1493", "#ff0000",
         "#ffd700", "#7fff00", "#ba55d3", "#00ff7f", "#00ffff", "#00bfff", "#0000ff", "#ff7f50", "#ff00ff", "#1e90ff",
         "#f0e68c", "#dda0dd", "#ff1493"]

    legend_elements = []
    for i, x_lorenz1 in enumerate(arr_lst):
        color = color_lst[i]
        label = label_lst[i]
        if scatter:
            ax.scatter(np.arange(x_lorenz1.size) / (x_lorenz1.size - 1), x_lorenz1,
                       marker='o', color=color, s=1)
        ax.plot(np.arange(x_lorenz1.size) / (x_lorenz1.size - 1), x_lorenz1, color=color, linewidth=0.5)

        legend_elements.append(Line2D([0], [0], marker='o', color=color, label=label, ms=1))

    ax.plot([0, 1], [0, 1], color='k')
    if title:
        ax.set_title(title, loc='left', fontdict={'size': 12, 'weight': 'bold'})

    ax.legend(handles=legend_elements, loc='upper left', title=legend_title)
    if annotation:
        ax.annotate(annotation, xy=(.2,.8))
    plt.xlim([.0, 1])
    plt.ylim([.0, 1])
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    plt.savefig(out_pth, dpi=600)
    plt.close()


def lorenz_curve_with_zoom_in(arr_lst, label_lst, color_lst, legend_title, out_pth, xlabel=None, ylabel=None, title=None,
                 scatter=False, annotation=None):
    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt
    import numpy as np

    plt.ioff()
    fig, axs = plt.subplots(ncols=2, figsize=cm2inch((43, 20)))

    if len(color_lst) < len(arr_lst):
        print('Color list too short, choosing other colors.')
        color_lst = ["#000000",
            "#191970", "#006400", "#bc8f8f", "#ff4500", "#00ff00", "#00ffff", "#0000ff", "#ffff54", "#ff1493",
            "#ff0000",
            "#ffd700", "#7fff00", "#ba55d3", "#00ff7f", "#00ffff", "#00bfff", "#0000ff", "#ff7f50", "#ff00ff",
            "#1e90ff",
            "#f0e68c", "#dda0dd", "#ff1493"]

    legend_elements = []
    for i, x_lorenz1 in enumerate(arr_lst):
        color = color_lst[i]
        label = label_lst[i]
        if scatter:
            axs[0].scatter(np.arange(x_lorenz1.size) / (x_lorenz1.size - 1), x_lorenz1,
                       marker='o', color=color, s=1)
            axs[1].scatter(np.arange(x_lorenz1.size) / (x_lorenz1.size - 1), x_lorenz1,
                           marker='o', color=color, s=1)
        axs[0].plot(np.arange(x_lorenz1.size) / (x_lorenz1.size - 1), x_lorenz1, color=color, linewidth=0.5)
        axs[1].plot(np.arange(x_lorenz1.size) / (x_lorenz1.size - 1), x_lorenz1, color=color, linewidth=0.5)

        legend_elements.append(Line2D([0], [0], marker='o', color=color, label=label, ms=1))

    axs[0].plot([0, 1], [0, 1], color='k')
    if title:
        axs[0].set_title(title, loc='left', fontdict={'size': 12, 'weight': 'bold'})

    axs[0].legend(handles=legend_elements, loc='upper left', title=legend_title)
    if annotation:
        axs[0].annotate(annotation, xy=(.6, .8))
    if xlabel:
        axs[0].set_ylabel(xlabel)
    if ylabel:
        axs[0].set_xlabel(ylabel)

    axs[1].set_title("Zoom in", loc='left', fontdict={'size': 12, 'weight': 'bold'})
    if xlabel:
        axs[1].set_ylabel(xlabel)
    if ylabel:
        axs[1].set_xlabel(ylabel)
    axs[1].set_xlim([.98, 1])
    axs[1].set_ylim([.4, 1])

    plt.savefig(out_pth, dpi=600)
    plt.close()


def boxplot_by_categories(df, category_col, value_col, out_pth, ylims=None, x_label=None, y_label=None, title=None,
                          order=None, showfliers=False, colour_dict=None):
    import matplotlib.pyplot as plt
    import seaborn as sns

    if not colour_dict:
        colour_dict = sns.color_palette()

    fig, ax = plt.subplots(figsize=cm2inch((20, 10)))
    ax.grid(visible=True, which="major", axis="y", zorder=0)
    ax.set_axisbelow(True)
    if order:
        sns.boxplot(x=category_col, y=value_col, data=df, ax=ax, linewidth=0.5, showfliers=showfliers, order=order,
                    flierprops=dict(markersize=3), palette=colour_dict)
    else:
        sns.boxplot(x=category_col, y=value_col, data=df, ax=ax, linewidth=0.5, showfliers=showfliers,
                    flierprops=dict(markersize=3), palette=colour_dict)
    if ylims:
        ax.set_ylim(ylims)

    if title:
        ax.set_title(title)

    if x_label:
        plt.xlabel(x_label)
    if y_label:
        plt.ylabel(y_label)

    plt.xticks(rotation=45)

    fig.tight_layout()
    plt.savefig(out_pth, dpi=300)
    plt.close()


def kernel_density_plot(df, col, out_pth, hue=None):
    import matplotlib.pyplot as plt
    import seaborn as sns

    print(f'Plotting kernel density of {col} a.\n Saving at {out_pth}')

    fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=cm2inch(10, 10))
    sns.kdeplot(data=df, x=col, ax=ax, legend=True, linestyle='-', linewidth=1,
                color='black', hue=hue) #log_scale=True, cut=0,
    ax.set_title(col, x=.05, y=.995, pad=-14, fontdict={'size': 10})
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # plt.sca(ax)
    # plt.xticks([3, 10, 100, 1000], ['3', '10', '100', '1000'])
    # ax.set_xticks([], minor=True)
    #
    # plt.sca(ax)
    # plt.xticks([.3, 1, 10, 100], ['0.3', '1', '10', '100'])
    # ax.set_xticks([], minor=True)

    fig.tight_layout()
    plt.savefig(out_pth, dpi=300)


def histogramm(df, col, out_pth, x_label=None, y_label=None):
    import matplotlib.pyplot as plt
    import seaborn as sns

    print(f'Plotting histogramm of {col}\n Saving at {out_pth}')

    fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=cm2inch(10, 10))
    x = sns.histplot(data=df, x=col)
    max_y_val = x.dataLim.bounds[-1]
    plt.axvline(df[col].mean(), 0, max_y_val, color='black', linewidth=0.5)
    plt.text(df[col].mean(), 0.8 * max_y_val, f"Mean: {round(df[col].mean(), 1)}", horizontalalignment='left',
             color='black')
    plt.axvline(df[col].median(), 0, max_y_val, color='orange', linewidth=0.5)
    plt.text(df[col].median(), 0.7 * max_y_val, f"Median: {round(df[col].median(), 1)}", horizontalalignment='left',
             color='orange')
    if x_label:
        plt.xlabel(x_label)
    if y_label:
        plt.ylabel(y_label)
    fig.tight_layout()
    plt.savefig(out_pth, dpi=600)


def histogramm_in_grid(df, cols, out_pth, nrow=None,  x_labels=None, y_label=None, titles=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    import seaborn as sns
    import math

    print(f'Plotting histogramms of {cols}\n Saving at {out_pth}')

    matplotlib.rcParams.update({'font.size': 16})

    if not nrow:
        nrow = 1
        y_size = 7
        ncol = len(cols)
    else:
        y_size = nrow * 2
        ncol = math.ceil(len(cols) / nrow)

    x_size = (ncol + 1) * 4

    fig, axs = plt.subplots(nrow, ncol, sharey=True, figsize=[x_size, y_size])

    for i, col in enumerate(cols):
        ix = np.unravel_index(i, axs.shape)

        x = sns.histplot(
            data=df,
            x=col,
            ax=axs[ix],
            facecolor="#a3cc91",
            edgecolor="none"
        )
        max_y_val = x.dataLim.bounds[-1]
        axs[ix].axvline(df[col].mean(), 0, max_y_val, color='black', linewidth=0.5)
        axs[ix].axvline(df[col].max(), 0, max_y_val, color='black', linewidth=0.5)
        bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
        kw = dict(bbox=bbox_props, zorder=2, va="center")
        axs[ix].annotate(f"Mean: {round(df[col].mean(), 2)}", (df[col].mean(), 0.8 * max_y_val), **kw)
        axs[ix].annotate(f"Max: {round(df[col].max(), 2)}", (df[col].max(), 0.8 * max_y_val), **kw)
        # axs[ix].axvline(df[col].median(), 0, max_y_val, color='orange', linewidth=0.5)
        # axs[ix].text(df[col].median(), 0.7 * max_y_val, f"Median: {round(df[col].median(), 1)}", horizontalalignment='left',
        #      color='orange')
        if not x_labels:
            x_label = col
        else:
            x_label = x_labels[i]
        if not y_label:
            y_label = "Count"

        axs[ix].set_ylabel(ylabel=y_label)
        axs[ix].set_xlabel(xlabel=x_label)
        axs[ix].set_xlim(df[col].min(), df[col].max() + 0.1*df[col].max())

        axs[ix].spines['right'].set_visible(False)
        axs[ix].spines['top'].set_visible(False)
        axs[ix].tick_params(labeltop=False, labelright=False)

        if titles:
            title = titles[i]
            axs[ix].set_title(title, size=16, loc="left")

    fig.tight_layout()
    plt.savefig(out_pth, dpi=600)


def plot_kde_plot_comparison(dfs, col, out_pth, labels=None):
    from matplotlib.lines import Line2D
    import matplotlib.pyplot as plt
    import seaborn as sns

    print(f'Plotting kernel density of {col} a.\n\tSaving at {out_pth}')

    fig, axs = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=cm2inch(7, 7))
    linestyles = ['-', '--', '-.', ':']
    legend_elements = [Line2D([0], [0], color='black', linestyle=linestyles[i], label=label) for i, label in enumerate(labels)]

    for i, df in enumerate(dfs):

        sns.kdeplot(data=df, x=col, ax=axs, cut=0, log_scale=False, legend=False, linestyle=linestyles[i],
                    linewidth=1, color='black')
        axs.set_ylabel("")
        axs.set_xlabel(col, size=8)
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)
        axs.legend(legend_elements, labels, frameon=False)

    fig.tight_layout()
    plt.savefig(out_pth, dpi=300)


def plot_map(shp, out_pth, col, vmin=None, vmax=None, shp2_pth=None, extremes=False, markersize=None):
    """
    Plots values of a variable in a map.
    :param shp: Shapefile with the variable as a column/field.
    :param out_pth: Path for map.
    :param col: Column to plot.
    :param vmin: Optional, if provided all values below will be set to vmin. Otherwise the 2 % Quantile will be used.
    :param vmax: Optional, if provided all values above will be set to vmax. Otherwise the 98 % Quantile will be used.
    :param shp2_pth: Optional, plot borders of a second shapefile, e.g. municipalities.
    :param extremes: Optional, highlight all area that are above mean + 2std
    :param markersize: Optional, only use for point geometries. Sets markersize.
    :return:
    """
    import matplotlib.pyplot as plt
    import geopandas as gpd

    if not vmin:
        vmin = shp[col].min()
        vmin = shp[col].quantile(q=0.02)

    if not vmax:
        vmax = shp[col].max()
        vmax = shp[col].quantile(q=0.98)

    print(f'Plotting {col} to {out_pth}')
    fig, ax = plt.subplots(1, 1, figsize=[6, 6])

    tr = shp[col].mean() + (2*shp[col].std())

    if markersize:
        shp.plot(
            column=col,
            ax=ax,
            legend=True,
            legend_kwds={'label': f"{col}", 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
            vmin=vmin,
            vmax=vmax,
            edgecolor='none',
            markersize=markersize
        )
    else:
        shp.plot(
            column=col,
            ax=ax,
            legend=True,
            legend_kwds={'label': f"{col}", 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
            vmin=vmin,
            vmax=vmax,
            edgecolor='none'
        )

    if extremes:
        shp.loc[shp[col] > tr].plot(
            edgecolor='red',
            facecolor="none",
            ax=ax,
            lw=0.5,
            zorder=2
        )

    if shp2_pth:
        shp2 = gpd.read_file(shp2_pth)
        shp2.plot(
            edgecolor='black',
            facecolor="none",
            ax=ax,
            lw=0.1,
            zorder=3
        )

    ax.axis("off")

    ax.margins(0)
    # ax.apply_aspect()
    # bbox = ax.get_window_extent().transformed(transform.inverted(fig.transFigure))
    # w, h = fig.get_size_inches()
    # fig.set_size_inches(w * bbox.width, h * bbox.height)

    plt.tight_layout()
    plt.savefig(out_pth)
    plt.close()


def plot_map_categorical(shp, out_pth, col, colour_dict, shp2_pth=None, markersize=None):
    """
    Plots values of a variable in a map.
    :param shp: Shapefile with the variable as a column/field.
    :param out_pth: Path for map.
    :param col: Column to plot.
    :param shp2_pth: Optional, plot borders of a second shapefile, e.g. municipalities.
    :param markersize: Optional, only use for point geometries. Sets markersize.
    :return:
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import geopandas as gpd
    import math

    shp["color"] = shp[col].map(colour_dict)
    custom_patches = [Patch(facecolor=colour_dict[v], label=v) for v in shp[col].unique()]

    print(f'Plotting {col} to {out_pth}')
    fig, ax = plt.subplots(1, 1, figsize=cm2inch(12, 12))

    if markersize:
        shp.plot(
            color=shp["color"],
            ax=ax,
            legend=False,
            edgecolor='none',
            markersize=markersize
        )
    else:
        shp.plot(
            color=shp["color"],
            ax=ax,
            legend=False,
            edgecolor='none'
        )
    ax.legend(handles=custom_patches, bbox_to_anchor=(1.1, .01), ncol=math.ceil(len(custom_patches)/2)) #,

    if shp2_pth:
        shp2 = gpd.read_file(shp2_pth)
        shp2.plot(
            edgecolor='black',
            facecolor="none",
            ax=ax,
            lw=0.1,
            zorder=3
        )

    ax.axis("off")

    ax.margins(0)
    # ax.apply_aspect()
    # bbox = ax.get_window_extent().transformed(transform.inverted(fig.transFigure))
    # w, h = fig.get_size_inches()
    # fig.set_size_inches(w * bbox.width, h * bbox.height)

    # plt.tight_layout()
    plt.savefig(out_pth)
    plt.close()


def legendgram(f, ax, y, breaks, pal, rescale=False):
    '''
    Add a histogram in a choropleth with colors aligned with map
    ...

    Arguments
    ---------
    f           : Figure
    ax          : AxesSubplot
    y           : ndarray/Series
                  Values to map
    breaks      : list
                  Sequence with breaks for each class (i.e. boundary values
                  for colors)
    rescale     : False/tuple
                  [Optional. Default=False] If a tuple, limits to set the X
                  axis of the histogram
    '''

    import numpy as np

    k = len(breaks)
    p = ax.get_position()
    histpos = [p.x0 + p.width * 0.01, \
               p.y0 + p.height * 0.1, \
               p.width * 0.27, \
               p.height * 0.2]
    histax = f.add_axes(histpos)
    N, bins, patches = histax.hist(y, bins=50, color='0.1')
    # ---
    pl = pal.get_mpl_colormap()
    bucket_breaks = [0] + [np.searchsorted(bins, i) for i in breaks]
    for c in range(k):
        for b in range(bucket_breaks[c], bucket_breaks[c + 1]):
            patches[b].set_facecolor(pl(c / k))
    # ---
    if rescale:
        histax.set_xlim(rescale)
    histax.set_frame_on(False)
    histax.get_yaxis().set_visible(False)
    histax.tick_params(labelsize=12)
    return None


def plot_maps_in_grid(shp, out_pth, cols, cmap=None, nrow=None, vmin=None, vmax=None, shp2_pth=None, labels=None, titles=False, markersize=None):

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    import geopandas as gpd
    import math

    print(f'Plotting {cols} to {out_pth}')

    matplotlib.rcParams.update({'font.size': 14})

    if not nrow:
        nrow = 1
        y_size = 6
        ncol = len(cols)
    else:
        y_size = nrow * 6
        ncol = math.ceil(len(cols)/nrow)

    x_size = (ncol + 1) * 4

    fig, axs = plt.subplots(nrow, ncol, figsize=[x_size, y_size])

    if not cmap:
        cmap = "viridis"

    for i, col in enumerate(cols):
        ix = np.unravel_index(i, axs.shape)

        if not vmin:
            vmin_use = shp[col].quantile(q=0.02)
        else:
            vmin_use = vmin

        if not vmax:
            vmax_use = shp[col].quantile(q=0.98)
        else:
            vmax_use = vmax

        if labels:
            label = labels[i]
        else:
            label = col

        if markersize:
            p = shp.plot(
                column=col,
                ax=axs[ix],
                legend=True,
                legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
                vmin=vmin_use,
                vmax=vmax_use,
                cmap=cmap,
                markersize=markersize
            )
        else:
            p = shp.plot(
                column=col,
                ax=axs[ix],
                legend=True,
                legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
                vmin=vmin_use,
                vmax=vmax_use,
                cmap=cmap
            )

        if titles:
            title = titles[i]
            axs[ix].set_title(title, size=16, x=0.01, y=0.9)
        axs[ix].axis("off")

        if shp2_pth:
            shp2 = gpd.read_file(shp2_pth)
            shp2.plot(edgecolor='black', facecolor="none", ax=axs[ix], lw=0.3, zorder=2)

    # fig.colorbar(np.array([vmin, vmax]), ax=axs[1], shrink=0.8, location='bottom')

    # cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    # sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    # # fake up the array of the scalar mappable. Urgh...
    # sm._A = []
    # fig.colorbar(sm, cax=cax)

    fig.tight_layout()
    plt.savefig(out_pth)
    plt.close()