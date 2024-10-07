# Author:
# github repository:

## ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import os
import time
import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib
import math
import seaborn as sns

import helper_functions
import plotting_lib
## ------------------------------------------ USER INPUT ------------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input:
ALK_IACS_INTERS_PTH = r"09_alkis_intersection_with_other_layers\alkis_iacs_inters.shp"
OWNERS_W_THRESH_LOC_AND_DIST_CLASSIFIED_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-dist+cleaned+loc+class.csv"
COMMUNITY_INFO_FROM_DAFNE_W_THRESH = r"10_owner_network_classification\community_infos_from_dafne_thr{0}.csv"
DISTRICT_SHP_PTH = r"00_data\vector\administrative\BB_districts.shp"

## Output:
OWNER_DF_FOR_PLOTTING = r"14_paper_figures\14_owners_w_parcels+comm_w_thr{0}.csv"
## Other output paths are defined in the functions below in this script.


## ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#

def prepare_df(owner_df_pth, alkis_iacs_intersection_pth, df_info_dafne_pth, threshold, out_pth):
    print("\tRead owner data")
    df_owners = pd.read_csv(owner_df_pth.format(threshold), sep=";")
    print("\tUnique owner classes:", df_owners["new_category"].unique())

    df_owners.loc[df_owners['mcomp_loc'].notna(), "mcomp_dist"] = df_owners.loc[
        df_owners['mcomp_loc'].notna()].apply(
        lambda row: helper_functions.wkt_point_distance(row.parcel_loc, row.mcomp_loc), axis=1)

    print("\tRead alkis-iacs data intersected")
    gdf = gpd.read_file(alkis_iacs_intersection_pth)
    gdf["area"] = gdf["geometry"].area
    gdf["area"] = gdf["area"] / 10000

    df_info = pd.read_csv(df_info_dafne_pth.format(threshold))

    print("\tCombine parcels with owner data")
    df = helper_functions.combine_parcels_with_owners(
        parcels=gdf[["BTNR", "OGC_FID", "area"]],
        owner_df=df_owners[["OGC_FID", "mother_company", f"community_{threshold}", "owner_merge", "level1", "level3",
                            "level_c_category", "new_category", "mcomp_dist", "distance",
                             "agric", 'city_mcomp', 'fstate_mcomp', "netw_max_dist"]])


    df["mcomp_dist"] = df["mcomp_dist"] / 1000
    df["distance"] = df["distance"] / 1000
    print(df["new_category"].unique())
    print("\t", len(df))

    df = pd.merge(df, df_info[[f"community_{threshold}", 'number_people', 'number_companies', 'num_owners_in_alkis']],
                  how="left", on=f"community_{threshold}")

    t = df.loc[df["new_category"].isna()].copy()
    ## Short fix for some missing cases, which for some reason occur. Fix later in  11_1
    helper_functions.print_red(f"Some entries with these level_c_categories are NAN. Setting to noagPR:\n"
                               f"{df.loc[(df['new_category'].isna()), 'level_c_category']}")
    df.loc[(df["new_category"].isna()), "new_category"] = 'noagPR'
    print("\tUnique owner classes:", df_owners["new_category"].unique())

    df.to_csv(out_pth, sep=";", index=False)


def fig_share_and_characteristics_owner_categories_post_revision(owner_df_pth, threshold, out_pth):
    print("Figure share and characteristics of owner categories")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    t = df.loc[df[cat_col].isna()].copy()
    df[cat_col] = df[cat_col].map({"PUBLIC": "Öffentliche Institutionen", "nCONETW": "Nicht landwirtschaftlich\ntätige Unternehmensnetzwerke",
                                   "aCONETW": "Landwirtschaftlich\ntätige Unternehmensnetzwerke", "NONPRO": "Andere Institutionen",
                                   "siCOMP": "Nicht landwirtschaftlich\ntätige Einzelunternehmen", "a_siCOMP": "Landwirtschaftlich\ntätige Einzelunternehmen",
                                   "noagPR": "Nicht landwirtschaftlich\ntätige Privatpersonen", "agriPR": "Landwirtschaftlich\ntätige Privatpersonen",
                                   "CHURCH": "Religiöse Institutionen"})
    # labels = ["PU", "nCN", "OTH", "nSC", "aSC", "nPR", "aPR", "CH"]
    labels = ["Landwirtschaftlich\ntätige Privatpersonen", "Nicht landwirtschaftlich\ntätige Privatpersonen",
              "Landwirtschaftlich\ntätige Einzelunternehmen", "Nicht landwirtschaftlich\ntätige Einzelunternehmen",
              "Landwirtschaftlich\ntätige Unternehmensnetzwerke", "Nicht landwirtschaftlich\ntätige Unternehmensnetzwerke",
              "Öffentliche Institutionen", "Andere Institutionen", "Religiöse Institutionen"]

    print("\tPrepare dfs for plotting")
    ## ger average plot size per owner and owner class
    df_own_parc = df.groupby(["mother_company", "OGC_FID"]).agg(
        area=pd.NamedAgg("area", "sum"),
        category=pd.NamedAgg(cat_col, "first"),
        mcomp_dist=pd.NamedAgg("mcomp_dist", "mean")
    ).reset_index()
    df_own_parc.rename(columns={"category": cat_col}, inplace=True)
    bins = [0, 10, 25, 50, 250, 12000]
    df_own_parc["dist_class"] = pd.cut(df_own_parc["mcomp_dist"],
                                          bins=bins,
                                          labels=["<10", "10 to\n<25", "25 to\n<50", "50 to\n<250", ">250"],
                                          )
    df_own_parc["dist_class"] = df_own_parc["dist_class"].cat.add_categories("unkown").fillna("unkown")

    df_dist_bar = df_own_parc.groupby([cat_col, "dist_class"]).agg(
        area=pd.NamedAgg("area", "sum")
    ).reset_index()
    df_dist_bar[cat_col] = pd.Categorical(df_dist_bar[cat_col], categories=labels, ordered=False)
    df_dist_share = df_dist_bar.groupby("dist_class").agg(
        area=pd.NamedAgg("area", "sum")
    ).reset_index()
    df_dist_share["share"] = round(df_dist_share["area"] / df_dist_share["area"].sum() * 100, 1)

    ## get number of owners and area per owner class
    df_own_agg = df_own_parc.groupby("mother_company").agg(
        area=pd.NamedAgg("area", "sum"),
        category=pd.NamedAgg(cat_col, "first"),
        num_plots=pd.NamedAgg("OGC_FID", "nunique"),
        avg_dist=pd.NamedAgg("mcomp_dist", "mean"),
        min_dist=pd.NamedAgg("mcomp_dist", "min"),
        max_dist=pd.NamedAgg("mcomp_dist", "max")
    ).reset_index()
    df_own_agg.rename(columns={"category": cat_col}, inplace=True)
    df_own_agg["range_dist"] = df_own_agg["max_dist"] - df_own_agg["min_dist"]

    df_own_agg["avg_dist_class"] = pd.cut(df_own_agg["avg_dist"],
                                          bins=[0, 10, 25, 50, 250, df_own_agg["avg_dist"].max()],
                                          labels=["<10", "10 to\n<25", "25 to\n<50", "50 to\n<250", ">250"]
                                          )
    df_own_agg["avg_dist_class"] = df_own_agg["avg_dist_class"].cat.add_categories("unkown").fillna("unkown")

    df_num = df_own_agg.groupby(cat_col).agg(
        number_owners=pd.NamedAgg("mother_company", "count"),
        area=pd.NamedAgg("area", "sum"),
        avg_area=pd.NamedAgg("area", "mean"),
        avg_dist=pd.NamedAgg("avg_dist", "mean"),
        avg_max_dist=pd.NamedAgg("max_dist", "mean"),
    ).reset_index()
    df_num["area_share"] = round((df_num["area"] / df_num["area"].sum()) * 100, 1)

    df_num2 = df_own_parc.groupby(cat_col).agg(
        number_owners=pd.NamedAgg("mother_company", "nunique"),
        area=pd.NamedAgg("area", "sum"),
        avg_plot_area=pd.NamedAgg("area", "mean"),
        avg_dist=pd.NamedAgg("mcomp_dist", "mean"),
        median_dist=pd.NamedAgg("mcomp_dist", "median"),
        num_plots=pd.NamedAgg("OGC_FID", "nunique")
    ).reset_index()
    df_num2["avg_no_plots"] = df_num2["num_plots"] / df_num2["number_owners"]
    df_num2["avg_area"] = df_num2["area"] / df_num2["number_owners"]
    df_num2.to_csv(out_pth[:-4] + '.csv', sep=";", index=False)

    df_num_dist = df_own_agg.groupby([cat_col, "avg_dist_class"]).agg(
        number_owners=pd.NamedAgg("mother_company", "count")
    ).reset_index()

    df_num_dist = pd.pivot(df_num_dist, index=cat_col, columns="avg_dist_class", values="number_owners")
    df_num_dist.loc['ALL'] = df_num_dist.sum()
    df_num_dist = round(df_num_dist.div(df_num_dist.sum(axis=1), axis=0) * 100, 0)
    df_num_dist = df_num_dist.reindex(labels)

    df_area_dist = df_own_agg.groupby([cat_col, "avg_dist_class"]).agg(
        area=pd.NamedAgg("area", "sum")
    ).reset_index()
    df_area_dist = pd.pivot(df_area_dist, index=cat_col, columns="avg_dist_class", values="area")
    df_area_dist.loc['ALL'] = df_area_dist.sum()
    df_area_dist = round(df_area_dist.div(df_area_dist.sum(axis=1), axis=0) * 100, 0)
    df_area_dist = df_area_dist.reindex(labels)

    ## PLOTTING
    print("\tPlotting")

    ############## OWNER MERGE ##############
    num_colors = len(labels)
    color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#c898f5', '#bcbd22', '#f6f739', '#7f7f7f', '#fcb97e']
    # colour_dict = dict(zip(labels, color_list[:num_colors]))
    labels = ["Agricultural\nprivate persons", "Non-agricultural\nprivate persons", "Agricultural\nsingle companies",
              "Non-agricultural\nsingle companies", "Agricultural\ncompany networks", "Non-agricultural\ncompany networks",
              "Public institutions", "Other institutions", "Religious institutions"]
    labels = ["Landwirtschaftlich\ntätige Privatpersonen", "Nicht landwirtschaftlich\ntätige Privatpersonen",
              "Nicht landwirtschaftlich\ntätige Einzelunternehmen", "Landwirtschaftlich\ntätige Einzelunternehmen",
              "Nicht landwirtschaftlich\ntätige Unternehmensnetzwerke",
              "Landwirtschaftlich\ntätige Unternehmensnetzwerke",
              "Öffentliche Institutionen", "Andere Institutionen", "Religiöse Institutionen"]

    colour_dict = {
        "Landwirtschaftlich\ntätige Unternehmensnetzwerke": '#f6f739',
        "Nicht landwirtschaftlich\ntätige Unternehmensnetzwerke": '#bcbd22',
        "Landwirtschaftlich\ntätige Einzelunternehmen": '#c898f5',
        "Nicht landwirtschaftlich\ntätige Einzelunternehmen": '#9467bd',
        "Landwirtschaftlich\ntätige Privatpersonen": '#ff7f0e',
        "Nicht landwirtschaftlich\ntätige Privatpersonen": '#fcb97e',
        "Öffentliche Institutionen": '#1f77b4',
        "Andere Institutionen": '#2ca02c',
        "Religiöse Institutionen": '#7f7f7f'
    }

    matplotlib.rcParams.update({'font.size': 11})

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=plotting_lib.cm2inch(24, 11))
    # gs = GridSpec(1, 2, figure=fig)
    for i, ax in enumerate(fig.axes):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_visible(False)
        # ax.spines['left'].set_visible(False)
        ax.tick_params(labeltop=False, labelright=False)

    ax2 = axs[0]
    ax4 = axs[1]
    ax2.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    ax4.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

    ## Area share
    df_num.sort_values(by="area_share", ascending=False, inplace=True)
    df_num[cat_col] = pd.Categorical(df_num[cat_col], categories=df_num[cat_col].tolist(), ordered=True)
    # ax2.set_title("a)", loc="left")
    ax2.grid(visible=True, which="both", axis="x", zorder=0)
    ax2.set_axisbelow(True)
    ax2.grid(visible=True, which="major", axis="y", zorder=0)
    sns.barplot(y=cat_col, x="area_share", data=df_num, ax=ax2, palette=colour_dict)
    ax2.set_xlabel("Anteil pro Kategorie [%]")
    ax2.xaxis.set_label_position('top')
    ax2.set_ylabel(None)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)

    ## Mean area per owner
    # ax4.set_title("b)", loc="left")
    ax4.grid(visible=True, which="both", axis="x", zorder=0)
    ax4.set_axisbelow(True)
    sns.barplot(y=cat_col, x="avg_area", data=df_num, ax=ax4, palette=colour_dict)
    ax4.set_xlabel("Mittlere Fläche pro Kategorie [ha]")
    ax4.xaxis.set_label_position('top')
    ax4.set_ylabel(None)
    ax4.get_yaxis().set_ticks([])
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_pth[:-4] + '_a.png')
    plt.close()



def table_largest_owner_examples(owner_df_pth, threshold, out_pth):
    print("Do something!")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    t = df.loc[df[cat_col].isna()].copy()
    df[cat_col] = df[cat_col].map({"PUBLIC": "PU", "nCONETW": "nCN", "aCONETW": "aCN", "NONPRO": "OTH", "siCOMP": "nSC", "a_siCOMP": "aSC",
                                           "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    # labels = ["PU", "nCN", "OTH", "nSC", "aSC", "nPR", "aPR", "CH"]
    labels = ["aPR", "nPR", "aSC", "nSC", "aCN", "nCN", "PU", "OTH", "RE"]

    print("\tPrepare dfs for plotting")
    ## ger average plot size per owner and owner class

    df_mcomp_parc = df.groupby(["mother_company", "OGC_FID"]).agg(
        area=pd.NamedAgg("area", "sum"),
        category=pd.NamedAgg(cat_col, "first"),
        number_people=pd.NamedAgg("number_people", "first"),
        number_companies=pd.NamedAgg("number_companies", "first"),
        num_owners_in_alkis=pd.NamedAgg("num_owners_in_alkis", "first"),
        num_hierarch_levels=pd.NamedAgg("netw_max_dist", "first"),
        mcomp_dist=pd.NamedAgg("mcomp_dist", "mean")
    ).reset_index()
    df_mcomp_parc.rename(columns={"category": cat_col}, inplace=True)

    ## get number of owners and area per owner class
    df_mcomp_agg = df_mcomp_parc.groupby("mother_company").agg(
        area=pd.NamedAgg("area", "sum"),
        category=pd.NamedAgg(cat_col, "first"),
        num_plots=pd.NamedAgg("OGC_FID", "nunique"),
        avg_dist=pd.NamedAgg("mcomp_dist", "mean"),
        min_dist=pd.NamedAgg("mcomp_dist", "min"),
        max_dist=pd.NamedAgg("mcomp_dist", "max"),
        number_people=pd.NamedAgg("number_people", "first"),
        number_companies=pd.NamedAgg("number_companies", "first"),
        num_owners_in_alkis=pd.NamedAgg("num_owners_in_alkis", "first"),
        num_hierarch_levels=pd.NamedAgg("num_hierarch_levels", "first")
    ).reset_index()
    df_mcomp_agg.rename(columns={"category": cat_col}, inplace=True)
    df_mcomp_agg["range_dist"] = df_mcomp_agg["max_dist"] - df_mcomp_agg["min_dist"]
    df_mcomp_agg["avg_dist_class"] = pd.cut(df_mcomp_agg["avg_dist"],
                                          bins=[0, 10, 25, 50, 250, df_mcomp_agg["avg_dist"].max()],
                                          labels=["<10", "10 to <25", "25 to <50", "50 to <250", ">250"])
    df_mcomp_agg["avg_dist_class"] = df_mcomp_agg["avg_dist_class"].cat.add_categories("unkown").fillna("unkown")
    df_mcomp_agg["share"] = round(df_mcomp_agg["area"] / df_mcomp_agg["area"].sum() * 100, 4)
    df_mcomp_agg = df_mcomp_agg[["new_category", "area", "share", "avg_dist_class", "avg_dist", "mother_company"]]
    df_mcomp_agg.sort_values(by="area", ascending=False, inplace=True)
    df_mcomp_agg["area"] = round(df_mcomp_agg["area"], 0)
    df_mcomp_agg["share"] = round(df_mcomp_agg["share"], 1)
    df_mcomp_agg["avg_dist"] = round(df_mcomp_agg["avg_dist"], 0)
    df_mcomp_agg.columns = ["Owner class", "Area [ha]", "Share land [%]", "Distance class", "Average distance", "GUO"]
    df_mcomp_agg = df_mcomp_agg.loc[df_mcomp_agg["GUO"] != "unbekannt"].copy()

    df_mcomp_agg = df_mcomp_agg[:100]
    df_mcomp_agg.index = range(1, len(df_mcomp_agg)+1)

    df_mcomp_agg.to_csv(out_pth, sep=";")


def fig_comparison_map_and_histograms_concentration_measures_post_revision(df_res_pth, threshold, district_shp_pth, out_pth, area_threshold=7000):

    print("Plot maps and histogramms of concentration measures.")
    grid_res = 4
    version_4km = 1
    gdf_grid = gpd.read_file(rf"00_data\vector\grids\square_grid_4km_v01_with_12km_POLYIDs_BB.gpkg")
    #rf"00_data\vector\grids\square_grid_{grid_res}km_v{version_4km:02d}_with_12km_POLYIDs.shp"

    df_res = pd.read_csv(df_res_pth.format(threshold))

    shp = pd.merge(gdf_grid, df_res, how='inner', left_on='POLYID', right_on='id_sp_unit')
    # shp = shp.loc[shp['gini_coeff'] != 0.0]
    shp = shp.loc[shp['total_area'] > area_threshold].copy()

    ## For QGIS projects
    # shp.to_file(r"14_paper_figures\vector\grid_with_values.gpkg", driver="GPKG")

    # shp["total_area"] = shp["total_area"] / 10000

    shp2 = gpd.read_file(district_shp_pth)

    matplotlib.rcParams.update({'font.size': 12})

    fig = plt.figure(constrained_layout=True, figsize=(12, 5.5))
    widths = [1, 1, 1]
    heights = [4, 1]
    spec = fig.add_gridspec(
        ncols=3, nrows=2, width_ratios=widths,
        height_ratios=heights)
    cmap = "crest"

    col = "cr1"
    label = "CR1 [%]"
    title = "a)"
    ax = fig.add_subplot(spec[0, 0])
    vmin_use = shp[col].min() #shp[col].quantile(q=0.02)
    vmax_use = shp[col].quantile(q=0.98) # shp[col].max() #41 #shp[col].quantile(q=0.999)
    shp.plot(
        column=col,
        ax=ax,
        legend=True,
        legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.5, 1.5), "pad": 0.01},
        vmin=vmin_use,
        vmax=vmax_use,
        cmap=cmap
    )
    # shp.loc[shp[col] > 40].plot(
    #     column=col,
    #     ax=ax,
    #     legend=False,
    #     edgecolor='orange',
    #     facecolor="none",
    #     lw=0.3, zorder=2
    # )
    # ax.set_title(title, size=16, x=0.01, y=0.9)
    ax.axis("off")
    shp2.plot(edgecolor='black', facecolor="none", ax=ax, lw=0.3, zorder=2)

    ax_sub = fig.add_subplot(spec[1, 0])
    x = sns.histplot(
        data=shp,
        x=col,
        ax=ax_sub,
        facecolor="#a3cc91",
        edgecolor="none"
    )
    max_y_val = x.dataLim.bounds[-1]
    ax_sub.axvline(shp[col].mean(), 0, max_y_val, color='black', linewidth=0.5)
    bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(bbox=bbox_props, zorder=2, va="center")
    ax_sub.annotate(f"Mean: {round(shp[col].mean(), 2)}", (shp[col].mean(), 0.8 * max_y_val), **kw)
    ax_sub.set_ylabel(ylabel="Anzahl\nSubregionen")
    ax_sub.set_ylim(0, 220)
    ax_sub.set_xlabel(None) #(xlabel=label)
    ax_sub.spines['right'].set_visible(False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.tick_params(labeltop=False, labelright=False)


    col = "hhi"
    label = "HHI"
    title = "c)"
    ax = fig.add_subplot(spec[0, 1])
    vmin_use = shp[col].min() #shp[col].quantile(q=0.02)
    vmax_use = shp[col].quantile(q=0.98) # shp[col].max() #1100
    shp.plot(
        column=col,
        ax=ax,
        legend=True,
        legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.5, 1.5), "pad": 0.01},
        vmin=vmin_use,
        vmax=vmax_use,
        cmap=cmap
    )
    # shp.loc[shp[col] > 1000].plot(
    #     column=col,
    #     ax=ax,
    #     legend=False,
    #     edgecolor='orange',
    #     facecolor="none",
    #     lw=0.3, zorder=2
    # )
    # ax.set_title(title, size=16, x=0.01, y=0.9)
    ax.axis("off")
    shp2.plot(edgecolor='black', facecolor="none", ax=ax, lw=0.3, zorder=2)

    ax_sub = fig.add_subplot(spec[1, 1])
    x = sns.histplot(
        data=shp,
        x=col,
        ax=ax_sub,
        facecolor="#a3cc91",
        edgecolor="none"
    )
    max_y_val = x.dataLim.bounds[-1]
    ax_sub.axvline(shp[col].mean(), 0, max_y_val, color='black', linewidth=0.5)
    bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(bbox=bbox_props, zorder=2, va="center")
    ax_sub.annotate(f"Mean: {round(shp[col].mean(), 2)}", (shp[col].mean(), 0.8 * max_y_val), **kw)
    # ax_sub.set_ylabel(ylabel="No. polygons")
    ax_sub.set_ylabel(ylabel=None)
    ax_sub.set_ylim(0,220)
    ax_sub.set_xlabel(None) #(xlabel=label)
    ax_sub.spines['right'].set_visible(False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.tick_params(labeltop=False, labelright=False)

    col = "gini_coeff"
    label = "Gini"
    title = "d)"
    ax = fig.add_subplot(spec[0, 2])
    vmin_use = shp[col].min() #shp[col].quantile(q=0.02)
    vmax_use = shp[col].quantile(q=0.98) #shp[col].max()
    shp.plot(
        column=col,
        ax=ax,
        legend=True,
        legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.5, 1.5), "pad": 0.01},
        vmin=vmin_use,
        vmax=vmax_use,
        cmap=cmap
    )
    # ax.set_title(title, size=16, x=0.01, y=0.9)
    ax.axis("off")
    shp2.plot(edgecolor='black', facecolor="none", ax=ax, lw=0.3, zorder=2)

    ax_sub = fig.add_subplot(spec[1, 2])
    x = sns.histplot(
        data=shp,
        x=col,
        ax=ax_sub,
        facecolor="#a3cc91",
        edgecolor="none"
    )
    max_y_val = x.dataLim.bounds[-1]
    ax_sub.axvline(shp[col].mean(), 0, max_y_val, color='black', linewidth=0.5)
    bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(bbox=bbox_props, zorder=2, va="center")
    ax_sub.annotate(f"Mean: {round(shp[col].mean(), 2)}", (shp[col].mean(), 0.8 * max_y_val), **kw)
    ax_sub.set_ylabel(ylabel=None)
    ax_sub.set_ylim(0,220)
    ax_sub.set_xlabel(None) #(xlabel=label)
    ax_sub.spines['right'].set_visible(False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.tick_params(labeltop=False, labelright=False)

    # fig.tight_layout()
    plt.savefig(out_pth)
    plt.close()


def calculate_radius_agricultural_owners(alkis_iacs_inters_pth, owners_pth, comm_col, owner_col, out_pth):
    print("Combine parcels with owner data")
    gdf = gpd.read_file(alkis_iacs_inters_pth)
    gdf["area"] = gdf["geometry"].area
    gdf['area'] = gdf['area'] / 10000
    df = helper_functions.read_table_to_df(owners_pth)

    df_alk = helper_functions.combine_parcels_with_owners(
        gdf[["OGC_FID", "EIGENTUEME", "BTNR", "area", "geometry"]],
        df[["OGC_FID", "mother_company", "owner_merge", "level_c_category", "agric", comm_col]])

    df_alk = gpd.GeoDataFrame(df_alk)
    df_alk["centroid"] = df_alk.geometry.centroid
    df_alk["x"] = df_alk["centroid"].x
    df_alk["y"] = df_alk["centroid"].y

    ## subset to agricultural owners
    df_alk = df_alk.loc[df_alk["agric"] == 1].copy()

    ## Get owner centroids
    owners = df_alk.groupby(owner_col)[["x", "y"]].mean().reset_index()
    owners["geometry"] = gpd.points_from_xy(owners['x'], owners['y'])

    buffer_lst = []

    for f, owner_id in enumerate(list(owners[owner_col])):
        print(f+1, len(owners))
        df_alk_sub = df_alk.loc[df_alk[owner_col] == owner_id].copy()
        minx, miny, maxx, maxy = df_alk_sub.geometry.total_bounds
        owner_centroid = owners.loc[owners[owner_col] == owner_id].copy()
        owner_centroid.crs = df_alk.crs

        def coord_lister(geom):
            if geom.geom_type == 'MultiPolygon':
                coords = []
                for x in geom.geoms:
                    coords += list(x.exterior.coords)
            else:
                coords = list(geom.exterior.coords)
            return (coords)

        coords = []
        for row in df_alk_sub.itertuples():
           coords += coord_lister(row.geometry)

        coords_df = gpd.points_from_xy([i[0] for i in coords], [i[1] for i in coords])
        buffer_size = max([owner_centroid.distance(point).iloc[0] for point in coords_df])

        owner = owner_centroid.copy()
        owner.geometry = owner.buffer(buffer_size)

        buffer_lst.append(buffer_size)

    owners["buffer_radius"] = buffer_lst

    owners.csr = df_alk.crs
    owners.to_file(out_pth, driver="GPKG")


## ------------------------------------------ RUN PROCESSES ---------------------------------------------------#
def main():
    s_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + s_time)
    os.chdir(WD)

    threshold = 50

    helper_functions.create_folder(r"14_paper_figures\figures")
    helper_functions.create_folder(r"14_paper_figures\tables")

    # prepare_df(
    #     owner_df_pth=OWNERS_W_THRESH_LOC_AND_DIST_CLASSIFIED_PTH,
    #     alkis_iacs_intersection_pth=ALK_IACS_INTERS_PTH,
    #     df_info_dafne_pth=COMMUNITY_INFO_FROM_DAFNE_W_THRESH,
    #     threshold=50,
    #     out_pth=OWNER_DF_FOR_PLOTTING.format(50)
    # )

    # fig_share_and_characteristics_owner_categories_post_revision(
    #     owner_df_pth=OWNER_DF_FOR_PLOTTING,
    #     threshold=50,
    #     out_pth=rf"16_other_figures\policy_brief\fig03_share_and_characteristics_owner_categories.png")
    #
    # ## 12km concentrations
    # fig_comparison_map_and_histograms_concentration_measures_post_revision(
    #     threshold=50,
    #     df_res_pth=rf"11_ownership_concentration\mw_grid_buffers\mw_conc_meas-mother_companies-comm_w_thr{threshold}-iacs_areas_12km.csv",
    #     district_shp_pth=DISTRICT_SHP_PTH,
    #     out_pth=r"16_other_figures\policy_brief\fig04_comparison_histograms_concentrations_measures_12km.png",
    #     area_threshold=200)

    ## Analyse radi of owners
    threshold = 50
    comm_col = f"community_{threshold}"
    calculate_radius_agricultural_owners(
        alkis_iacs_inters_pth=ALK_IACS_INTERS_PTH,
        owners_pth=OWNERS_W_THRESH_LOC_AND_DIST_CLASSIFIED_PTH.format(threshold),
        comm_col=comm_col,
        owner_col="owner_merge",
        out_pth=r"16_other_figures\owner_centroids_max_radius.gpkg")

    e_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + s_time)
    print("end: " + e_time)


if __name__ == '__main__':
    main()
