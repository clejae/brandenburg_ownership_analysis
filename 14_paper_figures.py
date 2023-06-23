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


def fig_share_and_characteristics_owner_categories(owner_df_pth, threshold, out_pth):
    print("Figure share and characteristics of owner categories")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    t = df.loc[df[cat_col].isna()].copy()
    df[cat_col] = df[cat_col].map({"PUBLIC": "PU", "nCONETW": "CN", "aCONETW": "AH", "NONPRO": "NP", "siCOMP": "nSC", "a_siCOMP": "aSC",
                                           "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    # labels = ["PU", "CN", "NP", "nSC", "aSC", "nPR", "aPR", "CH"]
    labels = ["aPR", "nPR", "aSC", "nSC", "AH", "CN", "PU", "NP", "RE"]

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
    colour_dict = {
        "aPR": '#f6f739',
        "nPR": '#bcbd22',
        "aSC": '#c898f5',
        "nSC": '#9467bd',
        "CN": '#ff7f0e',
        "AH": '#fcb97e',
        "PU": '#1f77b4',
        "NP": '#2ca02c',
        "RE": '#7f7f7f'
    }

    matplotlib.rcParams.update({'font.size': 11})

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=plotting_lib.cm2inch(24, 8))
    # gs = GridSpec(1, 2, figure=fig)
    for i, ax in enumerate(fig.axes):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_visible(False)
        # ax.spines['left'].set_visible(False)
        ax.tick_params(labeltop=False, labelright=False)

    ax2 = axs[0]
    ax4 = axs[1]

    ## Area share
    ax2.set_title("a)", loc="left")
    ## sort rows in manual order
    df_mapping = pd.DataFrame({cat_col: ["aPR", "CN", "AH", "nSC", "PU", "NP", "nPR", "aSC", "RE"],})
    sort_mapping = df_mapping.reset_index().set_index(cat_col)
    df_num['cat_num'] = df_num[cat_col].map(sort_mapping['index'])
    df_num.sort_values('cat_num', inplace=True)
    # df_num[cat_col] = pd.Categorical(df_num[cat_col], categories=["CN", "AH", "PU", "RE", "nPR", "aPR", "aSC", "nSC", "NP"], ordered=True)
    # df_num.sort_values(by="area_share", inplace=True)
    wedges, texts = ax2.pie(df_num["area"], colors=[colour_dict[key] for key in df_num[cat_col]],
                            radius=1, wedgeprops=dict(width=0.3, edgecolor='w', linewidth=.01), startangle=50)
    bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        text = f"{df_num['area_share'].iloc[i]} %"
        ax2.annotate(text, xy=(x, y), xytext=(1.1 * np.sign(x), 1.1 * y), horizontalalignment=horizontalalignment,
                     **kw)
    legend_elements = [Patch(facecolor=colour_dict[key], edgecolor=None, label=key) for key in colour_dict]
    ax2.legend(handles=legend_elements, bbox_to_anchor=(1.2, .8), ncol=2, title="Owner class")

    # no. plots per owner
    df_num2.sort_values(by="avg_area", ascending=False, inplace=True)
    df_num2[cat_col] = pd.Categorical(df_num2[cat_col], categories=df_num2[cat_col].tolist(), ordered=True)
    ax4.set_title("b)", loc="left")
    ax4.set_axisbelow(True)
    ax4.grid(visible=True, which="major", axis="y", zorder=0)
    # sns.scatterplot(x="avg_no_plots", y="avg_area", data=df_num2, hue=cat_col, palette=colour_dict, legend=False)
    sns.barplot(x=cat_col, y="avg_area", data=df_num2, ax=ax4, palette=colour_dict)
    ax4.set_ylabel("Mean area per owner [ha]")
    ax4.set_xlabel(None)
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_pth[:-4] + '_a.png')
    plt.close()

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=plotting_lib.cm2inch(25, 8))
    # gs = GridSpec(1, 2, figure=fig)
    for i, ax in enumerate(fig.axes):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_visible(False)
        # ax.spines['left'].set_visible(False)
        ax.tick_params(labeltop=False, labelright=False)

    ax1 = axs[0]
    ax5 = axs[1]

    ## Distance (x) + share area (y)
    ax1.set_title("c)", loc="left")
    ax1.grid(visible=True, which="major", axis="y", zorder=0)
    ax1.set_axisbelow(True)
    df_dist_share['area'] = df_dist_share['area'] / 1000
    sns.barplot(x="dist_class", y='area', data=df_dist_share,  ax=ax1, edgecolor='none', facecolor='#a3cc91')
    ax1.legend([], [], frameon=False)
    bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(bbox=bbox_props, zorder=1, va="center")
    for i, dist_class in enumerate(df_dist_share["dist_class"].tolist()):
        y = df_dist_share.loc[df_dist_share["dist_class"] == dist_class, "area"].iloc[0]
        text = f"{df_dist_share.loc[df_dist_share['dist_class'] == dist_class, 'share'].iloc[0]} %"
        ax1.annotate(text, (i, y+50), (i-0.3, y+50), **kw)
    ax1.set_ylabel("Area [1,000 ha]")
    ax1.set_xlabel("Distance to parcel [km]")
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.set_ylim(0, 800)

    ## share of area per group and distance class
    ax5.set_title("d)", loc="left")
    sns.heatmap(df_area_dist, annot=True, ax=ax5, cmap="crest", cbar=False)
    ax5.set_xlabel("Distance to parcel [km]")
    ax5.set_ylabel(None)
    ax5.tick_params(axis='y', labelrotation=0)
    ax5.tick_params(axis='x', labelrotation=0)
    for i in range(len(labels)+1):
        ax5.axhline(y=i, color='w', linewidth=1)

    # plt.subplots_adjust(hspace=0.1, wspace=0.4)
    plt.tight_layout()
    plt.savefig(out_pth[:-4] + '_b.png')
    plt.close()

    # ## Entwurf 2
    # df_plt = pd.concat([df_own_agg[[cat_col, "num_plots", "area", "max_dist"]], df_own_agg[[cat_col, "num_plots", "area", "max_dist"]]])
    # df_plt[cat_col] = df_own_agg[cat_col].tolist() + ["ALL" for x in df_own_agg[cat_col]]
    #
    # fig, axs = plt.subplots(nrows=2, ncols=3, figsize=plotting_lib.cm2inch(32, 16))
    # # area [%]
    # ix = (0, 0)
    # axs[ix].set_title("a)", loc="left")
    # axs[ix].pie(df_num["area"], colors=[colour_dict[key] for key in df_num[cat_col]],
    #              radius=1, wedgeprops=dict(width=0.3, edgecolor='w', linewidth=.01), autopct='%1.0f%%',
    #              startangle=90)
    #
    # legend_elements = [Patch(facecolor=colour_dict[key], edgecolor=None, label=key) for key in colour_dict]
    # axs[ix].legend(handles=legend_elements, bbox_to_anchor=(.7, .7), ncol=2, title="Owner class")
    #
    # # no. plots per owner
    # ix = (0, 1)
    # axs[ix].set_title("b)", loc="left")
    # sns.boxplot(x=cat_col, y="num_plots", data=df_plt, ax=axs[ix], order=labels,
    #                 showfliers=False, palette=colour_dict)
    # axs[ix].set_ylabel("No. parcels")
    # axs[ix].set_xlabel(None)
    # axs[ix].spines['right'].set_visible(False)
    # axs[ix].spines['top'].set_visible(False)
    #
    # # area per owner [ha]
    # ix = (0, 2)
    # axs[ix].set_title("c)", loc="left")
    # sns.boxplot(x=cat_col, y="area", data=df_plt, ax=axs[ix], order=labels,
    #                 showfliers=False, palette=colour_dict)
    # axs[ix].set_ylabel("Area [ha]")
    # axs[ix].set_xlabel(None)
    # axs[ix].spines['right'].set_visible(False)
    # axs[ix].spines['top'].set_visible(False)
    #
    # # no. owner [%] vs. distance
    # ix = (1, 0)
    # axs[ix].set_title("d)", loc="left")
    # sns.heatmap(df_num_dist, annot=True, ax=axs[ix], cmap="crest", cbar=False)
    # axs[ix].set_xlabel("Distance [km]")
    # axs[ix].set_ylabel(None)
    # axs[ix].tick_params(axis='y', labelrotation=0)
    # axs[ix].tick_params(axis='x', labelrotation=0)
    # for i in range(6):
    #     axs[ix].axhline(y=i, color='w', linewidth=1)
    #
    # # area [%] vs. distance
    # ix = (1, 1)
    # axs[ix].set_title("e)", loc="left")
    # sns.heatmap(df_area_dist, annot=True, ax=axs[ix], cmap="crest", cbar=False)
    # axs[ix].set_xlabel("Distance [km]")
    # axs[ix].set_ylabel(None)
    # axs[ix].tick_params(axis='y', labelrotation=0)
    # axs[ix].tick_params(axis='x', labelrotation=0)
    # for i in range(6):
    #     axs[ix].axhline(y=i, color='w', linewidth=1)
    #
    # # area [%] vs. distance
    # ix = (1, 2)
    # axs[ix].set_title("f)", loc="left")
    # sns.boxplot(x="max_dist", y=cat_col, data=df_plt, ax=axs[ix], order=labels,
    #                 showfliers=False, palette=colour_dict)
    # axs[ix].set_xlabel("Max. distance [km]")
    # axs[ix].spines['right'].set_visible(False)
    # axs[ix].spines['top'].set_visible(False)
    # axs[ix].set_ylabel(None)
    #
    # fig.tight_layout()
    # plt.savefig(rf"12_area_calculations\figures\fig_share_and_characteristics_owner_categories{descr}.png")
    # plt.close()


def fig_appendix_change_in_distances_aggregation_levels(owner_df_pth, threshold, out_pth):
    print("Plot appendix figure with change in distances between the aggregation levels.")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    t = df.loc[df[cat_col].isna()].copy()
    df[cat_col] = df[cat_col].map({"PUBLIC": "PU", "nCONETW": "CN", "aCONETW": "AH", "NONPRO": "NP", "siCOMP": "nSC", "a_siCOMP": "aSC",
                                           "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    # labels = ["PU", "CN", "NP", "nSC", "aSC", "nPR", "aPR", "CH"]
    labels = ["aPR", "nPR", "aSC", "nSC", "AH", "CN", "PU", "NP", "RE"]

    df = df.loc[df[cat_col].isin(["AH", "CN"])]


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

    df_own_parc = df.groupby(["owner_merge", "OGC_FID"]).agg(
        area=pd.NamedAgg("area", "sum"),
        category=pd.NamedAgg(cat_col, "first"),
        number_people=pd.NamedAgg("number_people", "first"),
        number_companies=pd.NamedAgg("number_companies", "first"),
        num_owners_in_alkis=pd.NamedAgg("num_owners_in_alkis", "first"),
        num_hierarch_levels=pd.NamedAgg("netw_max_dist", "first"),
        own_dist=pd.NamedAgg("distance", "mean")
    ).reset_index()
    df_own_parc.rename(columns={"category": cat_col}, inplace=True)

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
                                          labels=["<10", "10 to\n<25", "25 to\n<50", "50 to\n<250", ">250"]
                                          )
    df_mcomp_agg["avg_dist_class"] = df_mcomp_agg["avg_dist_class"].cat.add_categories("unkown").fillna("unkown")

    df_own_agg = df_own_parc.groupby("owner_merge").agg(
        area=pd.NamedAgg("area", "sum"),
        category=pd.NamedAgg(cat_col, "first"),
        num_plots=pd.NamedAgg("OGC_FID", "nunique"),
        avg_dist=pd.NamedAgg("own_dist", "mean"),
        min_dist=pd.NamedAgg("own_dist", "min"),
        max_dist=pd.NamedAgg("own_dist", "max")
    ).reset_index()
    df_own_agg.rename(columns={"category": cat_col}, inplace=True)
    df_own_agg["range_dist"] = df_own_agg["max_dist"] - df_own_agg["min_dist"]
    df_own_agg["avg_dist_class"] = pd.cut(df_own_agg["avg_dist"],
                                          bins=[0, 10, 25, 50, 250, df_own_agg["avg_dist"].max()],
                                          labels=["<10", "10 to\n<25", "25 to\n<50", "50 to\n<250", ">250"]
                                          )
    df_own_agg["avg_dist_class"] = df_own_agg["avg_dist_class"].cat.add_categories("unkown").fillna("unkown")

    df_area_dist_own = df_own_agg.groupby(["avg_dist_class"]).agg(
        area=pd.NamedAgg("area", "sum")
    ).reset_index()
    df_area_dist_own["category"] = "Aggr1"
    df_area_dist_own["area"] = round(df_area_dist_own["area"] / df_area_dist_own["area"].sum() * 100, 0)

    df_area_dist_mcomp = df_mcomp_agg.groupby(["avg_dist_class"]).agg(
        area=pd.NamedAgg("area", "sum")
    ).reset_index()
    df_area_dist_mcomp["category"] = "Aggr2"
    df_area_dist_mcomp["area"] = round(df_area_dist_mcomp["area"] / df_area_dist_mcomp["area"].sum() * 100, 0)

    df_area_dist = pd.concat([df_area_dist_own, df_area_dist_mcomp], axis=0)
    df_area_dist = pd.pivot(df_area_dist, index="category", columns="avg_dist_class", values="area")

    ## PLOTTING
    print("\tPlotting")

    matplotlib.rcParams.update({'font.size': 11})

    ############## OWNER MERGE ##############
    fig, ax = plt.subplots(figsize=plotting_lib.cm2inch(12.5, 8))
    ## share of area per group and distance class
    sns.heatmap(df_area_dist, ax = ax, annot=True,  cmap="crest", cbar=False)
    ax.set_xlabel("Distance to parcel [km]")
    ax.set_ylabel(None)
    ax.tick_params(axis='y', labelrotation=0)
    ax.tick_params(axis='x', labelrotation=0)
    for i in range(len(labels) + 1):
        ax.axhline(y=i, color='w', linewidth=1)

    # plt.subplots_adjust(hspace=0.1, wspace=0.4)
    plt.tight_layout()
    plt.savefig(out_pth)
    plt.close()

    print("done!")


def table_largest_owner_examples(owner_df_pth, threshold, out_pth):
    print("Do something!")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    t = df.loc[df[cat_col].isna()].copy()
    df[cat_col] = df[cat_col].map({"PUBLIC": "PU", "nCONETW": "CN", "aCONETW": "AH", "NONPRO": "NP", "siCOMP": "nSC", "a_siCOMP": "aSC",
                                           "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    # labels = ["PU", "CN", "NP", "nSC", "aSC", "nPR", "aPR", "CH"]
    labels = ["aPR", "nPR", "aSC", "nSC", "AH", "CN", "PU", "NP", "RE"]

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
                                          labels=["<10", "10 to <25", "25 to <50", "50 to <250", ">250"]
                                          )
    df_mcomp_agg["avg_dist_class"] = df_mcomp_agg["avg_dist_class"].cat.add_categories("unkown").fillna("unkown")
    df_mcomp_agg["share"] = round(df_mcomp_agg["area"] / df_mcomp_agg["area"].sum() * 100, 4)
    df_mcomp_agg = df_mcomp_agg[[ "new_category", "area", "share", "avg_dist_class", "mother_company"]]
    df_mcomp_agg.sort_values(by="area", ascending=False, inplace=True)
    df_mcomp_agg["area"] = round(df_mcomp_agg["area"], 0)
    df_mcomp_agg["share"] = round(df_mcomp_agg["share"], 1)
    df_mcomp_agg.columns = ["Owner class", "Area [ha]", "Share land [%]", "Distance class", "GUO"]
    df_mcomp_agg = df_mcomp_agg.loc[df_mcomp_agg["GUO"] != "unbekannt"].copy()

    df_mcomp_agg = df_mcomp_agg[:100]
    df_mcomp_agg.index = range(1, len(df_mcomp_agg)+1)

    df_mcomp_agg.to_csv(out_pth, sep=";")


def fig_appendix_characteristics_company_networks(owner_df_pth, threshold, out_pth):
    print("Plot appendix figure of the chararcteristics of the company networks.")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    t = df.loc[df[cat_col].isna()].copy()
    df[cat_col] = df[cat_col].map({"PUBLIC": "PU", "nCONETW": "CN", "aCONETW": "AH", "NONPRO": "NP", "siCOMP": "nSC", "a_siCOMP": "aSC",
                                           "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    # labels = ["PU", "CN", "NP", "nSC", "aSC", "nPR", "aPR", "CH"]
    labels = ["aPR", "nPR", "aSC", "nSC", "AH", "CN", "PU", "NP", "RE"]

    df = df.loc[df[cat_col].isin(["AH", "CN"])]

    ## ToDo: Compare mean area and distance per owner_merge and mother_company

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

    df_own_parc = df.groupby(["owner_merge", "OGC_FID"]).agg(
        area=pd.NamedAgg("area", "sum"),
        category=pd.NamedAgg(cat_col, "first"),
        number_people=pd.NamedAgg("number_people", "first"),
        number_companies=pd.NamedAgg("number_companies", "first"),
        num_owners_in_alkis=pd.NamedAgg("num_owners_in_alkis", "first"),
        num_hierarch_levels=pd.NamedAgg("netw_max_dist", "first"),
        own_dist=pd.NamedAgg("distance", "mean")
    ).reset_index()
    df_own_parc.rename(columns={"category": cat_col}, inplace=True)

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

    df_own_agg = df_own_parc.groupby("owner_merge").agg(
        area=pd.NamedAgg("area", "sum"),
        category=pd.NamedAgg(cat_col, "first"),
        num_plots=pd.NamedAgg("OGC_FID", "nunique"),
        avg_dist=pd.NamedAgg("own_dist", "mean"),
        min_dist=pd.NamedAgg("own_dist", "min"),
        max_dist=pd.NamedAgg("own_dist", "max")
    ).reset_index()
    df_own_agg.rename(columns={"category": cat_col}, inplace=True)
    df_own_agg["range_dist"] = df_own_agg["max_dist"] - df_own_agg["min_dist"]

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.precision', 1)

    descr_mcomps = df_mcomp_agg.agg(
        {
            "area": ["count", "min", "max", "median", "mean", "sum", "skew"],
            "avg_dist": ["count", "min", "max", "median", "mean", "sum", "skew"],
            "max_dist": ["count", "min", "max", "median", "mean", "sum", "skew"],
            "range_dist": ["count", "min", "max", "median", "mean", "sum", "skew"],
            "num_owners_in_alkis": ["count", "min", "max", "median", "mean", "sum", "skew"]
        }
    )
    descr_owners = df_own_agg.agg(
        {
            "area": ["count", "min", "max", "median", "mean", "sum", "skew"],
            "avg_dist": ["count", "min", "max", "median", "mean", "sum", "skew"],
            "max_dist": ["count", "min", "max", "median", "mean", "sum", "skew"],
            "range_dist": ["count", "min", "max", "median", "mean", "sum", "skew"]
        }
    )


    ## PLOTTING
    print("\tPlotting")

    matplotlib.rcParams.update({'font.size': 11})

    ############## OWNER MERGE ##############
    num_colors = len(labels)
    color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#c898f5', '#bcbd22', '#f6f739', '#7f7f7f', '#fcb97e']
    # colour_dict = dict(zip(labels, color_list[:num_colors]))
    colour_dict = {
        "aPR": '#f6f739',
        "nPR": '#bcbd22',
        "aSC": '#c898f5',
        "nSC": '#9467bd',
        "CN": '#ff7f0e',
        "AH": '#fcb97e',
        "PU": '#1f77b4',
        "NP": '#2ca02c',
        "RE": '#7f7f7f'
    }

    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=plotting_lib.cm2inch(25, 8))
    # gs = GridSpec(1, 2, figure=fig)
    for i, ax in enumerate(fig.axes):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_visible(False)
        # ax.spines['left'].set_visible(False)
        ax.tick_params(labeltop=False, labelright=False)

    ax1 = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]
    ax4 = axs[3]

    sns.boxplot(
        data=df_mcomp_agg,
        x=cat_col,
        y="num_owners_in_alkis",
        fliersize=1,
        palette=colour_dict,
        linewidth=.5,
        ax=ax1
    )
    ax1.grid(visible=True, which="major", axis="y", zorder=0)
    ax1.set_ylabel("Total no. owners in network")
    ax1.set_xlabel(None)
    sns.boxplot(
        data=df_mcomp_agg,
        x=cat_col,
        y="number_companies",
        fliersize=1,
        palette=colour_dict,
        linewidth=.5,
        ax=ax2
    )
    ax2.grid(visible=True, which="major", axis="y", zorder=0)
    ax2.set_ylabel("No. companies in network")
    ax2.set_xlabel(None)
    sns.boxplot(
        data=df_mcomp_agg,
        x=cat_col,
        y="number_people",
        fliersize=1,
        palette=colour_dict,
        linewidth=.5,
        ax=ax3
    )
    ax3.grid(visible=True, which="major", axis="y", zorder=0)
    ax3.set_ylabel("No. private people in network")
    ax3.set_xlabel(None)
    sns.boxplot(
        data=df_mcomp_agg,
        x=cat_col,
        y="num_hierarch_levels",
        fliersize=1,
        palette=colour_dict,
        linewidth=.5,
        ax=ax4
    )
    ax4.grid(visible=True, which="major", axis="y", zorder=0)
    ax4.set_ylabel("No. hierarchical levels of network")
    ax4.set_xlabel(None)
    fig.tight_layout()
    plt.savefig(out_pth)
    plt.close()

    # sns.pairplot(
    #     data=df_own_agg[[cat_col, "area", "number_people", "number_companies", "num_owners_in_alkis", "num_hierarch_levels"]],
    #     hue=cat_col,
    #     palette=colour_dict,
    #     plot_kws={"s": 3}
    # )
    # plt.show()

    print("done!")


def fig_share_and_location_largest_owners(df_res_pth, df_ct_pth, df_sh_pth, district_shp_pth, out_pth):

    print("Plot map and figure with share and location of the largest owner.")

    grid_res = 4
    version_4km = 1
    gdf_grid = gpd.read_file(rf"00_data\vector\grids\square_grid_{grid_res}km_v{version_4km:02d}_with_12km_POLYIDs.shp")

    df_res = pd.read_csv(df_res_pth)
    df_res["total_area"] = df_res["total_area"] / 10000

    df_ct = pd.read_csv(df_ct_pth)
    df_sh = pd.read_csv(df_sh_pth)
    sub_cols = list(df_sh.columns)
    sub_cols = [col for col in sub_cols if "100" in col]
    sub_cols.append("POLYID")
    df_sh = df_sh[sub_cols]

    df_res = pd.merge(df_res, df_ct, how="left", left_on="POLYID", right_on="id_sp_unit")
    df_res = pd.merge(df_res, df_sh, how="left", left_on="POLYID", right_on="POLYID")

    shp = pd.merge(gdf_grid, df_res, how='inner', left_on='POLYID', right_on='POLYID')
    shp = shp.loc[shp['gini_coeff'] != 0.0]

    shp["top_owner_cat"] = np.nan
    owner_classes = ["PUBLIC", "nCONETW", "aCONETW", "NONPRO", "siCOMP", "a_siCOMP", "noagPR", "agriPR", "CHURCH"]
    for owner_class in owner_classes:
        count_col = f"count_{owner_class}_top1"
        if count_col in list(shp.columns):
            shp.loc[shp[count_col] == 1, "top_owner_cat"] = owner_class

    shp["main_owner_cat"] = ""
    cols = [f"share_{owner_class}_sh100" for owner_class in owner_classes]
    shp["main_owner_cat"] = shp[cols].idxmax(axis=1)
    shp["main_owner_cat"] = shp["main_owner_cat"].apply(lambda x: x.split("_")[1])
    shp["main_owner_cat_share"] = shp[cols].max(axis=1)

    shp["top_owner_cat"] = shp["top_owner_cat"].map(
        {"PUBLIC": "PU", "nCONETW": "CN", "aCONETW": "AH", "NONPRO": "NP", "siCOMP": "nSC", "a_siCOMP": "aSC",
         "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    shp["top_owner_cat"] = pd.Categorical(shp["top_owner_cat"],
                                            categories=["CN", "AH", "PU", "RE", "nPR", "aPR", "aSC", "nSC", "NP"],
                                            ordered=True)

    shp["main_owner_cat"] = shp["main_owner_cat"].map(
        {"PUBLIC": "PU", "nCONETW": "CN", "aCONETW": "AH", "NONPRO": "NP", "siCOMP": "nSC", "a_siCOMP": "aSC",
         "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    shp["main_owner_cat"] = pd.Categorical(shp["main_owner_cat"],
                                          categories=["CN", "AH", "PU", "RE", "nPR", "aPR", "aSC", "nSC", "NP"],
                                          ordered=True)

    ## Map with main owner_class
    colour_dict = {
        "aPR": '#f6f739',
        "nPR": '#bcbd22',
        "aSC": '#c898f5',
        "nSC": '#9467bd',
        "CN": '#ff7f0e',
        "AH": '#fcb97e',
        "PU": '#1f77b4',
        "NP": '#2ca02c',
        "RE": '#7f7f7f'
    }

    col = "top_owner_cat"
    shp["color"] = shp[col].map(colour_dict)
    shp = shp.loc[shp["color"].notna()].copy()

    custom_patches = [Patch(facecolor=colour_dict[v], label=v) for v in ["CN", "AH", "nPR", "aPR", "aSC", "nSC", "PU", "RE", "NP"]]

    print("Plotting")

    matplotlib.rcParams.update({'font.size': 11})

    fig, axs = plt.subplots(1, 2, figsize=plotting_lib.cm2inch(25, 12))
    ax1 = axs[0]
    ax2 = axs[1]
    shp.plot(
        color=shp["color"],
        ax=ax1,
        legend=False,
        edgecolor='none'
    )
    ax1.legend(handles=custom_patches, bbox_to_anchor=(1.1, .01), ncol=math.ceil(len(custom_patches) / 2))  # ,

    shp2 = gpd.read_file(district_shp_pth)
    shp2.plot(
        edgecolor='black',
        facecolor="none",
        ax=ax1,
        lw=0.1,
        zorder=3
    )
    ax1.axis("off")
    ax1.margins(0)
    ax1.set_title("a)", loc="left")

    value_col = "cr1"
    category_col = "top_owner_cat"

    mean_x = shp[value_col].mean()
    ax2.grid(visible=True, which="major", axis="y", zorder=0)
    ax2.set_axisbelow(True)
    sns.kdeplot(x=value_col, hue=category_col, data=shp, ax=ax2, palette=colour_dict, legend=False)
    ax2.set_title("b)", loc="left")
    ax2.set_xlabel("CR1 [%]")
    ax2.set_ylabel(None)
    ax2.tick_params(axis='y', labelrotation=0)
    ax2.axvline(x=mean_x, linestyle="--", color="black")
    bbox_props = dict(boxstyle="round,pad=0.3", fc="w", ec="k", lw=0.72)
    kw = dict(bbox=bbox_props, zorder=2, va="center")
    ax2.annotate(f"Mean: {round(mean_x,1)} %", (mean_x, .04), (mean_x-0.5, .04),  **kw)

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.tick_params(labeltop=False, labelright=False)

    fig.tight_layout()
    plt.tight_layout()
    plt.savefig(out_pth, dpi=300)
    plt.close()


    #####
    col = "main_owner_cat"
    shp["color_main"] = shp[col].map(colour_dict)
    shp = shp.loc[shp["color_main"].notna()].copy()

    fig, axs = plt.subplots(1, 2, figsize=plotting_lib.cm2inch(25, 12))
    ax1 = axs[0]
    ax2 = axs[1]
    shp.plot(
        color=shp["color_main"],
        ax=ax1,
        legend=False,
        edgecolor='none'
    )
    ax1.legend(handles=custom_patches, bbox_to_anchor=(1.1, .01), ncol=math.ceil(len(custom_patches) / 2))  # ,

    shp2 = gpd.read_file(district_shp_pth)
    shp2.plot(
        edgecolor='black',
        facecolor="none",
        ax=ax1,
        lw=0.1,
        zorder=3
    )
    ax1.axis("off")
    ax1.margins(0)
    ax1.set_title("a)", loc="left")

    col = "main_owner_cat"
    shp["color_main"] = shp[col].map(colour_dict)
    shp = shp.loc[shp["color"].notna()].copy()

    shp.plot(
        column="main_owner_cat_share",
        ax=ax2,
        legend=True,
        legend_kwds={'label': "Share of main owner category", 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
        cmap="crest"
    )

    shp2 = gpd.read_file(district_shp_pth)
    shp2.plot(
        edgecolor='black',
        facecolor="none",
        ax=ax2,
        lw=0.1,
        zorder=3
    )
    ax2.axis("off")
    ax2.margins(0)
    ax2.set_title("b)", loc="left")

    fig.tight_layout()
    plt.tight_layout()
    plt.savefig(out_pth[:-4] + "MAIN_OWNER_CAT_TEST.png", dpi=300)
    plt.close()



def fig_comparison_concentrations_measures(df_res_pth, threshold, out_pth):
    grid_res = 4
    version_4km = 1
    gdf_grid = gpd.read_file(
        rf"00_data\vector\grids\square_grid_{grid_res}km_v{version_4km:02d}_with_12km_POLYIDs.shp")

    df_res = pd.read_csv(df_res_pth.format(threshold))
    df_res["total_area"] = df_res["total_area"] / 10000

    shp1 = pd.merge(gdf_grid, df_res, how='inner', left_on='POLYID', right_on='POLYID')
    shp1 = shp1.loc[shp1['gini_coeff'] != 0.0]

    plotting_lib.plot_maps_in_grid(
        shp=shp1,
        out_pth=out_pth,
        cols=["gini_coeff", "hhi", "palma_v2", "cr1", "cr3", "cr5"],
        titles=["a)", "b)", "c)", "d)", "e)", "f)"],
        labels=["Gini coefficient", "HHI", 'Palma' + helper_functions.get_sub("a"), "CR1 [%]", "CR3 [%]", "CR5 [%]"],
        nrow=2,
        shp2_pth=r"00_data\vector\administrative\BB_Landkreise.shp",
        cmap="crest"
    )


def fig_comparison_map_and_histograms_concentration_measures(df_res_pth, threshold, district_shp_pth, out_pth):

    print("Plot maps and histogramms of concentration measures.")
    grid_res = 4
    version_4km = 1
    gdf_grid = gpd.read_file(
        rf"00_data\vector\grids\square_grid_{grid_res}km_v{version_4km:02d}_with_12km_POLYIDs.shp")

    df_res = pd.read_csv(df_res_pth.format(threshold))
    df_res["total_area"] = df_res["total_area"] / 10000

    shp = pd.merge(gdf_grid, df_res, how='inner', left_on='POLYID', right_on='POLYID')
    shp = shp.loc[shp['gini_coeff'] != 0.0]

    shp2 = gpd.read_file(district_shp_pth)

    matplotlib.rcParams.update({'font.size': 12})

    fig = plt.figure(constrained_layout=True, figsize=(12, 11))
    widths = [1, 1, 1]
    heights = [4, 1, 4, 1]
    spec = fig.add_gridspec(
        ncols=3, nrows=4, width_ratios=widths,
        height_ratios=heights)
    cmap = "crest"

    col = "gini_coeff"
    label = "Gini"
    title = "a)"
    ax = fig.add_subplot(spec[0, 0])
    vmin_use = shp[col].quantile(q=0.02)
    vmax_use = shp[col].quantile(q=0.98)
    shp.plot(
        column=col,
        ax=ax,
        legend=True,
        legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
        vmin=vmin_use,
        vmax=vmax_use,
        cmap=cmap
    )
    ax.set_title(title, size=16, x=0.01, y=0.9)
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
    ax_sub.set_ylabel(ylabel="No. polygons")
    ax_sub.set_xlabel(xlabel=label)
    ax_sub.spines['right'].set_visible(False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.tick_params(labeltop=False, labelright=False)

    col = "hhi"
    label = "HHI"
    title = "b)"
    ax = fig.add_subplot(spec[0, 1])
    vmin_use = shp[col].quantile(q=0.02)
    vmax_use = shp[col].quantile(q=0.98)
    shp.plot(
        column=col,
        ax=ax,
        legend=True,
        legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
        vmin=vmin_use,
        vmax=vmax_use,
        cmap=cmap
    )
    ax.set_title(title, size=16, x=0.01, y=0.9)
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
    ax_sub.set_ylabel(ylabel="No. polygons")
    ax_sub.set_xlabel(xlabel=label)
    ax_sub.spines['right'].set_visible(False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.tick_params(labeltop=False, labelright=False)

    col = "palma_v2"
    label = 'Palma' + helper_functions.get_sub("a")
    title = "c)"
    ax = fig.add_subplot(spec[0, 2])
    vmin_use = shp[col].quantile(q=0.02)
    vmax_use = shp[col].quantile(q=0.98)
    shp.plot(
        column=col,
        ax=ax,
        legend=True,
        legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
        vmin=vmin_use,
        vmax=vmax_use,
        cmap=cmap
    )
    ax.set_title(title, size=16, x=0.01, y=0.9)
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
    ax_sub.set_ylabel(ylabel="No. polygons")
    ax_sub.set_xlabel(xlabel=label)
    ax_sub.spines['right'].set_visible(False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.tick_params(labeltop=False, labelright=False)

    col = "cr1"
    label = "CR1 [%]"
    title = "d)"
    ax = fig.add_subplot(spec[2, 0])
    vmin_use = shp[col].quantile(q=0.02)
    vmax_use = shp[col].quantile(q=0.98)
    shp.plot(
        column=col,
        ax=ax,
        legend=True,
        legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
        vmin=vmin_use,
        vmax=vmax_use,
        cmap=cmap
    )
    ax.set_title(title, size=16, x=0.01, y=0.9)
    ax.axis("off")
    shp2.plot(edgecolor='black', facecolor="none", ax=ax, lw=0.3, zorder=2)

    ax_sub = fig.add_subplot(spec[3, 0])
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
    ax_sub.set_ylabel(ylabel="No. polygons")
    ax_sub.set_xlabel(xlabel=label)
    ax_sub.spines['right'].set_visible(False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.tick_params(labeltop=False, labelright=False)

    col = "cr3"
    label = "CR3 [%]"
    title = "e)"
    ax = fig.add_subplot(spec[2, 1])
    vmin_use = shp[col].quantile(q=0.02)
    vmax_use = shp[col].quantile(q=0.98)
    shp.plot(
        column=col,
        ax=ax,
        legend=True,
        legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
        vmin=vmin_use,
        vmax=vmax_use,
        cmap=cmap
    )
    ax.set_title(title, size=16, x=0.01, y=0.9)
    ax.axis("off")
    shp2.plot(edgecolor='black', facecolor="none", ax=ax, lw=0.3, zorder=2)

    ax_sub = fig.add_subplot(spec[3, 1])
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
    ax_sub.set_ylabel(ylabel="No. polygons")
    ax_sub.set_xlabel(xlabel=label)
    ax_sub.spines['right'].set_visible(False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.tick_params(labeltop=False, labelright=False)

    col = "cr5"
    label = "CR5 [%]"
    title = "f)"
    ax = fig.add_subplot(spec[2, 2])
    vmin_use = shp[col].quantile(q=0.02)
    vmax_use = shp[col].quantile(q=0.98)
    shp.plot(
        column=col,
        ax=ax,
        legend=True,
        legend_kwds={'label': label, 'orientation': "horizontal", "fraction": 0.04, "anchor": (0.1, 1.5), "pad": 0.01},
        vmin=vmin_use,
        vmax=vmax_use,
        cmap=cmap
    )
    ax.set_title(title, size=16, x=0.01, y=0.9)
    ax.axis("off")
    shp2.plot(edgecolor='black', facecolor="none", ax=ax, lw=0.3, zorder=2)

    ax_sub = fig.add_subplot(spec[3, 2])
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
    ax_sub.set_ylabel(ylabel="No. polygons")
    ax_sub.set_xlabel(xlabel=label)
    ax_sub.spines['right'].set_visible(False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.tick_params(labeltop=False, labelright=False)

    # fig.tight_layout()
    plt.savefig(out_pth)
    plt.close()


def fig_histograms_change_concentration_measures(df_res_omerge_pth, df_res_mcomp_pth, out_pth):

    print("Plot histograms with change in concentration measures.")

    df_om = pd.read_csv(df_res_omerge_pth)
    df_comp = pd.read_csv(df_res_mcomp_pth)

    df_om.columns = [f"{col}_om" for col in list(df_om.columns)]
    df_comp.columns = [f"{col}_comp" for col in list(df_comp.columns)]

    df_comb = pd.merge(df_om, df_comp, how="left", left_on="POLYID_om", right_on="POLYID_comp")
    df_comb.drop(columns=["POLYID_comp"], inplace=True)
    df_comb.rename(columns={"POLYID_om": "id_sp_unit"}, inplace=True)

    cols = ["gini_coeff", "cr1", "cr3", "cr5", "lac", "palma_v1", "palma_v2", "hhi", "num_owners", "share_p100",
            "share_p95_99", "share_v19", "share_m50", "share_b40"]

    for col in cols:
        df_comb[f"{col}_diff"] = df_comb[f"{col}_comp"] - df_comb[f"{col}_om"]
        df_comb[f"{col}_incr"] = (df_comb[f"{col}_comp"] / df_comb[f"{col}_om"]) - 1

    palma_text = 'Diff. Palma' + helper_functions.get_sub("a")

    plotting_lib.histogramm_in_grid(
        df=df_comb,
        cols=["gini_coeff_diff", "hhi_diff", "palma_v2_diff", "cr1_diff", "cr3_diff", "cr5_diff"],
        nrow=2,
        titles=["a)", "b)", "c)", "d)", "e)", "f)"],
        x_labels=["Diff. Gini", "Diff. HHI", palma_text, "Diff. CR\u2081 [%]", "Diff. CR\u2083 [%]", "Diff. CR\u2083 [%]"],
        y_label="No. polygons",
        out_pth=out_pth,
    )


def table_change_in_concentration_measures(cm_omerge_grid_pth, cm_mcomp_grid_pth, cm_omerge_state_pth, cm_mcomp_state_pth, out_pth):
    cm_omerge_grid = pd.read_csv(cm_omerge_grid_pth)
    cm_mcomp_grid = pd.read_csv(cm_mcomp_grid_pth)
    cm_omerge_state = pd.read_csv(cm_omerge_state_pth)
    cm_mcomp_state = pd.read_csv(cm_mcomp_state_pth)

    cm_omerge_state = cm_omerge_state[["total_area", "num_owners", "gini_coeff",  "hhi", "palma_v2",
                                       "cr1", "cr3", "cr5", "share_p100", "share_p95_99",
                                       "share_v19", "share_m50", "share_b40"]] #"palma_v1", "rosenbluth_index",
    cm_mcomp_state = cm_mcomp_state[["total_area", "num_owners", "gini_coeff",  "hhi", "palma_v2",
                                       "cr1", "cr3", "cr5", "share_p100", "share_p95_99",
                                       "share_v19", "share_m50", "share_b40"]] #"palma_v1", "rosenbluth_index",

    cm_omerge_grid["total_area"] = cm_omerge_grid["total_area"] / 10000
    cm_mcomp_grid["total_area"] = cm_mcomp_grid["total_area"] / 10000

    state_comp = pd.concat([cm_omerge_state, cm_mcomp_state])
    state_comp = state_comp.round({"total_area": 0, "num_owners": 0, "gini_coeff": 2, "cr1": 1, "cr3": 1, "cr5": 1,
                                 "hhi": 0, "palma_v1": 1, "palma_v2": 1, "rosenbluth_index": 3, "share_p100": 1,
                                 "share_p95_99": 1, "share_v19": 1, "share_m50": 1, "share_b40": 1})
    state_comp = state_comp.T
    state_comp.columns = "State\nAggr1", "State\nAggr2"
    state_comp["State diff.\nAggr2-Aggr1"] = state_comp["State\nAggr2"] - state_comp["State\nAggr1"]

    sum_om_gr = cm_omerge_grid.agg(
        {
            "total_area": ["mean", "std"],
            "num_owners": ["mean", "std"],
            "gini_coeff": ["mean", "std"],
            "hhi": ["mean", "std"],
            "palma_v2": ["mean", "std"],
            "cr1": ["mean", "std"],
            "cr3": ["mean", "std"],
            "cr5": ["mean", "std"],
            # "palma_v1": ["mean", "std"],
            # "rosenbluth_index": ["mean", "std"],
            "share_p100": ["mean", "std"],
            "share_p95_99": ["mean", "std"],
            "share_v19": ["mean", "std"],
            "share_m50": ["mean", "std"],
            "share_b40": ["mean", "std"],
        }
    )

    sum_om_gr = sum_om_gr.round({"total_area": 0, "num_owners": 0, "gini_coeff": 2, "cr1": 1, "cr3": 1, "cr5": 1,
                                 "hhi": 0, "palma_v1": 1, "palma_v2": 1, "rosenbluth_index": 3, "share_p100": 1,
                                 "share_p95_99": 1, "share_v19": 1, "share_m50": 1, "share_b40": 1})
    sum_om_gr = sum_om_gr.T
    sum_om_gr["Grid mean\xb1std"] = [f"{row.mean} \xb1({row.std})" for row in sum_om_gr.itertuples()]
    sum_om_gr.columns = [f"{col}\nAggr1" for col in sum_om_gr]

    sum_mc_gr = cm_mcomp_grid.agg(
        {
            "total_area": ["mean", "std"],
            "num_owners": ["mean", "std"],
            "gini_coeff": ["mean", "std"],
            "hhi": ["mean", "std"],
            # "palma_v1": ["mean", "std"],
            "palma_v2": ["mean", "std"],
            "cr1": ["mean", "std"],
            "cr3": ["mean", "std"],
            "cr5": ["mean", "std"],
            # "rosenbluth_index": ["mean", "std"],
            "share_p100": ["mean", "std"],
            "share_p95_99": ["mean", "std"],
            "share_v19": ["mean", "std"],
            "share_m50": ["mean", "std"],
            "share_b40": ["mean", "std"],
        }
    )

    sum_mc_gr = sum_mc_gr.round({"total_area": 0, "num_owners": 0, "gini_coeff": 2, "cr1": 1, "cr3": 1, "cr5": 1,
                                 "hhi": 0, "palma_v1": 1, "palma_v2": 1, "rosenbluth_index": 3, "share_p100": 1,
                                 "share_p95_99": 1, "share_v19": 1, "share_m50": 1, "share_b40": 1})
    sum_mc_gr = sum_mc_gr.T
    sum_mc_gr["Grid mean\xb1std"] = [f"{row.mean} \xb1({row.std})" for row in sum_mc_gr.itertuples()]
    sum_mc_gr.columns = [f"{col}\nAggr2" for col in sum_mc_gr]

    grid_comp = sum_om_gr.join(sum_mc_gr)
    grid_comp["Grid diff.\nAggr2-Aggr1"] = grid_comp["mean\nAggr2"] - grid_comp["mean\nAggr1"]
    grid_comp = grid_comp[["Grid mean\xb1std\nAggr1", "Grid mean\xb1std\nAggr2", "Grid diff.\nAggr2-Aggr1"]]

    compl_comp = state_comp.join(grid_comp)
    compl_comp = compl_comp.round({"Grid diff.\nAggr2-Aggr1": 1, "State diff.\nAggr2-Aggr1":1})
    compl_comp.index = ["Agric. area [ha]", "No. owners", "Gini coefficient", "HHI", 'Palma' + helper_functions.get_sub("a"), "CR1", "CR3", "CR5", "P100", "P95-99", "P90-94", "P41-89", "P1-40"]
    compl_comp.to_csv(out_pth, sep=";")
    print("done!")


def fig_presentation_avg_area(owner_df_pth, threshold, out_pth):
    print("Do something!")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    df[cat_col] = df[cat_col].map({"PUBLIC": "Öffentliche Institutionen",
                                   "nCONETW": "Unternehmensnetzwerke",
                                   "aCONETW": "Landw. Unternehmensverbünde",
                                   "NONPRO": "Non-Profit, Vereine, Stiftungen",
                                   "siCOMP": "Andere Einzelunternehmen",
                                   "a_siCOMP": "Landw. Einzelunternehmen",
                                   "noagPR": "Andere Privatpersonen",
                                   "agriPR": "Landw. Privatpersonen",
                                   "CHURCH": "Religiöse Institutionen"})

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

    ## PLOTTING
    print("\tPlotting")

    colour_dict = {
        "Landw. Privatpersonen": '#f6f739',
        "Andere Privatpersonen": '#bcbd22',
        "Landw. Einzelunternehmen": '#c898f5',
        "Andere Einzelunternehmen": '#9467bd',
        "Unternehmensnetzwerke": '#ff7f0e',
        "Landw. Unternehmensverbünde": '#fcb97e',
        "Öffentliche Institutionen": '#1f77b4',
        "Non-Profit, Vereine, Stiftungen": '#2ca02c',
        "Religiöse Institutionen": '#7f7f7f'
    }

    fig, ax = plt.subplots(figsize=plotting_lib.cm2inch(16, 8))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(labeltop=False, labelright=False)

    ax4 = ax

    # no. plots per owner
    df_num2.sort_values(by="avg_area", ascending=False, inplace=True)
    df_num2[cat_col] = pd.Categorical(df_num2[cat_col], categories=df_num2[cat_col].tolist(), ordered=True)
    ax4.set_title("b)", loc="left")
    ax4.set_axisbelow(True)
    ax4.grid(visible=True, which="both", axis="x", zorder=0)
    sns.barplot(y=cat_col, x="avg_area", data=df_num2, ax=ax4, palette=colour_dict)
    ax4.set_xlabel("Mittlere Fläche pro Eigentümer*In [ha]")
    ax4.set_ylabel(None)
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_pth)
    plt.close()


def fig_stacked_lorenz_curve(owner_df_pth, threshold, out_pth):
    print("Do something!")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    df[cat_col] = df[cat_col].map({"PUBLIC": "Öffentliche Institutionen",
                                   "nCONETW": "Unternehmensnetzwerke",
                                   "aCONETW": "Landw. Unternehmensverbünde",
                                   "NONPRO": "Non-Profit, Vereine, Stiftungen",
                                   "siCOMP": "Andere Einzelunternehmen",
                                   "a_siCOMP": "Landw. Einzelunternehmen",
                                   "noagPR": "Andere Privatpersonen",
                                   "agriPR": "Landw. Privatpersonen",
                                   "CHURCH": "Religiöse Institutionen"})

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

    df_mcomp_agg.sort_values(by="area", inplace=True, ascending=False)
    df_mcomp_agg["percentiles"] = pd.qcut(df_mcomp_agg["area"], q=100, labels=range(1, 101))

    df_sorting = df_mcomp_agg.groupby([cat_col, "percentiles"]).sum().reset_index()
    df_sorting["share"] = (df_sorting["area"] / df_sorting["area"].sum()) * 100
    df_sorting['cumsum'] = df_sorting.groupby([cat_col])['area'].cumsum()
    df_sorting = df_sorting.pivot(index=cat_col, columns='percentiles', values='cumsum')
    df_sorting.sort_values(by=99, inplace=True, ascending=False)
    rename_dict = {col: f"{i+1:02d}_{col}" for i, col in enumerate(df_sorting.index)}

    df_mcomp_agg[cat_col] = df_mcomp_agg[cat_col].map(rename_dict)

    df_agg = df_mcomp_agg.groupby([cat_col, "percentiles"]).sum().reset_index()
    df_agg["share"] = (df_agg["area"] / df_agg["area"].sum()) * 100
    df_agg['cumsum'] = df_agg.groupby([cat_col])['share'].cumsum()
    df_agg["cumsumplt"] = df_agg.groupby(["percentiles"])["cumsum"].cumsum()

    plt_cols = df_agg[cat_col].unique()
    colour_dict = {
        "Landw. Privatpersonen": '#f6f739',
        "Andere Privatpersonen": '#bcbd22',
        "Landw. Einzelunternehmen": '#c898f5',
        "Andere Einzelunternehmen": '#9467bd',
        "Unternehmensnetzwerke": '#ff7f0e',
        "Landw. Unternehmensverbünde": '#fcb97e',
        "Öffentliche Institutionen": '#1f77b4',
        "Non-Profit, Vereine, Stiftungen": '#2ca02c',
        "Religiöse Institutionen": '#7f7f7f'
    }
    colour_dict = {rename_dict[key]: colour_dict[key] for key in colour_dict}
    ############### 1

    df_plt = df_agg.pivot(index='percentiles', columns=cat_col, values='cumsumplt')
    df_plt.reset_index(inplace=True)

    fig, ax = plt.subplots(figsize=plotting_lib.cm2inch(16, 8))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(labeltop=False, labelright=False)
    sns.lineplot(data=df_plt[plt_cols])

    ############### 2
    df_plt2 = df_agg.pivot(index='percentiles', columns=cat_col, values='cumsum')

    # create stacked bar chart
    plt.figure(figsize=(plotting_lib.cm2inch(16, 16)))
    y = [df_plt2[col].tolist() for col in plt_cols]
    colors = [colour_dict[cat] for cat in plt_cols]
    plt.stackplot(df_plt2.index, y, colors=colors)
    legend_elements = [Patch(facecolor=colour_dict[key], edgecolor=None, label=key.split('_')[1]) for key in colour_dict]
    plt.legend(handles=legend_elements, loc='upper left')
    plt.xlabel('Anteil an Eigentümern [%]')
    plt.ylabel('Anteil an Land [%]')
    plt.savefig(rf"12_area_calculations\figures\stacked_lorenz_curve.png")
    plt.close()

    #############
    plt.figure(figsize=(plotting_lib.cm2inch(16, 16)))
    y = [df_plt2[col].tolist() for col in plt_cols]
    colors = [colour_dict[cat] for cat in plt_cols]
    plt.stackplot(df_plt2.index, y, colors=colors)
    legend_elements = [Patch(facecolor=colour_dict[key], edgecolor=None, label=key.split('_')[1]) for key in
                       colour_dict]
    plt.legend(handles=legend_elements, loc='upper left')
    plt.xlabel('Anteil an Eigentümern [%]')
    plt.ylabel('Anteil an Land [%]')
    plt.xlim(90, 100)
    plt.savefig(rf"12_area_calculations\figures\stacked_lorenz_curve_zoom_in.png")
    plt.close()

    #############
    plt.figure(figsize=(plotting_lib.cm2inch(16, 16)))
    y = [df_plt2[col].tolist() for col in plt_cols]
    colors = [colour_dict[cat] for cat in plt_cols]
    plt.stackplot(df_plt2.index, y, colors=colors)
    legend_elements = [Patch(facecolor=colour_dict[key], edgecolor=None, label=key.split('_')[1]) for key in
                       colour_dict]
    plt.legend(handles=legend_elements, loc='upper left')
    plt.xlabel('Anteil an Eigentümern [%]')
    plt.ylabel('Anteil an Land [%]')
    plt.xlim(80, 100)
    plt.savefig(rf"12_area_calculations\figures\stacked_lorenz_curve_zoom_in_80-100.png")
    plt.close()

    #############
    plt.figure(figsize=(plotting_lib.cm2inch(16, 16)))
    y = [df_plt2[col].tolist() for col in plt_cols]
    colors = [colour_dict[cat] for cat in plt_cols]
    plt.stackplot(df_plt2.index, y, colors=colors)
    legend_elements = [Patch(facecolor=colour_dict[key], edgecolor=None, label=key.split('_')[1]) for key in
                       colour_dict]
    plt.legend(handles=legend_elements, loc='upper left')
    plt.xlabel('Anteil an Eigentümern [%]')
    plt.ylabel('Anteil an Land [%]')
    plt.xlim(70, 100)
    plt.savefig(rf"12_area_calculations\figures\stacked_lorenz_curve_zoom_in_70-100.png")
    plt.close()

    ############### 2
    colors = [colour_dict[cat] for cat in plt_cols]
    df_plt2.plot(kind='bar', stacked=True, color=colors)
    legend_elements = [Patch(facecolor=colour_dict[key], edgecolor=None, label=key.split('_')[1]) for key in colour_dict]
    plt.legend(handles=legend_elements, loc='upper left')
    plt.xlabel('Anteil an Eigentümern [%]')
    plt.ylabel('Anteil an Land [%]')
    plt.xticks([0, 20, 40, 60, 80, 100], [0, 20, 40, 60, 80, 100])
    plt.tight_layout()
    plt.savefig(rf"12_area_calculations\figures\stacked_shares_bars.png")
    plt.close()

    ############### 3
    df_plt3 = df_mcomp_agg.sort_values(by="area")
    df_plt3["share"] = (df_plt3["area"] / df_plt3["area"].sum()) * 100
    df_plt3['cumsum'] = df_plt3.groupby([cat_col])['share'].cumsum()
    df_plt3["x"] = range(len(df_plt3))
    df_plt3["cumsumplt"] = df_plt3["share"].cumsum()

    fig, ax = plt.subplots(figsize=plotting_lib.cm2inch(16, 8))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(labeltop=False, labelright=False)
    # for cat in plt_cols:
    #     df_sub = df_plt3.loc[df_plt2[cat_col] == cat].copy()
    #     ax.plot(df_sub["x"], df_sub["cumsumplt"])

def table_brandenburg_owners(owner_df_pth, threshold, out_pth):
    print("Create table with information on owners in Brandenburg")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    df[cat_col] = df[cat_col].map({"PUBLIC": "Öffentliche Institutionen",
                                   "nCONETW": "Unternehmensnetzwerke",
                                   "aCONETW": "Landw. Unternehmensverbünde",
                                   "NONPRO": "Non-Profit, Vereine, Stiftungen",
                                   "siCOMP": "Andere Einzelunternehmen",
                                   "a_siCOMP": "Landw. Einzelunternehmen",
                                   "noagPR": "Andere Privatpersonen",
                                   "agriPR": "Landw. Privatpersonen",
                                   "CHURCH": "Religiöse Institutionen"})

    df_agg = df.groupby(["fstate_mcomp", "new_category"]).agg(
        area=pd.NamedAgg("area", "sum")
    ).reset_index()

    df_nonagri_priv = df_agg.loc[df_agg["new_category"] == "Andere Privatpersonen"].copy()
    df_nonagri_priv["share"] = (df_nonagri_priv["area"] / df_nonagri_priv["area"].sum()) * 100

    df_nonagri_priv.to_csv(r"12_area_calculations\figures\table_non_agric_private_owners_orgins.csv", index=False)

    print("\tDone!")


def table_number_owners_and_area_in_company_networks(owner_df_pth, threshold, out_pth):
    print("Do something!")
    print("\tRead owner data")
    df = pd.read_csv(owner_df_pth.format(threshold), sep=";")

    cat_col = "new_category"
    print(df[cat_col].unique())
    t = df.loc[df[cat_col].isna()].copy()
    df[cat_col] = df[cat_col].map(
        {"PUBLIC": "PU", "nCONETW": "CN", "aCONETW": "AH", "NONPRO": "NP", "siCOMP": "nSC", "a_siCOMP": "aSC",
         "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    # labels = ["PU", "CN", "NP", "nSC", "aSC", "nPR", "aPR", "CH"]
    labels = ["aPR", "nPR", "aSC", "nSC", "AH", "CN", "PU", "NP", "RE"]

    df = df.loc[df[cat_col].isin(["CN", "AH"])].copy()

    df_agg_om = df.groupby("owner_merge").agg(
        distance=pd.NamedAgg("distance", "mean")
    ).reset_index()

    df_agg_cn = df.groupby(f"community_{threshold}").agg(
        distance=pd.NamedAgg("mcomp_dist", "mean")
    ).reset_index()

    out_dict = {
        "no. networks": len(df[f"community_{threshold}"].unique()),
        "no. owners aggr1": len(df["owner_merge"].unique()),
        "no. people": len(df.loc[df["level1"] == 1, "owner_merge"].unique()),
        "no. companies": len(df.loc[df["level1"] == 2, "owner_merge"].unique()),
        "area": df["area"].sum(),
        "mean distance aggr1": df_agg_om["distance"].mean(),
        "mean_distance aggr1": df_agg_cn["distance"].mean(),
        "distance_increase": df_agg_cn["distance"].mean() - df_agg_om["distance"].mean()
    }

    df_out = pd.DataFrame.from_dict(out_dict, orient="index").reset_index()
    df_out.columns = ["variable", "value"]

    df_out.to_csv(out_pth, index=False)


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
    #
    # fig_share_and_characteristics_owner_categories(
    #     owner_df_pth=OWNER_DF_FOR_PLOTTING,
    #     threshold=50,
    #     out_pth=rf"14_paper_figures\figures\fig_share_and_characteristics_owner_categories.png"
    # )

    # fig_appendix_change_in_distances_aggregation_levels(
    #     owner_df_pth=OWNER_DF_FOR_PLOTTING,
    #     threshold=50,
    #     out_pth=rf"14_paper_figures\figures\fig_appendix_change_in_distances_aggregation_levels.png"
    # )
    #
    # fig_appendix_characteristics_company_networks(
    #     owner_df_pth=OWNER_DF_FOR_PLOTTING,
    #     threshold=50,
    #     out_pth=rf"14_paper_figures\figures\fig_appendix_characteristics_company_networks.png"
    # )
    #
    # table_largest_owner_examples(
    #     owner_df_pth=OWNER_DF_FOR_PLOTTING,
    #     threshold=50,
    #     out_pth=rf"14_paper_figures\tables\table_largest_owner_examples.csv"
    # )

    table_number_owners_and_area_in_company_networks(
        owner_df_pth=OWNER_DF_FOR_PLOTTING,
        threshold=50,
        out_pth=rf"14_paper_figures\tables\table_number_owners_and_area_in_company_networks.csv"
    )

    # fig_comparison_map_and_histograms_concentration_measures(
    #     threshold=50,
    #     df_res_pth=rf"11_ownership_concentration\mw_grid\mw_mean_conc_meas-mother_companies-comm_w_thr{threshold}-iacs_areas.csv",
    #     district_shp_pth=DISTRICT_SHP_PTH,
    #     out_pth=r"14_paper_figures\figures\fig_comparison_histograms_concentrations_measures.png"
    # )
    #
    # fig_histograms_change_concentration_measures(
    #     df_res_omerge_pth=rf"11_ownership_concentration\mw_grid\mw_mean_conc_meas-owner_merge-iacs_areas.csv",
    #     df_res_mcomp_pth=rf"11_ownership_concentration\mw_grid\mw_mean_conc_meas-mother_companies-comm_w_thr{threshold}-iacs_areas.csv",
    #     out_pth=r"14_paper_figures\figures\fig_histograms_change_concentration_measures.png"
    # )
    #
    # fig_share_and_location_largest_owners(
    #     df_res_pth=rf"11_ownership_concentration\mw_grid\mw_mean_conc_meas-mother_companies-comm_w_thr{threshold}-iacs_areas.csv",
    #     df_ct_pth=rf"11_ownership_concentration\mw_grid\mw_counts_categories_in_topx-grid_4km_v01-mother_companies-comm_w_thr{threshold}-iacs_areas.csv",
    #     df_sh_pth=rf"11_ownership_concentration\mw_grid\mw_mean_share_counts_categories-mother_companies-comm_w_thr50-iacs_areas.csv",
    #     district_shp_pth=DISTRICT_SHP_PTH,
    #     out_pth=r"14_paper_figures\figures\fig_share_and_location_largest_owners.png"
    # )

    # table_change_in_concentration_measures(
    #     cm_omerge_grid_pth=rf"11_ownership_concentration\mw_grid\mw_mean_conc_meas-owner_merge-iacs_areas.csv",
    #     cm_mcomp_grid_pth=rf"11_ownership_concentration\mw_grid\mw_mean_conc_meas-mother_companies-comm_w_thr{threshold}-iacs_areas.csv",
    #     cm_omerge_state_pth=rf"11_ownership_concentration\state\state_conc_meas-owner_merge-iacs_areas.csv",
    #     cm_mcomp_state_pth=rf"11_ownership_concentration\state\state_conc_meas-mother_companies-comm_w_thr{threshold}-iacs_areas.csv",
    #     out_pth=rf"14_paper_figures\tables\table_change_in_concentration_measures.csv"
    # )
    #
    # fig_presentation_avg_area(
    #         owner_df_pth=OWNER_DF_FOR_PLOTTING,
    #         threshold=50,
    #         out_pth=rf"14_paper_figures\figures\fig_presentation_avg_area.png"
    # )

    ## only private non-agricultural people
    # fig_comparison_map_and_histograms_concentration_measures(
    #     threshold=50,
    #     df_res_pth=rf"11_ownership_concentration\mw_grid\mw_mean_conc_meas-mother_companies-comm_w_thr{threshold}_noagPR.csv",
    #     municip_shp_pth=MUNICIP_SHP_PTH,
    #     out_pth=r"14_paper_figures\figures\fig_comparison_histograms_concentrations_measures_only_non-agricultural_people.png"
    # )
    #
    # fig_stacked_lorenz_curve(
    #     owner_df_pth=OWNER_DF_FOR_PLOTTING,
    #     threshold=50,
    #     out_pth=rf"14_paper_figures\figures\stacked_lorenz_curve.png")
    #
    # table_brandenburg_owners(
    #     owner_df_pth=OWNER_DF_FOR_PLOTTING,
    #     threshold=50,
    #     out_pth=rf"14_paper_figures\tables\stacked_lorenz_curve.png")

    e_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + s_time)
    print("end: " + e_time)


if __name__ == '__main__':
    main()
