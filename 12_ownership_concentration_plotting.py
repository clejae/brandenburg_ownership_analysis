# Author:
# github repository:

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import pandas as pd
import geopandas as gpd
import os
import time
import json
import statistics
import numpy as np

## Project library
import conc_meas_lib
import plotting_lib
import helper_functions

# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input paths
GRID_FOLDER = r"00_data\vector\grids"
FIGURES_FOLDER_GRID_MW = r"12_ownership_concentration_plotting\mw_grid"
CONC_MEASURES_MW_GRID_VERSIONS_PTH = r"11_ownership_concentration\mw_grid\mw_mean_conc_meas-{0}.csv"
SHARE_CATEG_MW_GRID_COMBINED_PTH = r"11_ownership_concentration\mw_grid\mw_mean_share_counts_categories-{0}.csv"
COUNT_CATEG_MW_GRID_COMBINED_PTH = r"11_ownership_concentration\mw_grid\mw_mean_counts_categories_in_topx-{0}.csv"

MUNICIP_SHP_PTH = r"00_data\vector\administrative\BB_municipalities.shp"
# ------------------------------------------ PROCESSING ------------------------------------------#

def plot_concentration_measures_for_mw_grid(df_res_pth, threshold, categories, figures_folder,
                                            descr, municip_shp_pth, df_sh_pth=None, df_ct_pth=None, vector_out_pth=False):

    """

    Args:
        df_res_pth:
        threshold:
        descr:
        categories:
        figures_folder:
        municip_shp_pth:
        df_sh_pth:
        df_ct_pth:
        vector_out_pth:

    Returns:

    """

    helper_functions.create_folder(figures_folder)
    helper_functions.create_folder(f"{figures_folder}\maps")
    helper_functions.create_folder(f"{figures_folder}\graphs")

    grid_res = 4
    version_4km = 1
    gdf_grid = gpd.read_file(rf"00_data\vector\grids\square_grid_{grid_res}km_v{version_4km:02d}_with_12km_POLYIDs.shp")

    df_res = pd.read_csv(df_res_pth)
    df_res["total_area"] = df_res["total_area"] / 10000
    if df_sh_pth:
        df_sh = pd.read_csv(df_sh_pth)
        df_res = pd.merge(df_res, df_sh, how="left", on="POLYID")
    if df_ct_pth:
        df_ct = pd.read_csv(df_ct_pth)
        df_res = pd.merge(df_res, df_ct, how="left", on="POLYID")

    shp1 = pd.merge(gdf_grid, df_res, how='inner', left_on='POLYID', right_on='POLYID')
    shp1 = shp1.loc[shp1['gini_coeff'] != 0.0]

    if vector_out_pth:
        shp1.to_file(vector_out_pth)

    # shp1["main_owner_cat"] = ""
    # for owner_class in categories:
    #     shp1.loc[shp1[f"count_{owner_class}_top1"] == 1, "main_owner_cat"] = owner_class
    #
    # shp1["main_owner_cat"] = shp1["main_owner_cat"].map(
    #     {"PUBLIC": "PU", "nCONETW": "CN", "aCONETW": "AH", "NONPRO": "NP", "siCOMP": "nSC", "a_siCOMP": "aSC",
    #      "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    # shp1["main_owner_cat"] = pd.Categorical(shp1["main_owner_cat"],
    #                                         categories=["CN", "AH", "PU", "RE", "nPR", "aPR", "aSC", "nSC", "NP"],
    #                                         ordered=True)
    #
    # ## Map with main owner_class
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
    # col = "main_owner_cat"
    # out_folder = f"{figures_folder}\maps\main_owner_category"
    # helper_functions.create_folder(out_folder)
    # out_pth = rf"{out_folder}\mw_map_{col}_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
    # plotting_lib.plot_map_categorical(
    #     shp=shp1,
    #     out_pth=out_pth,
    #     col=col,
    #     colour_dict=colour_dict,
    #     shp2_pth=municip_shp_pth)

    # col = "cr1"
    # out_folder = f"{figures_folder}\graphs\main_owner_category_shares"
    # helper_functions.create_folder(out_folder)
    # out_pth = rf"{out_folder}\mw_boxplot_main_owner_cat_vs_{col}_mother_companies_w_thr{threshold}{descr}.png"
    # plotting_lib.boxplot_by_categories(
    #     df=shp1,
    #     category_col="main_owner_cat",
    #     value_col=col,
    #     showfliers=True,
    #     out_pth=out_pth,
    #     x_label=None,
    #     y_label="CR1",
    #     colour_dict=colour_dict
    # )

    # Maps of concentration measures
    for col in ["hhi", "rosenbluth_index", "cr1", "cr3", "cr5", "hhi", "gini_coeff", "palma_v1", "palma_v2"]:
        tr = shp1[col].mean() + shp1[col].std()
        out_folder = f"{figures_folder}\maps\concentration_measures"
        helper_functions.create_folder(out_folder)
        out_pth = rf"{out_folder}\mw_map_{col}_gr{round(tr,2)}-mean+std_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
        plotting_lib.plot_map(
            shp=shp1.loc[shp1[col] > tr],
            out_pth=out_pth,
            col=col,
            shp2_pth=municip_shp_pth)

        out_pth = rf"{out_folder}\mw_map_{col}_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
        plotting_lib.plot_map(
            shp=shp1,
            out_pth=out_pth,
            col=col,
            shp2_pth=municip_shp_pth,
            extremes=True
        )

    out_folder = f"{figures_folder}\maps\concentration_measures"
    helper_functions.create_folder(out_folder)
    out_pth = rf"{out_folder}\\mw_map_comp_conc_meas_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
    plotting_lib.plot_maps_in_grid(
        shp=shp1,
        out_pth=out_pth,
        cols=["gini_coeff", "hhi", "palma_v2", "cr1", "cr3", "cr5"],
        titles=["a)", "b)", "c)", "d)", "e)", "f)"],
        labels=["Gini coefficient", "HHI", "T100/B90", "CR1", "CR3", "CR5"],
        nrow=2,
        shp2_pth=municip_shp_pth
    )

    ## Maps of shares per category
    ## Histograms of share per category
    ## Scatterplots of shares per category vs palma indeces
    # for cat in categories:
    #     for share in [1, 10, 100]:
    #         col1 = f"share_{cat}_sh{share}"
    #         out_folder = f"{figures_folder}\maps\shares_per_owner_category"
    #         helper_functions.create_folder(out_folder)
    #         out_pth = rf"{out_folder}\mw_map_{col1}_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
    #         plotting_lib.plot_map(
    #             shp=shp1,
    #             out_pth=out_pth,
    #             col=col1,
    #             shp2_pth=municip_shp_pth)
    #
    #         out_folder = f"{figures_folder}\graphs\shares_per_owner_category_histograms"
    #         helper_functions.create_folder(out_folder)
    #         out_pth = f"{out_folder}\mw_histogramm_{col1}_mother_companies_w_thr{threshold}{descr}.png"
    #         plotting_lib.histogramm(
    #             df=shp1,
    #             col=col1,
    #             out_pth=out_pth,
    #             x_label=f"{col1} in 4km mov.wind. grid",
    #             y_label="Number of 4km grid cells")
    #
    #         for col2 in ["palma_v2", "palma_v1"]:
    #             if col1 in shp1.columns:
    #                 out_folder = f"{figures_folder}\graphs\shares_per_owner_category_vs_palmas\{col2}\share{share}"
    #                 helper_functions.create_folder(out_folder)
    #                 out_pth = rf"{out_folder}\mw_scatter_{col1}_vs_{col2}_mother_companies_w_thr{threshold}{descr}.png"
    #                 plotting_lib.scatterplot_two_columns(
    #                     df=shp1,
    #                     col1=col1,
    #                     col2=col2,
    #                     out_pth=out_pth,
    #                     xminmax=(0, .5),
    #                     x_label=col1,
    #                     y_label=f"{col2}")
    #
    #                 out_pth = rf"{out_folder}\mw_scatter_{col1}_vs_{col2}_mother_companies_w_thr{threshold}_xlog{descr}.png"
    #                 plotting_lib.scatterplot_two_columns(
    #                     df=shp1,
    #                     col1=col1,
    #                     col2=col2,
    #                     out_pth=out_pth,
    #                     y_label=f"{col2}",
    #                     x_log=True)
    #             else:
    #                 print(col1,  "not in df")

    ## Maps of occurences per category in top x owners
    ## Scatterplots of occurences per category in top x owners vs. concentration ratios
    for x in [1, 3, 5]:
        for cat in categories:
            col1 = f"count_{cat}_top{x}"
            out_folder = f"{figures_folder}\maps\occurences_per_category_in_top_x_owners"
            helper_functions.create_folder(out_folder)
            out_pth = rf"{out_folder}\mw_map_{col1}_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
            plotting_lib.plot_map(
                shp=shp1,
                out_pth=out_pth,
                col=col1,
                shp2_pth=municip_shp_pth)

            col2 = f"cr{x}"
            if col1 in shp1.columns:
                out_folder = f"{figures_folder}\graphs\occurences_per_category_in_top_x_owners\{col2}"
                helper_functions.create_folder(out_folder)
                out_pth = rf"{out_folder}\mw_scatter_{col1}_vs_{col2}_mother_companies_w_thr{threshold}{descr}.png"
                plotting_lib.scatterplot_two_columns(
                    df=shp1,
                    col1=col1,
                    col2=col2,
                    out_pth=out_pth,
                    x_label=col1,
                    y_label=f"{col2}")
            else:
                print(col1, "not in df")

    ## Histograms of concentration measures
    cols = ["cr1", "cr3", "cr5", "hhi", "total_area", "num_owners", "lac", "palma_v1", "palma_v2", "gini_coeff",
            "share_p100", "share_p95_99", "share_v19", "share_m50", "share_b40"]
    for col in cols:
        out_folder = f"{figures_folder}\graphs\histograms_of_concentration_measures"
        helper_functions.create_folder(out_folder)
        out_pth = f"{out_folder}\mw_histogramm_{col}_mother_companies_w_thr{threshold}{descr}.png"
        plotting_lib.histogramm(
            df=shp1,
            col=col,
            out_pth=out_pth,
            x_label=f"{col} in 4km mov.wind. grid",
            y_label="Number of 4km grid cells")

    ## Scatterplots of concentration measures vs total area
    cols = ["gini_coeff", "cr1", "cr3", "cr5", "hhi", "rosenbluth_index", "palma_v1"]
    for col in cols:
        out_folder = f"{figures_folder}\graphs\scatterplots_of_concentration_measures_vs_total_area"
        helper_functions.create_folder(out_folder)
        out_pth = f"{out_folder}\mw_scatter_{col}_vs_total_area_mcomp_w_thr{threshold}{descr}.png"
        plotting_lib.scatterplot_two_columns(
            df=shp1,
            col1=col,
            col2="total_area",
            out_pth=out_pth,
            x_label=f"{col} - moving window",
            y_label="mean agric. area in 12km mw [ha]"
        )



def explain_concentration_with_regression(df, y_col, xs_col):
    import statsmodels.api as sm

    X = np.array(df[xs_col])
    y = np.array(df[y_col]).T

    X2 = sm.add_constant(X)
    est = sm.OLS(y, X2)
    est2 = est.fit()
    xnames = ["const"] + xs_col
    print(est2.summary(xname=(xnames)))


def explore_concentration_at_mw_grid(pth_conc_owner_merge, pth_conc_mother_company, threshold, out_tables_folder):

    df_om = pd.read_csv(pth_conc_owner_merge)
    df_comp = pd.read_csv(pth_conc_mother_company)

    df_om.columns = [f"{col}_om" for col in list(df_om.columns)]
    df_comp.columns = [f"{col}_comp" for col in list(df_comp.columns)]

    df_comb = pd.merge(df_om, df_comp, how="left", left_on="POLYID_om", right_on="POLYID_comp")
    df_comb.drop(columns=["POLYID_comp"], inplace=True)
    df_comb.rename(columns={"POLYID_om": "id_sp_unit"}, inplace=True)

    cols = ["gini_coeff", "cr1", "cr3", "cr5", "hhi", "lac", "palma_v1", "palma_v2", "num_owners", "share_p100", "share_p95_99", "share_v19",
            "share_m50", "share_b40"]

    for col in cols:
        df_comb[f"{col}_diff"] = df_comb[f"{col}_comp"] - df_comb[f"{col}_om"]
        df_comb[f"{col}_incr"] = (df_comb[f"{col}_comp"] / df_comb[f"{col}_om"]) - 1

    out_pth = fr"{out_tables_folder}\conc_meas-mw_grid-diff_mcomp_w_thr{threshold}-owner_merge.csv"
    df_comb.to_csv(out_pth, index=False)

    df_descr = df_comb.describe(percentiles=[.05, .25, .50, .75, .95])
    out_pth = fr"{out_tables_folder}\descr-mw_grid-diff_mcomp_w_thr{threshold}-owner_merge.csv"
    df_descr.to_csv(out_pth, index=False)


def plot_concentration_differences_mw(df_diff_pth, shp_pth, figures_folder, threshold):

    df = pd.read_csv(df_diff_pth, dtype={'id_sp_unit': str})

    shp_municip = gpd.read_file(shp_pth)
    shp = pd.merge(shp_municip, df, how='left', left_on='POLYID', right_on='id_sp_unit')

    cols = ["gini_coeff", "hhi", "cr1", "cr3", "cr5", "palma_v1", "palma_v2"]
    for col in cols:
        plotting_lib.plot_map(
            shp=shp,
            out_pth=rf"{figures_folder}\maps\mw\change_maps_difference\map_difference_{col}_mw_grid_mother_companies_w_thr{threshold}.png",
            col=f"{col}_diff")

    ## Plot difference maps
    # cols = ["gini_coeff", "cr1", "cr3", "cr5", "lac", "share_p100", "share_p95_99", "share_v19",
    #         "share_m50", "share_b40"]
    # for col in cols:
    #     plotting_lib.plot_map(
    #         shp=shp,
    #         out_pth=rf"{figures_folder}\maps\mw\change_maps_difference\map_difference_{col}_mw_grid_mother_companies_w_thr{threshold}.png",
    #         col=f"{col}_diff",
    #         cmap=plt.cm.Reds,
    #         centralize=True)

    # cols = ["num_owners"]
    # for col in cols:
    #     plotting_lib.plot_map(
    #         shp=shp,
    #         out_pth=rf"{figures_folder}\maps\mw\change_maps_difference\map_difference_{col}_mw_grid_mother_companies_w_thr{threshold}.png",
    #         col=f"{col}_diff",
    #         cmap=plt.cm.Blues_r,
    #         centralize=False,
    #         use_norm=False
    #     )


    ## Plot histogramm
    # df_hist = shp.copy()
    # plotting_lib.histogramm(
    #     df=df_hist,
    #     col="palma_v1_diff",
    #     out_pth=rf"{figures_folder}\graphs\histogramm_difference_palma_v1_diff_mw_grid_mother_companies_w_thr{threshold}.png",
    #     bins=50,
    #     x_label="Difference in Palma index",
    #     y_label="Number of grid cells"
    # )



def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)

    os.chdir(WD)

    ## GRIDS
    threshold = 50
    descr = f"mother_companies-comm_w_thr{threshold}-iacs_areas"
    plot_concentration_measures_for_mw_grid(
        df_res_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH.format(descr),
        threshold=threshold,
        categories=["PUBLIC", "aCONETW", "noagPR", "NONPRO", "a_siCOMP", "agriPR", "CHURCH", "nCONETW"],
        figures_folder=FIGURES_FOLDER_GRID_MW,
        municip_shp_pth=MUNICIP_SHP_PTH,
        descr=f"_{descr}",
        # df_sh_pth=SHARE_CATEG_MW_GRID_COMBINED_PTH.format(descr),
        df_ct_pth=COUNT_CATEG_MW_GRID_COMBINED_PTH.format(descr),
        vector_out_pth=False)

    # for owner_class in ["noagPR"]:#["PUBLIC", "nCONETW", "aCONETW", "NONPRO", "siCOMP", "a_siCOMP", "noagPR", "agriPR", "CHURCH"]:
    #     plotting_wrapper_for_mw_grids(
    #         threshold=50,
    #         df_res_pth=CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold, f"{owner_class}"),
    #         # df_sh_pth=SHARE_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH,
    #         # df_ct_pth=COUNT_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH.format(4, 1, threshold),
    #         categories=["PUBLIC", "CONETW", "noagPR", "NONPRO", "siCOMP", "agriPR", "CHURCH", "COOPER", "BVVG"],
    #         descr=owner_class
    #     )
    #
    # ## Plotting
    # # categories = ['1_1_1', '1_9_1', '2_9_1', '2_1_4', '2_1_2', '5_2_2', '2_1_5', '5_2_1', '3_1_1',
    # #               '2_9_2', '5_2_5', '3_2_2', '5_2_4', '4_1_1', '5_3_1', '2_1_8', '2_2_2', '2_1_7',
    # #               '1_2_1', '5_9_1', '5_1_1', '1_2_3', '2_1_1', '90', '5_2_3', '1_2_2', '3_2_3']
    # # categories = ['2_2_1', '2_1_3', '3_9_1', '2_1_6', '1_2_4', '3_1_2', '1_1_2']
    #
    # vector_out_pth = rf"13_ownership_concentration\vector\grids\{os.path.basename(CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold))[:-3]}shp"
    #
    # plotting_wrapper_for_mw_grids(
    #     threshold=50,
    #     df_res_pth=CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold),
    #     # df_sh_pth=SHARE_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH,
    #     df_ct_pth=COUNT_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH.format(4, 1, threshold),
    #     categories=["PUBLIC", "CONETW", "noagPR", "NONPRO", "siCOMP", "agriPR", "CHURCH", "COOPER", "BVVG"],
    #     descr="",
    #     vector_out_pth=vector_out_pth
    # )
    #
    #
    # plotting_wrapper_for_mw_grids(
    #     threshold=50,
    #     df_res_pth=CONC_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH,
    #     # df_sh_pth=SHARE_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH,
    #     # df_ct_pth=COUNT_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH,
    #     categories=["CONETW", "noagPR", "NONPRO", "siCOMP", "agriPR", "COOPER"],
    #     descr="_privland"
    # )
    #
    # ## Explain concentration measures
    # df2 = pd.read_csv(SHARE_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold))
    # df1 = pd.read_csv(CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold))
    # df = pd.merge(df1, df2, how="left", on="POLYID")
    #
    # xs_col = ["share_PUBLIC_sh100", "share_CONETW_sh100", "share_noagPR_sh100", "share_NONPRO_sh100",
    #           "share_siCOMP_sh100", "share_agriPR_sh100", "share_CHURCH_sh100", "share_COOPER_sh100",
    #           "share_BVVG_sh100",
    #           "count_PUBLIC_sh100", "count_CONETW_sh100", "count_noagPR_sh100", "count_NONPRO_sh100",
    #           "count_siCOMP_sh100", "count_agriPR_sh100", "count_CHURCH_sh100", "count_COOPER_sh100",
    #           "count_BVVG_sh100"]
    # xs_col = ["count_PUBLIC_sh100", "count_CONETW_sh100", "count_noagPR_sh100", "count_NONPRO_sh100",
    #           "count_siCOMP_sh100", "count_agriPR_sh100", "count_CHURCH_sh100", "count_COOPER_sh100",
    #           "count_BVVG_sh100"]
    # xs_col = ["share_PUBLIC_sh100", "share_CONETW_sh100", "share_noagPR_sh100", "share_NONPRO_sh100",
    #           "share_siCOMP_sh100", "share_agriPR_sh100", "share_CHURCH_sh100", "share_COOPER_sh100",
    #           "share_BVVG_sh100"]
    # xs_col = ["share_CONETW_sh100", "share_noagPR_sh100", "share_NONPRO_sh100",
    #           "share_siCOMP_sh100", "share_agriPR_sh100", "share_COOPER_sh100"]
    # plot_corr(df[xs_col])
    # y_col = "palma_v2"
    # tr = df[y_col].mean() + df[y_col].std()
    # explain_concentration_with_regression(
    #     df=df.loc[df[y_col] > tr],
    #     y_col=y_col,
    #     xs_col=xs_col)
    #
    # xs_col = ["share_CONETW_sh100", "share_noagPR_sh100", "share_siCOMP_sh100", "share_agriPR_sh100",
    #           "share_COOPER_sh100"]
    # y_col = "gini_coeff"
    # tr = df[y_col].mean() + df[y_col].std()
    # explain_concentration_with_regression(
    #     df=df.loc[df[y_col] > tr],
    #     y_col=y_col,
    #     xs_col=xs_col)
    #
    # xs_col = ["share_CONETW_sh100", "share_noagPR_sh100", "share_siCOMP_sh100", "share_agriPR_sh100",
    #           "share_COOPER_sh100"]
    # y_col = "hhi"
    # tr = df[y_col].mean() + df[y_col].std()
    # explain_concentration_with_regression(
    #     df=df.loc[df[y_col] > tr],
    #     y_col=y_col,
    #     xs_col=xs_col)
    #
    # explore_concentration_at_mw_grid(pth_conc_owner_merge=CONC_OWNER_MERGE_GRID_MW_COMBINED_PTH,
    #                                  pth_conc_mother_company=CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(
    #                                      threshold),
    #                                  threshold=threshold,
    #                                  out_tables_folder=TABLES_FOLDER_GRID_MW)
    #
    # out_pth_diff_vals = fr"{TABLES_FOLDER_GRID_MW}\conc_meas-mw_grid-diff_mcomp_w_thr{threshold}-owner_merge.csv"
    # shp_pth = rf"00_data\vector\grids\square_grid_4km_v01_with_12km_POLYIDs.shp"
    # plot_concentration_differences_mw(df_diff_pth=out_pth_diff_vals,
    #                                   shp_pth=shp_pth,
    #                                   figures_folder=FIGURES_FOLDER_GRID_MW,
    #                                   threshold=threshold)


    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)


if __name__ == '__main__':
    main()