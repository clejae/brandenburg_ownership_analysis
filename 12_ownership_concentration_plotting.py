# Clemens Jänicke
# github Repo: https://github.com/clejae

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
ALKIS_IACS_GRID_PTH = r"10_alkis_intersection_with_other_layers\alkis_iacs_grid_{0}km_v{1:02d}_inters.shp"
OWNERS_W_THRESH_PTH = r"11_community_classification\11_owners_stretched+comm_w_thr{0}-dist+cleaned+loc+class.csv" #r"11_community_classification\11_owners_stretched+comm_w_thr{0}+loc.csv"

## Output paths
GRID_FOLDER = r"00_data\vector\grids"
TABLES_FOLDER_GRID_MW = r"13_ownership_concentration\tables\grids"
FIGURES_FOLDER_GRID_MW = r"13_ownership_concentration\figures\grids"


CONC_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH = r"13_ownership_concentration\tables\grids\mw_conc_meas-grid_{0}km_v{1:02d}-mother_companies-comm_w_thr{2}.csv"
SHARE_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH = r"13_ownership_concentration\tables\grids\mw_share_counts_categories-grid_{0}km_v{1:02d}-mother_companies-comm_w_thr{2}.csv"
COUNT_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH = r"13_ownership_concentration\tables\grids\mw_counts_categories_in_topx-grid_{0}km_v{1:02d}-mother_companies-comm_w_thr{2}.csv"
CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_VERSION_PTH = r"13_ownership_concentration\tables\grids\mw_conc_meas-grid_{0}km_v{1:02d}-mother_companies-comm_w_thr{2}_{3}.csv"

CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_conc_meas-mother_companies-comm_w_thr{0}.csv"
SHARE_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_share_counts_categories-mother_companies-comm_w_thr{0}.csv"
COUNT_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_counts_categories_in_topx-mother_companies-comm_w_thr{0}.csv"
CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_conc_meas-mother_companies-comm_w_thr{0}_{1}.csv"

CONC_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH = r"13_ownership_concentration\tables\grids\mw_conc_meas_privland-grid_{0}km_v{1:02d}-mother_companies-comm_w_thr{2}.csv"
SHARE_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH = r"13_ownership_concentration\tables\grids\mw_share_counts_categories_privland-grid_{0}km_v{1:02d}-mother_companies-comm_w_thr{2}.csv"
COUNT_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH = r"13_ownership_concentration\tables\grids\mw_counts_categories_in_topx_privland-grid_{0}km_v{1:02d}-mother_companies-comm_w_thr{2}.csv"

CONC_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_conc_meas_privland-mother_companies-comm_w_thr{0}.csv"
SHARE_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_share_counts_categories_privland-mother_companies-comm_w_thr{0}.csv"
COUNT_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_counts_categories_in_topx_privland-mother_companies-comm_w_thr{0}.csv"

CONC_OWNER_MERGE_GRID_MW_VERSION_PTH = r"13_ownership_concentration\tables\grids\mw_conc_meas-grid_{0}km_v{1:02d}-owner_merge.csv"
CONC_OWNER_MERGE_GRID_MW_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_conc_meas-owner_merge.csv"
CONC_BTNR_GRID_MW_VERSION_PTH = r"13_ownership_concentration\tables\grids\mw_conc_meas-grid_{0}km_v{1:02d}-BTNR.csv"
CONC_BTNR_GRID_MW_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_conc_meas-BTNR.csv"


#### Processing
def prepare_moving_window_calculation():
    """
    Prepares a dictionary that helps to assign the 4km polygons to the 12km moving windows.
    :return:
    """
    WD = r"C:\Users\IAMO\Documents\work_data\chapter1\ALKIS"
    ALKIS_IACS_GRID_PTH = r"10_alkis_intersection_with_other_layers\alkis_iacs_grid_{0}km_v{1:02d}_inters.shp"

    import json
    import geopandas as gpd
    import os

    os.chdir(WD)

    ## Open alkis iacs grid intersection
    gdf_alk = gpd.read_file(ALKIS_IACS_GRID_PTH.format(4, 1))

    ## open 4km grid that holds information on 12km polyids
    gdf = gpd.read_file(r"00_data\vector\grids\square_grid_4km_v01_with_12km_POLYIDs.shp")

    ## drop 4km polygons that do not intersect with 9 12km polygons
    gdf.dropna(inplace=True)

    ## only use 4km polygons that actually intersect with alkis
    ## (before: get all polyids of the 4km grid that intersect with the alkis data)
    polys_in_alkis = gdf_alk["POLYID"].tolist()
    gdf = gdf.loc[gdf["POLYID"].isin(polys_in_alkis)].copy()

    ## create a dictionary that translates each 4km polygon to a 12km polygon
    grid_v_dict = {}

    ## Loop over all 9 columns that belong to the 9 versions of the 12km polygon grid
    for version in range(1, 10):
        col = f"v{version:02d}_POLYID"

        ## Count the number of 4km polygons in each 12km polygon
        ## and only use the 12km polygons that cover 9 4km polygons that intersect with alkis data
        def count_polyids(poly_group):
            return len(poly_group)
        counts = gdf[[col, "POLYID"]].groupby(col).apply(count_polyids).reset_index()
        counts.columns = [col, "num_4km_polys"]
        counts = counts.loc[counts["num_4km_polys"] >= 9].copy()
        gdf_sub = gdf.loc[gdf[col].isin(counts[col].tolist())].copy()

        ## create a subdictionary translating each 4km polygon to a 12km polygon and pass to out dictionary
        polyids = gdf_sub["POLYID"].tolist()
        polyids_v = gdf_sub[col].tolist()

        sub_dict = {polyid: polyids_v[i] for i, polyid in enumerate(polyids)}
        grid_v_dict[f"v{version:02d}"] = sub_dict

    ## write out
    out_pth = rf"13_ownership_concentration\4km_polyid_to_12km_polyids.json"
    with open(out_pth, 'w') as fp:
        json.dump(grid_v_dict, fp, indent=4)


def wrapper_for_mw_grids(threshold, owner_merge=False, btnr=False):
    """
    Calculates the concentration measures for the 9 versions of the moving windows at the level of the 4km polygons.
    :return:
    """

    print("Calculate concentration measures with moving window over 4km grid.")

    ## Open alkis iacs 4km-grid intersection
    print("\tOpen ALKIS IACS 4km-grid intersection")
    grid_res = 4
    version_4km = 1
    gdf = gpd.read_file(ALKIS_IACS_GRID_PTH.format(grid_res, version_4km))
    gdf["area"] = gdf["geometry"].area

    ## Read owner data
    print(f"\tread owner data with community derived with {threshold}% threshold")
    df_wt = pd.read_csv(OWNERS_W_THRESH_PTH.format(threshold), sep=";")
    df_wt.loc[(df_wt["new_category"].isna()), "new_category"] = 'noagPR'

    ## read 4km to 12km polyid dictionary
    with open(rf"13_ownership_concentration\4km_polyid_to_12km_polyids.json") as json_file:
        grid_v_dict = json.load(json_file)

    ## loop over 12km polyid grid columns
    for version in range(1, 10):
        print("\n\t12km version", version)

        col = f"v{version:02d}_POLYID"
        v_dict = grid_v_dict[f"v{version:02d}"]

        ## translate 4km polyid to 12km polyid
        gdf[col] = gdf["POLYID"].map(v_dict)

        ## drop all parcels that have not 12km polyid assigned
        gdf_curr = gdf.dropna(subset=[col])

        ####################################################################################################################
        print(f"\tCalculate concentration measures with name of mother companies - with {threshold}% threshold")
        print("\tCombine parcels with owner data")
        df_comb_wt = helper_functions.combine_parcels_with_owners(
            gdf_curr[[col, "OGC_FID", "BTNR", "area"]],
            df_wt[["OGC_FID", "mother_company", "owner_merge", "level_c_category", "new_category", "new_category_ext"]])

    #     print("\tOn all land.")
    #     print("\tCalculate concentration measures")
    #     out_pth = CONC_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH.format(grid_res, version, threshold)
    #     print("\tWrite results to", out_pth)
        # conc_meas_lib.calculate_concentration_measures_from_df(
        #     df=df_comb_wt,
        #     target_unit_id_col=col,
        #     area_col='area',
        #     owner_col='mother_company',
        #     out_pth=out_pth)

        # print("\tCalculate shares of categories")
        # out_pth = SHARE_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH.format(grid_res, version, threshold)
        # print("\tWrite results to", out_pth)
        # conc_meas_lib.get_share_of_owner_categories_per_spatial_unit(
        #     df=df_comb_wt,
        #     target_unit_id_col=col,
        #     area_col='area',
        #     owner_col='mother_company',
        #     category_col="new_category",
        #     # category_sub=["1_1_1", "2_9_1"],
        #     out_pth=out_pth)

        # print("\tCount no. of categories in top owners")
        # out_pth = COUNT_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH.format(grid_res, version, threshold)
        # print("\tWrite results to", out_pth)
        # conc_meas_lib.get_count_of_owner_categories_per_spatial_unit(
        #     df=df_comb_wt,
        #     target_unit_id_col=col,
        #     area_col='area',
        #     owner_col='mother_company',
        #     category_col="new_category",
        #     # category_sub=["1_1_1", "2_9_1"],
        #     out_pth=out_pth
        # )

        print("\tLeave a class out AND exclusively for that class.")
        for owner_class in ["noagPR"]: #list(df_comb_wt["new_category"].unique()):
            print("\tLeave", owner_class, "out")
            df_sub = df_comb_wt.loc[df_comb_wt["new_category"] != owner_class].copy()
            conc_meas_lib.calculate_concentration_measures_from_df(
                df=df_sub,
                target_unit_id_col=col,
                area_col='area',
                owner_col='mother_company',
                out_pth=CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_VERSION_PTH.format(grid_res, version, threshold, f"out_{owner_class}")
            )

            print("\tOnly", owner_class)
            df_sub = df_comb_wt.loc[df_comb_wt["new_category"] == owner_class].copy()
            conc_meas_lib.calculate_concentration_measures_from_df(
                df=df_sub,
                target_unit_id_col=col,
                area_col='area',
                owner_col='mother_company',
                out_pth=CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_VERSION_PTH.format(grid_res, version, threshold, owner_class)
            )

        # print("\tOnly on private land.")
        # df_comb_wt_priv = df_comb_wt.loc[~df_comb_wt["level_c_category"].isin(["4_1_1", "4_9_1", "5_1_1", "5_2_1", "5_2_2", "5_2_3", "5_2_4", "5_2_5", "5_2_6", "5_3_1", "5_9_1"])].copy()
        # print("\tCalculate concentration measures")
        # out_pth = CONC_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH.format(grid_res, version, threshold)
        # print("\tWrite results to", out_pth)
        # conc_meas_lib.calculate_concentration_measures_from_df(
        #     df=df_comb_wt_priv,
        #     target_unit_id_col=col,
        #     area_col='area',
        #     owner_col='mother_company',
        #     out_pth=out_pth)

        # print("\tCalculate shares of categories")
        # out_pth = SHARE_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH.format(grid_res, version, threshold)
        # print("\tWrite results to", out_pth)
        # conc_meas_lib.get_share_of_owner_categories_per_spatial_unit(
        #     df=df_comb_wt_priv,
        #     target_unit_id_col=col,
        #     area_col='area',
        #     owner_col='mother_company',
        #     category_col="new_category",
        #     #category_sub=["1_1_1", "2_9_1"],
        #     out_pth=out_pth)

        # print("\tCalculate counts of categories in top owners")
        # out_pth = COUNT_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH.format(grid_res, version, threshold)
        # print("\tWrite results to", out_pth)
        # conc_meas_lib.get_count_of_owner_categories_per_spatial_unit(
        #     df=df_comb_wt_priv,
        #     target_unit_id_col=col,
        #     area_col='area',
        #     owner_col='mother_company',
        #     category_col="level_c_category",
        #     # category_sub=["1_1_1", "2_9_1"],
        #     out_pth=out_pth
        # )

        # ####################################################################################################################
        # if owner_merge:
        #     print("\tCalculate concentration measures with OWNER MERGE")
        #     print("\tUse df from owner with thresh")
        #
        #     print("\tCalculate concentration measures")
        #     out_pth = CONC_OWNER_MERGE_GRID_MW_VERSION_PTH.format(grid_res, version)
        #     print("\tWrite results to", out_pth)
        #     conc_meas_lib.calculate_concentration_measures_from_df(df=df_comb_wt,
        #                                              target_unit_id_col=col,
        #                                              area_col='area',
        #                                              owner_col='owner_merge',
        #                                              out_pth=out_pth)
        #
        # ####################################################################################################################
        # if btnr:
        #     print("\tCalculate concentration measures with BTNR")
        #
        #     print("\tCalculate concentration measures")
        #     out_pth = CONC_BTNR_GRID_MW_VERSION_PTH.format(grid_res, version)
        #     print("\tWrite results to", out_pth)
        #     conc_meas_lib.calculate_concentration_measures_from_df(df=df_comb_wt,
        #                                              target_unit_id_col=col,
        #                                              area_col='area',
        #                                              owner_col='BTNR',
        #                                              out_pth=out_pth)

    # print("\tCount no. of categories in top owners")
    # df_comb_wt = helper_functions.combine_parcels_with_owners(
    #     gdf[["POLYID", "OGC_FID", "BTNR", "area"]],
    #     df_wt[["OGC_FID", "mother_company", "owner_merge", "level_c_category", "new_category", "new_category_ext"]])
    #
    # out_pth = COUNT_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH.format(4, 1, threshold)
    # print("\tWrite results to", out_pth)
    # conc_meas_lib.get_count_of_owner_categories_per_spatial_unit(
    #     df=df_comb_wt,
    #     target_unit_id_col="POLYID",
    #     area_col='area',
    #     owner_col='mother_company',
    #     category_col="new_category",
    #     topx_lst=[1, 2, 3, 4, 5],
    #     out_pth=out_pth
    # )


def combine_mw_grid_results(conc_of_verison_pth, threshold, cols, out_pth, descr=None):
    grid_res = 4
    version_4km = 1

    ## Calculate means
    gdf_grid = gpd.read_file(rf"00_data\vector\grids\square_grid_{grid_res}km_v{version_4km:02d}_with_12km_POLYIDs.shp")

    val_dict = {col: {polyid_4km: [] for polyid_4km in gdf_grid["POLYID"]} for col in cols}

    for version in range(1, 10):
        print(version)

        id_col = f"v{version:02d}_POLYID"
        if not descr:
            pth = conc_of_verison_pth.format(grid_res, version, threshold)
        else:
            pth = conc_of_verison_pth.format(grid_res, version, threshold, descr)

        df = pd.read_csv(pth, sep=',')

        df_comb = pd.merge(gdf_grid[["POLYID", id_col]], df, how="inner", left_on=id_col, right_on="id_sp_unit")

        for col in cols:
            print(col)
            if col in df_comb:

                df_comb_curr = df_comb.dropna(subset=[col])
                polyids = df_comb_curr["POLYID"].tolist()
                values = df_comb_curr[col].tolist()

                for i, polyid in enumerate(polyids):
                    value = values[i]

                    val_dict[col][polyid].append(value)

    val_dict = {col: {polyid: val_dict[col][polyid] for polyid in val_dict[col] if len(val_dict[col][polyid]) >= 3} for col in val_dict}

    res_dict = {col: {polyid: statistics.mean(val_dict[col][polyid]) for polyid in val_dict[col]} for col in cols}

    df_res = pd.DataFrame.from_dict(res_dict, orient="columns").reset_index()
    df_res.rename(columns={"index": "POLYID"}, inplace=True)

    df_res.to_csv(out_pth, sep=',', index=False)


def plotting_wrapper_for_mw_grids(df_res_pth, threshold, descr, categories, df_sh_pth=None, df_ct_pth=None, vector_out_pth=False):
    grid_res = 4
    version_4km = 1
    gdf_grid = gpd.read_file(rf"00_data\vector\grids\square_grid_{grid_res}km_v{version_4km:02d}_with_12km_POLYIDs.shp")

    df_res = pd.read_csv(df_res_pth.format(threshold))
    df_res["total_area"] = df_res["total_area"] / 10000
    if df_sh_pth:
        df_sh = pd.read_csv(df_sh_pth.format(threshold))
        df_res = pd.merge(df_res, df_sh, how="left", on="POLYID")
    if df_ct_pth:
        df_ct = pd.read_csv(df_ct_pth.format(threshold))
        df_res = pd.merge(df_res, df_ct, how="left", left_on="POLYID", right_on="id_sp_unit")

    shp1 = pd.merge(gdf_grid, df_res, how='inner', left_on='POLYID', right_on='POLYID')
    shp1 = shp1.loc[shp1['gini_coeff'] != 0.0]

    if vector_out_pth:
        shp1.to_file(vector_out_pth)

    # shp1["main_owner_cat"] = ""
    # for owner_class in ["PUBLIC", "nCONETW", "aCONETW", "NONPRO", "siCOMP", "a_siCOMP","noagPR", "agriPR", "CHURCH"]:
    #     shp1.loc[shp1[f"count_{owner_class}_top1"] == 1, "main_owner_cat"] = owner_class
    #
    #
    # shp1["main_owner_cat"] = shp1["main_owner_cat"].map(
    #     {"PUBLIC": "PU", "nCONETW": "CN", "aCONETW": "AH", "NONPRO": "NP", "siCOMP": "nSC", "a_siCOMP": "aSC",
    #      "noagPR": "nPR", "agriPR": "aPR", "CHURCH": "RE"})
    # shp1["main_owner_cat"] = pd.Categorical(shp1["main_owner_cat"],
    #                                         categories=["CN", "AH", "PU", "RE", "nPR", "aPR", "aSC", "nSC", "NP"],
    #                                         ordered=True)

    # ## Map with main owner_class
    # colour_dict = {
    #     "aPR": '#f6f739',
    #     "nPR": '#bcbd22',
    #     "aSC": '#c898f5',
    #     "nSC": '#9467bd',
    #     "CN": '#ff7f0e',
    #     "AH": '#fcb97e',
    #     "PU": '#1f77b4',
    #     "NP": '#2ca02c',
    #     "RE": '#7f7f7f'
    # }
    # col = "main_owner_cat"
    # out_pth = rf"{FIGURES_FOLDER_GRID_MW}\maps\mw\mw_map_{col}_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
    # plotting_lib.plot_map_categorical(
    #     shp=shp1,
    #     out_pth=out_pth,
    #     col=col,
    #     colour_dict=colour_dict,
    #     shp2_pth=r"00_data\vector\administrative\BB_Landkreise.shp")
    #
    # col = "cr1"
    # out_pth = rf"{FIGURES_FOLDER_GRID_MW}\graphs\mw\mw_boxplot_main_owner_cat_vs_{col}_mother_companies_w_thr{threshold}{descr}.png"
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
    # for col in ["hhi", "rosenbluth_index", "cr1", "cr3", "cr5", "hhi", "gini_coeff", "palma_v1", "palma_v2"]:
    #     # tr = shp1[col].mean() + shp1[col].std()
    #     # out_pth = rf"{FIGURES_FOLDER_GRID}\maps\mw\mw_map_{col}_gr{round(tr,2)}-mean+std_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
    #     # plotting_lib.plot_map(
    #     #     shp=shp1.loc[shp1[col] > tr],
    #     #     out_pth=out_pth,
    #     #     col=col,
    #     #     shp2_pth=r"00_data\vector\administrative\BB_Landkreise.shp")
    #
        # out_pth = rf"{FIGURES_FOLDER_GRID_MW}\maps\mw\mw_map_{col}_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
        # plotting_lib.plot_map(
        #     shp=shp1,
        #     out_pth=out_pth,
        #     col=col,
        #     shp2_pth=r"00_data\vector\administrative\BB_Landkreise.shp",
        #     extremes=True
        # )

    # out_pth = rf"{FIGURES_FOLDER_GRID_MW}\maps\mw\mw_map_comp_conc_meas_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
    # plotting_lib.plot_maps_in_grid(
    #     shp=shp1,
    #     out_pth=out_pth,
    #     cols=["gini_coeff", "hhi", "palma_v2", "cr1", "cr3", "cr5"],
    #     titles=["a)", "b)", "c)", "d)", "e)", "f)"],
    #     labels=["Gini coefficient", "HHI", "T100/B90", "CR1", "CR3", "CR5"],
    #     nrow=2,
    #     shp2_pth=r"00_data\vector\administrative\BB_Landkreise.shp"
    # )

    ## Maps of shares per category
    ## Histograms of share per category
    ## Scatterplots of shares per category vs palma indeces
    # for cat in categories:
    #     for share in [1, 10, 100]:
    #         col1 = f"share_{cat}_sh{share}"
    #         out_pth = rf"{FIGURES_FOLDER_GRID_MW}\maps\mw\occ_cats\mw_map_{col1}_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
    #         plotting_lib.plot_map(
    #             shp=shp1,
    #             out_pth=out_pth,
    #             col=col1, shp2_pth=r"00_data\vector\administrative\BB_Landkreise.shp")
    #         out_pth = f"{FIGURES_FOLDER_GRID_MW}\graphs\mw\mw_histogramm_{col1}_mother_companies_w_thr{threshold}{descr}.png"
    #         plotting_lib.histogramm(
    #             df=shp1,
    #             col=col1,
    #             out_pth=out_pth,
    #             x_label=f"{col1} in 4km mov.wind. grid",
    #             y_label="Number of 4km grid cells")
    #
    #         for col2 in ["palma_v2", "palma_v1"]:
    #             if col1 in shp1.columns:
    #                 out_pth = rf"{FIGURES_FOLDER_GRID_MW}\graphs\mw\influence_share_cats\{col2}\share{share}\mw_scatter_{col1}_vs_{col2}_mother_companies_w_thr{threshold}{descr}.png"
    #                 plotting_lib.scatterplot_two_columns(
    #                     df=shp1,
    #                     col1=col1,
    #                     col2=col2,
    #                     out_pth=out_pth,
    #                     xminmax=(0, .5),
    #                     x_label=col1,
    #                     y_label=f"{col2}")
    #
    #                 out_pth = rf"{FIGURES_FOLDER_GRID_MW}\graphs\mw\influence_share_cats\{col2}\share{share}\mw_scatter_{col1}_vs_{col2}_mother_companies_w_thr{threshold}_xlog{descr}.png"
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
    # for x in [1, 3, 5]:
    #     for cat in categories:
    #         col1 = f"count_{cat}_top{x}"
    #         out_pth=rf"{FIGURES_FOLDER_GRID_MW}\maps\mw\occ_cats\mw_map_{col1}_grid_mean_mother_companies_w_thr{threshold}{descr}.png"
    #         plotting_lib.plot_map(
    #             shp=shp1,
    #             out_pth=out_pth,
    #             col=col1, shp2_pth=r"00_data\vector\administrative\BB_Landkreise.shp")
    #
    #         col2 = f"cr{x}"
    #         if col1 in shp1.columns:
    #             out_pth = rf"{FIGURES_FOLDER_GRID_MW}\graphs\mw\influence_share_cats\{col2}\mw_scatter_{col1}_vs_{col2}_mother_companies_w_thr{threshold}{descr}.png"
    #             plotting_lib.scatterplot_two_columns(df=shp1, col1=col1, col2=col2, out_pth=out_pth,
    #                                     x_label=col1, y_label=f"{col2}")
    #         else:
    #             print(col1, "not in df")

    ## Histograms of concentration measures
    # cols = ["cr1", "cr3", "cr5", "hhi", "total_area", "num_owners", "lac", "palma_v1", "palma_v2", "gini_coeff",
    #         "share_p100", "share_p95_99", "share_v19", "share_m50", "share_b40"]
    # for col in cols:
    #     out_pth = f"{FIGURES_FOLDER_GRID_MW}\graphs\mw\mw_histogramm_{col}_mother_companies_w_thr{threshold}{descr}.png"
    #     plotting_lib.histogramm(
    #         df=shp1,
    #         col=col,
    #         out_pth=out_pth,
    #         x_label=f"{col} in 4km mov.wind. grid",
    #         y_label="Number of 4km grid cells")

    ## Scatterplots of concentration measures vs total area
    # cols = ["gini_coeff", "cr1", "cr3", "cr5", "hhi", "rosenbluth_index", "palma_v1"]
    # for col in cols:
    #     out_pth = f"{FIGURES_FOLDER_GRID_MW}\graphs\mw\mw_scatter_{col}_vs_total_area_mother_companies_w_thr{threshold}{descr}.png"
    #     plotting_lib.scatterplot_two_columns(
    #         df=shp1,
    #         col1=col,
    #         col2="total_area",
    #         out_pth=out_pth,
    #         x_label=f"{col} - moving window",
    #         y_label="mean agric. area in 12km mw [ha]"
    #     )



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

    # prepare_moving_window_calculation()
    # moving_window_concentration_measure_calculation(
    #     threshold=threshold,
    #     owner_merge=True
    # )

    ## Calculate mean of moving window calculations of concentration measures
    # cols = ["gini_coeff", "cr1", "cr3", "cr5", "total_area", "num_owners"]
    cols = ["gini_coeff", "cr1", "cr3", "cr5", "hhi", "total_area", "num_owners", "lac", "palma_v1", "palma_v2",
            "rosenbluth_index", "share_p100", "share_p95_99", "share_v19", "share_m50", "share_b40"]

    # combine_mw_grid_results(
    #     conc_of_verison_pth=CONC_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH,
    #     threshold=threshold,
    #     cols=cols,
    #     out_pth=CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold))
    # combine_mw_grid_results(
    #     conc_of_verison_pth=CONC_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH,
    #     threshold=threshold,
    #     cols=cols,
    #     out_pth=CONC_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH.format(threshold))
    # combine_mw_grid_results(
    #     conc_of_verison_pth=CONC_OWNER_MERGE_GRID_MW_VERSION_PTH,
    #     threshold=threshold,
    #     cols=cols,
    #     out_pth=CONC_OWNER_MERGE_GRID_MW_COMBINED_PTH)
    # combine_mw_grid_results(
    #     conc_of_verison_pth=CONC_BTNR_GRID_MW_VERSION_PTH,
    #     threshold=threshold,
    #     cols=cols,
    #     out_pth=CONC_BTNR_GRID_MW_COMBINED_PTH)

    ## Calculate mean of moving window calculations of share per category calculations
    # cols = pd.read_csv(SHARE_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH.format(4, 1, threshold))
    # cols = list(cols.columns)
    # cols.remove("id_sp_unit")
    # combine_mw_grid_results(
    #     conc_of_verison_pth=SHARE_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH,
    #     threshold=threshold,
    #     cols=cols,
    #     out_pth=SHARE_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH.format(threshold))

    # cols = pd.read_csv(COUNT_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH.format(4, 1, threshold))
    # cols = list(cols.columns)
    # cols.remove("id_sp_unit")
    # # f1 = lambda x, y: f"count_{x}_top{y}"
    # # cols = [f1(cat, topx) for cat in ['1_1_1', '2_9_1'] for topx in [1, 3, 5]]
    # combine_mw_grid_results(
    #     conc_of_verison_pth=COUNT_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_VERSION_PTH,
    #     threshold=threshold,
    #     cols=cols,
    #     out_pth=COUNT_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH.format(threshold))

    # for owner_class in ["noagPR"]:#["PUBLIC", "nCONETW", "aCONETW", "NONPRO", "siCOMP", "a_siCOMP", "noagPR", "agriPR", "CHURCH"]:
    #     combine_mw_grid_results(
    #         conc_of_verison_pth=CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_VERSION_PTH,
    #         threshold=threshold,
    #         cols=cols,
    #         out_pth=CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold, f"{owner_class}"),
    #         descr=f"{owner_class}"
    #     )
    #
    #     combine_mw_grid_results(
    #         conc_of_verison_pth=CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_VERSION_PTH,
    #         threshold=threshold,
    #         cols=cols,
    #         out_pth=CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold, f"out_{owner_class}"),
    #         descr=f"out_{owner_class}"
    #     )
    #
    #     plotting_wrapper_for_mw_grids(
    #         threshold=50,
    #         df_res_pth=CONC_MOTHER_COMP_SUBGROUPS_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold, f"{owner_class}"),
    #         # df_sh_pth=SHARE_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH,
    #         # df_ct_pth=COUNT_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH.format(4, 1, threshold),
    #         categories=["PUBLIC", "CONETW", "noagPR", "NONPRO", "siCOMP", "agriPR", "CHURCH", "COOPER", "BVVG"],
    #         descr=owner_class
    #     )

    ## Plotting
    # categories = ['1_1_1', '1_9_1', '2_9_1', '2_1_4', '2_1_2', '5_2_2', '2_1_5', '5_2_1', '3_1_1',
    #               '2_9_2', '5_2_5', '3_2_2', '5_2_4', '4_1_1', '5_3_1', '2_1_8', '2_2_2', '2_1_7',
    #               '1_2_1', '5_9_1', '5_1_1', '1_2_3', '2_1_1', '90', '5_2_3', '1_2_2', '3_2_3']
    # categories = ['2_2_1', '2_1_3', '3_9_1', '2_1_6', '1_2_4', '3_1_2', '1_1_2']

    vector_out_pth = rf"13_ownership_concentration\vector\grids\{os.path.basename(CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold))[:-3]}shp"

    plotting_wrapper_for_mw_grids(
        threshold=50,
        df_res_pth=CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold),
        # df_sh_pth=SHARE_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH,
        df_ct_pth=COUNT_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_VERSION_PTH.format(4, 1, threshold),
        categories=["PUBLIC", "CONETW", "noagPR", "NONPRO", "siCOMP", "agriPR", "CHURCH", "COOPER", "BVVG"],
        descr="",
        vector_out_pth=vector_out_pth
    )


    plotting_wrapper_for_mw_grids(
        threshold=50,
        df_res_pth=CONC_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH,
        # df_sh_pth=SHARE_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH,
        # df_ct_pth=COUNT_CATEG_MOTHER_COMP_W_THRESH_PRIVLAND_GRID_MW_COMBINED_PTH,
        categories=["CONETW", "noagPR", "NONPRO", "siCOMP", "agriPR", "COOPER"],
        descr="_privland"
    )

    ## Explain concentration measures
    df2 = pd.read_csv(SHARE_CATEG_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold))
    df1 = pd.read_csv(CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(threshold))
    df = pd.merge(df1, df2, how="left", on="POLYID")

    xs_col = ["share_PUBLIC_sh100", "share_CONETW_sh100", "share_noagPR_sh100", "share_NONPRO_sh100",
              "share_siCOMP_sh100", "share_agriPR_sh100", "share_CHURCH_sh100", "share_COOPER_sh100",
              "share_BVVG_sh100",
              "count_PUBLIC_sh100", "count_CONETW_sh100", "count_noagPR_sh100", "count_NONPRO_sh100",
              "count_siCOMP_sh100", "count_agriPR_sh100", "count_CHURCH_sh100", "count_COOPER_sh100",
              "count_BVVG_sh100"]
    xs_col = ["count_PUBLIC_sh100", "count_CONETW_sh100", "count_noagPR_sh100", "count_NONPRO_sh100",
              "count_siCOMP_sh100", "count_agriPR_sh100", "count_CHURCH_sh100", "count_COOPER_sh100",
              "count_BVVG_sh100"]
    xs_col = ["share_PUBLIC_sh100", "share_CONETW_sh100", "share_noagPR_sh100", "share_NONPRO_sh100",
              "share_siCOMP_sh100", "share_agriPR_sh100", "share_CHURCH_sh100", "share_COOPER_sh100",
              "share_BVVG_sh100"]
    xs_col = ["share_CONETW_sh100", "share_noagPR_sh100", "share_NONPRO_sh100",
              "share_siCOMP_sh100", "share_agriPR_sh100", "share_COOPER_sh100"]
    plot_corr(df[xs_col])
    y_col = "palma_v2"
    tr = df[y_col].mean() + df[y_col].std()
    explain_concentration_with_regression(
        df=df.loc[df[y_col] > tr],
        y_col=y_col,
        xs_col=xs_col)

    xs_col = ["share_CONETW_sh100", "share_noagPR_sh100", "share_siCOMP_sh100", "share_agriPR_sh100",
              "share_COOPER_sh100"]
    y_col = "gini_coeff"
    tr = df[y_col].mean() + df[y_col].std()
    explain_concentration_with_regression(
        df=df.loc[df[y_col] > tr],
        y_col=y_col,
        xs_col=xs_col)

    xs_col = ["share_CONETW_sh100", "share_noagPR_sh100", "share_siCOMP_sh100", "share_agriPR_sh100",
              "share_COOPER_sh100"]
    y_col = "hhi"
    tr = df[y_col].mean() + df[y_col].std()
    explain_concentration_with_regression(
        df=df.loc[df[y_col] > tr],
        y_col=y_col,
        xs_col=xs_col)

    explore_concentration_at_mw_grid(pth_conc_owner_merge=CONC_OWNER_MERGE_GRID_MW_COMBINED_PTH,
                                     pth_conc_mother_company=CONC_MOTHER_COMP_W_THRESH_GRID_MW_COMBINED_PTH.format(
                                         threshold),
                                     threshold=threshold,
                                     out_tables_folder=TABLES_FOLDER_GRID_MW)

    out_pth_diff_vals = fr"{TABLES_FOLDER_GRID_MW}\conc_meas-mw_grid-diff_mcomp_w_thr{threshold}-owner_merge.csv"
    shp_pth = rf"00_data\vector\grids\square_grid_4km_v01_with_12km_POLYIDs.shp"
    plot_concentration_differences_mw(df_diff_pth=out_pth_diff_vals,
                                      shp_pth=shp_pth,
                                      figures_folder=FIGURES_FOLDER_GRID_MW,
                                      threshold=threshold)


    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)


if __name__ == '__main__':
    main()