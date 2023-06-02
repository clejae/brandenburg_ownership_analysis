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
GRID_4km_WITH_12KM_IDS_PTH = r"00_data\vector\grids\square_grid_4km_v01_with_12km_POLYIDs.shp"

## Output
DICT_4KM_TO_12KM_IDS = rf"13_ownership_concentration\4km_polyid_to_12km_polyids.json"
GRID_FOLDER = r"00_data\vector\grids"
TABLES_FOLDER_GRID_MW = r"13_ownership_concentration\grids"

CONC_MEASURES_MW_GRID_VERSIONS_PTH = r"13_ownership_concentration\tables\grids\mw_conc_meas-grid_{0}km_v{1:02d}-{2}.csv"
SHARE_CATEG_MW_GRID_VERSIONS_PTH = r"13_ownership_concentration\tables\grids\mw_share_counts_categories-grid_{0}km_v{1:02d}-{2}.csv"
COUNT_CATEG_MW_GRID_VERSIONS_PTH = r"13_ownership_concentration\tables\grids\mw_counts_categories_in_topx-grid_{0}km_v{1:02d}-{2}.csv"

CONC_MEASURES_MW_GRID_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_conc_meas-{0}.csv"
SHARE_CATEG_MW_GRID_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_share_counts_categories-{0}.csv"
COUNT_CATEG_MW_GRID_COMBINED_PTH = r"13_ownership_concentration\tables\grids\mw_mean_counts_categories_in_topx-{0}.csv"



# ------------------------------------------ PROCESSING ------------------------------------------#

def prepare_moving_window_calculation(gdf_alk, grid_4km_with_12km_ids_pth, out_pth):
    """
    Prepares a dictionary that helps to assign the 4km polygons to the 12km moving windows.
    Args:
        alkis_grid_iacs_pth: Path to ALKIS-Grid-IACS intersection/union shapefile.
        grid_4km_with_12km_ids_pth: Path to 4x4km grid shapefile with 12km IDs
        out_pth: Output path to dictionary (json) translating 4km IDs to 12km IDs.

    Returns:

    """

    ## open 4km grid that holds information on 12km polyids
    gdf = gpd.read_file(grid_4km_with_12km_ids_pth)

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
    with open(out_pth, 'w') as fp:
        json.dump(grid_v_dict, fp, indent=4)


def combine_mw_grid_results(grid_4km_with_12km_ids_pth, conc_of_grid_version_pth, threshold, cols, out_pth, descr=None):
    grid_res = 4
    version_4km = 1

    ## Calculate means
    gdf_grid = gpd.read_file(grid_4km_with_12km_ids_pth)

    val_dict = {col: {polyid_4km: [] for polyid_4km in gdf_grid["POLYID"]} for col in cols}

    for version in range(1, 10):
        print(version)

        id_col = f"v{version:02d}_POLYID"
        if not descr:
            pth = conc_of_grid_version_pth.format(grid_res, version, threshold)
        else:
            pth = conc_of_grid_version_pth.format(grid_res, version, threshold, descr)

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


def moving_window_concentration_measure_calculation(parcels, owner_df, dict_4km_to_12km_ids_pth, threshold,
                                                    target_grid_res, descr, owner_col,
                                                    conc_measures_at_mw_grid_versions_mcomp_pth,
                                                    conc_measures_at_mw_grid_combined_mcomp_pth,
                                                    grid_4km_with_12km_ids_pth,
                                                    conc_cols,
                                                    shares_grid_versions_pth=None,
                                                    shares_grid_combined_pth=None,
                                                    counts_grid_versions_pth=None,
                                                    counts_grod_combined_pth=None
                                                    ):
    """
    Calculates the concentration measures for the 9 versions of the moving windows at the level of the 4km polygons.
    :return:
    """

    print("Calculate concentration measures with moving window over 4km grid.")

    ## read 4km to 12km polyid dictionary
    with open(dict_4km_to_12km_ids_pth) as json_file:
        grid_v_dict = json.load(json_file)

    ## loop over 12km polyid grid columns
    for version in range(1, 10):
        print("\n\t12km version", version)

        col = f"v{version:02d}_POLYID"
        v_dict = grid_v_dict[f"v{version:02d}"]

        ## translate 4km polyid to 12km polyid
        parcels[col] = parcels["POLYID"].map(v_dict)

        ## drop all parcels that have not 12km polyid assigned
        gdf_curr = parcels.dropna(subset=[col])

        ####################################################################################################################
        print(f"\tCalculate concentration measures with name of mother companies - with {threshold}% threshold")
        print("\tCombine parcels with owner data")
        df_comb_wt = helper_functions.combine_parcels_with_owners(
            gdf_curr[[col, "OGC_FID", "BTNR", "area"]],
            owner_df[["OGC_FID", "mother_company", "owner_merge", "level_c_category", "new_category", "new_category_ext"]])

        print("\tCalculation")
        out_pth = conc_measures_at_mw_grid_versions_mcomp_pth.format(target_grid_res, version, descr)
        print("\tWrite results to", out_pth)
        conc_meas_lib.calculate_concentration_measures_from_df(
            df=df_comb_wt,
            target_unit_id_col=col,
            area_col='area',
            owner_col=owner_col,
            out_pth=out_pth)

        if shares_grid_versions_pth:
            print("\tCalculate shares of categories")
            out_pth = shares_grid_versions_pth.format(target_grid_res, version, descr)
            print("\tWrite results to", shares_grid_versions_pth)
            conc_meas_lib.get_share_of_owner_categories_per_spatial_unit(
                df=df_comb_wt,
                target_unit_id_col=col,
                area_col='area',
                owner_col='mother_company',
                category_col=owner_col,
                out_pth=shares_grid_versions_pth)

        if counts_grid_versions_pth:
            print("\tCount no. of categories in top owners")
            out_pth = counts_grid_versions_pth.format(target_grid_res, version, descr)
            print("\tWrite results to", out_pth)
            conc_meas_lib.get_count_of_owner_categories_per_spatial_unit(
                df=df_comb_wt,
                target_unit_id_col=col,
                area_col='area',
                owner_col=owner_col,
                category_col="new_category",
                out_pth=out_pth
            )

    combine_mw_grid_results(
        conc_of_grid_version_pth=conc_measures_at_mw_grid_versions_mcomp_pth,
        grid_4km_with_12km_ids_pth=grid_4km_with_12km_ids_pth,
        threshold=threshold,
        cols=conc_cols,
        out_pth=conc_measures_at_mw_grid_combined_mcomp_pth.format(descr))

    if shares_grid_versions_pth:
        cols = pd.read_csv(shares_grid_versions_pth.format(4, 1, descr))
        cols = list(cols.columns)
        cols.remove("id_sp_unit")
        combine_mw_grid_results(
            conc_of_grid_version_pth=shares_grid_versions_pth,
            grid_4km_with_12km_ids_pth=grid_4km_with_12km_ids_pth,
            threshold=threshold,
            cols=cols,
            out_pth=shares_grid_combined_pth.format(descr))

    if counts_grid_versions_pth:
        cols = pd.read_csv(counts_grid_versions_pth.format(4, 1, descr))
        cols = list(cols.columns)
        cols.remove("id_sp_unit")
        # # # f1 = lambda x, y: f"count_{x}_top{y}"
        # # # cols = [f1(cat, topx) for cat in ['1_1_1', '2_9_1'] for topx in [1, 3, 5]]
        combine_mw_grid_results(
            conc_of_grid_version_pth=counts_grid_versions_pth,
            grid_4km_with_12km_ids_pth=grid_4km_with_12km_ids_pth,
            threshold=threshold,
            cols=cols,
            out_pth=counts_grod_combined_pth.format(descr))


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)

    os.chdir(WD)
    threshold = 50

    #################################### MOVING WINDOW CONCENTRATION CALCULATION ####################################
    helper_functions.create_folder(TABLES_FOLDER_GRID_MW)

    ## Open alkis iacs 4km-grid intersection
    print("\tOpen ALKIS IACS 4km-grid intersection")
    grid_res = 4
    version_4km = 1
    parcels = gpd.read_file(ALKIS_IACS_GRID_PTH.format(grid_res, version_4km))
    parcels["area"] = parcels["geometry"].area

    prepare_moving_window_calculation(
        gdf_alk=parcels.loc[parcels["BTNR"].notna()],
        grid_4km_with_12km_ids_pth=GRID_4km_WITH_12KM_IDS_PTH,
        out_pth=DICT_4KM_TO_12KM_IDS
    )

    ## Read owner data
    print(f"\tRead owner data with community derived with {threshold}% threshold")
    owner_df = pd.read_csv(OWNERS_W_THRESH_PTH.format(threshold), sep=";")
    owner_df.loc[(owner_df["new_category"].isna()), "new_category"] = 'noagPR'

    ## On all IACS areas for company networks
    cols = ["gini_coeff", "cr1", "cr3", "cr5", "hhi", "total_area", "num_owners", "lac", "palma_v1", "palma_v2",
            "rosenbluth_index", "share_p100", "share_p95_99", "share_v19", "share_m50", "share_b40"]

    moving_window_concentration_measure_calculation(
        parcels=parcels.loc[parcels["BTNR"].notna()],
        owner_df=owner_df,
        dict_4km_to_12km_ids_pth=DICT_4KM_TO_12KM_IDS,
        owner_col='mother_company',
        threshold=threshold,
        target_grid_res=grid_res,
        descr=f"mother_companies-comm_w_thr{threshold}-iacs_areas",
        conc_measures_at_mw_grid_versions_mcomp_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH,
        conc_measures_at_mw_grid_combined_mcomp_pth=CONC_MEASURES_MW_GRID_COMBINED_PTH,
        grid_4km_with_12km_ids_pth=GRID_4km_WITH_12KM_IDS_PTH,
        conc_cols=cols,
        counts_grid_versions_pth=COUNT_CATEG_MW_GRID_VERSIONS_PTH,
        counts_grod_combined_pth=COUNT_CATEG_MW_GRID_COMBINED_PTH)

    ## On all IACS areas for single_owners
    cols = ["gini_coeff", "cr1", "cr3", "cr5", "hhi", "total_area", "num_owners", "lac", "palma_v1", "palma_v2",
            "rosenbluth_index", "share_p100", "share_p95_99", "share_v19", "share_m50", "share_b40"]

    moving_window_concentration_measure_calculation(
        parcels=parcels.loc[parcels["BTNR"].notna()],
        owner_df=owner_df,
        dict_4km_to_12km_ids_pth=DICT_4KM_TO_12KM_IDS,
        owner_col='owner_merge',
        threshold=threshold,
        target_grid_res=grid_res,
        descr="owner_merge-iacs_areas",
        conc_measures_at_mw_grid_versions_mcomp_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH,
        conc_measures_at_mw_grid_combined_mcomp_pth=CONC_MEASURES_MW_GRID_COMBINED_PTH,
        grid_4km_with_12km_ids_pth=GRID_4km_WITH_12KM_IDS_PTH,
        conc_cols=cols)

    # ## On private land
    # private_categories = ["4_1_1", "4_9_1", "5_1_1", "5_2_1", "5_2_2", "5_2_3", "5_2_4", "5_2_5", "5_2_6", "5_3_1",
    #                       "5_9_1"]
    # moving_window_concentration_measure_calculation(
    #     parcels=parcels.loc[parcels["BTNR"].nona()],
    #     owner_df=owner_df.loc[~owner_df["level_c_category"].isin(private_categories)],
    #     dict_4km_to_12km_ids_pth=DICT_4KM_TO_12KM_IDS,
    #     owner_col='mother_company',
    #     threshold=threshold,
    #     target_grid_res=grid_res,
    #     descr=f"mother_companies-comm_w_thr{threshold}-private_land",
    #     conc_measures_at_mw_grid_versions_mcomp_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH,
    #     conc_measures_at_mw_grid_combined_mcomp_pth=CONC_MEASURES_MW_GRID_COMBINED_PTH,
    #     grid_4km_with_12km_ids_pth=GRID_4km_WITH_12KM_IDS_PTH,
    #     conc_cols=cols)
    #
    # ## On all IACS areas for company networks but only for certain owner classes
    # owner_classes = ["PUBLIC", "nCONETW", "aCONETW", "NONPRO", "siCOMP", "a_siCOMP", "noagPR", "agriPR", "CHURCH"]
    #
    # for owner_class in owner_classes:
    #     moving_window_concentration_measure_calculation(
    #         parcels=parcels.loc[parcels["BTNR"].nona()],
    #         owner_df=owner_df.loc[owner_df["new_category"] == owner_class],
    #         dict_4km_to_12km_ids_pth=DICT_4KM_TO_12KM_IDS,
    #         owner_col='mother_company',
    #         threshold=threshold,
    #         target_grid_res=grid_res,
    #         descr=f"mother_companies-comm_w_thr{threshold}-only_{owner_class}",
    #         conc_measures_at_mw_grid_versions_mcomp_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH,
    #         conc_measures_at_mw_grid_combined_mcomp_pth=CONC_MEASURES_MW_GRID_COMBINED_PTH,
    #         grid_4km_with_12km_ids_pth=GRID_4km_WITH_12KM_IDS_PTH,
    #         conc_cols=cols)
    #
    #     moving_window_concentration_measure_calculation(
    #         parcels=parcels.loc[parcels["BTNR"].nona()],
    #         owner_df=owner_df.loc[owner_df["new_category"] != owner_class],
    #         dict_4km_to_12km_ids_pth=DICT_4KM_TO_12KM_IDS,
    #         owner_col='mother_company',
    #         threshold=threshold,
    #         target_grid_res=grid_res,
    #         descr=f"mother_companies-comm_w_thr{threshold}-not_{owner_class}",
    #         conc_measures_at_mw_grid_versions_mcomp_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH,
    #         conc_measures_at_mw_grid_combined_mcomp_pth=CONC_MEASURES_MW_GRID_COMBINED_PTH,
    #         grid_4km_with_12km_ids_pth=GRID_4km_WITH_12KM_IDS_PTH,
    #         conc_cols=cols)

    #################################### STATE LEVEL ####################################

    #################################### STATE LEVEL ####################################

    #################################### FARM BUFFERS ####################################

    #################################### MUNICIPALITIES ####################################


    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)


if __name__ == '__main__':
    main()