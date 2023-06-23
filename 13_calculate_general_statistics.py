# Author:
# github repository:

## ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import time
import os
import geopandas as gpd
import numpy as np
import pandas as pd
import os
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter

## Project library
import helper_functions
import plotting_lib
## ------------------------------------------ USER INPUT ------------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input
## A) IACS, original ALKIS, IACS-ALKIS intersection
IACS_PTH = r"00_data\vector\IACS\IACS_BB_2020.shp"
ALKIS_ORIG_PTH = r"00_data\vector\ALKIS\v_eigentuemer_bb_reduced.shp"
ALK_IACS_INTERS_PTH = r"09_alkis_intersection_with_other_layers\alkis_iacs_inters.shp"
BB_SHP = r"00_data\vector\administrative\BB_municipalities.shp"
## B) ALKIS with community information
OWNERS_W_THRESH_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-dist+cleaned+loc+class.csv"
COMMUNITY_INFO_DICT_W_THRESH = r"08_network_analysis\owners+communities_thr{0}\08_alkis_owners_comm_thr{0}_info_dict.json"

## Output
## A) Area calculation on original data
IACS_STATISTICS_OUT_PTH = r"13_general_statistics\13_area_information_iacs.xlsx"
ALKIS_OUT_PTH = r"13_general_statistics\13_area_information_alkis+comm_thr{0}_iacs_intersection.xlsx"
ALKIS_ORIG_OUT_PTH = r"13_general_statistics\13_area_information_alkis_original.csv"

## B) Area calculation per owner or community
AREA_PER_OWN_PTH = rf"13_general_statistics\ALKIS_IACS_intersection\area_per_owner_merge_alkis_iacs.csv"

AREA_PER_OWN_PTH_EXCL_PUBL_CHURCH_PTH = rf"13_general_statistics\ALKIS_IACS_intersection\area_per_owner_merge_alkis_iacs_excl_publ_church.csv"
AREA_PER_COMM_W_THRESH = r"13_general_statistics\ALKIS_IACS_intersection\area_per_community_w_thr{0}_alkis_iacs.csv"
AREA_PER_COMM_W_THRESH_EXCL_PUBL_CHURCH_PTH = r"13_general_statistics\ALKIS_IACS_intersection\area_per_community_w_thr{0}_alkis_iacs_excl_publ_church.csv"
# C) Parcel distribution of selected owners
OUT_SHP_FOLDER_W_THRESH = r"13_general_statistics\community_shapes\w_thr{0}"
OUT_FIG_FOLDER_W_THRESH = r"13_general_statistics\figures\community_distribution\w_thr{0}"


## ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#
def calculate_area_per_category(df, category_col, area_col):
    """
    Calculate the area per category in the given DataFrame.

    Args:
        df (DataFrame): Input DataFrame.
        category_col (str): Column name representing the category.
        area_col (str): Column name representing the area.

    Returns:
        DataFrame: DataFrame with calculated area per category.
    """
    print(f"Start calculating area per {category_col}")

    # Calculate area per category
    df_summ = df[[category_col, area_col]].groupby(category_col).sum().reset_index()

    # Insert total row
    total_val = df_summ['area'].sum()
    total = df_summ.apply(np.sum)
    total[category_col] = 'Total'
    df_summ = df_summ.append(pd.DataFrame(total.values, index=total.keys()).T, ignore_index=True)
    df_summ.sort_values(by="area", ascending=False, inplace=True)

    # Calculate shares of categories
    df_summ["share"] = (df_summ["area"] / total_val) * 100

    return df_summ


def calculate_iacs_area_statistics_wrapper(iacs_pth, iacs_statistics_out_pth):
    """
    Generate area statistics for IACS data.

    Args:
        iacs_pth (str): File path to the IACS data.
        iacs_statistics_out_pth (str): File path for the output Excel file.

    Returns:
        None
    """
    print("Generate area statistics for IACS data")

    # Read IACS data and calculate area
    shp_iacs = gpd.read_file(iacs_pth)  # Read IACS data
    shp_iacs['area'] = shp_iacs['geometry'].area  # Calculate area of each geometry
    shp_iacs['area'] = shp_iacs['area'] / 10000  # Convert area to hectares

    # Calculate area per TF_TYP category
    df_tftyp = calculate_area_per_category(
        df=shp_iacs,
        category_col="TF_TYP",
        area_col="area"
    )

    # Calculate area per ID_KTYP category
    df_ktyp = calculate_area_per_category(
        df=shp_iacs,
        category_col="ID_KTYP",
        area_col="area"
    )

    # Calculate area per Oeko category
    df_oeko = calculate_area_per_category(
        df=shp_iacs,
        category_col="Oeko",
        area_col="area"
    )

    # Write the dataframes to separate sheets in an Excel file
    writer = pd.ExcelWriter(iacs_statistics_out_pth)
    for c, df in enumerate([df_tftyp, df_ktyp, df_oeko]):
        name = ['IACS_TFTYP', 'IACS_KTYP', 'IACS_OEKO'][c]
        df.to_excel(writer, sheet_name=name, index=False)
    writer.save()

    print("Done.")


def calculate_alkis_area_statistics_wrapper(alkis_iacs_inters_pth, owners_pth, out_pth):
    """
    Generate area statistics for ALKIS data.

    Args:
        alkis_iacs_inters_pth (str): File path to the ALKIS-IACS intersections data.
        owners_pth (str): File path to the owners data.
        out_pth (str): File path for the output Excel file.

    Returns:
        None
    """
    print("Generate area statistics for ALKIS data")
    print("Combine parcels with owner data")

    # Read ALKIS-IACS intersections data and calculate area
    gdf = gpd.read_file(alkis_iacs_inters_pth)  # Read ALKIS-IACS intersections data
    gdf["area"] = gdf["geometry"].area  # Calculate area of each geometry
    gdf['area'] = gdf['area'] / 10000  # Convert area to hectares

    # Read owners data
    df = helper_functions.read_table_to_df(owners_pth)

    # Combine parcels with owner data
    df_comb = helper_functions.combine_parcels_with_owners(
        gdf[["OGC_FID", "area"]],
        df[["OGC_FID", "mother_company", "owner_merge"]]
    )

    # Calculate area per owner merge category
    df_omerge = calculate_area_per_category(
        df=df_comb,
        category_col="owner_merge",
        area_col="area"
    )

    # Calculate area per mother company category
    df_mcomp = calculate_area_per_category(
        df=df_comb,
        category_col="mother_company",
        area_col="area"
    )

    # Write the dataframes to separate sheets in an Excel file
    writer = pd.ExcelWriter(out_pth)
    for c, df in enumerate([df_omerge, df_mcomp]):
        name = ['ALKIS_OWNER_MERGE', 'ALKIS_MOTHER_COMPANY'][c]
        df.to_excel(writer, sheet_name=name, index=False)
    writer.save()

    print("Done.")


def create_overall_stats_per_level(input_df, level_col, area_col, fid_col, owners_col):
    """
    Calculate overall statistics per level in the input DataFrame.

    Args:
        input_df (DataFrame): Input DataFrame containing ALKIS data.
        level_col (str): Column name representing the levels.
        area_col (str): Column name representing the area.
        fid_col (str): Column name representing the feature ID.
        owners_col (str): Column name representing the owners.

    Returns:
        DataFrame: DataFrame with overall statistics per level.
    """
    out_df = []

    # Get unique levels in the input DataFrame
    uni_levels = input_df[level_col].unique()

    # Calculate statistics for each level
    for lev in uni_levels:
        sub = input_df[input_df[level_col] == lev].copy()  # Create a subset for the specific level

        # Calculate the number of unique owners
        num_owners = len(sub[owners_col].unique())

        # Calculate the mean area per owner
        mean_area = round(sub[[area_col, owners_col]].groupby(owners_col).sum().mean().iloc[0], 1)

        # Calculate the mean number of land parcels per owner
        mean_num = sub[[fid_col, owners_col]].groupby(owners_col).count().mean().iloc[0]

        # Calculate the total area for the level
        tot_area = round(sub[area_col].sum(), 1)

        # Append the statistics to the output DataFrame
        out_df.append([lev, num_owners, mean_area, mean_num, tot_area])

    # Convert the output list to a DataFrame
    out_df = pd.DataFrame(out_df)
    out_df.columns = ['Category', 'Number of owners', 'Mean area per owner [ha]', 'Mean number of land parcels per owner', 'Total [ha]']
    out_df = out_df.sort_values(by=['Category'])

    return out_df


def calculate_original_alkis_statistics_wrapper(alkis_pth, stats_out_pth):
    """
    Calculate statistics about original ALKIS data and owner entries.

    Args:
        alkis_pth (str): File path to the ALKIS data.
        stats_out_pth (str): File path for the output statistics CSV.

    Returns:
        None
    """
    print("Generate statistics about original ALKIS data and owner entries.")
    print("Read input data.")
    gdf = gpd.read_file(alkis_pth)  # Read ALKIS data from the file path

    # Calculate area for each geometry and convert to hectares.
    print("Calculate area for each geometry.")
    gdf['area'] = gdf['geometry'].area  # Calculate area of each geometry
    gdf['area'] = gdf['area'] / 10000  # Convert area to hectares

    # Set level column to "1" for all rows.
    print("Set level column to '1' for all rows.")
    gdf['level_overall'] = "1"

    # Calculate overall statistics.
    print("Calculate overall statistics.")
    df_stats = create_overall_stats_per_level(
        input_df=gdf,  # Input data frame
        level_col="level_overall",  # Level column
        area_col='area',  # Area column
        fid_col='OGC_FID',  # Feature ID column
        owners_col="EIGENTUEME"  # Owners column
    )

    # Write out the calculated statistics to a CSV file.
    print("Write out the calculated statistics to a CSV file.")
    df_stats.to_csv(stats_out_pth)


def calculate_area_per_community(alkis_iacs_inters_pth, owners_pth, comm_info_dict_pth, area_per_community_out_pth,
                                 comm_col, area_per_owner_out_pth=None, excl_categories_dict=None):
    print("Calculate areas.")
    print("Combine parcels with owner data")
    gdf = gpd.read_file(alkis_iacs_inters_pth)
    gdf["area"] = gdf["geometry"].area
    gdf['area'] = gdf['area'] / 10000
    df = pd.read_csv(owners_pth, sep=';', dtype={comm_col: str})

    df_alk = helper_functions.combine_parcels_with_owners(gdf[["OGC_FID", "EIGENTUEME", "BTNR", "area"]],
                                          df[["OGC_FID", "mother_company", "owner_merge", "level_c_category", "new_category", "new_category_ext", "level3", "agric", "fstateofowner", comm_col]])

    ## Temporary:
    t = df_alk.loc[df_alk["level_c_category"] == '2_9_1'].copy()
    t2 = t.drop_duplicates(subset="owner_merge")
    t3 = t2.groupby("level3").agg({"mother_company": "count"}).reset_index()
    print("Number of different company networks:", len(t["mother_company"].unique()))
    print("Number of different owners:", len(t2["owner_merge"].unique()))
    print(t3)

    ## Exclude certain categories
    if excl_categories_dict:
        column_name = excl_categories_dict["column_name"]
        excl_categs = excl_categories_dict["categories_to_exclude"]
        df_alk = df_alk.loc[~df_alk[column_name].isin(excl_categs)]

    ## Open dictionary with information about mother company names and names and numbers of subcompanies per community
    with open(comm_info_dict_pth) as json_file:
        comm_info_dict = json.load(json_file)
    comm_dict = comm_info_dict["comm_dict"]
    num_subcomps_dict = comm_info_dict["num_subcomps_dict"]
    name_subcomps_dict = comm_info_dict["name_subcomps_dict"]

    ## Reformat keys from string to integer
    # comm_dict = {int(float(key)): comm_dict[key] for key in comm_dict}
    # num_subcomps_dict = {int(float(key)): num_subcomps_dict[key] for key in num_subcomps_dict}
    # name_subcomps_dict = {int(float(key)): name_subcomps_dict[key] for key in name_subcomps_dict}

    ## Get number of original entries for "Unbekannt"
    df_unknown = df_alk.loc[df_alk["mother_company"] == 'unbekannt'].copy()
    num_unkown_owners = len(df_unknown["EIGENTUEME"].unique())
    print("number of original entries that are now classified into 'unbekannt':", num_unkown_owners)
    # print("Examples:")
    # for i in t["EIGENTUEME"].unique()[:10]:
    #     print(i)

    ## Get the area per community, add number and names of all owners from ALKIS
    def get_all_unique_attributes_in_list(group):
        out = list(set(list(group)))
        return out

    df_area = df_alk.groupby(by=[comm_col]).agg(
        area=pd.NamedAgg(column="area", aggfunc="sum"),
        num_farms_covered=pd.NamedAgg(column="BTNR", aggfunc=get_all_unique_attributes_in_list),
        number_agric_parcels=pd.NamedAgg(column="agric", aggfunc="sum"),
        all_fstates_owners=pd.NamedAgg(column="fstateofowner", aggfunc=get_all_unique_attributes_in_list),
        level_c_category=pd.NamedAgg(column="level_c_category", aggfunc=get_all_unique_attributes_in_list),
        new_category=pd.NamedAgg(column="new_category", aggfunc=get_all_unique_attributes_in_list),
        new_category_ext=pd.NamedAgg(column="new_category_ext", aggfunc=get_all_unique_attributes_in_list)
    ).reset_index()

    df_area["agric_owners"] = np.where(df_area["number_agric_parcels"] > 0, 1, 0)
    df_area["num_farms_covered"] = [len(lst) for lst in df_area["num_farms_covered"]]
    df_area["level_c_category"] = [lst[0] for lst in df_area["level_c_category"]]
    df_area["new_category"] = [lst[0] for lst in df_area["new_category"]]
    df_area["new_category_ext"] = [lst[0] for lst in df_area["new_category_ext"]]
    df_area["all_fstates_owners"] = [[i for i in lst if type(i) == str] for lst in df_area["all_fstates_owners"]]
    df_area["all_fstates_owners"] = ['_'.join(lst) for lst in df_area["all_fstates_owners"]]

    df_area["mother_company"] = df_area[comm_col].map(comm_dict)
    df_area["num_comp_network"] = df_area[comm_col].map(num_subcomps_dict)
    df_area["names_comps_network"] = df_area[comm_col].map(name_subcomps_dict)
    df_area.sort_values(by="area", ascending=False, inplace=True)
    df_area.index = range(len(df_area))
    total_area = df_alk['area'].sum()
    df_area['share'] = df_area['area'] / total_area * 100
    df_area['cum_share'] = df_area['share'].cumsum()
    df_area['area'] = round(df_area['area'], 5)
    df_area["place"] = df_area.index + 1

    dir_name = os.path.dirname(area_per_community_out_pth)
    helper_functions.create_folder(dir_name)

    df_area.to_csv(area_per_community_out_pth, index=False)

    ## --> calculate area of all companies in lindhorst family network compare to theos analyse
    # t = df_area.loc[df_area[comm_col] == 172].copy()
    # comps = t["names_comps_network"].iloc[0]
    # comps = comps.split('\n')
    # t2 = df_alk.loc[df_alk["owner_clean"].isin(comps)].copy()
    # t_area = t2[["area", "owner_clean"]].groupby("owner_clean").sum().reset_index()
    # t_area["area"] = round(t_area["area"] / 10000, 2)
    # t_area.to_csv(r"C:\Users\IAMO\Documents\work_data\chapter1\ALKIS\09_network_analysis\09_owners_aggregated_ALKIS_IACS_GEM\lindhorst_network_areas.csv")

    if area_per_owner_out_pth:

        dir_name = os.path.dirname(area_per_owner_out_pth)
        helper_functions.create_folder(dir_name)

        df_area2 = df_alk.groupby(by=["owner_merge"]).agg(
            area=pd.NamedAgg(column="area", aggfunc="sum"),
            num_farms_covered=pd.NamedAgg(column="BTNR", aggfunc=get_all_unique_attributes_in_list),
            number_agric_owners=pd.NamedAgg(column="agric", aggfunc="sum"),
            all_fstates_owners=pd.NamedAgg(column="fstateofowner", aggfunc=get_all_unique_attributes_in_list),
        ).reset_index()

        df_area2["num_farms_covered"] = [len(lst) for lst in df_area2["num_farms_covered"]]
        df_area2["all_fstates_owners"] = [[i for i in lst if type(i) == str] for lst in df_area2["all_fstates_owners"]]
        df_area2["all_fstates_owners"] = ['_'.join(lst) for lst in df_area2["all_fstates_owners"]]
        df_area2.sort_values(by="area", inplace=True, ascending=False)
        df_area.index = range(len(df_area))
        df_area2['share'] = df_area2['area'] / total_area * 100
        df_area2['cum_share'] = df_area2['share'].cumsum()
        df_area2['area'] = round(df_area2['area'], 3)
        df_area2["place"] = df_area2.index + 1
        df_area2.to_csv(area_per_owner_out_pth, index=False)


def plot_distribution_of_parcels_per_owner(alkis_iacs_inters_pth, owners_pth, bb_shp_pth, comm_col, community_lst, out_shp_folder, out_fig_folder):

    print("Combine parcels with owner data")
    gdf = gpd.read_file(alkis_iacs_inters_pth)
    gdf["area"] = gdf["geometry"].area
    gdf['area'] = gdf['area'] / 10000
    df = helper_functions.read_table_to_df(owners_pth)

    df_alk = helper_functions.combine_parcels_with_owners(
        gdf[["OGC_FID", "EIGENTUEME", "BTNR", "area", "geometry"]],
        df[["OGC_FID", "mother_company", "owner_merge", "level_c_category", comm_col]])

    df_alk = gpd.GeoDataFrame(df_alk)

    for i, comm in enumerate(community_lst):
        sub = df_alk.loc[df_alk[comm_col] == comm].copy()
        helper_functions.create_folder(out_shp_folder)
        out_shp_pth = f"{out_shp_folder}/community_{comm}.shp"
        if not os.path.exists(out_shp_pth):
            sub.to_file(out_shp_pth)

        mother_company = sub["mother_company"].iloc[0]

        shp_bb = gpd.read_file(bb_shp_pth)

        helper_functions.create_folder(out_fig_folder)
        out_fig_pth = f"{out_fig_folder}/{i:02}_community_{comm}.png"

        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=plotting_lib.cm2inch(12, 12))
        shp_bb.plot(edgecolor='none', facecolor='#bebebe', ax=axs, zorder=0)
        sub.plot(edgecolor='none', facecolor='#0000FF', ax=axs, zorder=1)
        fig.suptitle(f'{comm}_{mother_company}', fontsize=9)

        plt.savefig(out_fig_pth, dpi=300)
        plt.close()


## ------------------------------------------ RUN PROCESSES ---------------------------------------------------#
def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    # ## Calulate general statistics from IACS
    # calculate_iacs_area_statistics_wrapper(
    #     iacs_pth=IACS_PTH,
    #     iacs_statistics_out_pth=IACS_STATISTICS_OUT_PTH
    # )
    #
    # ## Calulate general statistics from original ALKIS
    # calculate_original_alkis_statistics_wrapper(
    #     alkis_pth=ALKIS_ORIG_PTH,
    #     stats_out_pth=ALKIS_ORIG_OUT_PTH)


    for threshold in [50]:
        print(threshold)
        comm_col = f"community_{threshold}"

        # ## Calulate general statistics from ALKIS
        # calculate_alkis_area_statistics_wrapper(
        #     alkis_iacs_inters_pth=ALK_IACS_INTERS_PTH,
        #     owners_pth=OWNERS_W_THRESH_PTH.format(threshold),
        #     out_pth=ALKIS_OUT_PTH.format(threshold))

        # Calculate area
        print("#####################################################\n "
              "Community IDs WITH threshold. +  ALKIS intersected with IACS data\n "
              "#####################################################")

        calculate_area_per_community(
            alkis_iacs_inters_pth=ALK_IACS_INTERS_PTH,
            owners_pth=OWNERS_W_THRESH_PTH.format(threshold),
            comm_info_dict_pth=COMMUNITY_INFO_DICT_W_THRESH.format(threshold),
            area_per_community_out_pth=AREA_PER_COMM_W_THRESH.format(threshold),
            comm_col=comm_col)
            # area_per_owner_out_pth=AREA_PER_OWN_PTH.format(threshold))

        print("... only looking at private land (excluding public and church land).")
        calculate_area_per_community(
            alkis_iacs_inters_pth=ALK_IACS_INTERS_PTH,
            owners_pth=OWNERS_W_THRESH_PTH.format(threshold),
            comm_info_dict_pth=COMMUNITY_INFO_DICT_W_THRESH.format(threshold),
            area_per_community_out_pth=AREA_PER_COMM_W_THRESH_EXCL_PUBL_CHURCH_PTH.format(threshold),
            comm_col=comm_col,
            area_per_owner_out_pth=AREA_PER_OWN_PTH_EXCL_PUBL_CHURCH_PTH.format(threshold),
            excl_categories_dict={"column_name": "level_c_category",
                                  "categories_to_exclude": ['4_1_1', '4_9_1', '5_1_1', '5_2_1',
                                                            '5_2_2', '5_2_3', '5_2_4', '5_2_5',
                                                            '5_3_1', '5_9_1']})

        df = pd.read_csv(AREA_PER_COMM_W_THRESH.format(threshold))
        community_lst = list(df[f"community_{threshold}"][:30])
        # community_lst = ['0_0']
        plot_distribution_of_parcels_per_owner(
            alkis_iacs_inters_pth=ALK_IACS_INTERS_PTH,
            owners_pth=OWNERS_W_THRESH_PTH.format(threshold),
            bb_shp_pth=BB_SHP,
            comm_col=comm_col,
            community_lst=community_lst,
            out_shp_folder=OUT_SHP_FOLDER_W_THRESH.format(threshold),
            out_fig_folder=OUT_FIG_FOLDER_W_THRESH.format(threshold))


    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)



if __name__ == '__main__':
    main()
