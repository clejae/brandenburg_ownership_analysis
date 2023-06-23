# Author:
# github repository:

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import os
import time
import geopandas as gpd
import pandas as pd
import re
import warnings

## Project functions
import helper_functions

# ------------------------------------------ USER VARIABLES ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input
ORIGINAL_ALKIS_SHP = r"00_data\vector\ALKIS\v_eigentuemer_bb_pgsql2shp.shp"

## Output (the output folder will be created automatically!)
REDUCED_SHP_PTH = r"00_data\vector\ALKIS\v_eigentuemer_bb_reduced.shp"
MANUAL_ASSIGNMENT_PTH = r"01_clean_owner_strings\owners_for_manual_assignment.csv"
OWNER_DF_PTH = r"01_clean_owner_strings\01_owners_and_addresses_separated.csv"
# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#

def separate_geometries_from_owner_information(orig_alkis_shp, reduced_shp_out_pth):
    """
    In the following steps, the ownership information will be separated from the geometries. However, we will need them
    later on but do not want to carry and load always unnecessary information. Therefore, we will save the geoemtries
    with the parcel IDs as as reduced shapefile.

    :param orig_alkis_shp: Path to original ALKIS shapefile, as was extracted from the database.
    :param reduced_shp_out_pth: Output path for reduced
    :return:
    """

    ## Load data
    print("Load data")
    orig_shp = gpd.read_file(orig_alkis_shp)

    ## Reduced to necessary information
    print("Reduce data.")
    reduced_shp = orig_shp[["OGC_FID", "EIGENTUEME", "AMTLFLSFL", "geometry"]].copy()
    reduced_shp["area"] = reduced_shp["geometry"].area
    for col in reduced_shp.columns:
        print(col, reduced_shp[col].dtype)
    reduced_shp["EIGENTUEME"] = reduced_shp["EIGENTUEME"].astype(str)
    reduced_shp = reduced_shp[["OGC_FID", "EIGENTUEME", "AMTLFLSFL", "area", "geometry"]]
    # reduced_shp.to_crs(25833)

    ## Write out
    print("Write out.")
    ## ToDo check why this is not working.
    reduced_shp.to_file(reduced_shp_out_pth)


def separate_owner_names_and_addresses(reduced_shp_pth, out_csv_manual_assignment_pth, owner_df_out_pth):
    """
    Cleans owner strings and brings them into uniform state.
    :param reduced_shp_pth: Path to reduced shapefile.
    :param out_csv_manual_assignment_pth: Output path for cases that need manual separation.
    :param owner_df_out_pth: Output path for ownership information.
    :return:
    """
    ## read land parcels, extract parcel ID and owner string
    df = gpd.read_file(reduced_shp_pth)
    df = pd.DataFrame(df)
    df = df[['OGC_FID', 'EIGENTUEME', 'AMTLFLSFL', 'area']]

    ## to lowercase, replace empty strings, remove commas at end of string
    df['owner_string'] = df['EIGENTUEME'].str.lower()
    df.loc[df['EIGENTUEME'].isnull(), 'owner_string'] = 'unbekannt'
    df['owner_string'] = df['owner_string'].str.rstrip(',')

    ## Get some characteristics of owner names
    ## Count breaks (indicates number of owners)
    ## Count asterisks (indicates number of owners)
    df['break_count'] = df['owner_string'].str.count('\n')
    df['asterisk_count'] = df['owner_string'].str.count('\*')

    ## clean owner string
    df['owner_string'] = df['owner_string'].apply(
        helper_functions.replace_characters_in_string, char_lst=['\n'], replace_string=',')
    df['owner_string'] = df['owner_string'].apply(
        helper_functions.replace_characters_in_string, char_lst=['  '], replace_string=' ')
    df['owner_string'] = df['owner_string'].apply(
        helper_functions.replace_characters_in_string, char_lst=['(mehrere)'], replace_string='')
    df['owner_string'] = df['owner_string'].apply(
        helper_functions.replace_characters_in_string, char_lst=['mehrere'], replace_string='')
    df['owner_string'] = df['owner_string'].str.strip(',')

    ## remove "geb. xxx" (born as xxx) strings
    df['owner_string'] = df['owner_string'].apply(
        helper_functions.remove_search_term_and_followup_in_substrings, search_term='geb.'
    )

    ## count commas
    df['comma_count'] = df['owner_string'].str.count(',')

    ## remove adresses
    df['owner_names'] = df['owner_string'].apply(helper_functions.remove_address)

    ## separate single owners from multiple owners
    df['own_num'] = 0
    df.loc[(df['break_count'] == 0) & (df['asterisk_count'] < 2) & (df['comma_count'] < 5), 'own_num'] = 1

    ## in a very few cases the owner_string is NaN, because the whole string is only one
    ## MAYBE INSTEAD OF isna() use == ''
    df.loc[df['owner_string'].isna(), 'owner_names'] = 'unbekannt'
    df.loc[df['owner_string'].isna(), 'owner_string'] = 'unbekannt'

    ## in some cases the owner names are mistaken by addresses. After checking, we assume that there is pretty certain
    ## only one owner entity and we can assume that the first part of the owner string is the owner name
    df.loc[(df['own_num'] == 1) & (df['owner_names'] == ''), 'owner_names'] = df.apply(
        lambda row: row.owner_string.split(',')[0], axis=1)

    ## in a
    df.loc[(df['own_num'] == 0), 'owner_names'] = df['EIGENTUEME'].apply(helper_functions.identify_owners, return_count=False)
    df.loc[(df['own_num'] == 0), 'own_num'] = df['EIGENTUEME'].apply(helper_functions.identify_owners, return_count=True)

    ## OPTIONAL: after exporting to csv and reading in again, some cases are NaNs,
    ## because their names are confused with streets. All these cases seem to be only single owners,
    ## thus the first occurence is their name
    ## MAYBE INSTEAD OF isna() use == ''
    df.loc[(df['own_num'] == 1) & (df['owner_names'].isna()), 'owner_names'] = df.apply(
        lambda row: row.owner_string.split(',')[0], axis=1)

    # replace owner strings for special cases, that let the find_addresses_in_string function crash
    df.loc[df[
               'owner_names'] == 'die gesamtheit der beteiligten', 'owner_string'] = 'die gesamtheit der beteiligten, unkown 1, 00000 unkown'
    df.loc[df['owner_names'] == 'wirwich, inge, * 1933-01-01 | +17.11.1979', 'owner_string'] = 'wirwich, inge, * 1933-01-01'

    df['addresses'] = df.apply(lambda row: helper_functions.find_addresses_in_string(row.owner_names, row.owner_string), axis=1)

    ## Get owners for which no addresses could be assigned and write to file for manual address identification
    df_miss = df.loc[df['addresses'] == '']
    df_miss = df_miss[['EIGENTUEME', 'owner_names', 'addresses']].copy()
    df_miss.drop_duplicates(subset='EIGENTUEME', inplace=True)
    if not df_miss.empty:
        warnings.warn(f"There are entries that need a manual separation of names and addresses."
                      f" Separate names from addresses before running the next steps in {out_csv_manual_assignment_pth}."
                      f" Save as {out_csv_manual_assignment_pth[:-4]}_corrected.csv.")
        df_miss.to_csv(out_csv_manual_assignment_pth, sep=';', index=False, encoding='ISO 8859-15')

    ## Write current result df to disc
    df.to_csv(owner_df_out_pth, sep=';', index=False)


def assign_missing_addresses_to_owners(owner_df_pth, csv_manual_assignment_pth):
    """
    Assigns missing addresses to owners.
    :param owner_df_pth: Path to owner df. Will also be used to save the corrected version.
    :param csv_manual_assignment_pth: Path to csv with separated owner names and addresses.
    :return:
    """

    if not os.path.exists(csv_manual_assignment_pth):
        warnings.warn(f"{csv_manual_assignment_pth} does not exist! Please separate names and addresses manually.")
        return

    df = pd.read_csv(owner_df_pth, sep=';')
    df_man = pd.read_csv(csv_manual_assignment_pth, sep=';', encoding='ISO 8859-15')

    for i, eig in enumerate(df_man['EIGENTUEME']):
        owner_names = eig.lower()
        owner_names = owner_names.strip(',')
        owner_names = owner_names.strip()
        df.loc[df['EIGENTUEME'] == eig, 'owner_names'] = owner_names
        address = df_man['addresses'].iloc[i]
        address = address.lower()
        address = address.strip(',')
        address = address.strip()
        print(i, owner_names)
        print(address + '\n')
        df.loc[df['EIGENTUEME'] == eig, 'addresses'] = address

    df.to_csv(owner_df_pth, sep=';', index=False)


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    separate_geometries_from_owner_information(
        orig_alkis_shp=ORIGINAL_ALKIS_SHP,
        reduced_shp_out_pth=REDUCED_SHP_PTH
    )
    separate_owner_names_and_addresses(
        reduced_shp_pth=REDUCED_SHP_PTH,
        out_csv_manual_assignment_pth=MANUAL_ASSIGNMENT_PTH,
        owner_df_out_pth=OWNER_DF_PTH
    )
    assign_missing_addresses_to_owners(
        owner_df_pth=OWNER_DF_PTH,
        csv_manual_assignment_pth=MANUAL_ASSIGNMENT_PTH[:-4] + '_corrected.csv'
    )

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)


if __name__ == '__main__':
    main()

