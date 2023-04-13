# Clemens Jänicke
# github Repo: https://github.com/clejae

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import pandas as pd
import os
import time
import math
from collections import defaultdict
import json
from collections import Counter
from fuzzywuzzy import process
import pandas as pd
from joblib import delayed, Parallel

## Project library
import helper_functions
# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input paths
ALKIS_PTH = r"05_georeference_addresses\05_owners_stretched_addresses_geolocated_distances.csv"
PTH_CLEANING_MANUAL = r"07_owner_name_cleaning\cleaning_table_manuelle_namensbereinigung.xlsx"
PTH_ALKIS_TO_DAFNE_NAMES = r"07_owner_name_cleaning\__name-connections_clean.csv" ## Alkis names to Dafne name
PTH_CLEANING_CHURCHES = r"07_owner_name_cleaning\cleaning_table_churches.xlsx"
PTH_CLEANING_NONPROF = r"07_owner_name_cleaning\cleaning_table_non-profit.xlsx"
PTH_CLEANING_PUBLIC = r"07_owner_name_cleaning\cleaning_table_public.xlsx"
PTH_CLEANING_PEOPLE = r"07_owner_name_cleaning\cleaning_dict_private_people.json"
PTH_CLEANING_PEOPLE_MANUAL = r"07_owner_name_cleaning\cleaning_dict_private_people_manual_assignment.json"
PTH_CLEANING_REST = r"07_owner_name_cleaning\cleaning_table_rest.xlsx"

AGRICULTURAL_OWNERS_PTH = r"00_data\tables\matching_results_alkis_farmsubs_futtermittel.csv" # This table was created by student collaborator
DAFNE_SEARCH_RESULTS_PTH = r"06_prepare_dafne_search\combined_search_results\Unternehmensdaten.xlsx"
PTH_SUBSIDIES_2019 = r"00_data\tables\subsidies\de_2019.csv"
PTH_SUBSIDIES_2020 = r"00_data\tables\subsidies\de_2020.csv"

## Output (the output folder will be created automatically!)
PTH_ALKIS_CLEAN_PRELIMINARY = r"07_owner_name_cleaning\07_owners_stretched_pre.csv"
PTH_ALKIS_CLEAN = r"07_owner_name_cleaning\07_owners_stretched.csv"

# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#


def manually_clean_owners(alkis_pth, manual_cleaning_pth, alkis_clean_pth):
    """
    Clean the ALKIS data with help of a table where manual adjustments were made to names or addresses.
    Args:
        alkis_pth: Path to dataframe with owners.
        manual_cleaning_pth: Path to dataframe with manual adjustments.
        alkis_clean_pth: Output path to dataframe with cleaned names.

    Returns:

    """

    print(f"##################################################\n"
          f"Correct owner names with df {manual_cleaning_pth}\n"
          f"##################################################\n")

    print("\tRead data.")
    df_alk = pd.read_csv(alkis_pth, sep=';')
    df_clean = pd.read_excel(manual_cleaning_pth)

    print("\tLength of input table:", len(df_alk))

    ## Replace NaNs with None
    for col in ["clean_address", "level1", "level2", "level3", "category"]:
        df_clean.loc[df_clean[col].isna(), col] = None

    ## Create dictionaries from table for mapping
    clean_name_dict = {}
    clean_addr_dict = {}
    clean_cate_dict = {}
    clean_lev1_dict = {}
    clean_lev2_dict = {}
    clean_lev3_dict = {}
    for row in df_clean.itertuples():
        clean_name_dict[row.owner_clean_old] = row.owner_clean_new
        if row.clean_address != None and type(row.clean_address) != float:
            clean_addr_dict[row.owner_clean_old] = row.clean_address
        if row.level2 != None and type(row.level2) != float:
            clean_lev1_dict[row.owner_clean_old] = row.level1
            clean_lev2_dict[row.owner_clean_old] = row.level2
            clean_lev3_dict[row.owner_clean_old] = row.level3
            clean_cate_dict[row.owner_clean_old] = row.category

    ## Create two lists that help to locate the entries in the ALKIS df that need updating
    ## As not all entries in the cleaning df have updated information on the levels and categories we need two lists
    sub1 = df_clean.loc[(df_clean["level2"].notna()), "owner_clean_old"].tolist()
    sub2 = df_clean.loc[(df_clean["level2"].isna()), "owner_clean_old"].tolist()

    ## For all instances where owner clean is in sub1 (i.e. level2 is categorized)
    ## update category, level1, level2, level3, owner_clean based on the new values
    ## that are assigned to the owner_clean column in the cleaning-df
    # sub = list(df_clean["owner_clean_old"].unique())
    df_alk.loc[df_alk["owner_clean"].isin(sub1), "category"] = df_alk.loc[
        df_alk["owner_clean"].isin(sub1), "owner_clean"].replace(clean_cate_dict)
    df_alk.loc[df_alk["owner_clean"].isin(sub1), "level1"] = df_alk.loc[
        df_alk["owner_clean"].isin(sub1), "owner_clean"].replace(clean_lev1_dict)
    df_alk.loc[df_alk["owner_clean"].isin(sub1), "level2"] = df_alk.loc[
        df_alk["owner_clean"].isin(sub1), "owner_clean"].replace(clean_lev2_dict)
    df_alk.loc[df_alk["owner_clean"].isin(sub1), "level3"] = df_alk.loc[
        df_alk["owner_clean"].isin(sub1), "owner_clean"].replace(clean_lev3_dict)
    df_alk.loc[df_alk["owner_clean"].isin(sub1), "owner_clean"] = df_alk.loc[
        df_alk["owner_clean"].isin(sub1), "owner_clean"].replace(clean_name_dict)

    ## For all instances where owner clean is in sub2 (i.e. level2 is missing)
    ## update clean_address, owner_clean based on the new values
    ## that are assigned to the owner_clean column in the cleaning-df
    df_alk.loc[df_alk["owner_clean"].isin(sub2), "clean_address"] = df_alk.loc[
        df_alk["owner_clean"].isin(sub2), "owner_clean"].replace(clean_addr_dict)
    df_alk.loc[df_alk["owner_clean"].isin(sub2), "owner_clean"] = df_alk.loc[
        df_alk["owner_clean"].isin(sub2), "owner_clean"].replace(clean_name_dict)

    t = df_alk.loc[df_alk["owner_clean"].isin(df_clean["owner_clean_new"].tolist())].copy()

    print(f"\tWrite output to {alkis_clean_pth}")
    df_alk.to_csv(alkis_clean_pth, sep=';', index=False)


def clean_company_names(alkis_clean_pth, alkis_to_dafne_names_pth):
    """
    Clean company names, i.e. replace ALKIS names with DAFNE names. The "ALKIS-to-DAFNE-names" dataframe needs to be
    manually created from the intermediate output of DAFNE.
    Args:
        alkis_clean_pth: Path to dataframe with owners that were already manually cleaned. Will also be used for output.
        alkis_to_dafne_names_pth: Path to dataframe with ALKIS names and DAFNE names.

    Returns:

    """
    print("Clean company names, i.e. replace ALKIS names with DAFNE names.")

    ## Read data
    print("\tRead data.")
    owner_df = pd.read_csv(alkis_clean_pth, sep=';')
    df_conn = pd.read_csv(alkis_to_dafne_names_pth)

    print("\tLength of input table:", len(owner_df))

    ## Dictionary assigns bvd company names to alkis names
    ## if there are no bvd names use the old alkis names to avoid NaNs while mapping
    bvd_dict = {}
    for row in df_conn.itertuples():
        alk_name = row.search_name
        bvd_name = row.comp_name
        uni_name = row.uni_name_not_found
        if not type(bvd_name) == float:
            bvd_dict[alk_name] = bvd_name
        elif not type(uni_name) == float:
            bvd_dict[alk_name] = uni_name
        else:
            bvd_dict[alk_name] = alk_name

    ## print information before cleaning
    num_comp = len(owner_df.loc[owner_df["level1"] == 2, "owner_merge"].unique())
    print("Number of companies before unification:", num_comp)

    ## map bvd company names to alkis company names
    owner_df.loc[owner_df["level1"] == 2, "owner_clean"] = owner_df.loc[owner_df["level1"] == 2, "owner_clean"].map(bvd_dict)

    ## save old owner merge
    owner_df["owner_merge_old"] = owner_df["owner_merge"]

    ## Update new owner merge. use owner clean for it and clean the string
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_clean"].str.lower()
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('mit sitz in ', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('ä', 'ae', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('ö', 'oe', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('ü', 'ue', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('ß', 'ss', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('sitz in ', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace(' mit sitz ', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace(' sitz ', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace(' in ', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace(' ', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('-', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace(',', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('&', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('+', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('`', '', regex=False)
    owner_df.loc[owner_df["level1"] == 2, "owner_merge"] = owner_df.loc[
        owner_df["level1"] == 2, "owner_merge"].str.replace('\n', '', regex=False)

    ## print information after cleaning
    num_comp = len(owner_df.loc[owner_df["level1"] == 2, "owner_clean"].unique())
    print("\tNumber of companies after unification:", num_comp)

    #Number of companies before unification: 5379
    #Number of companies after unification: 4658

    # ## check if two companies have the same address
    # ## --> cleaning is not worth the effort, there is only a hand full of cases
    # comp_df = owner_df.loc[owner_df["level1"] == 2].copy()
    # comp_df.drop_duplicates(subset="owner_clean", inplace=True)
    #
    # ## to compare with number of companies
    # num_addr = len(comp_df["clean_address"].unique())
    #
    # ## get number of different companies per address
    # counts = Counter(comp_df["clean_address"])
    #
    # ## keep only those entries where the same addresses are assigned to multiple companies, not unkown, and with
    # ## a full address
    # ## there are several problems here: the full address columns is not 100% correct, sometimes the address is wrongly assigned
    # counts = [key for key in counts if counts[key] > 1]
    # counts.remove("unbekannt")
    # comp_df = comp_df.loc[comp_df["clean_address"].isin(counts)].copy()
    # comp_df = comp_df.loc[~comp_df["full_address"].isin([6])].copy()
    #
    # ## keep only those entries which could not be found in DAFNE, as the others are probably corrected by the Dafne alg
    # comps = comp_df["owner_clean"].tolist()
    # comps = [comp for comp in comps if comp.islower()]
    # addr = comp_df.loc[comp_df["owner_clean"].isin(comps), "clean_address"].copy()
    # comp_df = comp_df.loc[comp_df["clean_address"].isin(addr)].copy()
    #
    # ## look at the remaining cases in debug mode
    # comp_df.sort_values(by="clean_address", inplace=True)
    print(f"\tWrite output to {alkis_clean_pth}")
    owner_df.to_csv(alkis_clean_pth, sep=';', index=False)



def clean_churchs_nonprof_public_names(alkis_pth_name, clean_table_name):
    """
    Cleans names of public, religious, non-profit and other organizations with help of a manually created table with
    old names, clean names and clean addresses of all the institutions.
    Args:
        alkis_pth_name: Path to dataframe with owners that were already manually cleaned. Will also be used for output.
        clean_table_name: Path to manually created table with clean names.

    Returns:

    """
    print(f"Correct owner names with df {clean_table_name}")

    df_clean = helper_functions.read_table_to_df(clean_table_name) #pd.read_excel(clean_table_name)
    df_alk = helper_functions.read_table_to_df(alkis_pth_name)  # pd.read_csv(alkis_pth_name, sep=';')

    print("\tLength of input table:", len(df_alk))

    ## Replace NaNs with None
    for col in ["clean_address"]:
        df_clean.loc[df_clean[col].isna(), col] = None

    ## Create dictionaries from table for mapping
    clean_name_dict = {}
    clean_addr_dict = {}
    for row in df_clean.itertuples():
        clean_name_dict[row.owner_clean_old] = row.owner_clean_new
        if row.clean_address != None and type(row.clean_address) != float:
            clean_addr_dict[row.owner_clean_old] = row.clean_address

    ## Create two lists that help to locate the entries in the ALKIS df that need updating
    ## As not all entries in the cleaning df have updated information on the levels and categories we need two lists
    sub1 = df_clean.loc[(df_clean["clean_address"].notna()), "owner_clean_old"].tolist()
    sub2 = df_clean.loc[(df_clean["owner_clean_new"].notna()), "owner_clean_old"].tolist()

    ## For all instances where owner clean is in sub1 (i.e. level2 is categorized)
    ## update owner_clean and clean_address based on the new values
    ## that are assigned to the owner_clean column in the cleaning-df
    ## Important: first update the address based on the old name, then update the old name
    df_alk.loc[df_alk["clean_address"].isin(sub1), "clean_address"] = df_alk.loc[
        df_alk["clean_address"].isin(sub1), "owner_clean"].replace(clean_addr_dict)
    df_alk.loc[df_alk["owner_clean"].isin(sub2), "owner_clean"] = df_alk.loc[
        df_alk["owner_clean"].isin(sub2), "owner_clean"].replace(clean_name_dict)

    ## For manual checking in debug mode:
    ## create subset of alkis-df that can be found in cleaning-df
    t = df_alk.loc[df_alk["owner_clean"].isin(df_clean["owner_clean_new"].tolist())].copy()

    ## Write to disc
    print(f"\tWrite output to {alkis_pth_name}")
    df_alk.to_csv(alkis_pth_name, sep=';', index=False)


def create_cleaning_dictionary_private_persons(df_alk, clean_dict_pth=None):
    """
    Creates automatically a dictionary with unified addresses and names.
    Args:
        df_alk: Path to dataframe with owners that were already manually cleaned.
        clean_dict_pth: Output path to dictionary with cleaned entries.

    Returns:

    """
    print("Create cleaning dict of private persons.")
    ## Create cleaning dictionary
    ## Get subset of people with birthdate, indicated by *
    sub = df_alk.loc[df_alk["category"] == 1].copy()
    sub = sub.loc[sub["owner_clean"].str.count("\*") > 0].copy()
    sub.drop_duplicates(subset='owner_clean', inplace=True)

    ## derive birthdate, family name and surnames
    sub["birthdate"] = sub.apply(lambda row: row.owner_clean.split(',')[-1], axis=1)
    sub = sub.loc[sub["birthdate"].str.count("-") > 1].copy()
    sub["name"] = sub.apply(lambda row: row.owner_clean.split('*')[0], axis=1)
    sub["comma_count"] = sub["name"].str.count(',')
    sub = sub[sub["comma_count"] > 0].copy()
    sub["family_name"] = sub.apply(lambda row: row.owner_clean.split(',')[0], axis=1)
    sub["surname"] = sub.apply(lambda row: row.owner_clean.split(',')[1], axis=1)
    sub["surname1"] = sub.apply(lambda row: row.owner_clean.split(' ')[1], axis=1)

    ## create columne of family_name + birthdate combination
    sub["comb"] = sub["family_name"] + sub["birthdate"]

    ## Count occurence of this combination. If a combination occurs several times,
    ## this indicates the possible occurrence of a name with different spellings
    counts = Counter(sub["comb"].tolist())
    threshold = 1
    lst = []
    for k, v in counts.items():
        if v > threshold:
            lst.append(k)
    df_work = sub.loc[sub["comb"].isin(lst)].copy()
    df_work.reset_index(inplace=True)
    df_work.drop(columns=['index'], inplace=True)
    # df_done = sub.loc[~sub["comb"].isin(lst)].copy()

    ## Open a json file to write the results in
    if clean_dict_pth:
        with open(clean_dict_pth, 'w', encoding='ISO-8859-1') as file:
            file.write("{\n")

    correction_dict = {}
    print("\tLoop over entries")
    for row in df_work.itertuples():
        ## Get characteristics of current person
        curr_person = row.owner_clean
        family_name = row.family_name
        surname = row.surname1
        if '-' in surname:
            surname = surname.split('-')[0]
        birthdate = row.birthdate

        ## if current person is already listed in dictionary, then skip
        if curr_person in correction_dict:
            continue

        ## create a subset of all people that have same family name, birthdate and a similar surname
        sub2 = df_work.loc[df_work["family_name"].str.count(family_name) > 0].copy()
        sub2 = sub2.loc[sub2["surname"].str.count(surname) > 0].copy()
        sub2 = sub2.loc[sub2["birthdate"].str.count(birthdate) > 0].copy()

        ## remove current person from subset, if df is still not empty, then add to dictionary
        sub3 = sub2.loc[sub2["owner_clean"] != curr_person].copy()
        if not sub3.empty:
            persons = sub3["owner_clean"].tolist()
            ## create another subset with possible addresses
            ## if possible addresses fulfill certain quality requirement as indicated by code 5,6,7
            ## then use the best address. If not use any
            sub4 = sub2.loc[sub2["full_address"].isin([5, 6, 7])]
            if not sub4.empty:
                max_f = sub4["full_address"].max()
                clean_address = sub4.loc[sub4["full_address"] == max_f, "clean_address"].iloc[0]
            else:
                clean_address = sub2["clean_address"].iloc[0]

            ## for all persons left, add the information to dictionary
            for person in persons:
                if clean_dict_pth:
                    with open(clean_dict_pth, "a", encoding='ISO-8859-1') as file:
                        # file.write(f'\t"{person}": "{curr_person}",\n')
                        file.write(f'\t"{person}":' + '{' + f'"clean_name": "{curr_person}", "clean_address": "{clean_address}"'+ '},' + '\n')
                sub_dict = {}
                sub_dict["clean_name"] = curr_person
                sub_dict["clean_address"] = clean_address
                correction_dict[person] = sub_dict

    if clean_dict_pth:
        with open(clean_dict_pth, "a", encoding='ISO-8859-1') as file:
            file.write("}")

    return correction_dict


def create_aristrocracy_dict(alkis_pth, out_pth):
    """
    Creates a dictionary of aristrocratic names that can be used for manual cleaning.
    Args:
        alkis_pth: Path to dataframe with owners.
        out_pth: Output path to dictionary that needs some manualy adjsutment.

    Returns:

    """
    print("Create a dictionary of aristrocratic names that can be used for manual cleaning.")
    ## Open Data
    print("\tRead data.")
    df_alk = pd.read_csv(alkis_pth, sep=";")

    df_alk["conn_left"] = df_alk["owner_clean"]

    ## For private owner create merge name with pattern: "surname familyname_location" to merge with dafne entries
    df_alk.loc[(df_alk["level3"] == "1_1_1") & (df_alk["conn_left"].str.count(',') > 0), "conn_left"] = \
        df_alk.loc[(df_alk["level3"] == "1_1_1") & (df_alk["conn_left"].str.count(',') > 0), "conn_left"].apply(
            lambda row: row.split(',')[1] + ' ' + row.split(',')[0])

    df_alk["city"] = df_alk["clean_address"].apply(helper_functions.identify_city)
    df_alk.loc[(df_alk["level3"] == '1_1_1') &
               (df_alk["owner_clean"].str.count(',') > 0), "conn_left"] = df_alk.loc[(df_alk["level1"] == 1) &
                                                                                     (df_alk["owner_clean"].str.count(
                                                                                         ',') > 0), "conn_left"] + '_' + \
                                                                          df_alk.loc[(df_alk["level1"] == 1) &
                                                                                     (df_alk["owner_clean"].str.count(
                                                                                         ',') > 0), "city"]

    df_alk["conn_left"] = df_alk["conn_left"].str.lstrip(' ')
    df_alk["conn_left"] = df_alk["conn_left"].str.rstrip(' ')

    df_alk.loc[df_alk["level3"] == '1_1_1'].copy()
    df_alk.drop_duplicates(subset=["owner_clean"], inplace=True)

    t = df_alk.loc[
        df_alk["conn_left"].str.contains(" von ") | df_alk["conn_left"].str.contains(" van ") | df_alk[
            "conn_left"].str.contains("von ") | df_alk["conn_left"].str.contains("van ")
        ].copy()
    t.drop_duplicates(subset="conn_left", inplace=True)
    t["birthdate"] = t["owner_clean"].apply(lambda x: x.split("*")[1][:11] if "*" in x else '')
    t["birthdate"] = t["birthdate"].str.strip()
    t["search_term"] = t["conn_left"] + '_' + t["birthdate"]
    t = t.loc[t["level3"] == '1_1_1'].copy()

    print("\tNumber of people to search:", len(t))

    out_dict = {}
    sterm_lst = t["search_term"].tolist()
    for i, search_term in enumerate(sterm_lst):
        print("\t", i, search_term)
        vorname = search_term.split(" ")[0]
        nachname = search_term.split("_")[0].split(' ')[-1]
        ort = search_term.split("_")[1]
        bdate = search_term.split("_")[-1]
        t2 = df_alk.loc[
            (df_alk["owner_clean"].str.contains(vorname, regex=False) &
            df_alk["owner_clean"].str.contains(nachname, regex=False) &
            df_alk["owner_clean"].str.contains(bdate, regex=False)) |
            (df_alk["owner_clean"].str.contains(vorname, regex=False) &
             df_alk["owner_clean"].str.contains(nachname, regex=False) &
             df_alk["clean_address"].str.contains(ort, regex=False))
            ].copy()
        t2.drop_duplicates(subset="owner_clean", inplace=True)
        if len(t2) > 1:
            print("\tNumber matches:", len(t2) - 1)
            t3 = t2.loc[
                (t2["owner_clean"].str.contains("von") &
                 t2["clean_address"] != "unbekannt")
                |(t2["owner_clean"].str.contains("van")&
                 t2["clean_address"] != "unbekannt")
                ].copy()
            clean_name = t3["owner_clean"].iloc[0]
            clean_address = t3["clean_address"].iloc[0]
            t4 = t2[t2["owner_clean"] != clean_name]
            lst = t4["owner_clean"].tolist()
            for item in lst:
                out_dict[item] = {"clean_name": clean_name, "clean_address": clean_address}

            lst = t2["owner_clean"].tolist()
            t = t.loc[~t["owner_clean"].isin(lst)].copy()

    with open(out_pth, 'w', encoding='ISO-8859-1') as file:
        json.dump(out_dict, file, indent=4)

    helper_functions.print_red("\tClean out dict manually and rename to manual cleaning dict for private people!!")


def clean_private_people(alkis_pth_name, clean_dict_pth, clean_dict_pth_manual):
    """
    Cleans names of private people. Requires an input dictionary with uncleaned names as keys with subdictionaries as
    items. These carry "clean_name" and "clean_address" as keys with the and manually adjusted names and addresses as items.
    Args:
        alkis_pth_name: Path to dataframe with owners that were already manually cleaned. Will also be used for output.
        clean_dict_pth: Path to dictionary created during automatic cleaning. It will be created automatically.
        clean_dict_pth_manual: Path to dictionary with manualy adjustments for private people.

    Returns:

    """

    print(f"Correct private owner names based on family names, parts of surnames and birthdates")
    print("\tRead data.")
    df_alk = helper_functions.read_table_to_df(alkis_pth_name)

    print("\tLength of input table:", len(df_alk))

    if not os.path.exists(clean_dict_pth):
        clean_dict_pre = create_cleaning_dictionary_private_persons(df_alk=df_alk, clean_dict_pth=clean_dict_pth)
    else:
        file = open(clean_dict_pth)
        clean_dict_pre = json.load(file)

    ## Add manual assigned people to cleaning dict
    file = open(clean_dict_pth_manual)
    clean_dict_manual = json.load(file)

    pop_lst = []
    for key_manual in clean_dict_manual:
        manual_val = clean_dict_manual[key_manual]["clean_name"]
        manual_adr = clean_dict_manual[key_manual]["clean_address"]
        if manual_val in clean_dict_pre:
            print("\tClean name from manual_clean_dict is key (wrong name) in aut_clean_dict. Remove:", manual_val)
            clean_dict_pre.pop(manual_val)
        if key_manual in clean_dict_pre:
            print("\tWrong name from manual_clean_dict is key (wrong name) in aut_clean_dict. Unify.", key_manual)
            clean_dict_pre[key_manual]["clean_name"] = manual_val
            clean_dict_pre[key_manual]["clean_address"] = manual_adr

    for key_manual in clean_dict_manual:
        manual_val = clean_dict_manual[key_manual]["clean_name"]
        manual_adr = clean_dict_manual[key_manual]["clean_address"]
        for key_pre in clean_dict_pre:
            if clean_dict_pre[key_pre]["clean_name"] == key_manual:
                print("\tWrong name from manual_clean_dict is value (correct name) in aut_clean_dict. Unify.", key_manual)
                clean_dict_pre[key_pre]["clean_name"] = manual_val
                clean_dict_pre[key_pre]["clean_address"] = manual_adr

    clean_dict = clean_dict_pre | clean_dict_manual

    # for key in clean_dict_manual:
    #     if key not in clean_dict:
    #         clean_dict[key] = clean_dict_manual[key]
    # t = list(v for v in clean_dict if 'matthes' in v)

    ## check case that made problems: matthes reiner
    clean_owner_names = []
    wrong_owner_names = []
    for key in clean_dict:
        clean_owner_names.append(clean_dict[key]["clean_name"])
        wrong_owner_names.append(key)
    clean_owner_names = set(list(clean_owner_names))
    wrong_owner_names = set(list(wrong_owner_names))

    ## get a subset of all entries with wrong names
    df_wrong1 = df_alk.loc[df_alk["owner_clean"].isin(wrong_owner_names)].copy()
    print("\tNumber of entries with wrong names", len(df_wrong1))

    ## key is the name that needs correction
    clean_name_dict = {key: clean_dict[key]["clean_name"] for key in clean_dict}
    clean_addr_dict = {key: clean_dict[key]["clean_address"] for key in clean_dict}

    clean_addr_dict = {}
    for key in clean_dict:
        ## key is wrong owner name
        new_address = clean_dict[key]["clean_address"]
        old_address = df_wrong1.loc[df_wrong1["owner_clean"] == key, "clean_address"].iloc[0]

        if (old_address != None) and type(old_address) != float:
            clean_addr_dict[old_address] = new_address

    ## first clean addresses based on wrong names
    df_alk.loc[df_alk["owner_clean"].isin(wrong_owner_names), "clean_address"] = df_alk.loc[
        df_alk["owner_clean"].isin(wrong_owner_names), "clean_address"].replace(clean_addr_dict)
    ## then clean names
    df_alk.loc[df_alk["owner_clean"].isin(wrong_owner_names), "owner_clean"] = df_alk.loc[
        df_alk["owner_clean"].isin(wrong_owner_names), "owner_clean"].replace(clean_name_dict)

    ## get again a subset of all entries with wrong names - theoretically this should be empty now
    df_wrong2 = df_alk.loc[df_alk["owner_clean"].isin(wrong_owner_names)].copy()
    if df_wrong2.empty:
        print("\tAll entries cleaned.")
    else:
        print("\tNot all entries cleaned.")
        df_wrong2.drop_duplicates(subset="owner_clean", inplace=True)

    ## get a subset with corrected names
    df_correct = df_alk.loc[df_alk["owner_clean"].isin(clean_owner_names)]

    ## Write to disc
    t = df_alk.loc[df_alk["owner_clean"].str.count("matthes") > 0].copy()
    t2 = t.loc[t["owner_names"].str.count("rainer") > 0].copy()
    t3 = df_alk.loc[df_alk["owner_names"].str.count("rainer") > 0].copy()
    t4 = t.loc[t["owner_names"].str.count("reiner") > 0].copy()
    t5 = df_alk.loc[df_alk["clean_address"].str.count("\*") > 0].copy()

    ## remove from private people names: graf, graefin, baronin, baron, freiherr, freifrau, herzog, herzogin
    def remove_words_from_family_names(text, search_term):
        subs = text.split(',')
        # get family name
        fam_name = subs[0]
        # get rest
        rest = subs[1:]
        rest = [i.strip() for i in rest]
        rest = ', '.join(rest)
        # separate wores
        fam_name_lst = fam_name.split(' ')
        # only remove a word if there are more then one word as a family name
        if len(fam_name_lst) > 1:
            if search_term in fam_name_lst:
                fam_name_lst.remove(search_term)
                fam_name_lst = [i.strip() for i in fam_name_lst]
                fam_name_new = ' '.join(fam_name_lst)
            else:
                fam_name_lst = [i.strip() for i in fam_name_lst]
                fam_name_new = ' '.join(fam_name_lst)
        else:
            fam_name_new = ' '.join(fam_name_lst)
        out_text = f"{fam_name_new}, {rest}"
        return out_text

    for word in ["van", "von", "graf", "graefin", "baronin", "baron", "freiherr", "freifrau", "herzog", "herzogin"]:
        df_alk.loc[df_alk["level3"] == '1_1_1', 'owner_clean'] = df_alk.loc[df_alk["level3"] == '1_1_1', 'owner_clean'].apply(remove_words_from_family_names, search_term=word)

    print(f"\tWrite output to {alkis_pth_name}")
    df_alk.to_csv(alkis_pth_name, sep=';', index=False)


def clean_gbrs(alkis_pth_name):
    """
    Quick intermediate cleaning of GbRs (Civil law partnerships).
    Args:
        alkis_pth_name: Path to dataframe with owners that were already manually cleaned. Will also be used for output.

    Returns:

    """
    print("Clean quickly GbRs")
    df_alk = pd.read_csv(alkis_pth_name, sep=';')

    print("\tLength of input table:", len(df_alk))

    ## Clean GBRs
    df_alk.loc[df_alk["level3"] == '1_2_1', 'owner_clean'] = df_alk.loc[
        df_alk["level3"] == '1_2_1', 'owner_clean'].str.replace('bestehend aus', '', regex=False)
    df_alk.loc[df_alk["level3"] == '1_2_1', 'owner_clean'] = df_alk.loc[
        df_alk["level3"] == '1_2_1', 'owner_clean'].str.replace('bestehen aus', '', regex=False)

    t = df_alk.loc[df_alk["level3"] == '1_2_1'].copy()
    t.drop_duplicates(subset="owner_clean", inplace=True)
    t.sort_values(by="owner_clean", inplace=True)

    print(f"\tWrite output to {alkis_pth_name}")
    df_alk.to_csv(alkis_pth_name, sep=';', index=False)

def create_new_unique_identifier_and_unify_addresses(alkis_pth_name):
    """
    Based on the new cleaned names create a new unique identifier and unify the addresses for all unique entries.
    Args:
        alkis_pth_name: Path to dataframe with owners that were already manually cleaned. Will also be used for output.

    Returns:

    """
    print("Create new unique identifier and unify addresses.")
    ## Open data
    print("\tRead data.")
    df_alk = pd.read_csv(alkis_pth_name, sep=';')

    print("\tLength of input table:", len(df_alk))

    df_alk["owner_clean"] = df_alk["owner_clean"].str.strip()
    df_alk["owner_clean"] = df_alk["owner_clean"].str.strip(',')

    ## Create new owner merge column
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_clean"].str.lower()
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_merge"].str.replace("-", "", regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_merge"].str.replace(" ", "", regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_merge"].str.replace(".", "", regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_merge"].str.replace("ä", "ae", regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_merge"].str.replace("ö", "oe", regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_merge"].str.replace("ü", "ue", regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_merge"].str.replace("ß", "ss", regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_merge"].str.replace("!", "", regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', "owner_merge"] = df_alk.loc[
        df_alk["level3"] != '1_1_1', "owner_merge"].str.replace("-", "", regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[
        df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace('mit sitz in ', '', regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[
        df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace('sitz in ', '', regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[
        df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace(' mit sitz ', '', regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace(' sitz ', '', regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[
        df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace(' in ', '', regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[
        df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace(',', '', regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[
        df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace('&', '', regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[
        df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace('+', '', regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[
        df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace('`', '', regex=False)
    df_alk.loc[df_alk["level3"] != '1_1_1', 'owner_merge'] = df_alk.loc[
        df_alk["level3"] != '1_1_1', 'owner_merge'].str.replace('\n', '', regex=False)

    def count_addresses(name_group):

        uni_addresses = list(name_group["clean_address"].unique())
        count = len(uni_addresses)

        return count

    ## Reduce to all unique entries based on combination of owner_clean and clean address
    ## Goal: Unify addresses for all unique owner_clean names
    df_red = df_alk.drop_duplicates(subset=["owner_clean", "clean_address"])
    clean_owner_address_counts = df_red[["owner_clean", "clean_address"]].groupby("owner_clean").apply(
        count_addresses).reset_index()
    clean_owner_address_counts.columns = ["clean_owner", "num_addresses"]
    clean_owner_address_counts = clean_owner_address_counts.loc[clean_owner_address_counts["num_addresses"] > 1]

    ## Get entries with multiple addresses
    multis = clean_owner_address_counts["clean_owner"].tolist()

    ## Remove GBRs (Why again?)
    for oc in ["gbr", "gbr bestehend aus"]:
        if oc in multis:
            print(f"\tremoved '{oc}' from multis")
            multis.remove(oc)

    ## Do the magic. Unify addresses
    if len(multis) > 0:
        df_red = df_red.loc[df_red["owner_clean"].isin(multis)].copy()
        df_red.sort_values(by="owner_clean", inplace=True)

        def find_best_address(name_group):

            max_val = name_group["full_address"].max()

            selection = name_group.loc[name_group["full_address"] == max_val].copy()

            if len(selection) == 1:
                return list(selection["clean_address"])[0]
            else:
                for val in ["google_api", "nominatim", "osm_level1", "fuzzy_matching", "osm_level2"]:
                    selection2 = selection.loc[selection["geocoding"] == val].copy()
                    if len(selection2) > 0:
                        return list(selection2["clean_address"])[0]

            return "unbekannt"

        df_addresses = df_red[["owner_clean", "clean_address", "full_address", "geocoding"]].groupby("owner_clean").apply(
            find_best_address).reset_index()
        df_addresses.columns = ["owner_clean", "clean_address"]

        addresses_dict = {row.owner_clean: row.clean_address for row in df_addresses.itertuples()}

        df_alk.loc[df_alk["owner_clean"].isin(multis), "clean_address"] = df_alk.loc[
            df_alk["owner_clean"].isin(multis), "owner_clean"].map(addresses_dict)

    ## Identify all non-private entries where there are multiple clean owner entries for a single owner merge.
    ## Goal: to identify cases where more cleaning needs to be done
    ## But first clean one manually:
    df_alk.loc[
        df_alk[
            "owner_clean"] == "gbr mit dem namen,  gbr ulrich, jochen und joerg graebert grundstuecksgemeinschaft",
        "owner_clean"] = "gbr mit dem namen gbr ulrich, jochen und joerg graebert grundstuecksgemeinschaft"

    def count_names(name_group):

        uni_names = list(name_group["owner_clean"].unique())
        count = len(uni_names)

        return count

    df_red = df_alk.loc[df_alk["level3"] != '1_1_1'].copy()
    df_red.drop_duplicates(subset=["owner_merge", "owner_clean"], inplace=True)
    owner_merge_name_counts = df_red[["owner_merge", "owner_clean"]].groupby("owner_merge").apply(
        count_names).reset_index()
    owner_merge_name_counts.columns = ["owner_merge", "num_names"]
    owner_merge_name_counts = owner_merge_name_counts.loc[owner_merge_name_counts["num_names"] > 1]
    multis = owner_merge_name_counts["owner_merge"].tolist()

    for oc in ["gbr", "gbrbestehendaus"]:
        if oc in multis:
            print(f"removed '{oc}' from multis")
            multis.remove(oc)
    df_red = df_red.loc[df_red["owner_merge"].isin(multis)].copy()
    df_red.sort_values(by="owner_merge", inplace=True)

    ## In cases where there is the word "unbekannt" in column owner_merge  (likely due to unkown address),
    ## replace the owner_merge with owner clean
    df_alk.loc[df_alk["owner_merge"].str.count("unbekannt") > 0, "owner_merge"] = \
        df_alk.loc[df_alk["owner_merge"].str.count("unbekannt") > 0, "owner_clean"]

    if not df_red.empty:
        df_red[["level1", "owner_clean", "owner_merge", "clean_address"]].to_csv(fr"07_owner_name_cleaning\all_levels_cleaning.csv", sep=",", index=False)
    else:
        print("\tWork table empty, writing table out. Length of final table:", len(df_alk))
        df_alk.to_csv(alkis_pth_name, sep=';', index=False)

    t = df_alk.loc[df_alk["level3"] == '1_2_1'].copy()
    t.drop_duplicates(subset="owner_clean", inplace=True)
    t.sort_values(by="owner_clean", inplace=True)

    uni_levels = df_red["level1"].unique()


def fuzzy_matching_alkis_data(alkis_path, matching_data, matched_owners_pth):
    """
    Matches ALKIS names with subsidy recipients. Written by Tilman Schmitz. Not checked whether it runs.
    Args:
        alkis_path: Path to list of ALKIS owners (very old version).
        matching_data: Path to dataframe created with R-script (07_data_preparations_alkis_farmsubsidy.R)
        matched_owners_pth: Output path.

    Returns:

    """

    def highest_match_fuzzy_wuzzy(word_postcode_tuple, acc_value_1, acc_value_2):
        # get information from input df
        name_subs = word_postcode_tuple[0]
        postcode = word_postcode_tuple[1]
        original_names = word_postcode_tuple[2]
        owner_cat = word_postcode_tuple[3]
        origin = word_postcode_tuple[4]

        # subset alkis list based on postcode and transform it into list
        new_alkis = alkis.loc[alkis['postcode'] == postcode]
        new_alkis_names_list = new_alkis.new_names.tolist()

        # get fuzzy matching result with highest accuracy
        highest = process.extractOne(name_subs, new_alkis_names_list)

        if owner_cat == 1:

            try:
                if highest[1] > acc_value_1:
                    match = list(highest)
                    match.insert(0, name_subs)
                    match.extend((original_names, owner_cat, origin, postcode))
                    return match
            except:
                print(highest, "is none")
        else:

            try:
                if highest[1] > acc_value_2:
                    match = list(highest)
                    match.insert(0, name_subs)
                    match.extend((original_names, owner_cat, origin, postcode))
                    return match
            except:
                print(highest, "is none")

    # load in data
    alkis = pd.read_csv(alkis_path, sep=';')
    # subs_2019_recipient_list = pd.read_csv(subs_2019_recipient_pth)
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!! HERE IS SOMETHING MISSING WHICH I CAN'T RECONSTRUCT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subs_2019_recipient_list = []

    matching_data_list = pd.read_csv(matching_data, sep=';', dtype={"postcode": object})
    matching_data_list = matching_data_list[subs_2019_recipient_list['owner_cat'] == 1]

    # subs as list
    subs_list = list(zip(matching_data_list['new_names'], matching_data_list['postcode'],
                         matching_data_list['original_names'], matching_data_list['owner_cat'],
                         matching_data_list['origin']))

    # apply function
    res = Parallel(n_jobs=15, verbose=20)(
        delayed(highest_match_fuzzy_wuzzy)(word_postcode_tuple, 85, 70) for word_postcode_tuple in subs_list)

    # delete none elements from list
    res = list(filter(None.__ne__, res))
    # define col names
    def_columns = ["input_words", "alkis", "accuracy", "original_names", "owner_cat", "origin", "postcode"]

    df = pd.DataFrame(res, columns=def_columns)

    # save as csv
    df.to_csv(
        matched_owners_pth,
        sep=";", na_rep="NA",
        float_format='%.3f',
        index=True, encoding='utf-8-sig')


def combine_owners_with_agricultural_information(alkis_pth_name, agri_information_subsidies, alkis_to_dafne_names,
                                                 branch_information_dafne, alkis_pth_out):

    """
    Adds the additional information whether owners are agricultural active and from which economic branches
    the companies come to the owner dataframe.
    Args:
        alkis_pth_name: Path to cleaned owner dataframe.
        agri_information_subsidies: Path to dataframe of subsidy recipients.
        alkis_to_dafne_names: Path to dataframe with ALKIS names and DAFNE names.
        branch_information_dafne: Path to dataframe with company names and economic branches.
        alkis_pth_out: Output path to dataframe with cleaned names and additional information.

    Returns:

    """
    print("Attach information if owner is agriculturally active.")

    ## Read ALKIS data, agricultural companies from "farmsubsidies.org" and "Futtermittelbetriebe",
    ## and alkis to bvd dict
    print("\tRead input")
    owner_df = pd.read_csv(alkis_pth_name, sep=';')

    print("\tLength of input table:", len(owner_df))

    df_subs_futt = pd.read_csv(agri_information_subsidies, sep=';')
    df_conn = pd.read_csv(alkis_to_dafne_names)
    df_dafne_res = pd.read_excel(branch_information_dafne)
    cols = ['Name des Unternehmens', 'Tätigkeitsbeschreibung', 'Tatsächliche Tätigkeit',
    'WZ 2008 - Haupttätigkeit - Code',
    'WZ 2008 - Haupttätigkeit - Beschreibung',
    'WZ 2008 - Nebentätigkeit - Code',
    'WZ 2008 Nebentätigkeit - Beschreibung', 'BvD ID Nummer', ]

    df_dafne = pd.merge(df_conn, df_dafne_res[cols], how="left", left_on="bvd_id", right_on="BvD ID Nummer")

    # ## Read already matched companies
    # lst = glob.glob(rf"{WD}\{DAFNE_SEARCH_PTH}\*results_v1.xlsx")
    # df_lst = [pd.read_excel(pth, sheet_name="Ergebnisse") for pth in lst]
    # for i, df in enumerate(df_lst):
    #     df["ID"] = i
    # df_dafne = pd.concat(df_lst, axis=0)
    # df_dafne.drop(columns='Unnamed: 0', inplace=True)
    #
    # lst = glob.glob(rf"{WD}\{DAFNE_SEARCH_PTH}\*matches_v1.xlsx")
    # df_lst = [pd.read_excel(pth) for pth in lst]
    # for i, df in enumerate(df_lst):
    #     df["ID"] = i
    # df_matches = pd.concat(df_lst, axis=0)
    # df_matches.drop_duplicates(subset='Unternehmensname', inplace=True)
    # df_dafne = pd.merge(df_dafne, df_matches, how='left', left_on="BvD ID Nummer", right_on="Gematchte BvD ID")

    occup_dict_prep = {"haupt_codes": df_dafne['WZ 2008 - Haupttätigkeit - Code'].tolist(),
                       "haupt_descr": df_dafne['WZ 2008 - Haupttätigkeit - Beschreibung'].tolist(),
                       "neben_codes": df_dafne['WZ 2008 - Nebentätigkeit - Code'].tolist(),
                       "neben_descr": df_dafne['WZ 2008 Nebentätigkeit - Beschreibung'].tolist()}

    occup_dict_prep["neben_codes"] = [x for x in occup_dict_prep["neben_codes"] if str(x) != 'nan']
    occup_dict_prep["neben_descr"] = [x for x in occup_dict_prep["neben_descr"] if str(x) != 'nan']

    occup_dict_prep["neben_codes"] = [item.split('\n') for item in occup_dict_prep["neben_codes"]]
    occup_dict_prep["neben_descr"] = [item.split('\n') for item in occup_dict_prep["neben_descr"]]

    ## Flatten list of lists
    occup_dict_prep["neben_codes"] = [item for sublist in occup_dict_prep["neben_codes"] for item in sublist]
    occup_dict_prep["neben_descr"] = [item for sublist in occup_dict_prep["neben_descr"] for item in sublist]

    occup_dict = {"codes": occup_dict_prep["haupt_codes"] + occup_dict_prep["neben_codes"],
                  "descr": occup_dict_prep["haupt_descr"] + occup_dict_prep["neben_descr"]}

    occup_df = pd.DataFrame(occup_dict)
    occup_df.drop_duplicates(subset=["codes", "descr"], inplace=True)
    occup_df.dropna(inplace=True)
    occup_df['agric'] = 0
    occup_df['descr'] = occup_df['descr'].str.lower()
    occup_df.loc[occup_df['codes'].str.count('A') > 0, 'agric'] = 1
    occup_df.loc[occup_df['descr'].str.count('fischerei') > 0, 'agric'] = 0
    occup_df.loc[occup_df['descr'].str.count('folz') > 0, 'agric'] = 0
    occup_df.loc[occup_df['descr'].str.count('forst') > 0, 'agric'] = 0
    occup_df.loc[occup_df['descr'].str.count('baum') > 0, 'agric'] = 0
    occup_df.loc[occup_df['descr'].str.count('aquakultur') > 0, 'agric'] = 0
    occup_df.loc[occup_df['descr'].str.count('holz') > 0, 'agric'] = 0
    occup_df.drop_duplicates(subset="codes", inplace=True)
    # print(len(occup_df['codes'].unique()))

    df_dafne = pd.merge(df_dafne, occup_df[["codes", "agric"]], how='left', left_on='WZ 2008 - Haupttätigkeit - Code',
                        right_on="codes")
    df_dafne.rename(columns={'agric': 'agric_haupt'}, inplace=True)

    df_dafne['agric_neben'] = 0
    for code in occup_df.loc[occup_df["agric"] == 1, 'codes']:
        print(f"\tSet all useres with {code} to agricultural.")
        df_dafne.loc[df_dafne['WZ 2008 - Nebentätigkeit - Code'].str.count(code) > 0, "agric_neben"] = 1

    df_dafne.loc[df_dafne['agric_haupt'].isna(), 'agric_haupt'] = 0
    df_dafne.loc[df_dafne['agric_neben'].isna(), 'agric_neben'] = 0
    df_dafne.loc[(df_dafne['agric_neben'] + df_dafne['agric_haupt']) > 0, "agric"] = 1
    df_dafne.loc[df_dafne['agric'].isna(), 'agric'] = 0

    df_dafne.rename(columns={"comp_name": "Unternehmensname"}, inplace=True)
    df_dafne.drop_duplicates(subset='Unternehmensname', inplace=True)
    # df_dafne.drop(columns=["ID_x", "Nationale ID", "Stadt", "Land", "ID_y", "codes"], inplace=True)
    # df_dafne.to_csv(DAFNE_AGRIC_PTH, sep=';', index=False)

    df1 = df_dafne[['Unternehmensname', 'agric']].copy()
    df1.rename(columns={'Unternehmensname': 'owner_clean'}, inplace=True)
    df1 = df1.loc[df1['agric'] == 1].copy()
    dict1 = defaultdict(lambda: 0)
    for row in df1.itertuples():
        dict1[row.owner_clean] = row.agric

    df2 = df_subs_futt[['alkis', 'correct']].copy()
    df2.rename(columns={'correct': 'agric', 'alkis': 'owner_merge'}, inplace=True)
    dict2 = defaultdict(lambda: 0)
    for row in df2.itertuples():
        dict2[row.owner_merge] = row.agric

    owner_df['agric'] = 0
    owner_df['agric'] = owner_df['owner_merge_old'].map(dict2)
    print(len(owner_df.loc[owner_df['agric'] == 1]))
    owner_df.loc[owner_df['agric'] == 0, 'agric'] = owner_df.loc[owner_df['agric'] == 0, 'owner_clean'].map(dict1)
    owner_df.loc[owner_df['agric'] == 0, 'agric'] = owner_df.loc[owner_df['agric'] == 0, 'owner_names'].map(dict1)
    print(len(owner_df.loc[owner_df['agric'] == 1]))

    # t = owner_df.loc[owner_df["level1"] == 2].copy()
    # t = t[["owner_clean", "agric"]].copy()
    # t.drop_duplicates(subset="owner_clean", inplace=True)

    print(f"\tWrite output to {alkis_pth_out}")
    if 'address' in owner_df.columns:
        owner_df.drop(columns='address', inplace=True)
    owner_df.to_csv(alkis_pth_out, sep=';', index=False)

    df_out = df_dafne[["search_name", "bvd_name", "bvd_id", 'WZ 2008 - Haupttätigkeit - Code', 'WZ 2008 - Nebentätigkeit - Code', "agric"]]
    df_out.columns = ["search_name", "bvd_name", "bvd_id", 'main_branch', 'side_branch', "agric"]
    df_out.to_csv(r"07_owner_name_cleaning\companies_with_branches_and_agric.csv")

    # comp_df = owner_df.loc[owner_df["level1"] == 2].copy()
    # comp_df.to_csv(r"C:\Users\IAMO\Documents\work_data\chapter1\ALKIS\09_network_analysis\temp_company_df.csv", index=False)


def control_subsidy_recipients(alkis_clean_pth, pth_subsidies_2019, pth_subsidies_2020):
    """
    Controls if all subsidy recipients were matched. Tries to match again with different algorithm.
    Args:
        alkis_clean_pth: Path to dataframe with cleaned names and additional information. Will also be used for output.
        pth_subsidies_2019: Path to dataframe with subsidy recipients from 2019.
        pth_subsidies_2020: Path to dataframe with subsidy recipients from 2020.
    Returns:

    """
#
    print("Control if all subsidy recipients were matched.")
    ## read ALKIS
    print("\tRead data")
    owner_df = pd.read_csv(alkis_clean_pth, sep=';')
    print("\tLength of input table:", len(owner_df))

    count_agric_owners_pre = len(owner_df.loc[owner_df["agric"] == 1, "owner_merge"].unique())

    def check_additional_matches_with_subsidy_data(subsidy_df_pth, owner_df):

        ## get entries that are not yet agricultural owners
        non_agri_all = owner_df.loc[owner_df["agric"] == 0].copy()

        ## only look at privat persons and drop duplicates
        non_agri_df = non_agri_all.loc[owner_df["category"] == 1].copy()
        non_agri_df.drop_duplicates(subset="owner_clean",
                                    inplace=True)  # don't use owner_merge as multiple person can be combined under owner_merge

        ## get unique attributes of the people
        non_agri_df["postcode"] = non_agri_df["clean_address"].apply(helper_functions.identify_plz)
        non_agri_df["city"] = non_agri_df["clean_address"].apply(helper_functions.identify_city)
        non_agri_df["family_name"] = non_agri_df["owner_clean"].apply(lambda x: x.split(',')[0])
        non_agri_df["personal_name"] = non_agri_df["owner_clean"].apply(
            lambda x: x.split(',')[1] if ',' in x else 'unkown')
        non_agri_df["personal_name"] = non_agri_df["personal_name"].str.strip()
        ## create a second personal name for people with multiple surnames
        ## (as they are sometimes only listed with their first name in the subsidy data)
        non_agri_df["personal_name2"] = non_agri_df["personal_name"].apply(lambda x: x.split(' ')[0] if ' ' in x else x)

        ## create unique ID2 for combination with subsidy data
        non_agri_df["pcode_fname"] = non_agri_df["postcode"] + '_' + non_agri_df["family_name"]
        non_agri_df["pcode_fname_pname"] = non_agri_df["postcode"] + '_' + non_agri_df["family_name"] + '_' + \
                                           non_agri_df["personal_name"]
        non_agri_df["pcode_fname_pname2"] = non_agri_df["postcode"] + '_' + non_agri_df["family_name"] + '_' + \
                                            non_agri_df[
                                                "personal_name2"]
        non_agri_df["city_fname"] = non_agri_df["city"] + '_' + non_agri_df["family_name"]

        ## Open subsidy data and drop duplicates
        subs_df = pd.read_csv(subsidy_df_pth)
        subs_df.drop_duplicates(subset=["recipient_id"], inplace=True)

        ## clean name and address column to match to ALKIS format
        cols = ["recipient_name", "recipient_location"]
        for col in cols:
            subs_df[col] = subs_df[col].str.lower()
            subs_df[col] = subs_df[col].str.replace("ä", "ae")
            subs_df[col] = subs_df[col].str.replace("ö", "oe")
            subs_df[col] = subs_df[col].str.replace("ü", "ue")
            subs_df[col] = subs_df[col].str.replace("ß", "ss")

        ## Clean postcode
        subs_df["recipient_postcode"] = subs_df["recipient_postcode"].apply(lambda x: x.split('-')[1])

        ## Only use persons (identified by 1 comma in their names
        subs_df["comma_count"] = subs_df["recipient_name"].str.count(",")
        subs_persons = subs_df.loc[subs_df["comma_count"] == 1].copy()
        subs_other = subs_df.loc[subs_df["comma_count"] != 1].copy()

        ## get sur- and familiyname
        subs_persons["family_name"] = subs_persons["recipient_name"].apply(lambda x: x.split(',')[0])
        subs_persons["personal_name"] = subs_persons["recipient_name"].apply(lambda x: x.split(',')[1])
        subs_persons["personal_name"] = subs_persons["personal_name"].str.strip()

        ## split entries with enumeration of person
        ## create separate entries and attach again to common df
        subs_persons_unds = subs_persons.loc[subs_persons["personal_name"].str.count(" und") > 0].copy()
        subs_persons = subs_persons.loc[subs_persons["personal_name"].str.count(" und") == 0].copy()
        subs_persons_unds1 = subs_persons_unds.copy()
        subs_persons_unds2 = subs_persons_unds.copy()
        subs_persons_unds1["personal_name"] = subs_persons_unds1["personal_name"].apply(lambda x: x.split(" und ")[0])
        subs_persons_unds2["personal_name"] = subs_persons_unds2["personal_name"].apply(lambda x: x.split(" und ")[1] if ' und ' in x else x)
        subs_persons_unds = pd.concat([subs_persons_unds1, subs_persons_unds2], axis=0)
        subs_persons_unds.sort_values(by="recipient_name", inplace=True)
        subs_persons = pd.concat([subs_persons, subs_persons_unds], axis=0)

        ## create unique ID for combination with ALKIS data
        subs_persons["pcode_fname"] = subs_persons["recipient_postcode"] + '_' + subs_persons["family_name"]
        subs_persons["pcode_fname_pname"] = subs_persons["recipient_postcode"] + '_' + subs_persons["family_name"] + '_' + subs_persons[
            "personal_name"]
        subs_persons["city_fname"] = subs_persons["recipient_location"] + '_' + subs_persons["family_name"]

        ## get list of unique IDs (simple) from ALKIS and subset the df with persons to only these IDs
        unique_postcodes = list(non_agri_df["pcode_fname"].unique())
        subs_persons = subs_persons.loc[subs_persons["pcode_fname"].isin(unique_postcodes)].copy()

        ## merge subsidy person data with ALKIS data with unique ID (advanced) to identify agricultural owners
        df_merge = pd.merge(subs_persons, non_agri_df, how="inner", on="pcode_fname_pname")
        print("\tFirst round:", len(df_merge))

        ## identify ALKIS owners that aren't agricultural and try again with second unique ID (advanced)
        matches = list(df_merge["pcode_fname_pname"].unique())
        df_miss = non_agri_df.loc[~non_agri_df["pcode_fname_pname"].isin(matches)].copy()
        df_merge2 = pd.merge(subs_persons, df_miss, how="inner", left_on="pcode_fname_pname", right_on="pcode_fname_pname2")

        df_merge = pd.concat([df_merge, df_merge2], axis=0)
        df_merge.index = range(len(df_merge))
        print("\tSecond round:", len(df_merge))

        # ## identify ALKIS owners that aren't agricultural and try again with third unique ID (advanced)
        # matches = list(df_merge["pcode_fname_pname2"].unique())
        # df_miss = non_agri_df.loc[~non_agri_df["pcode_fname_pname2"].isin(matches)].copy()
        # df_merge3 = pd.merge(subs_persons, df_miss, how="inner", on="city_fname")
        #
        # df_merge = pd.concat([df_merge, df_merge3], axis=0)
        # df_merge.index = range(len(df_merge))
        # print("Third round:", len(df_merge))

        df_merge["agric"] = 1

        ########################################
        non_agri_oth = non_agri_all.loc[owner_df["level3"] != '1_1_1'].copy()
        non_agri_oth.drop_duplicates(subset="owner_clean",
                                     inplace=True)  # don't use owner_merge as multiple person can be combined under owner_merge

        ## create unique ID for combination with subsidy data
        non_agri_oth["postcode"] = non_agri_oth["clean_address"].apply(helper_functions.identify_plz)
        non_agri_oth["pcode_name"] = non_agri_oth["postcode"] + '_' + non_agri_oth["owner_clean"]

        unique_postcodes = list(non_agri_oth["postcode"].unique())

        subs_other = subs_other.dropna(subset=["recipient_name"]).copy()
        subs_other["recipient_name"] = subs_other["recipient_name"].str.replace("gdbr", "gbr")

        ## create unique ID2 for combination with ALKIS data
        subs_other = subs_other.loc[subs_other["recipient_name"].str.count(" gbr") > 0].copy()
        subs_other["pcode_name"] = subs_other["recipient_postcode"] + '_' + subs_other["recipient_name"]

        ## merge subsidy person data with ALKIS data with unique ID (advanced) to identify agricultural owners
        oth_merge = pd.merge(subs_other, non_agri_oth, how="inner", on="pcode_name")

        ## identify ALKIS owners that aren't agricultural and try again with second unique ID (advanced)
        matches_oth = list(oth_merge["pcode_name"].unique())
        oth_miss = non_agri_oth.loc[~non_agri_oth["pcode_name"].isin(matches_oth)].copy()

        subs_oth = subs_other.loc[subs_other["recipient_postcode"].isin(unique_postcodes)]

        oth_miss["agric"] = 0

        for row in subs_oth.itertuples():
            pcode = row.recipient_postcode
            name = row.recipient_name
            name_lst = name.split(' ')
            name_lst = [s.replace('.', '') for s in name_lst]
            name_lst = [s.replace('/', '') for s in name_lst]
            name_lst = [s for s in name_lst if len(s) > 2]

            for s in ["gbr", "&", "+", "", "und", "-", "mbh", "gmbh"]:
                if s in name_lst:
                    name_lst.remove(s)

            if len(name_lst) == 0:
                continue

            base = r'^{}'
            expr = '(?=.*{})'
            search_string = base.format(''.join(expr.format(w) for w in name_lst))

            oth_miss.loc[(oth_miss["owner_clean"].str.contains(search_string)) & (
                    oth_miss["postcode"] == pcode), "agric"] = 1

            # matches = oth_miss.loc[(oth_miss["owner_clean"].str.contains(search_string)) & (
            #             oth_miss["postcode"] == pcode), "owner_clean"].tolist()
            # if len(matches) > 0:
            #     print("\n", name, f"Matchwords:  {search_string}", matches)
            # else:
            #     name1 = name_lst[0]
            #     if len(name1) > 2:
            #         matches = oth_miss.loc[(oth_miss["owner_clean"].str.contains(name1)) & (
            #                 oth_miss["postcode"] == pcode), "owner_clean"].tolist()
            #         if len(matches) > 0:
            #             print_term = f"\n{name} Matchword:  {name1} Not-working Matchstring {search_string}"
            #             print_red(print_term)
            #             print(matches)

            # oth_miss.loc[(oth_miss["owner_clean"].str.contains(search_string)) & (
            #         oth_miss["postcode"] == pcode), "owner_clean"] = name

        oth_merge2 = oth_miss.loc[oth_miss["agric"] == 1]


        oth_merge = oth_merge[["owner_merge", "agric"]].copy()
        oth_merge2 = oth_merge2[["owner_merge", "agric"]].copy()
        df_merge = df_merge[["owner_merge", "agric"]].copy()

        df_merge = pd.concat([df_merge, oth_merge, oth_merge2], axis=0)

        return df_merge

    df_merge1 = check_additional_matches_with_subsidy_data(pth_subsidies_2019, owner_df)
    df_merge1["origin"] = "farmsubsidy_2019"
    df_merge2 = check_additional_matches_with_subsidy_data(pth_subsidies_2020, owner_df)
    df_merge2["origin"] = "farmsubsidy_2020"
    df_merge = pd.concat([df_merge1, df_merge2], axis=0)
    df_merge.drop_duplicates(subset=["owner_merge"], inplace=True)
    df_merge.to_csv(r"07_owner_name_cleaning\07_agricultural_owners_subsidies.csv", sep=';', index=False)

    owner_merges = list(df_merge["owner_merge"].unique())
    owner_df.loc[owner_df["owner_merge"].isin(owner_merges), "agric"] = 1
    count_agric_owners_post = len(owner_df.loc[owner_df["agric"] == 1, "owner_merge"].unique())
    print(f"\tNumber of agricultural owners before correction: {count_agric_owners_pre}")
    print(f"\tNumber of agricultural owners after correction:  {count_agric_owners_post}")

    print("\tLength of output table:", len(owner_df))
    print(f"\tWrite output to {alkis_clean_pth}")
    owner_df.to_csv(alkis_clean_pth, sep=';', index=False)


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    #################################### Cleaning procedure ####################################
    create_aristrocracy_dict(
        alkis_pth=ALKIS_PTH,
        out_pth=r"07_owner_name_cleaning\cleaning_dict_aristrocracy.json"
    )
    ## Some manual adjustments
    manually_clean_owners(
        alkis_pth=ALKIS_PTH,
        manual_cleaning_pth=PTH_CLEANING_MANUAL,
        alkis_clean_pth=PTH_ALKIS_CLEAN_PRELIMINARY
    )
    ## DAFNE company names
    clean_company_names(
        alkis_clean_pth=PTH_ALKIS_CLEAN_PRELIMINARY,
        alkis_to_dafne_names_pth=PTH_ALKIS_TO_DAFNE_NAMES
    )
    ## Manual adjustments for churches
    clean_churchs_nonprof_public_names(
        alkis_pth_name=PTH_ALKIS_CLEAN_PRELIMINARY,
        clean_table_name=PTH_CLEANING_CHURCHES
    )
    ## Manual adjustments for non-profit etc.
    clean_churchs_nonprof_public_names(
        alkis_pth_name=PTH_ALKIS_CLEAN_PRELIMINARY,
        clean_table_name=PTH_CLEANING_NONPROF
    )
    ## Manual adjustments for public
    clean_churchs_nonprof_public_names(
        alkis_pth_name=PTH_ALKIS_CLEAN_PRELIMINARY,
        clean_table_name=PTH_CLEANING_PUBLIC
    )
    ## Manual and automatic adjustments for private people
    clean_private_people(
        alkis_pth_name=PTH_ALKIS_CLEAN_PRELIMINARY,
        clean_dict_pth=PTH_CLEANING_PEOPLE,
        clean_dict_pth_manual=PTH_CLEANING_PEOPLE_MANUAL
    )
    ## Manual adjustments for others
    clean_churchs_nonprof_public_names(
        alkis_pth_name=PTH_ALKIS_CLEAN_PRELIMINARY,
        clean_table_name=PTH_CLEANING_REST
    )
    ## Quick cleaning of GbRs
    clean_gbrs(alkis_pth_name=PTH_ALKIS_CLEAN_PRELIMINARY
               )
    ## Creation of new unique identifier and unifying of addresses
    create_new_unique_identifier_and_unify_addresses(
        alkis_pth_name=PTH_ALKIS_CLEAN_PRELIMINARY
    )

    #################################### Add supplementary information ####################################

    ## 1. Match owners with data frames on recipients of agricultural subsidies
    ## This part is only in here for documentation reasons. It was done in a very early stage of the project and
    ## includes some steps that could not be reconstructed.
    ## Written by Tilman Schmitz.
    # # # helper_functions.print_red("!!!! You have to run 07_data_preparations_alkis_farmsubsidy.R in R before proceeding !!!!")
    # # #
    # # # fuzzy_matching_alkis_data(
    # # #     alkis_path=r"07_owner_name_cleaning\string_match_alkis_farmer_information\alkis_2020_unique_owner_clean.csv",
    # # #     matching_data=r"07_owner_name_cleaning\string_match_alkis_farmer_information\binded_output_all_dfs.csv",
    # # #     matched_owners_pth=r"07_owner_name_cleaning\result_matching_alkis_binded_only_private.csv"
    # # # )

    combine_owners_with_agricultural_information(
        alkis_pth_name=PTH_ALKIS_CLEAN_PRELIMINARY,
        agri_information_subsidies=AGRICULTURAL_OWNERS_PTH,
        branch_information_dafne=DAFNE_SEARCH_RESULTS_PTH,
        alkis_to_dafne_names=PTH_ALKIS_TO_DAFNE_NAMES,
        alkis_pth_out=PTH_ALKIS_CLEAN
    )

    control_subsidy_recipients(
        alkis_clean_pth=PTH_ALKIS_CLEAN,
        pth_subsidies_2019=PTH_SUBSIDIES_2019,
        pth_subsidies_2020=PTH_SUBSIDIES_2020
    )

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)



if __name__ == '__main__':
    main()
