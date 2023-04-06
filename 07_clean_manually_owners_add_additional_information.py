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

## Project library
import helper_functions
# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

PTH_ALKIS = r"05_calculate_geolocation\05_owners_stretched_addresses_geolocated_distances.csv"
PTH_CLEANING_MANUAL = r"07_owner_name_cleaning\cleaning_table_manuelle_namensbereinigung.xlsx"
PTH_CLEANING_CHURCHES = r"07_owner_name_cleaning\cleaning_table_churches.xlsx"
PTH_CLEANING_NONPROF = r"07_owner_name_cleaning\cleaning_table_non-profit.xlsx"
PTH_CLEANING_PUBLIC = r"07_owner_name_cleaning\cleaning_table_public.xlsx"
PTH_CLEANING_PEOPLE = r"07_owner_name_cleaning\cleaning_dict_private_people.json"
PTH_CLEANING_PEOPLE_MANUAL = r"07_owner_name_cleaning\cleaning_dict_private_people_manual_assignment.json"
PTH_CLEANING_REST = r"07_owner_name_cleaning\cleaning_table_rest.xlsx"
CONNECTION_PTH = r"07_owner_name_cleaning\__name-connections_clean.csv" ## Alkis names to Dafne name

AGRI_COMP_PTH = r"07_owner_name_cleaning\matching_results_alkis_farmsubs_futtermittel.csv"
DAFNE_PTH = r"08_network_analysis\Unternehmensdaten.xlsx"

## Path for final output
PTH_ALKIS_CLEAN = r"07_owner_name_cleaning\07_owners_stretched_pre.csv"

# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#

def manually_clean_owners():

    print(f"##################################################\n"
          f"Correct owner names with df {PTH_CLEANING_MANUAL}\n"
          f"##################################################\n")

    df_alk = pd.read_csv(PTH_ALKIS, sep=';')
    df_clean = pd.read_excel(PTH_CLEANING_MANUAL)

    print("Length of table:", len(df_alk))

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

    df_alk.to_csv(PTH_ALKIS_CLEAN, sep=';', index=False)


def clean_company_names():
    print("##################################################\n"
          "Clean company names, i.e. replace ALKIS names with DAFNE names.\n"
          "##################################################\n")

    ## Read data
    os.chdir(WD)
    owner_df = pd.read_csv(PTH_ALKIS_CLEAN, sep=';')
    df_conn = pd.read_csv(CONNECTION_PTH)

    print("Length of table:", len(owner_df))

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
    print("Number of companies after unification:", num_comp)

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

    owner_df.to_csv(PTH_ALKIS_CLEAN, sep=';', index=False)



def clean_churchs_nonprof_public_names(alkis_pth_name, clean_table_name, alkis_clean_out_pth):
    print(f"##################################################\n"
          f"Correct owner names with df {clean_table_name}\n"
          f"##################################################\n")

    df_clean = read_table_to_df(clean_table_name) #pd.read_excel(clean_table_name)
    df_alk = read_table_to_df(alkis_pth_name)  # pd.read_csv(alkis_pth_name, sep=';')

    print("Length of table:", len(df_alk))

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
    df_alk.to_csv(alkis_clean_out_pth, sep=';', index=False)



def create_cleaning_dictionary_private_persons(df_alk, clean_dict_pth=None):
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

    ## count occurence of this combination, if a combination occurs multiple times that it indicates
    ## the possible occurence of a name with different writings
    counts = Counter(sub["comb"].tolist())
    threshold = 1
    lst = []
    for k, v in counts.items():
        if v > threshold:
            lst.append(k)
    df_work = sub.loc[sub["comb"].isin(lst)].copy()
    df_work.reset_index(inplace=True)
    df_work.drop(columns=['index'], inplace=True)
    df_done = sub.loc[~sub["comb"].isin(lst)].copy()

    ## open a json file to write the results in
    if clean_dict_pth:
        with open(clean_dict_pth, 'w', encoding='ISO-8859-1') as file:
            file.write("{\n")

    correction_dict = {}
    print("loop")
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
            ## if possible addresses fullfill certain quality requirement as indicated by code 5,6,7
            ## then use the best address. If not use any
            sub4 = sub2.loc[sub2["full_address"].isin([5, 6, 7])]
            if not sub4.empty:
                max_f = sub4["full_address"].max()
                clean_address = sub4.loc[sub4["full_address"] == max_f, "clean_address"].iloc[0]
            else:
                clean_address = sub2["clean_address"].iloc[0]

            ## for all persons left add information to dictionary
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

def identify_city(text):
    import re
    text = text.lower()

    ## identify street names and postal codes + city name
    lst = text.split(',')
    city = ''
    for l, sub in enumerate(lst):
        ## look for a postcal code + city name, if there hasn't been found any yet
        if city == '':
            ## look  for 5-digit number and a subsequent word
            if re.search(r'\d{5} \b\w+\b', sub):
                city = re.findall(r'\D{1,100}', sub)
                city = ' '.join(city)
        if city == '':
            ## if there was no 5-digit number, try a 4-digit number + word (e.g. Switzerland has 4-digit postal codes)
            if re.search(r'\d{4} \b\w+\b', sub):
                city = re.findall(r'\D{1,100}', sub)
                city = ' '.join(city)

    city = city.strip()

    return city

def get_adlige_dict():
    df_alk = pd.read_csv(PTH_ALKIS, sep=";")

    df_alk["conn_left"] = df_alk["owner_clean"]

    ## For private owner create merge name with pattern: "surname familyname_location" to merge with dafne entries
    df_alk.loc[(df_alk["level3"] == "1_1_1") & (df_alk["conn_left"].str.count(',') > 0), "conn_left"] = \
        df_alk.loc[(df_alk["level3"] == "1_1_1") & (df_alk["conn_left"].str.count(',') > 0), "conn_left"].apply(
            lambda row: row.split(',')[1] + ' ' + row.split(',')[0])

    df_alk["city"] = df_alk["clean_address"].apply(identify_city)
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

    print("Number of people to search:", len(t))

    out_dict = {}
    sterm_lst = t["search_term"].tolist()
    for i, search_term in enumerate(sterm_lst):
        print(i, search_term)
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
            print("Number matches:", len(t2) - 1)
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

    with open(r"07_owner_name_cleaning\cleaning_dict_aristrocracy.json", 'w', encoding='ISO-8859-1') as file:
        json.dump(out_dict, file, indent=4)

    print("Clean out dict manually and rename!!")

    # print("done!")
    # l1 = [out_dict[key]["clean_name"] for key in out_dict]
    # l1 = list(set(l1))
    # l2 = [key for key in out_dict]
    # l2 = list(set(l2))
    # l = l1 + l2




def clean_private_people(alkis_pth_name, clean_dict_pth, clean_dict_pth_manual, alkis_clean_out_pth):

    os.chdir(WD)

    print(f"##################################################\n"
          f"Correct private owner names based on family names, parts of surnames and birthdates\n"
          f"##################################################\n")
    df_alk = read_table_to_df(alkis_pth_name)

    print("Length of table:", len(df_alk))

    if not os.path.exists(clean_dict_pth):
        clean_dict_pre = create_cleaning_dictionary_private_persons(df_alk, clean_dict_pth=clean_dict_pth)
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
            print("Clean name from manual_clean_dict is key (wrong name) in aut_clean_dict. Remove:", manual_val)
            clean_dict_pre.pop(manual_val)
        if key_manual in clean_dict_pre:
            print("Wrong name from manual_clean_dict is key (wrong name) in aut_clean_dict. Unify.", key_manual)
            clean_dict_pre[key_manual]["clean_name"] = manual_val
            clean_dict_pre[key_manual]["clean_address"] = manual_adr

    for key_manual in clean_dict_manual:
        manual_val = clean_dict_manual[key_manual]["clean_name"]
        manual_adr = clean_dict_manual[key_manual]["clean_address"]
        for key_pre in clean_dict_pre:
            if clean_dict_pre[key_pre]["clean_name"] == key_manual:
                print("Wrong name from manual_clean_dict is value (correct name) in aut_clean_dict. Unify.", key_manual)
                clean_dict_pre[key_pre]["clean_name"] = manual_val
                clean_dict_pre[key_pre]["clean_address"] = manual_adr

    clean_dict = clean_dict_pre | clean_dict_manual

    # for key in clean_dict_manual:
    #     if key not in clean_dict:
    #         clean_dict[key] = clean_dict_manual[key]
    # t = list(v for v in clean_dict if 'matthes' in v)

    ## check matthes reiner
    clean_owner_names = []
    wrong_owner_names = []
    for key in clean_dict:
        clean_owner_names.append(clean_dict[key]["clean_name"])
        wrong_owner_names.append(key)
    clean_owner_names = set(list(clean_owner_names))
    wrong_owner_names = set(list(wrong_owner_names))

    ## get a subset of all entries with wrong names
    df_wrong1 = df_alk.loc[df_alk["owner_clean"].isin(wrong_owner_names)].copy()
    print("Number of entries with wrong names", len(df_wrong1))

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
        print("All entries cleaned.")
    else:
        print("Not all entries cleaned.")
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

    df_alk.to_csv(alkis_clean_out_pth, sep=';', index=False)


def clean_gbrs(alkis_pth_name, alkis_clean_out_pth):
    df_alk = pd.read_csv(alkis_pth_name, sep=';')

    print("Length of table:", len(df_alk))

    ## Clean GBRs
    df_alk.loc[df_alk["level3"] == '1_2_1', 'owner_clean'] = df_alk.loc[
        df_alk["level3"] == '1_2_1', 'owner_clean'].str.replace('bestehend aus', '', regex=False)
    df_alk.loc[df_alk["level3"] == '1_2_1', 'owner_clean'] = df_alk.loc[
        df_alk["level3"] == '1_2_1', 'owner_clean'].str.replace('bestehen aus', '', regex=False)

    t = df_alk.loc[df_alk["level3"] == '1_2_1'].copy()
    t.drop_duplicates(subset="owner_clean", inplace=True)
    t.sort_values(by="owner_clean", inplace=True)

    df_alk.to_csv(alkis_clean_out_pth, sep=';', index=False)

def clean_addresses_and_owner_merge_names(alkis_pth_name, alkis_clean_out_pth):
    ## Open data
    df_alk = pd.read_csv(alkis_pth_name, sep=';')

    print("Length of table:", len(df_alk))

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
            print(f"removed '{oc}' from multis")
            multis.remove(oc)

    ## To the magic. Unify addresses
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

    ## Get all occurences of non privat persons entries where for one owner merge there are multiple owner clean entries
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

    ## In cases where in owner merge there is the word "unbekannt" (likely due to unkown address),
    ## replace the owner_merge with owner clean
    df_alk.loc[df_alk["owner_merge"].str.count("unbekannt") > 0, "owner_merge"] = \
        df_alk.loc[df_alk["owner_merge"].str.count("unbekannt") > 0, "owner_clean"]

    if not df_red.empty:
        df_red[["level1", "owner_clean", "owner_merge", "clean_address"]].to_csv(fr"07_owner_name_cleaning\all_levels_cleaning.csv", sep=",", index=False)
    else:
        print("table empty, writing table out. Length of final table:", len(df_alk))
        df_alk.to_csv(alkis_clean_out_pth, sep=';', index=False)

    t = df_alk.loc[df_alk["level3"] == '1_2_1'].copy()
    t.drop_duplicates(subset="owner_clean", inplace=True)
    t.sort_values(by="owner_clean", inplace=True)

    uni_levels = df_red["level1"].unique()




def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)


    # get_adlige_dict()

    # manually_clean_owners()
    # clean_company_names()
    # clean_churchs_nonprof_public_names(alkis_pth_name=PTH_ALKIS_CLEAN,
    #                                    clean_table_name=PTH_CLEANING_CHURCHES,
    #                                    alkis_clean_out_pth=PTH_ALKIS_CLEAN)
    # clean_churchs_nonprof_public_names(alkis_pth_name=PTH_ALKIS_CLEAN,
    #                                    clean_table_name=PTH_CLEANING_NONPROF,
    #                                    alkis_clean_out_pth=PTH_ALKIS_CLEAN)
    # clean_churchs_nonprof_public_names(alkis_pth_name=PTH_ALKIS_CLEAN,
    #                                    clean_table_name=PTH_CLEANING_PUBLIC,
    #                                    alkis_clean_out_pth=PTH_ALKIS_CLEAN)
    # clean_private_people(alkis_pth_name=PTH_ALKIS_CLEAN,
    #                      clean_dict_pth=PTH_CLEANING_PEOPLE,
    #                      clean_dict_pth_manual=PTH_CLEANING_PEOPLE_MANUAL,
    #                      alkis_clean_out_pth=PTH_ALKIS_CLEAN)
    #
    # clean_churchs_nonprof_public_names(alkis_pth_name=PTH_ALKIS_CLEAN,
    #                                    clean_table_name=PTH_CLEANING_REST,
    #                                    alkis_clean_out_pth=PTH_ALKIS_CLEAN)
    #
    # clean_gbrs(alkis_pth_name=PTH_ALKIS_CLEAN,
    #            alkis_clean_out_pth=PTH_ALKIS_CLEAN)

    clean_addresses_and_owner_merge_names(alkis_pth_name=PTH_ALKIS_CLEAN,
                                          alkis_clean_out_pth=PTH_ALKIS_CLEAN)


    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)



if __name__ == '__main__':
    main()
