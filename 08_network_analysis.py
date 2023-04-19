# Clemens Jänicke
# github Repo: https://github.com/clejae

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#

import time
import math
import pandas as pd
import os
import json
import glob
import networkx as nx
from community import community_louvain
import matplotlib.pyplot as plt

## Project library
import helper_functions
import plotting_lib

# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input paths
DAFNE_PTH = r"06_prepare_dafne_search\combined_search_results\Unternehmensdaten.xlsx"
ALKIS_PTH = r"07_owner_name_cleaning\07_owners_stretched.csv"
CORRECTED_MATCHES_PTH = r"08_network_analysis\network_private_persons_matches_in_alkis_corrected.json"
FOLDER_DAFNE_SECOND_SEARCH = r"08_network_analysis\second_search" #Will be create automatically
# DONE_COMP_WITH_BRANCHES = r"07_owner_name_cleaning\companies_with_branches_and_agric.csv"

## Output paths
COL_RENAMING = r"08_network_analysis\renaming_columns_first_dafne_search.json"
NETW_LONG_PTH = r"08_network_analysis\network_connections_long.csv"
NETW_PTH = r"08_network_analysis\network_connections.csv"
POSSIBLE_MATCHES_FOLDER = r"08_network_analysis"
DAFNE_SEARCH_PTH = r"08_network_analysis\dafne_search_second_round.csv"
NETW_COMM_PTH = r"08_network_analysis\network_connections_with_community.csv"
COMP_COMMUNITY_MCOMP_DICT = r"08_network_analysis\owners+communities_thr{0}\08_comp_comm_mcomp_thr{0}_dict.json"
COMMUNITY_INFO_DICT = r"08_network_analysis\owners+communities_thr{0}\08_alkis_owners_comm_thr{0}_info_dict.json"
OWNERS_STRETCHED_COMM = r"08_network_analysis\owners+communities_thr{0}\08_owners_stretched+comm_thr{0}.csv"
COMMUNITY_MAX_DIST_DICT_PTH = r"08_network_analysis\community{0}_max_dist_dict.json"

# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#


def check_share_is_number(share_var):
    """
    Checks whether the provided share string from DAFNE is a number.
    Args:
        share_var: Share string.

    Returns:

    """
    share = str(share_var)

    ## GP - general partner - komplementär
    ## WO - Wholly Owned - Totalbesitz >= 98.-%
    ## MO - Majority Owned- Mehrheitsbeteiligung > 50.00%
    ## T
    ## VE
    ## FME (only occurs in TGS2_dshare)
    ## ADV (only occurs in TGS2_dshare)
    ## NG - Negligible - Vernachlässigbar <= 0.01%

    def has_numbers(input_string):

        check = any(char.isdigit() for char in input_string)

        return check

    if has_numbers(share):
        share = share.replace(',', '.')
        if '>' in share:
            share = share.replace('>', '')
            share = float(share) + 0.09999
        share = float(share)
    elif share in ['n.v.', '-', "GP", "T", "FME", "ADV", "VE", "NG", 'nan']:
        share = math.nan
    elif share == "WO":
        share = 98.09999
    elif share == "MO":
        share = 50.09999
    elif share == "NG":
        share = 0.009

    return share


def company_information_to_long_data_frame(dafne_search_results, col_renaming_json_pth, network_data_long_df_pth):
    """
    Reformats the data frame with the DAFNE search results into a long data frame for easier processing later on.
    Args:
        dafne_search_results: Path to dataframe with DAFNE search results.
        col_renaming_json_pth: Output path to json with dictionary of columns renaming.
        network_data_long_df_pth: Output path to long DAFNE dataframe.

    Returns:

    """
    print("Reformat company information to long data frame.")

    print("\tRead data.")
    df = pd.read_excel(dafne_search_results)

    col_names = {'Unnamed: 0': "index",
                 'Name des Unternehmens': "comp_name",
                 'Ort / Gemeinde': "comp_munic",
                 'Vollständiger Name': "comp_name_full",
                 'Postleitzahl': "post_code",
                 'Ort': "comp_loc",
                 'Handelsregisternummer': "hrn",
                 'Tätigkeitsbeschreibung': "occup",
                 'Tatsächliche Tätigkeit': "real_occ",
                 'MAN - Name (*)': "MAN_name",
                 'MAN - Geburtsjahr (*)': "MAN_year",
                 'MAN - Ort (*)': "MAN_loc",
                 'MAN - Funktion (*)': "MAN_func",
                 'Funktion (*)': 'func',
                 'DM\nVollständiger Name (*)': "DM_name",
                 'DM\nAktuell oder ehemalig (*)': "DM_current",
                 'DM\nAlter (*)': "DM_age",
                 'DM\nAdresse (*)': "DM_addr",
                 'DM\nArt der Position (*)': "DM_posit",
                 'ATE - Unternehmen/Person (*)': "ATE_type",
                 'ATE - Name (*)': "ATE_name",
                 'ATE - Funktion (*)': "ATE_func",
                 'ATE - Prozent gehalten\n% (*)': "ATE_share",
                 'SHH - Geschäftsführender Gesellschafter (*)': "ATE_man",
                 'MAN - Strasse und Hausnummer (*)': "MAN_street",
                 'ENT - Name (*)': "ENT_name",
                 'ENT - Strasse und Hausnummer (*)': "ENT_street",
                 'ENT - Ort (*)': "ENT_loc",
                 'ENT - Geburtsjahr (*)': "ENT_year",
                 'COS - Name (*)': "COS_name",
                 'COS - Strasse und Hausnummer (*)': "COS_street",
                 'COS - Ort (*)': "COS_loc",
                 'COS - Postleitzahl (*)': "COS_pcode",
                 'ENT - Postleitzahl (*)': "ENT_pcode",
                 'MAN - Postleitzahl (*)': "MAN_pcode",
                 'TGS - Name (*)': "TGS_name",
                 'TGS - Strasse und Hausnummer (*)': "TGS_street",
                 'TGS - Ort (*)': "TGS_loc",
                 'TGS - Postleitzahl (*)': "TGS_pcode",
                 'TGS - Prozente gehalten\n% (*)': "TGS_share",
                 'UFR - Name (*)': "UFR_name",
                 'UFR - Strasse und Hausnummer (*)': "UFR_street",
                 'UFR - Postleitzahl (*)': "UFR_pcode",
                 'UFR - Ort (*)': "UFR_loc",
                 'Anzahl Unternehmen in der Konzerngruppe (*)': "num_comp",
                 'Gesellschafter-\nName (*)': "shold_name",
                 'Gesellschafter - Direkt % (*)': "shold_dshare",
                 'Gesellschafter - Gesamt % (*)': "shold_tshare",
                 'Gesellschafter - Ort (*)': "shold_loc",
                 'Gesellschafter - Postleitzahl (*)': "shold_pcode",
                 'ISH - Name (*)': "ISH_Name",
                 'ISH - Direkte % (*)': "ISH_dshare",
                 'ISH - Gesamt % (*)': "ISH_tshare",
                 'ISH - Ort (*)': "ISH_loc",
                 'ISH - Postleitzahl (*)': "ISH_pcode",
                 'Gesellschafter - Auch Manager (*)': "shold_man",
                 'Gesellschafter - Auch Manager (*).1': "shold_man2",
                 'Globale KM - Name (*)': "KMg_name",
                 'Globale  KM - Direkt% (*)': "KMg_dshare",
                 'Globale KM - Gesamt% (*)': "KMg_tshare",
                 'Gesellschafter - Auch Manager (*).2': "KMg_man",
                 'Globale KM - Ort (*)': "KMg_loc",
                 'Globale KM - Postleitzahl (*)': "KMg_pcode",
                 'Nationale KM - Name (*)': "KMn_name",
                 'Nationale KM - Direkt% (*)': "KMn_dshare",
                 'Nationale KM - Gesamt (*)': "KMn_tshare",
                 'Gesellschafter - Auch Manager (*).3': "KMn_man",
                 'Nationale KM - Ort (*)': "KMn_loc",
                 'Nationale KM - Postleitzahl (*)': "KM_pcode",
                 'CSH - Name (*)': "CSH_name",
                 'CSH - Direkt % (*)': "CSH_dshare",
                 'CSH - Gesamt % (*)': "CSH_tshare",
                 'Gesellschafter - Auch Manager (*).4': "CSH_man",
                 'CSH - Ort (*)': "CSH_loc",
                 'CSH - Postleitzahl (*)': "CSH_pcode",
                 'Anzahl der dokumentierten Tochtergesellschaften (*)': "num_subcomp",
                 'Textinformation Tochtergesellschaften (*)': "subcomp_desc",
                 'Tochterges. - Name (*)': "TGS2_name",
                 'Tochterges. - Direkt %  (*)': "TGS2_dshare",
                 'Tochterges. - Gesamt (*)': "TGS2_tshare",
                 'Tochtergesellschaft - Ort (*)': "TGS2_loc",
                 'Tochtergesellschaft - Postleitzahl (*)': "TGS2_pcode",
                 'BvD ID Nummer': "bvd_id"}

    with open(file=col_renaming_json_pth, mode="w") as outfile:
        json.dump(col_names, outfile, indent=4)

    df.rename(columns=col_names, inplace=True)

    # df = df[1071:1072]
    cols = ["comp_name", "CSH_name", "KMn_name", "KMg_name", "TGS_name", "TGS2_name", "MAN_name", "ATE_name", "DM_name"]

    ## Replace Umlaute & All names to lower case letters only
    for col in cols:
        df[col] = df[col].str.lower()
        df[col] = df[col].str.replace("ä", "ae")
        df[col] = df[col].str.replace("ö", "oe")
        df[col] = df[col].str.replace("ü", "ue")
        df[col] = df[col].str.replace("ß", "ss")
        df[col] = df[col].str.replace("  ", " ")

    all_conns = {"company": [],
                 "connection": [],
                 "attribute": [],
                 "distance": [],
                 "share": [],
                 "ismanager": [],
                 "person": [],
                 "function": [],
                 "location": []}

    for row in df.itertuples():
        comp_name = row.comp_name

        ## Get chain of shareholders
        cshs = row.CSH_name
        cshs_shares = row.CSH_dshare
        cshs_mans = row.CSH_man

        ## If there is a chain of shareholder, then add all connections
        if type(cshs) == str:
            cshs = cshs.split("\n")
            num_cshs = len(cshs)
            if num_cshs > 1:
                cshs_shares = cshs_shares.split("\n")
                cshs_mans = cshs_mans.split("\n")
            else:
                ## Put single share into a list to be able to retrieve it with an index later on
                cshs_shares = [cshs_shares]
                cshs_mans = [cshs_mans]

            ## go the ladder up (reverse order), the first name is the current mother company, the last is the global owner
            for i in reversed(range(num_cshs)):
                csh = cshs[i]

                ## check if its a person
                csh_split = csh.split(' ')
                if csh_split[0] in ["mr", "mrs"]:
                    person = 'Person'
                    csh = ' '.join(csh_split[1:])
                else:
                    person = 'Unternehmen'

                share = cshs_shares[i]
                check_share_is_number(share)
                ismanager = cshs_mans[i]
                loc = ''

                all_conns["company"].append(comp_name)
                all_conns["connection"].append(csh)
                all_conns["attribute"].append("CSH")
                all_conns["distance"].append(num_cshs - i)
                all_conns["share"].append(share)
                all_conns["ismanager"].append(ismanager)
                all_conns["person"].append(person)
                all_conns["function"].append("mother company")
                all_conns["location"].append(loc)

        ## Get national parent company and add it to the list if its not current company itself
        pcompn_name = row.KMn_name
        pcompn_man = row.KMn_man
        share = row.KMn_tshare

        if type(pcompn_name) == str:
            ## check if its a person
            pcompn_name_split = pcompn_name.split(' ')
            if pcompn_name_split[0] in ["mr", "mrs"]:
                person = 'Person'
                pcompn_name = ' '.join(pcompn_name_split[1:])
            else:
                person = 'Unternehmen'

            ## if it is not the current company itself then add
            if comp_name != pcompn_name:
                ismanager = pcompn_man
                share = check_share_is_number(share)
                loc = ''

                all_conns["company"].append(comp_name)
                all_conns["connection"].append(pcompn_name)
                all_conns["attribute"].append("KMn")
                all_conns["distance"].append(1)
                all_conns["share"].append(share)
                all_conns["ismanager"].append(ismanager)
                all_conns["person"].append(person)
                all_conns["function"].append("mother company")
                all_conns["location"].append(loc)

        ## Get global parent company and add it to the list if its not the national parent company
        pcompg_name = row.KMg_name
        pcompg_man = row.KMg_man
        share = row.KMg_tshare

        if type(pcompg_name) == str:
            ## check if its a person
            pcompg_name_split = pcompg_name.split(' ')
            if pcompg_name_split[0] in ["mr", "mrs"]:
                person = 'Person'
                pcompg_name = ' '.join(pcompg_name_split[1:])
            else:
                person = 'Unternehmen'

            ## if it is not the current company itself then add
            if comp_name != pcompg_name:
                ismanager = pcompg_man
                share = check_share_is_number(share)
                loc =''

                all_conns["company"].append(comp_name)
                all_conns["connection"].append(pcompg_name)
                all_conns["attribute"].append("KMg")
                all_conns["distance"].append(None)
                all_conns["share"].append(share)
                all_conns["ismanager"].append(ismanager)
                all_conns["person"].append(person)
                all_conns["function"].append("mother company")
                all_conns["location"].append(loc)

        ## Get all shareholder, also the ones with minor shares, and add them to the list
        shareholders = row.ATE_name
        sh_shares = row.ATE_share
        sh_funcs = row.ATE_func
        sh_mans = row.ATE_man
        sh_persons = row.ATE_type

        if type(shareholders) == str:
            shareholders = shareholders.split("\n")
            num_shs = len(shareholders)
            if num_shs > 1:
                sh_shares = sh_shares.split("\n")
                sh_funcs = sh_funcs.split("\n")
                sh_mans = sh_mans.split("\n")
                sh_persons = sh_persons.split("\n")
            else:
                sh_shares = [sh_shares]
                sh_funcs = [sh_funcs]
                sh_mans = [sh_mans]
                sh_persons = [sh_persons]

            for i, shs in enumerate(shareholders):
                ismanager = sh_mans[i]
                share = sh_shares[i]
                share = check_share_is_number(share)
                func = sh_funcs[i]
                person = sh_persons[i]
                loc = ''

                all_conns["company"].append(comp_name)
                all_conns["connection"].append(shs)
                all_conns["attribute"].append("ATE")
                all_conns["distance"].append(1)
                all_conns["share"].append(share)
                all_conns["ismanager"].append(ismanager)
                all_conns["person"].append(person)
                all_conns["function"].append(func)
                all_conns["location"].append(loc)

        ## Get all subsidaries and add them to the list
        subsidiaries = row.TGS_name
        subs_shares = row.TGS_share

        if type(subsidiaries) == str:
            subsidiaries = subsidiaries.split("\n")
            num_subs = len(subsidiaries)
            if num_subs > 1:
                subs_shares = subs_shares.split("\n")
            else:
                subs_shares = [subs_shares]

            for i, subs in enumerate(subsidiaries):
                ismanager = 'unkown'
                share = subs_shares[i]
                share = check_share_is_number(share)
                func = "subsidary"
                person = "Unternehmen"
                loc = ''

                all_conns["company"].append(comp_name)
                all_conns["connection"].append(subs)
                all_conns["attribute"].append("TGS")
                all_conns["distance"].append(-1)
                all_conns["share"].append(share)
                all_conns["ismanager"].append(ismanager)
                all_conns["person"].append(person)
                all_conns["function"].append(func)
                all_conns["location"].append(loc)

        ## Get all subsidaries from second columns and add them to the list, if they are not already in there
        subsidiaries = row.TGS2_name
        subs_shares = row.TGS2_dshare

        if type(subsidiaries) == str:
            subsidiaries = subsidiaries.split("\n")
            num_subs = len(subsidiaries)
            if num_subs > 1:
                subs_shares = subs_shares.split("\n")
            else:
                subs_shares = [subs_shares]

            for i, subs in enumerate(subsidiaries):
                ismanager = 'unkown'
                share = subs_shares[i]
                share = check_share_is_number(share)
                func = "subsidary"
                person = "Unternehmen"
                loc = ''


                all_conns["company"].append(comp_name)
                all_conns["connection"].append(subs)
                all_conns["attribute"].append("TGS2")
                all_conns["distance"].append(-1)
                all_conns["share"].append(share)
                all_conns["ismanager"].append(ismanager)
                all_conns["person"].append(person)
                all_conns["function"].append(func)
                all_conns["location"].append(loc)

        ## Get all managers  and add them to the list
        managers = row.MAN_name
        mans_funcs = row.MAN_func
        mans_loc = row.MAN_loc

        managers2 = row.DM_name

        if type(managers) == str:
            managers = managers.split("\n")
            managers2 = managers2.split("\n")
            managers2_clean = [item.replace("hr. ", "") for item in managers2]
            managers2_clean = [item.replace("frau ", "") for item in managers2_clean]
            managers2_clean = [item.replace("dr. ", "") for item in managers2_clean]
            managers2_clean = [item.replace("dipl.-kfm. ", "") for item in managers2_clean]
            managers2_clean = [item.replace("dipl.-volkswirt ", "") for item in managers2_clean]
            managers2_clean = [item.replace("naturwissenschaften ", "") for item in managers2_clean]
            num_mans = len(managers)
            if num_mans > 1:
                # mans_bdates = mans_bdates.split("\n")
                mans_funcs = mans_funcs.split("\n")
                mans_loc = mans_loc.split("\n")
            else:
                # mans_bdates = [mans_bdates]
                mans_funcs = [mans_funcs]
                mans_loc = [mans_loc]

            for i, manager in enumerate(managers):
                ismanager = mans_funcs[i]
                share = 0.0
                func = mans_funcs[i]

                if manager in managers2_clean:
                    name_ind = managers2_clean.index(manager)
                    name2 = managers2[name_ind]

                    if ("hr." in name2) or \
                            ("frau" in name2) or \
                            ("dr." in name2) or \
                            ("dipl.-kfm." in name2) or \
                            ("dipl.-volkswirt" in name2) or \
                            ("naturwissenschaften" in name2):
                        person = 'Person'
                    else:
                        person = "Unternehmen"
                else:
                    if ("mbh" in name2) or \
                            ("gesellschaft" in name2) :
                        person = 'Unternehmen'
                    else:
                        person = "Person"

                # bdate = mans_bdates[i]
                try:
                    loc = mans_loc[i]
                except:
                    loc = ''

                all_conns["company"].append(comp_name)
                all_conns["connection"].append(manager)
                all_conns["attribute"].append("MAN")
                all_conns["distance"].append(1)
                all_conns["share"].append(share)
                all_conns["ismanager"].append(ismanager)
                all_conns["person"].append(person)
                all_conns["function"].append(func)
                all_conns["location"].append(loc)


    df_out = pd.DataFrame(all_conns)

    ## Unify ismanager column
    print("\t", list(df_out["ismanager"].unique()))
    uni_dict = {"Not a manager": "no",
                "nein": "no",
                "unkown": "unkown",
                "ja": "yes",
                "Current manager": "yes",
                "Previous manager": "no",
                "Geschäftsführer": "yes",
                "Prokurist": "no",
                "Vorstand": "no",
                "Aufsichtsrat": "no",
                "Komplementär": "yes",
                "Previous manager": "no",
                "Verwaltungsrat": "no",
                "Liquidator": "no",
                "persönlich haftender Gesellschafter": "no",
                "Geschäftsführender Direktor": "yes"}
    df_out["ismanager"] = df_out["ismanager"].map(uni_dict)

    row_lst = []
    conn_lst = []
    comp_lst = df_out["company"].unique()

    for company in comp_lst:
        df_sub = df_out.loc[df_out["company"] == company].copy()

        ##
        ## Clean person names
        ##

        df_comp = df_sub.loc[df_sub["person"] == 'Unternehmen'].copy()
        df_pers = df_sub.loc[df_sub["person"] == 'Person'].copy()

        df_pers["surname"] = df_pers["connection"].apply(lambda row: row.split(" ")[0])
        df_pers["familyname"] = df_pers["connection"].apply(lambda row: row.split(" ")[-1])

        ## Unify all names that where a person occurs twice - once with two family names (Müller-Schulze) and
        ## another time with one (Müller)
        df_pers["twofname"] = 0
        df_pers.loc[df_pers["familyname"].str.contains('-'), "twofname"] = 1
        df_pers["familyname2"] = df_pers["familyname"].apply(lambda row: row.split("-")[0])
        df_pers["name2"] = df_pers["surname"] + ' ' + df_pers["familyname2"]

        name2s = df_pers.loc[(df_pers["twofname"] == 1), "name2"].unique()
        name1s = df_pers["connection"].unique()
        for name in name2s:
            if name in name1s:
                df_pers.loc[df_pers["name2"] == name, "connection"] = name

        ## Unify all names where a person occurs twice - 1) surname familyname 2) familyname surname
        df_pers["name_rev"] = df_pers["familyname2"] + ' ' + df_pers["surname"]
        names_rev = df_pers["name_rev"].unique()
        names_norm = df_pers["connection"].unique()
        for name in names_rev:
            if name in names_norm:
                print("\tThis name occurs as surname-familyname and familyname-surname:", name)

        name_correction = {
            "koch ralf": "ralf koch",
            "gordon michael": "michael gordon",
            "schmidt andreas": "andreas schmidt",
            "schwarz oliver": "oliver schwarz",
            "schielicke thomas": "thomas schielicke",
            "dussmann peter": "peter dussmann",
            "gretschuskin hans": "hans gretschuskin",
            "schmidt walter": "walter schmidt",
            "braun michael": "michael braun",
            "ryan john": "john ryan",
            "john peter": "peter john",
            "lars anders": "anders lars",
            "gorodiscer mark": "mark gorodiscer",
            "biermann detlev": "detlev biermann",
            "hammerl robert": "robert hammerl",
            "huesmann robin": "robin huesmann",
            "sassenhagen wolfgang": "wolfgang sassenhagen",
            "thomas marco":"marco thomas",
            "prettl rolf": "rolf prettl",
            "bakiasi muharrem": "muharrem bakiasi",
            "kellner albrecht": "albrecht kellner",
            "witzel otto-friedrich": "otto-friedrich witzel",
            "gondert stephan": "stephan gondert",
            "gerlinger tobias": "tobias gerlinger",
            "lenz ulrich": "ulrich lenz",
            "hegyesi norbert": "norbert hegyesi"
        }

        uni_names = df_pers["connection"].unique()
        for name in uni_names:
            if name not in name_correction:
                name_correction[name] = name

        df_pers["connection"] = df_pers["connection"].map(name_correction)
        df_pers = df_pers[['company', 'connection', 'attribute', 'distance', 'share', 'ismanager', 'person', 'function', 'location']]
        df_sub = pd.concat([df_comp, df_pers])
        df_sub["connection"] = df_sub["connection"].str.lstrip(' ')
        df_sub["connection"] = df_sub["connection"].str.rstrip(' ')

        ##
        ## Clean all connections
        ##

        attrs = list(df_sub["attribute"].unique())

        df_comp = pd.DataFrame(data=None, columns=df_sub.columns)
        mother_comp = ''

        if 'CSH' in attrs:
            df_comp = df_sub.loc[df_sub["attribute"] == 'CSH'].copy()
            df_lst = df_comp.to_numpy().tolist()
            row_lst += df_lst
            mother_comp = df_comp.loc[df_comp["distance"] == df_comp["distance"].max(), "connection"].iloc[0]
            for x in df_comp.itertuples():
                conn_name = x.connection
                conn_lst.append((company, conn_name))

        elif "KMn" in attrs:
            df_comp = df_sub.loc[df_sub["attribute"] == 'KMn'].copy()
            mother_comp = df_comp["connection"].iloc[0]
            if mother_comp != company:
                df_lst = df_comp.to_numpy().tolist()
                row_lst += df_lst
                conn_lst.append((company, mother_comp))

        if "KMg" in attrs:
            df_kmg = df_sub.loc[df_sub["attribute"] == 'KMg'].copy()
            mother_compg = df_kmg["connection"].iloc[0]
            if mother_compg != mother_comp:
                if mother_comp != company:
                    df_lst = df_comp.to_numpy().tolist()
                    row_lst += df_lst
                    conn_lst.append((company, mother_comp))


        ## if shareholder not already listed, then append
        if "ATE" in attrs:
            df_ate = df_sub.loc[df_sub["attribute"] == 'ATE'].copy()
            for row in df_ate.itertuples():
                name = row.connection
                share = row.share
                row = list(row)[1:]
                ## Append if not listed
                ## Else check if information can be updated
                if (company, name) not in conn_lst:
                    row_lst.append(row)
                    conn_lst.append((company, name))
                else:
                    index = conn_lst.index((company, name))
                    share2 = row_lst[index][4]
                    share2 = check_share_is_number(share2)
                    if share > share2:
                        row_lst[index][4] = share
                    ismanager = row_lst[index][5]
                    if ismanager == "yes":
                        row_lst[index][5] = "yes"

        ## if subsidary not already listed, then append
        if "TGS" in attrs:
            df_tgs = df_sub.loc[df_sub["attribute"] == 'TGS'].copy()
            for row in df_tgs.itertuples():
                name = row.connection
                share = row.share
                ismanager = row.ismanager
                row = list(row)[1:]
                ## Append if not listed
                ## Else check if information can be updated
                if (company, name) not in conn_lst:
                    row_lst.append(row)
                    conn_lst.append((company, name))
                else:
                    index = conn_lst.index((company, name))
                    share2 = row_lst[index][4]
                    share2 = check_share_is_number(share2)
                    if share > share2:
                        row_lst[index][4] = share
                    if ismanager == "yes":
                        row_lst[index][5] = "yes"

        ## if subsidary (from second source) not already listed, then append
        if "TGS2" in attrs:
            df_tgs = df_sub.loc[df_sub["attribute"] == 'TGS2'].copy()
            for row in df_tgs.itertuples():
                name = row.connection
                share = row.share
                ismanager = row.ismanager
                row = list(row)[1:]
                ## Append if not listed
                ## Else check if information can be updated
                if (company, name) not in conn_lst:
                    row_lst.append(row)
                    conn_lst.append((company, name))
                else:
                    index = conn_lst.index((company, name))
                    share2 = row_lst[index][4]
                    share2 = check_share_is_number(share2)
                    if share > share2:
                        row_lst[index][4] = share
                    if ismanager == "yes":
                        row_lst[index][5] = "yes"

        ## if manager not already listed, then append
        if "MAN" in attrs:
            df_man = df_sub.loc[df_sub["attribute"] == 'MAN'].copy()
            for row in df_man.itertuples():
                name = row.connection
                share = row.share
                ismanager = row.ismanager
                location = row.location
                row = list(row)[1:]
                ## Append if not listed
                ## Else check if information can be updated
                if (company, name) not in conn_lst:
                    row_lst.append(row)
                    conn_lst.append((company, name))
                else:
                    index = conn_lst.index((company, name))
                    share2 = row_lst[index][4]
                    share2 = check_share_is_number(share2)
                    if share > share2:
                        row_lst[index][4] = share
                    if ismanager == "yes":
                        row_lst[index][5] = "yes"
                    if location != '':
                        row_lst[index][8] = location

    df_out_clean = pd.DataFrame(row_lst)
    df_out_clean.columns = df_sub.columns
    df_out_clean["share"] = df_out_clean["share"].apply(check_share_is_number)
    df_out_clean.to_csv(network_data_long_df_pth, sep=',', index=False)


def identify_possible_alkis_matches_for_manual_assignment(dafne_search_results, network_data_long_df_pth, alkis_pth,
                                                          possible_matches_folder):
    """
    Identifies possible entries in DAFNE that are worth a manual look (i.e. an assignment of the locations to the
    owner names to match the DAFNE entries with the ALKIS entries).
    Args:
        dafne_search_results: Path to dataframe with DAFNE search results.
        network_data_long_df_pth: Path to long DAFNE dataframe.
        alkis_pth: Path to dataframe with owners.
        possible_matches_folder: Output folder to store the json files for manual assignment in.

    Returns: Returns different dictionaries in possible_matches_folder. Each dictionary holds ALKIS land
    owners in different size classes. This enables to decide which land owners might be worth to look into.

    """

    print("Identify possible DAFNE-ALKIS matches that need some manual correction.")

    ## Open and prepare dafne data
    print("\tRead data.")
    df_dafne = pd.read_excel(dafne_search_results)
    df = pd.read_csv(network_data_long_df_pth)
    df_alkis = pd.read_csv(alkis_pth, sep=';')

    print("\tPrepare data")
    col_names = {'Unnamed: 0': "index",
                 'Name des Unternehmens': "comp_name",
                 'Ort / Gemeinde': "comp_munic",
                 'Vollständiger Name': "comp_name_full",
                 'Postleitzahl': "post_code",
                 'Ort': "comp_loc",
                 'Handelsregisternummer': "hrn",
                 'Tätigkeitsbeschreibung': "occup"}

    df_dafne.rename(columns=col_names, inplace=True)
    cols = ["comp_name", "comp_munic", "comp_name_full", "post_code", "comp_loc"]

    ## Replace Umlaute & All names to lower case letters only
    for col in cols:
        df_dafne.loc[df_dafne[col].isna(), col] = 'unbekannt'
        df_dafne[col] = df_dafne[col].astype(str)
        df_dafne[col] = df_dafne[col].str.lower()
        df_dafne[col] = df_dafne[col].str.replace("ä", "ae")
        df_dafne[col] = df_dafne[col].str.replace("ö", "oe")
        df_dafne[col] = df_dafne[col].str.replace("ü", "ue")
        df_dafne[col] = df_dafne[col].str.replace("ß", "ss")
        df_dafne[col] = df_dafne[col].str.replace("  ", " ")

    ## Open Network data, subset to relevant people, derive certain characteristics
    sub = df.loc[df['person'] == 'Person'].copy()
    sub = sub.loc[(sub["ismanager"] == 'yes') | (sub["share"] > 0)].copy()
    sub["family_name"] = sub.apply(lambda row: row.connection.split(" ")[-1], axis=1)
    sub["surname"] = sub.apply(lambda row: row.connection.split(" ")[0], axis=1)
    sub = sub.loc[sub["location"].isna()]
    sub.drop_duplicates(subset=["connection"], inplace=True)
    unique_lst = list(sub["family_name"].unique())

    ## Subset ALKIS to private people, calculate area per owner and address, derive important characteristics
    t = df_alkis.loc[df_alkis["owner_clean"].str.count("matthes") > 0].copy()
    t2 = t.loc[t["owner_names"].str.count("rainer") > 0].copy()
    t3 = df_alkis.loc[df_alkis["owner_names"].str.count("rainer") > 0].copy()
    t4 = t.loc[t["owner_names"].str.count("reiner") > 0].copy()
    sub_alkis = df_alkis.loc[df_alkis["owner_clean"].str.count("\*") > 0].copy()
    sub_alkis = sub_alkis.loc[sub_alkis["category"] == 1].copy()
    sub_alkis = sub_alkis[["owner_clean", "clean_address", "area"]].groupby(
        ["owner_clean", "clean_address"]).sum().reset_index()
    sub_alkis["area"] = sub_alkis["area"] / 10000
    # sub = sub.loc[sub['area'] > 50].copy()
    sub_alkis["birthdate"] = sub_alkis.apply(lambda row: row.owner_clean.split(',')[-1], axis=1)
    sub_alkis["name"] = sub_alkis.apply(lambda row: row.owner_clean.split('*')[0], axis=1)
    sub_alkis["comma_count"] = sub_alkis["name"].str.count(',')
    sub_alkis = sub_alkis[sub_alkis["comma_count"] > 0].copy()
    sub_alkis["family_name"] = sub_alkis.apply(lambda row: row.owner_clean.split(',')[0], axis=1)
    # t = sub_alkis.loc[sub_alkis["owner_clean"].str.count("olearius") > 0].copy()
    sub_alkis["surname"] = sub_alkis.apply(lambda row: row.owner_clean.split(',')[1], axis=1)

    ## Subset ALKIS with family name list from DAFNE data
    sub_alkis = sub_alkis.loc[sub_alkis["family_name"].isin(unique_lst)]
    sub_alkis.drop_duplicates(subset=["owner_clean"], inplace=True)
    print(len(sub_alkis))

    ## 1. calculate area per owner
    ## 2. attach family_name and name and address to owners with areas
    ## 3. subset to owners with more than x ha

    out_dict = {}

    print("\tLoop over private persons in DAFNE network data.")
    for row in sub.itertuples():
        curr_person = row.connection
        company = row.company
        family_name = row.family_name
        surname = row.surname
        if '-' in surname:
            surname = surname.split('-')[0]

        if family_name == 'matthes':
            print("\t")

        df_curr = sub_alkis.loc[sub_alkis["family_name"].isin([family_name])].copy()
        df_curr = df_curr.loc[df_curr["surname"].str.count(surname) > 0].copy()
        # df_curr = df_curr.loc[df_curr["area"] > 100].copy()

        if df_curr.empty:
            continue

        max_area = df_curr['area'].max()
        comp_munic = df_dafne.loc[df_dafne["comp_name"] == company, "comp_munic"].iloc[0]
        post_code = df_dafne.loc[df_dafne["comp_name"] == company, "post_code"].iloc[0]
        comp_loc = df_dafne.loc[df_dafne["comp_name"] == company, "comp_loc"].iloc[0]

        sub_dict = {}
        sub_dict["company"] = company
        sub_dict["company_postal_coda"] = post_code
        sub_dict["company_city"] = comp_loc
        sub_dict["company_gemeinde"] = comp_munic
        sub_dict["max_area"] = max_area
        sub_dict["new_name"] = curr_person

        i = 1
        for row2 in df_curr.itertuples():
            sub_dict[f"match{i}"] = f"{row2.owner_clean}_{row2.clean_address}"
            i += 1

        out_dict[curr_person] = sub_dict

    print("\tTotal possible matches:", len(out_dict))

    sub_dict = {key: out_dict[key] for key in out_dict if out_dict[key]["max_area"] > 100}
    with open(fr"{possible_matches_folder}\network_private_persons_matches_in_alkis_100-max.json", 'w') as fp:
        json.dump(sub_dict, fp, indent=4)
    print("\tTotal matches with area > 100 ha:", len(sub_dict))

    sub_dict = {key: out_dict[key] for key in out_dict if
                (out_dict[key]["max_area"] > 50) and (out_dict[key]["max_area"] <= 100)}
    with open(fr"{possible_matches_folder}\network_private_persons_matches_in_alkis_50-100.json", 'w') as fp:
        json.dump(sub_dict, fp, indent=4)
    print("\tTotal matches with area > 50 ha:", len(sub_dict))

    sub_dict = {key: out_dict[key] for key in out_dict if
                out_dict[key]["max_area"] > 25 and (out_dict[key]["max_area"] <= 50)}
    with open(fr"{possible_matches_folder}\network_private_persons_matches_in_alkis_25-50.json", 'w') as fp:
        json.dump(sub_dict, fp, indent=4)
    print("\tTotal matches with area > 25 ha:", len(sub_dict))

    sub_dict = {key: out_dict[key] for key in out_dict if
                out_dict[key]["max_area"] > 10 and (out_dict[key]["max_area"] <= 25)}
    with open(fr"{possible_matches_folder}\network_private_persons_matches_in_alkis_10-15.json", 'w') as fp:
        json.dump(sub_dict, fp, indent=4)
    print("\tTotal matches with area > 10 ha:", len(sub_dict))

    print("\tWrite out.")
    with open(fr"{possible_matches_folder}\network_private_persons_matches_in_alkis.json", 'w') as fp:
        json.dump(out_dict, fp, indent=4)


def clean_long_network_table(network_data_long_df_pth, prepared_network_data_pth,
                             corrected_dafne_alkis_matches_json=None):
    """
    Reformats the long dataframe into a more clean version.
    Args:
        network_data_long_df_pth: Path to long DAFNE dataframe
        corrected_dafne_alkis_matches_json: A dictionary with old DAFNE names and the corresponding version of
        the ALKIS names + locations in form "surname familynam_location"
        prepared_network_data_pth: Output path to clean network dataframe.

    Returns:

    """
    ## ToDo location is added here to all the connections
    print("Clean long dataframe version.")
    print("\tRead data.")
    df = pd.read_csv(network_data_long_df_pth)

    print("\tReformat")
    df.loc[df["location"].isna(), "location"] = ""
    cols = ["location"]
    for col in cols:
        df[col] = df[col].str.lower()
        df[col] = df[col].str.replace("ä", "ae")
        df[col] = df[col].str.replace("ö", "oe")
        df[col] = df[col].str.replace("ü", "ue")
        df[col] = df[col].str.replace("ß", "ss")
        df[col] = df[col].str.lstrip(' ')
        df[col] = df[col].str.rstrip(' ')
        df.loc[df[col] != '', col] = df.loc[df[col] != '', col].apply(helper_functions.clean_city_text, search_terms=['bei', 'ot'])
        # df[col] = df[col].apply(clean_city_text, search_terms=['bei', 'ot'])
        df[col] = df[col].str.replace("/", "")
        df[col] = df[col].str.replace("Aarhus C", "Aarhus")

    ## Reformat to old format for further analysis.
    comp_lst = df["company"].unique()
    all_conns = {"conn_left": [], "conn_right": [], "share": [], "ismanager": [], "function": [], "attribute": [],
                 "isperson": []}

    for company in comp_lst:
        df_sub = df.loc[df["company"] == company].copy()
        attributes = df_sub["attribute"].tolist()
        if 'CSH' in attributes:
            df_csh = df_sub.loc[df_sub["attribute"] == 'CSH'].copy()
            max_dist = df_csh["distance"].max()
            for dist in reversed(range(int(max_dist))):
                if dist != 0:
                    conn_left = df_csh.loc[df_csh["distance"] == dist, "connection"].iloc[0]
                else:
                    conn_left = df_csh.loc[df_csh["distance"] == dist + 1, "company"].iloc[0]
                conn_right = df_csh.loc[df_csh["distance"] == dist + 1, "connection"].iloc[0]
                share = df_csh.loc[df_csh["distance"] == dist + 1, "share"].iloc[0]
                ismanager = df_csh.loc[df_csh["distance"] == dist + 1, "ismanager"].iloc[0]
                function = df_csh.loc[df_csh["distance"] == dist + 1, "function"].iloc[0]
                attr = df_csh.loc[df_csh["distance"] == dist + 1, "attribute"].iloc[0]
                location = df_csh.loc[df_csh["distance"] == dist + 1, "location"].iloc[0]
                person = df_csh.loc[df_csh["distance"] == dist + 1, "person"].iloc[0]

                if location != '':
                    conn_right = f"{conn_right}_{location}"

                all_conns["conn_left"].append(conn_left)
                all_conns["conn_right"].append(conn_right)
                all_conns["share"].append(share)
                all_conns["ismanager"].append(ismanager)
                all_conns["function"].append(function)
                all_conns["attribute"].append(attr)
                all_conns["isperson"].append(person)

        df_sh = df_sub.loc[(df_sub["attribute"] != "CSH") & (df_sub["distance"] > 0)].copy()
        for row in df_sh.itertuples():
            conn_left = row.company
            conn_right = row.connection
            share = row.share
            ismanager = row.ismanager
            function = row.function
            attr = row.attribute
            location = row.location
            person = row.person

            if location != '':
                conn_right = f"{conn_right}_{location}"

            all_conns["conn_left"].append(conn_left)
            all_conns["conn_right"].append(conn_right)
            all_conns["share"].append(share)
            all_conns["ismanager"].append(ismanager)
            all_conns["function"].append(function)
            all_conns["attribute"].append(attr)
            all_conns["isperson"].append(person)

        df_subs = df_sub.loc[(df_sub["attribute"] != "CSH") & (df_sub["distance"] < 0)].copy()
        for row in df_subs.itertuples():
            conn_left = row.connection
            conn_right = row.company
            share = row.share
            ismanager = row.ismanager
            function = row.function
            attr = row.attribute
            location = row.location
            person = row.person

            if location != '':
                conn_right = f"{conn_right}_{location}"

            all_conns["conn_left"].append(conn_left)
            all_conns["conn_right"].append(conn_right)
            all_conns["share"].append(share)
            all_conns["ismanager"].append(ismanager)
            all_conns["function"].append(function)
            all_conns["attribute"].append(attr)
            all_conns["isperson"].append(person)

    df_netw = pd.DataFrame(all_conns)

    ## correct names of important land owners that occur in DAFNE, but without an address
    if corrected_dafne_alkis_matches_json:
        print("\t Correct old DAFNE names with the corresponding version of the ALKIS names_locations")
        file = open(corrected_dafne_alkis_matches_json)
        corr_dict = json.load(file)
        clean_dict = {}

        for key in corr_dict:
            clean_dict[key] = corr_dict[key]["new_name"]

        df_netw["conn_right"] = df_netw["conn_right"].replace(clean_dict)

    df_netw.to_csv(prepared_network_data_pth, index=False)


def network_connections_to_list_for_dafne_search(prepared_network_data_pth, dafne_search_results, col_renaming_json_pth,
                                                 dafne_search_list_second_round, folder_second_search):
    """
    Create a list of company names that are part of the company networks that were identified.
    The goal is to use the classification of these companies to characterize the company networks.
    Args:
        prepared_network_data_pth: Path to clean network dataframe.
        dafne_search_results: Path to dataframe with DAFNE search results.
        col_renaming_json_pth: Path to json with dictionary of columns renaming.
        dafne_search_list_second_round: Output path to data frame with company names for second search round.
        folder_second_search: Output folder to lists with company names for second search.

    Returns:

    """

    print("Prepare second round of DAFNE search.")

    ## Get list of all companies
    print("\tRead data.")
    df = pd.read_csv(prepared_network_data_pth)
    df_old = pd.read_excel(dafne_search_results)

    print("\tPrepare.")
    df_left = df[["conn_left"]].copy()
    df_right = pd.DataFrame(df.loc[df["isperson"] == "Unternehmen", "conn_right"])
    df_right.rename(columns={"conn_right": "conn_left"}, inplace=True)

    df_out = pd.concat([df_left, df_right], axis=0)
    df_out.drop_duplicates(subset="conn_left", inplace=True)

    ## Drop companies that were already searched
    with open(col_renaming_json_pth) as json_file:
        col_names = json.load(json_file)
    df_old.rename(columns=col_names, inplace=True)
    cols = ["comp_name", "CSH_name", "KMn_name", "KMg_name", "TGS_name", "TGS2_name", "MAN_name", "ATE_name", "DM_name"]

    ## Replace Umlaute & All names to lower case letters only
    for col in cols:
        df_old[col] = df_old[col].str.lower()
        df_old[col] = df_old[col].str.replace("ä", "ae")
        df_old[col] = df_old[col].str.replace("ö", "oe")
        df_old[col] = df_old[col].str.replace("ü", "ue")
        df_old[col] = df_old[col].str.replace("ß", "ss")
        df_old[col] = df_old[col].str.replace("  ", " ")

    done_lst = df_old["comp_name"].tolist()

    print(len(done_lst), len(df_out))
    df_done = df_out.loc[df_out["conn_left"].isin(done_lst)].copy()
    df_out = df_out.loc[~df_out["conn_left"].isin(done_lst)].copy()
    df_out.index = range(len(df_out))
    print(len(df_out) + len(done_lst))

    df_out.to_csv(dafne_search_list_second_round)

    # x = 0
    # num_rows = len(df_out)
    # num_lists = math.ceil(num_rows / 1000)
    # s = 0
    # for i in range(1, num_lists + 1):
    #     e = i * 1000
    #     df_sub = df_out[s:e]
    #     s += 1000
    #     out_pth = rf"{OUT_FOLDER}\batch_search_{i + x:02d}.txt"
    #     df_sub[["conn_left"]].to_csv(out_pth, sep=';', header=None, index=None)

    helper_functions.create_folder(folder_second_search)
    x = 14
    num_rows = len(df_done)
    num_lists = math.ceil(num_rows / 1000)
    s = 0
    for i in range(1, num_lists + 1):
        e = i * 1000
        df_sub = df_done[s:e]
        s += 1000
        out_pth = rf"{folder_second_search}\batch_search_{i + x:02d}.txt"
        df_sub[["conn_left"]].to_csv(out_pth, sep=';', header=None, index=None)

    print("Done.")


def add_economic_branch_to_network_companies(folder_second_search, df_companies_with_branches_and_locations_pth):
    """
    Adds the economic branch and the location of a company as provided in DAFNE to the company names.
    Args:
        folder_second_search: Folder with search results of second round in DAFNE.
        df_companies_with_branches_and_locations_pth:  Output path to data frame with company names, branches and locations.

    Returns:

    """

    print("Add economic branch and location to company names.")
    ## read table with already searched companies
    # df_dafne_res = pd.read_csv(DONE_COMP_WITH_BRANCHES) # ["search_name", "bvd_name", "bvd_id", 'main_branch', 'side_branch', "agric"]
    # df_dafne_res.rename(columns={"search_name": "company_name", "bvd_name": "company_name_bvd"}, inplace=True)

    ## read dafne search results
    search_lst = glob.glob(fr"{folder_second_search}\batch_search_*.txt")
    result_lst = glob.glob(fr"{folder_second_search}\batch_search_result_*.xlsx")
    meta_lst = glob.glob(fr"{folder_second_search}\batch_search_*_result_metadata.xlsx")

    search_names = [pd.read_csv(txt, sep=';', header=None) for txt in search_lst]
    result_dfs = [pd.read_excel(pth, sheet_name="Ergebnisse", dtype={"Postleitzahl": str}) for pth in result_lst]
    meta_dfs = [pd.read_excel(pth, sheet_name="Page 1") for pth in meta_lst]

    search_df = pd.concat(search_names, axis=0)
    result_df = pd.concat(result_dfs, axis=0)
    meta_df = pd.concat(meta_dfs, axis=0)

    ## Clean tables
    search_df.columns = ["company_name"]
    result_df.columns = ["index", "company_name_bvd", "main_side_branch", "main_side_branch_descr", "location", "postcode", "fstate",
                          "country", "main_branch", "main_branch_descr", "side_branch", "side_branch_descr", "bvd_id", "self_descr"]
    meta_df.columns = ["company_name", "national_id", "city_empty", "country_empty", "score", "bvd_id", "company_name_bvd"]


    ## Merge results and search names
    met_cols = ["company_name", "bvd_id"]
    res_cols = ["company_name_bvd", "bvd_id", "location", "fstate", "country", "postcode", "main_branch", "side_branch", "self_descr"]
    merge_df = pd.merge(meta_df[met_cols], result_df[res_cols], how="left", on="bvd_id")

    for col in res_cols:
        merge_df.loc[merge_df[col].isna(), col] = 'unkown'

    merge_df = pd.merge(merge_df, search_df, how="left", on="company_name")

    ## clean classification
    merge_df["main_branch"] = merge_df["main_branch"].apply(lambda x: ''.join([i for i in x if not i.isdigit()]))
    merge_df["side_branch"] = merge_df["side_branch"].apply(lambda x: ''.join([i for i in x if not i.isdigit()]))
    merge_df["agric"] = 0
    merge_df.loc[(merge_df["main_branch"].str.contains("A")) | (merge_df["side_branch"].str.contains("A")), "agric"] = 1

    # for col in ["main_branch", "side_branch"]:
    #     df_dafne_res.loc[df_dafne_res[col].isna(), col] = 'unkown'
    #
    # df_dafne_res["main_branch"] = df_dafne_res["main_branch"].apply(lambda x: ''.join([i for i in x if not i.isdigit()]))
    # df_dafne_res["side_branch"] = df_dafne_res["side_branch"].apply(lambda x: ''.join([i for i in x if not i.isdigit()]))
    #
    # ## create combined table only with branches and agric
    # cols = ["company_name", "bvd_id", "main_branch", "side_branch", "agric"]
    # df_branches = pd.concat([merge_df[cols], df_dafne_res[cols]], axis=0)

    ## create table with location info for newly researched companies
    # df_loc_new = merge_df[["company_name", "bvd_id", ]]
    df_branches = merge_df[["company_name", "bvd_id", "main_branch", "side_branch", "agric",
                            "location", "fstate", "postcode", "self_descr"]].copy()

    df_branches["international"] = 0
    international_code_words = ["sp.", "s.p.a.", "ltd.", "limited", "s.a.", " se ", "sp.z o.o.", "b.v.", "bhd.", "srl",
                                "ltd", "s.à r.l.", " sa ", " bv ", " sl ", "inc.", "s.r.l.", " inc "]

    df_branches.loc[df_branches["main_branch"] == "unkown", "international"] = \
        df_branches.loc[df_branches["main_branch"] == "unkown", "company_name"].apply(
        helper_functions.check_occ_of_words_v2, search_terms=international_code_words, return_code=1)

    df_branches.to_csv(df_companies_with_branches_and_locations_pth)

    print("Done!")


def identify_communities(df, thresh=None):

    """

    :param df: df or group of rows from df (from group_by) with left and right connections of a company network and
    the columns "share" for the shares that the right connection holds on the left connections and "ismanager"
    indicating wether right connection is the manager of left connection
    :param thresh: Treshhold which shares to include (all greater and equal).
    :return: df with connections.
    """

    if thresh:
        ## exclude all connections that are below share but include managers
        df1 = df.loc[df['share'] >= thresh].copy()
        df2 = df.loc[(df["share"] < thresh) & (df["ismanager"] == 'yes')].copy()
        df3 = df.loc[(df["share"] < thresh) & (df["function"] == 'Vorstand')].copy()
        # df2 = df2.loc[(df2["share"] > 0)].copy()
        df = pd.concat([df1, df2, df3], axis=0)
        df.sort_index(inplace=True)

    node_names = list(set(df["conn_left"].unique().tolist() + df["conn_right"].unique().tolist()))

    ## create edges from df
    edges = [(row.conn_left, row.conn_right) for row in df.itertuples()]

    def get_communities(node_names, edges):

        """
        This function separates sub communities from several communites
        :param node_names:
        :param edges:
        :return:
        """

        community_dict = {}

        sub_comm = 0
        for node_name in node_names:
            if node_name in community_dict:
                continue

            ## put in list that will dynamically be appended
            nnames = [node_name]

            for nname in nnames:
                ## look for all edges in which the current names occurs
                edges_lst = [edge for edge in edges if nname in edge]
                ## get all names that occur in the all these edges if they are nor already in nnames
                nname_lst = [nname for edge in edges_lst for nname in edge if nname not in nnames]

                ## append nnames with newly found names.
                ## This will be done until no new name is added
                nnames += nname_lst

            ## assign a number to each name
            for nname in nnames:
                if nname not in community_dict:
                    community_dict[nname] = sub_comm

            sub_comm += 1

        return community_dict

    partition = get_communities(node_names, edges)

    df[f"community_{thresh}"] = df["conn_left"].map(partition)
    df[f"community_{thresh}"] = df[f"community_{thresh}"].astype(int).astype(str)

    return df


def network_analysis_advanced(netw_pth, out_pth):
    """
    Derive the communities/networks from the table with all companies, managers etc. and their connections.
    :param netw_pth: Path to table with connections (conn_left, conn_right, shares)
    :param out_pth: Output path.
    :return:
    """
    print("Advanced network analysis.")

    ## All this was found at https://programminghistorian.org/en/lessons/exploring-and-analyzing-network-data-with-python
    os.chdir(WD)

    ## load data
    df = pd.read_csv(netw_pth)

    print("\tNumber rows of original df:", len(df))

    node_names = list(set(df["conn_left"].unique().tolist() + df["conn_right"].unique().tolist()))

    ## create edges from df
    edges = [(row.conn_left, row.conn_right) for row in df.itertuples()]

    G = nx.Graph()
    G.add_nodes_from(node_names)
    G.add_edges_from(edges)

    ## Add additional information
    share_dict = {}

    for row in df.itertuples():
        share_dict[(row.conn_left, row.conn_right)] = row.share

    nx.set_edge_attributes(G, share_dict, 'shares')

    ## Calculate communities in whole Graph
    partition = community_louvain.best_partition(G)
    nx.set_node_attributes(G, partition, 'company_cluster')

    # df["community_0"] = df["conn_left"].map(partition)
    # df[f"community_0"] = df[f"community_0"].astype(int).astype(str)

    df = identify_communities(df=df, thresh=0)

    df.drop_duplicates(subset=["conn_left", "conn_right"], inplace=True)

    print("\tUnify names of people within each community.")

    def unify_names(comm_group):

        people = comm_group.loc[comm_group["isperson"] == "Person"].copy()
        people.drop_duplicates(subset="conn_right", inplace=True)

        names = people["conn_right"].to_list()

        names_wo_loc = [name for name in names if '_' not in name]
        name_w_loc = [name for name in names if '_' in name]

        cleaning_dict = {name_wo: name_w for name_wo in names_wo_loc for name_w in name_w_loc if name_wo in name_w}

        comm_group.loc[comm_group["conn_right"].isin(cleaning_dict), "conn_right"] = comm_group.loc[
            comm_group["conn_right"].isin(cleaning_dict), "conn_right"].map(cleaning_dict)

        # if len(cleaning_dict) > 0:
        #     print(cleaning_dict)

        return comm_group

    df = df.groupby("community_0").apply(unify_names).reset_index()
    df.drop(columns=["index"], inplace=True)
    df.drop_duplicates(subset=["conn_left", "conn_right", "share"], inplace=True)

    print("\tNumber rows of df with communities, thresh=0:", len(df))

    df50 = df.groupby("community_0").apply(identify_communities, 50).reset_index(drop=True)
    print("\tNumber rows of df with communities, thresh=50:", len(df50))

    df25 = df.groupby("community_0").apply(identify_communities, 25).reset_index(drop=True)
    print("\tNumber rows of df with communities, thresh=25:", len(df25))

    df10 = df.groupby("community_0").apply(identify_communities, 10).reset_index(drop=True)
    print("\tNumber rows of df with communities, thresh=10:", len(df10))

    df = pd.merge(df, df50[["conn_left", "conn_right", "community_50"]], how="left", on=["conn_left", "conn_right"])
    df.loc[df["community_50"].isna(), "community_50"] = '99'
    df[f"community_50"] = df[["community_0", f"community_50"]].agg('_'.join, axis=1)
    df.loc[df["community_50"].str.count("_99") > 0, "community_50"] = None

    df = pd.merge(df, df25[["conn_left", "conn_right", "community_25"]], how="left", on=["conn_left", "conn_right"])
    df.loc[df["community_25"].isna(), "community_25"] = '99'
    df[f"community_25"] = df[["community_0", f"community_25"]].agg('_'.join, axis=1)
    df.loc[df["community_25"].str.count("_99") > 0, "community_25"] = None

    df = pd.merge(df, df10[["conn_left", "conn_right", "community_10"]], how="left", on=["conn_left", "conn_right"])
    df.loc[df["community_10"].isna(), "community_10"] = '99'
    df[f"community_10"] = df[["community_0", f"community_10"]].agg('_'.join, axis=1)
    df.loc[df["community_10"].str.count("_99") > 0, "community_10"] = None

    print("\tNumber rows of df for output:", len(df))

    print(f"\tWrite df to {out_pth}.")
    df.to_csv(out_pth, sep=',', index=False)

    # sub_nodes = set([key for key in partition if partition[key] == 0])
    # subgraph = G.subgraph(sub_nodes)
    #
    # get_information_from_network(subgraph)
    # time.sleep(1)
    # draw_network(subgraph)


def create_comp_community_mcomp_dictionary(net_comm_pth, comm_col, out_pth):
    """
    Create a company-community-mother company dictionary
    Args:
        net_comm_pth: Dataframe with network data and community numbers.
        comm_col: Column name with community IDs.
        out_pth: Output path for dictionary providing the mother companies and the community ids for each company.

    Returns:

    """
    print(f"Create a company-community-mother company dictionary with column {comm_col}")
    df_net = pd.read_csv(net_comm_pth, dtype={comm_col: str})
    df_net = df_net.loc[df_net[comm_col].notna()].copy()

    ## From left and right connections in network table create a table that assigns each network part to community
    df_net1 = df_net[["conn_left", comm_col]].copy()
    df_net2 = df_net[["conn_right", comm_col]].copy()
    df_net2.rename(columns={"conn_right": "conn_left"}, inplace=True)
    df_net_red = pd.concat([df_net1, df_net2], axis=0)
    df_net_red.drop_duplicates(subset=["conn_left", comm_col], inplace=True)

    ind = list(df_net_red.columns).index(comm_col) + 1
    comp_comm_dict = {row.conn_left: row[ind] for row in df_net_red.itertuples()}

    ## Find global owner per community number and save to dictionary
    comm_ids = df_net[comm_col].unique()
    comm_mcomp_dict = {}
    for comm_id in comm_ids:
        df_curr_net = df_net.loc[(df_net[comm_col] == comm_id)].copy()
        df_sub = df_curr_net.loc[(df_curr_net["function"].isin(["mother company"]))].copy()
        mother_companies = df_sub["conn_right"].unique()

        ## If a mother company in a community doesn't occur in the left connection columm (all lefts are "inferior")
        ## then it's likely the mother company of the whole community
        for mcomp in mother_companies:
            mcomp_clean = mcomp.replace(' ', '')
            conn_lefts = df_sub["conn_left"].to_list()
            conn_lefts = [item.replace(' ', '') for item in conn_lefts]
            if mcomp_clean not in conn_lefts:
                if comm_id in comm_mcomp_dict:
                    comm_mcomp_dict[comm_id].append(mcomp)
                    # print("ID:", comm_id, comm_mcomp_dict[comm_id])
                else:
                    comm_mcomp_dict[comm_id] = [mcomp]

        ## If there is no mother company, then get all companies that are not labeled as subsidaries
        ## from that list get companies/company with highest share and assign them/it as the mother company
        if len(mother_companies) == 0:
            df_sub = df_curr_net.loc[~df_curr_net["function"].isin(["subsidary"])].copy()
            df_sub = df_sub.loc[df_sub["share"] == df_sub["share"].max()].copy()
            if len(df_sub) == 1:
                mcomp = df_sub["conn_right"].iloc[0]
                comm_mcomp_dict[comm_id] = [mcomp]
            elif len(df_sub) > 1:
                mcomps = df_sub["conn_right"].tolist()
                comm_mcomp_dict[comm_id] = mcomps
            else:
                mcomps = list(df_curr_net["conn_right"].unique())
                comm_mcomp_dict[comm_id] = mcomps
                # print("ID:", comm_id, comm_mcomp_dict[comm_id])

    ## In a few cases there could not be identified a mother company
    ## In this case take the first company that occurs on the right connection side
    for comm_id in comm_ids:
        if comm_id not in comm_mcomp_dict:
            df_curr_net = df_net.loc[(df_net[comm_col] == comm_id)].copy()
            df_sub = df_curr_net.loc[(df_curr_net["function"].isin(["mother company"]))].copy()
            mother_companies = list(df_sub["conn_right"].unique())
            comm_mcomp_dict[comm_id] = mother_companies
            if len(mother_companies) > 1:
                comm_mcomp_dict[comm_id] = [mother_companies[0]]
                # print(comm_id, mother_companies)

    ## As in some cases there are multiple mother companies, they need to be transformed into a string
    for i in comm_mcomp_dict:
        comm_mcomp_dict[i] = ' AND '.join(list(set(comm_mcomp_dict[i])))

    out_dict = {"company_to_community": comp_comm_dict,
                "community_to_mother_company": comm_mcomp_dict}

    out_folder = os.path.dirname(out_pth)
    helper_functions.create_folder(out_folder)

    with open(out_pth, 'w') as fp:
        json.dump(out_dict, fp, indent=4)


def get_maximum_hierarchical_distance_in_network(net_comm_pth, comm_col, out_pth):
    """
    Derives the number of hierarchical levels per community.
    Args:
        net_comm_pth: Dataframe with network data and community numbers.
        comm_col: Column name with community IDs.
        out_pth: Output path to dictionary (json) with no. of hierarchical levels per community.

    Returns:

    """

    df_net = pd.read_csv(net_comm_pth, dtype={comm_col: str})
    df_net = df_net.loc[df_net[comm_col].notna()].copy()

    max_dist_dict = {}
    for comm in list(df_net[comm_col].unique()):
        df_sub = df_net.loc[df_net[comm_col] == comm].copy()
        df_sub.drop_duplicates(subset="conn_right", inplace=True)
        conn_rights = df_sub["conn_right"].tolist()
        conn_lefts = list(set(df_sub["conn_left"].tolist()))
        max_dist = 1
        for cl in conn_lefts:
            if cl in conn_rights:
                max_dist += 1
        max_dist_dict[comm] = max_dist

    with open(out_pth, 'w') as fp:
        json.dump(max_dist_dict, fp, indent=4)




def add_community_number_to_alkis_data(alkis_pth, comp_comm_mcomp_dict_pth, comm_col, out_pth,
                                       community_info_dict_pth, max_dist_dict_pth):
    """
    Adds the community ID to the ALKIS data, creates a dictionary with the information on subcompanies,
    number of subcompanies per community.
    Args:
        alkis_pth: Path to owner data frame.
        comp_comm_mcomp_dict_pth: Path to dictionary with company to community ID/mother company.
        comm_col: Column with community ID
        out_pth: Output path to ALKIS data frame with community information.
        community_info_dict_pth: Output path to dictionary with information on the communities.
        max_dist_dict_pth: Path to dictionary with no. of hierarchical levels per community.

    Returns:

    """
    print(f"Add community number to ALKIS data with column {comm_col}.")

    os.chdir(WD)

    ## Open data
    print("\tRead data.")
    df_alk = pd.read_csv(alkis_pth, sep=';')

    print("\tLength of table:", len(df_alk))

    with open(comp_comm_mcomp_dict_pth) as json_file:
        comp_comm_mcomp_dict = json.load(json_file)

    with open(max_dist_dict_pth) as json_file:
        max_dist_dict = json.load(json_file)

    comm_dict = comp_comm_mcomp_dict["community_to_mother_company"]
    comp_comm_dict = comp_comm_mcomp_dict["company_to_community"]
    df_comp_comm = pd.DataFrame.from_dict(comp_comm_dict, orient='index').reset_index()
    df_comp_comm.columns = ["conn_left", comm_col]

    ###############################################################
    ## Add mother company names and community numbers to ALKIS data
    ###############################################################

    ## Clean owner name in ALKIS data
    cols = ["owner_clean"]
    for col in cols:
        df_alk[col] = df_alk[col].str.lower()
        df_alk[col] = df_alk[col].str.replace("ä", "ae")
        df_alk[col] = df_alk[col].str.replace("ö", "oe")
        df_alk[col] = df_alk[col].str.replace("ü", "ue")
        df_alk[col] = df_alk[col].str.replace("ß", "ss")

    df_alk["conn_left"] = df_alk["owner_clean"]

    ## For private owner create merge name with pattern: "surname familyname_location" to merge with dafne entries
    df_alk.loc[(df_alk["level3"] == "1_1_1") & (df_alk["conn_left"].str.count(',') > 0), "conn_left"] = \
        df_alk.loc[(df_alk["level3"] == "1_1_1") & (df_alk["conn_left"].str.count(',') > 0), "conn_left"].apply(
            lambda row: row.split(',')[1] + ' ' + row.split(',')[0])

    df_alk["city"] = df_alk["clean_address"].apply(helper_functions.identify_city)
    df_alk.loc[(df_alk["level3"] == '1_1_1') &
               (df_alk["owner_clean"].str.count(',') > 0), "conn_left"] = df_alk.loc[(df_alk["level1"] == 1) &
               (df_alk["owner_clean"].str.count(',') > 0), "conn_left"] + '_' + df_alk.loc[(df_alk["level1"] == 1) &
               (df_alk["owner_clean"].str.count(',') > 0), "city"]

    df_alk["conn_left"] = df_alk["conn_left"].str.lstrip(' ')
    df_alk["conn_left"] = df_alk["conn_left"].str.rstrip(' ')

    ## Get a list of all DAFNE entries that are aristrocracy
    t = df_comp_comm.loc[
        df_comp_comm["conn_left"].str.contains(" von ") | df_comp_comm["conn_left"].str.contains(" van ")].copy()
    t = t.loc[t["conn_left"].str.contains("_")].copy()
    search_list = list(t["conn_left"].unique())

    ## create a reduced version of ALKIS for faster processing
    df_red = df_alk.drop_duplicates(subset="conn_left")
    df_red = df_red.loc[df_red["level3"] == '1_1_1'].copy()

    ## write a function that looks if the first name and the first part of the family name are in another string
    def func(x, search_lst):
        out = [item if ((item.split(" ")[0] in x) & (item.split("_")[0].split(' ')[-1] in x)) else False for item in search_lst]
        out = [item for item in out if item != False]
        if len(out) > 0:
            out = out[0]
        else:
            out = None
        return out

    ## Option 1 to check if any entry in ALKIS contains first name and first part of family name of any of the identified
    ## aristrocratic persons
    # df_red["select"] = df_red["owner_clean"].apply(lambda x: 1 if any((item.split(" ")[0] in x) & (item.split("_")[0].split(' ')[-1] in x) for item in search_list) else 0)
    # df_red2 = df_red.loc[df_red["select"] == 1].copy()

    ## Option 2 to do the same but more efficient (could even be improved)
    df_red["search_term"] = [func(x, search_list) for x in df_red["owner_clean"].values]
    df_red2 = df_red.loc[df_red["search_term"] != None].copy()
    df_red2 = df_red2[["owner_names", "owner_clean", "clean_address", "city", "conn_left", "search_term"]].copy()
    df_red2 = df_red2.loc[df_red2["search_term"].isin(search_list)].copy()
    df_red2.to_csv(r"08_network_analysis\dafne_cleaning_dict_aristrocracy.csv", sep=";", index=False)

    ##
    pth = r"08_network_analysis\dafne_cleaning_dict_aristrocracy_edited.csv"
    if not os.path.exists(pth):
        helper_functions.print_red("\tYou have to clean the identified aristrocacy entries manually and save them as dafne_cleaning_dict_aristrocracy_edited.csv!")
        exit()
    else:
        df_clean = pd.read_csv(pth, sep=";")
        clean_dict = {row.search_term: row.conn_left for row in df_clean.itertuples()}

        df_comp_comm.loc[df_comp_comm["conn_left"].isin(clean_dict.keys()), "conn_left"] = df_comp_comm.loc[
            df_comp_comm["conn_left"].isin(clean_dict.keys()), "conn_left"].map(clean_dict)


    ## Option 3, very slow
    # for owner_name in list(t["conn_left"].unique()):
    #     #     first_term = owner_name.split(" ")[0]
    #     #     second_term = owner_name.split("_")[0].split(' ')[-1]
    #     #     t2 = df_red.loc[
    #     #         df_red["owner_clean"].str.contains(first_term, regex=False) & df_red["owner_clean"].str.contains(second_term, regex=False)].copy()
    #     #     t3 = df_red.loc[df_red["owner_clean"].str.contains(f"(?=.*{second_term})(?=.*{first_term})", regex=True)]
    #     #     if not t2.empty:
    #     #         print("Test")
    #     #         print(f"{owner_name};{list(t2['conn_left'])};{list(t2['owner_clean'])};{list(t2['clean_address'])}")

    ## remove double entries from df_comp_comm
    print("\tDAFNE entries before dropping duplicates:", len(df_comp_comm))
    df_comp_comm.drop_duplicates(subset=["conn_left", comm_col], inplace=True)
    print("\tDAFNE entries after dropping duplicates:", len(df_comp_comm))

    ## Get all network connections that are in several communities
    ## This happens due to different writings prior to this script and its cleaning
    from collections import Counter
    counts = Counter(df_comp_comm["conn_left"])
    counts = [key for key in counts if counts[key] > 1]

    ## I assume that these connections actually build a bridge between the different communities they are in
    ## Thus I correct the community
    for val in counts:
        comms = df_comp_comm.loc[df_comp_comm["conn_left"] == val, comm_col].tolist()
        comm1 = comms[0]
        comms.remove(comm1)
        df_comp_comm.loc[df_comp_comm[comm_col].isin(comms), comm_col] = comm1
        ## The dictionary created in create_comp_community_mcomp_dictionary() is only
        ## used prior to this and does not need any correction
        ## ToDo clean aristrocratic DAFNE entries before network analysis (08_2_network_analysis.py)
        ## repeat all this afterwards.

    print("\tDAFNE entries before dropping duplicates:", len(df_comp_comm))
    df_comp_comm.drop_duplicates(subset=["conn_left", comm_col], inplace=True)
    print("\tDAFNE entries after dropping duplicates:", len(df_comp_comm))


    ## Merge ALKIS and DAFNE
    df_alk = pd.merge(df_alk, df_comp_comm, how="left", on="conn_left")
    print("\tLength of table after merging with DAFNE:", len(df_alk))

    ## Check for completeness
    # t = df_alk.loc[df_alk[comm_col].notna()].copy()
    # t2 = t.loc[t["conn_left"].str.count('_') > 0].copy()
    # t31 = df_alk.loc[df_alk[comm_col].isna()].copy()
    # t3 = t31.loc[t31["conn_left"].str.count('_') == 0].copy()
    # t3 = t3.loc[t3["level1"] == 2].copy()

    ## It is possible that two different owners with same family name and address were assigned to different communities
    ## In this case take the community with more parcels and assign it to the other
    df_comms = df_alk.loc[df_alk[comm_col].notna()].copy()
    df_comms.drop_duplicates(subset=['owner_merge', comm_col], inplace=True)
    print("\t1. Leschke has comm:", 'leschkeberlinerstr5604916herzberg' in df_comms["owner_merge"].tolist())

    print("\tNumber of unique owner_merge community combinations before cleaning:", len(df_comms))
    multi_count = Counter(df_comms["owner_merge"])
    multis = [key for key in multi_count if multi_count[key] > 1]
    df_comms_multis = df_alk.loc[df_alk["owner_merge"].isin(multis)].copy()
    print("\tNumber of owner_merge names that were assigned to different communities:", len(multis))
    print("\t2. Leschke has comm:", 'leschkeberlinerstr5604916herzberg' in df_comms["owner_merge"].tolist())
    print(df_alk.loc[df_alk["owner_merge"] == 'leschkeberlinerstr5604916herzberg', comm_col].unique())

    def get_community_with_most_occurences(name_group):

        counts = Counter(name_group[comm_col])
        max_key = max(counts, key=counts.get)

        return max_key

    om_to_comm = df_comms_multis[["owner_merge", comm_col]].groupby("owner_merge").apply(
        get_community_with_most_occurences).reset_index()
    om_to_comm.columns = ["owner_merge", comm_col]
    ind = list(om_to_comm.columns).index(comm_col) + 1
    om_to_comm_dict_multis = {row.owner_merge: row[ind] for row in om_to_comm.itertuples()}
    df_alk.loc[df_alk["owner_merge"].isin(multis), comm_col] = df_alk.loc[
        df_alk["owner_merge"].isin(multis), "owner_merge"].map(om_to_comm_dict_multis)
    print("\t", df_alk.loc[df_alk["owner_merge"] == 'leschkeberlinerstr5604916herzberg', comm_col].unique())

    print("\tNumber of unique owner_merge community combinations after cleaning:", len(df_alk.loc[df_alk[comm_col].notna()].drop_duplicates(subset=['owner_merge', comm_col])))
    print("\tLength of table after cleaning step 1:", len(df_alk))

    ## It is possible that from different owners with same family name and address, one was assigned
    ## to a community and the other(s) not. Assign the community to the others
    ## e.g. in case of networks WITH threshold:
    ## andreas ebel_beindersheim is in same community as rinderzucht lanz - lenzen ag,
    ## matthias ebel_beindersheim lives in same address as andreas (i.e. he is family) BUT is not part of the community
    ## He should be assigned to the same community and assigned to the same company
    ## owner_merge == ebelsiemensstr867259beindersheim

    df_comms = df_alk.loc[df_alk[comm_col].notna()].copy()
    print("\t3. Leschke has comm:", 'leschkeberlinerstr5604916herzberg' in df_comms["owner_merge"].tolist())
    df_comms.drop_duplicates(subset=['owner_merge', comm_col], inplace=True)
    print("\t4. Leschke has comm:", 'leschkeberlinerstr5604916herzberg' in df_comms["owner_merge"].tolist())


    print("\tNumber of unique communities in alkis data:", len(df_comms))
    ind = list(df_comms.columns).index(comm_col) + 1
    om_to_comm_dict = {row.owner_merge: row[ind] for row in df_comms.itertuples()}
    oms = [key for key in om_to_comm_dict]
    print("\tNumber of owner_merge names with a community number:", len(oms))
    print("\tLength of ALKIS data with a community number before correction:",
          len(df_alk.loc[df_alk[comm_col].notna()]))
    df_alk.loc[df_alk["owner_merge"].isin(oms), comm_col] = df_alk.loc[
        df_alk["owner_merge"].isin(oms), "owner_merge"].map(om_to_comm_dict)
    print("\tLength of ALKIS data with a community number after correction:",
          len(df_alk.loc[df_alk[comm_col].notna()]))
    print("Length of table after cleaning step 2:", len(df_alk))

    ## Set owner clean merge as mother company name, but for cases that already have a community number
    ## use the mother company name from the network data
    print("\tAssign mother company name to all community numbers")
    df_alk["mother_company"] = df_alk["owner_merge"]
    df_alk.loc[df_alk[comm_col].notna(), "mother_company"] = df_alk.loc[
        df_alk[comm_col].notna(), comm_col].map(comm_dict)
    print("\tLength of ALKIS data with a mother company name after assignment:",
          len(df_alk.loc[df_alk["mother_company"].notna()]))

    ## For all owner that dont belong to a community assign a community number
    print("\tAssign community number to all owners without a number yet.")
    no_comm = list(df_alk.loc[df_alk[comm_col].isna(), "owner_merge"].unique())
    no_comm_to_comm_dict = {name: i + 9999 for i, name in enumerate(no_comm)}
    df_alk.loc[df_alk[comm_col].isna(), comm_col] = df_alk.loc[df_alk[comm_col].isna(), "owner_merge"].map(no_comm_to_comm_dict)
    print("\tLength of ALKIS data with a community number after assignment:",
           len(df_alk.loc[df_alk[comm_col].notna()]))
    print("\tLength of table after cleaning step 3:", len(df_alk))

    ## get all cases that have the same owner merge name but different mother company names
    ## Create a combined mother company name
    ## Assign combined mother company name to the identified cases
    df_uni_om = df_alk.drop_duplicates(subset=["owner_merge", "mother_company"])
    multi_count = Counter(df_uni_om["owner_merge"])
    multis = [key for key in multi_count if multi_count[key] > 1]
    df_uni_om = df_uni_om.loc[df_uni_om["owner_merge"].isin(multis)].copy()

    if not df_uni_om.empty:

        def create_combined_mcomp_name(name_group):

            uni_names = list(name_group["mother_company"].unique())
            mcomp = ' AND '.join(uni_names)

            return mcomp

        om_to_mcomp = df_uni_om[["owner_merge", "mother_company"]].groupby("owner_merge").apply(create_combined_mcomp_name).reset_index()
        om_to_mcomp.columns = ["owner_merge", "mother_company"]

        om_to_mcomp_dict = {row.owner_merge: row.mother_company for row in om_to_mcomp.itertuples()}

        df_alk.loc[df_alk["owner_merge"].isin(multis), "mother_company"] = df_alk.loc[
                df_alk["owner_merge"].isin(multis), "owner_merge"].map(om_to_mcomp_dict)

    print("\tLength of table after cleaning step 4:", len(df_alk))

    ## identify the cases where different owners with the same owner merge were assigned to different communities
    ## they did not meet the criteria for belonging to the same comm.
    ## But in the process above they were assigned the same mother company, so we need to assign the to the same community
    df_uni_om = df_alk.drop_duplicates(subset=["owner_merge", comm_col])
    multi_count = Counter(df_uni_om["owner_merge"])
    multis = [key for key in multi_count if multi_count[key] > 1]

    df_uni_om = df_uni_om.loc[df_uni_om["owner_merge"].isin(multis)].copy()

    if not df_uni_om.empty:
        helper_functions.print_red("\tWARNING there are cases where two different owner from the same onwer_merge were assigned to different communities")

    ## Create a dictionary that holds for each community number the mother company name
    df_comms_fin = df_alk.drop_duplicates(subset=[comm_col, "mother_company"])

    if len(df_comms_fin[comm_col].unique()) > len(df_comms_fin["mother_company"].unique()):
        helper_functions.print_red("\nWARNING: There are more communities than mother companies\n")

        print("\tNumber of communities:", len(df_comms_fin[comm_col].unique()))
        print("\tNumber of mother companies:", len(df_comms_fin["mother_company"].unique()))

        multi_communities = Counter(df_comms_fin[comm_col])
        multi_communities = {comm: multi_communities[comm] for comm in multi_communities if multi_communities[comm] > 1}
        len(multi_communities)

        multi_mcomps = Counter(df_comms_fin["mother_company"])
        multi_mcomps = {mcomp: multi_mcomps[mcomp] for mcomp in multi_mcomps if multi_mcomps[mcomp] > 1}
        len(multi_mcomps)

    print("\tWrite ALKIS dataframe out.")
    alk_to_comm = df_alk.drop_duplicates(subset=["owner_merge", comm_col])
    alk_to_comm = alk_to_comm.set_index('owner_merge').to_dict()[comm_col]

    ## Add no. of hierachical levels to ALKIS data
    df_alk["netw_max_dist"] = df_alk[comm_col].map(max_dist_dict)
    df_alk.loc[df_alk["netw_max_dist"].isna(), "netw_max_dist"] = 0

    print("\tLength of table at writing out:", len(df_alk))
    df_alk.to_csv(out_pth, sep=';', index=False)
    print("\tWriting done.")

    print("\tCreate an information dictionary holding the names and the numbers of all sub companies/owners for each community")
    df_comms = df_alk.drop_duplicates(subset=[comm_col, "conn_left", "owner_merge", "mother_company"])
    print(len(df_comms["conn_left"]), len(df_comms))
    ## Get number of subcompanies per community
    comms_with_multi_subs = Counter(df_comms[comm_col])
    num_subcomps_dict = {comm: comms_with_multi_subs[comm] for comm in comms_with_multi_subs} #if comms_with_multi_subs[comm] > 1

    ## Print all communities with more than x (e.g. 10) subcompanies that occur in ALKIS
    # n = {comm: comms_with_multi_subs[comm] for comm in comms_with_multi_subs if comms_with_multi_subs[comm] > 10}
    # if len(n) > 0:
    #     print("Community")
    #     for comm in n:
    #         print("\t", comm, df_comms.loc[df_comms[comm_col] == comm, "owner_merge"].tolist())

    comm_dict_df = df_comms.loc[df_comms[comm_col].isin(num_subcomps_dict.keys())].copy()
    ind = list(comm_dict_df.columns).index(comm_col) + 1
    comm_dict = {row[ind]: row.mother_company for row in comm_dict_df.itertuples()}

    print("\tCreate subdictionary with names of subcompanies per community")
    comm_ids = list(comm_dict.keys())
    df_comm_ids = df_comms.loc[df_comms[comm_col].isin(comm_ids)]
    def get_names_of_subcomps(comm_group, name_col):
        names = list(comm_group[name_col].unique())
        names = '\n'.join(names)
        return names

    df_names_sc = df_comm_ids.groupby(comm_col).apply(get_names_of_subcomps, name_col="conn_left").reset_index()
    df_names_sc.columns = [comm_col, "names"]
    ind = list(df_names_sc.columns).index(comm_col) + 1
    name_subcomps_dict = {row[ind]: row.names for row in df_names_sc.itertuples()}

    print("\tCreate output dictionary including number and names per community")
    comm_dict = {str(key): comm_dict[key] for key in comm_dict}
    num_subcomps_dict = {str(key): num_subcomps_dict[key] for key in num_subcomps_dict}
    name_subcomps_dict = {str(key): name_subcomps_dict[key] for key in name_subcomps_dict}

    comm_info_dict = {"comm_dict": comm_dict,
                      "num_subcomps_dict": num_subcomps_dict,
                      "name_subcomps_dict": name_subcomps_dict}

    print("\tWrite json out")
    with open(community_info_dict_pth, 'w') as fp:
        json.dump(comm_info_dict, fp, indent=4)

    print("\tDone!")


def plot_network_graph_of_community_from_df(df, community_column, community_id, out_pth):
    """
    Plots the network graph for a community.
    Args:
        df: Network dataframe with companies and community IDs.
        community_column: Column name for community IDs.
        community_id: ID of community for which the plot should be generated.
        out_pth: Output path for figure.

    Returns:

    """
    df = df.loc[df[community_column] == community_id].copy()

    function_to_value = {"Aktionär": 1,
                         "Aufsichtsrat": 2,
                         "Geschäftsführender Direktor": 3,
                         "Geschäftsführer": 4,
                         "Gesellschafter": 5,
                         "Kommanditaktionär": 6,
                         "Kommanditist": 7,
                         "Komplementär": 8,
                         "Liquidator": 9,
                         "mother company": 10,
                         "persönlich haftender Gesellschafter": 11,
                         "Prokurist": 12,
                         "subsidary": 13,
                         "Verwaltungsrat": 14,
                         "Vorstand": 15}

    function_to_color = {"Aktionär": "#ffeda0", #gelb
                         "Aufsichtsrat": "#bdbdbd", #grau
                         "Geschäftsführender Direktor": "#a1d99b", #grün
                         "Geschäftsführer": "#a1d99b", #grün
                         "Gesellschafter": "#feb24c", #orange
                         "Kommanditaktionär": "#ffeda0", #gelb
                         "Kommanditist": "#feb24c", #orange
                         "Komplementär": "#feb24c", #orange
                         "Liquidator": "#bdbdbd", #grau
                         "mother company": "#f03b20", #darkorange
                         "persönlich haftender Gesellschafter": "#feb24c", #orange
                         "Prokurist": "#bdbdbd", #grau
                         "subsidary": "#f03b20", #darkorange
                         "Verwaltungsrat": "#bdbdbd", #grau
                         "Vorstand": "#bdbdbd"} #grau

    # df["value"] = pd.Categorical(df[color_column])
    # df['value'].cat.codes
    df["value"] = df["function"].map(function_to_value)
    df["color"] = df["function"].map(function_to_color)
    df.rename(columns={"share": "weight"}, inplace=True)

    df_len = len(df)

    # Build your graph
    G = nx.from_pandas_edgelist(df, 'conn_left', 'conn_right', ["weight", "function", "value", "color"])

    # # Custom the nodes:
    # fig, ax = plt.subplots()
    # nx.draw(G, ax=ax, with_labels=True, node_color='skyblue', node_size=250, edge_color=df['value'].cat.codes,
    #         width=3.0, edge_cmap=plt.cm.Set2, arrows=True, pos=nx.fruchterman_reingold_layout(G), font_size=8)
    # ax.set_xlim([1.32 * x for x in ax.get_xlim()])
    # ax.set_ylim([1.32 * y for y in ax.get_ylim()])
    # plt.legend()
    # # plt.margins(x=0.4)
    # # fig.tight_layout()
    # fig.savefig(out_pth)
    # plt.close()

    side_len = 11.5 * math.exp(0.0050 * df_len)

    fig, ax = plt.subplots(figsize=plotting_lib.cm2inch(side_len, side_len))

    edges = [(u, v) for (u, v, d) in G.edges(data=True)]
    colors = [d["color"] for (u, v, d) in G.edges(data=True)]
    # elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] >= 50]
    # emedium = [(u, v) for (u, v, d) in G.edges(data=True) if 25 <= d["weight"] < 50]
    # esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] < 25]

    pos = nx.spring_layout(G, seed=1, k=1, iterations=2)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=250, node_color="#e6550d")
    nx.draw_networkx_edges(G, pos, ax=ax, edgelist=edges, edge_color=colors, width=2, arrows=True)
    # nx.draw_networkx_edges(
    #     G, pos, ax=ax, edgelist=elarge, edge_color="#f03b20", width=2, arrows=True)
    # nx.draw_networkx_edges(
    #     G, pos, ax=ax, edgelist=emedium, edge_color="#feb24c", width=2, arrows=True)
    # nx.draw_networkx_edges(
    #     G, pos, ax=ax, edgelist=esmall, edge_color="#ffeda0", width=2, arrows=True)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=8, font_family="sans-serif") # node labels
    edge_labels = nx.get_edge_attributes(G, "weight") # edge weight labels
    nx.draw_networkx_edge_labels(G, pos, edge_labels, ax=ax, font_size=8)

    ax = plt.gca()
    ax.margins(0.08)
    plt.axis("off")
    plt.tight_layout()
    # plt.legend()
    # plt.show()
    out_pth = out_pth.replace('"', '')
    plt.savefig(out_pth)


def get_community_id_of_company(df, name, comm_column):
    """
    Provides the community ID for a company name.
    Args:
        df: Network dataframe with companies and community IDs
        name: Company name.
        comm_column: Column name for community IDs.

    Returns:

    """

    if name in df["conn_left"].tolist():
        comms = df.loc[df["conn_left"] == name, comm_column].to_list()
        return comms
    elif name in df["conn_right"].tolist():
        comms = df.loc[df["conn_right"] == name, comm_column].to_list()
        return comms
    else:
        print(f"{name} not found in df.")


def get_common_owners(net_comm_pth, out_pth):
    """
    Get global ultimate owners/investors that invest in multiple communities but are below the threshold
    Args:
        net_comm_pth: Network dataframe with companies and community IDs
        out_pth: Output path to dataframe.

    Returns:

    """

    df_netw = pd.read_csv(net_comm_pth)

    df1 = df_netw[["conn_left", "community_0", "community_50", "community_25"]]
    df1.columns = ["conn", "community_0", "community_50", "community_25"]
    df2 = df_netw[["conn_right", "community_0", "community_50", "community_25"]]
    df2.columns = ["conn", "community_0", "community_50", "community_25"]

    df = pd.concat([df1, df2], axis=0)

    def count_unique_values_in_colum(group, column_name):

        unique_communities = list(group[column_name].unique())
        unique_communities = [comm for comm in unique_communities if comm[-3:] != '_99']
        count = len(unique_communities)
        return count

    def return_unique_values_of_column(group, column_name):
        values = list(group[column_name].unique())
        values = [str(val) for val in values]
        values = [val for val in values if val[-3:] != '_99']
        values = '|'.join(values)
        return values

    df_count = df.groupby("conn").apply(count_unique_values_in_colum, "community_50").reset_index()
    df_count.columns = ["name", "count_of_communities"]
    df_values = df.groupby("conn").apply(return_unique_values_of_column, "community_50").reset_index()
    df_values.columns = ["name", "communites"]

    df_out = pd.merge(df_count, df_values, on="name", how="left")
    df_out.sort_values(by="count_of_communities", ascending=False, inplace=True)
    df_out.to_csv(out_pth)


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    #################################### Cleaning network data ####################################
    # company_information_to_long_data_frame(
    #     dafne_search_results=DAFNE_PTH,
    #     col_renaming_json_pth=COL_RENAMING,
    #     network_data_long_df_pth=NETW_LONG_PTH
    # )
    #
    # identify_possible_alkis_matches_for_manual_assignment(
    #     dafne_search_results=DAFNE_PTH,
    #     network_data_long_df_pth=NETW_LONG_PTH,
    #     alkis_pth=ALKIS_PTH,
    #     possible_matches_folder=POSSIBLE_MATCHES_FOLDER
    # )
    #
    # clean_long_network_table(
    #     network_data_long_df_pth=NETW_LONG_PTH,
    #     corrected_dafne_alkis_matches_json=CORRECTED_MATCHES_PTH,
    #     prepared_network_data_pth=NETW_PTH
    # )
    #
    # network_connections_to_list_for_dafne_search(
    #     prepared_network_data_pth=NETW_PTH,
    #     dafne_search_results=DAFNE_PTH,
    #     col_renaming_json_pth=COL_RENAMING,
    #     dafne_search_list_second_round=DAFNE_SEARCH_PTH,
    #     folder_second_search=FOLDER_DAFNE_SECOND_SEARCH
    # )
    #
    # add_economic_branch_to_network_companies(
    #     folder_second_search=FOLDER_DAFNE_SECOND_SEARCH,
    #     df_companies_with_branches_and_locations_pth=fr"08_network_analysis\all_companies_branches_and_locations.csv"
    # )
    #
    # #################################### Network analysis ####################################
    # network_analysis_advanced(
    #     netw_pth=NETW_PTH,
    #     out_pth=NETW_COMM_PTH
    # )
    ## Add the community IDs to ALKIS owners
    for threshold in [0, 25, 50]: # 10 could also be used.

        print(f"Add the community IDs to ALKIS owners with threshold {threshold}\n ")

        create_comp_community_mcomp_dictionary(
            net_comm_pth=NETW_COMM_PTH,
            comm_col=f"community_{threshold}",
            out_pth=COMP_COMMUNITY_MCOMP_DICT.format(threshold))

        get_maximum_hierarchical_distance_in_network(
            net_comm_pth=NETW_COMM_PTH,
            comm_col=f"community_{threshold}",
            out_pth=COMMUNITY_MAX_DIST_DICT_PTH.format(threshold)
        )

        add_community_number_to_alkis_data(
            alkis_pth=ALKIS_PTH,
            comp_comm_mcomp_dict_pth=COMP_COMMUNITY_MCOMP_DICT.format(threshold),
            comm_col=f"community_{threshold}",
            community_info_dict_pth=COMMUNITY_INFO_DICT.format(threshold),
            max_dist_dict_pth=COMMUNITY_MAX_DIST_DICT_PTH.format(threshold),
            out_pth=OWNERS_STRETCHED_COMM.format(threshold)
        )

    #################################### Visualization of network areas and network compostions ########################
    df = pd.read_csv(NETW_COMM_PTH)
    helper_functions.create_folder(r"08_network_analysis\figures")

    ## Examples
    name ='landwirtschaftliche vermoegensverwaltungsgesellschaft seelow mbh'

    ## Get community numbers for company name and plot
    comm_numbers = get_community_id_of_company(
        df=df,
        name=name,
        comm_column="community_0"
    )
    for comm_number in comm_numbers:
        plot_network_graph_of_community_from_df(
            df=df,
            community_column="community_0",
            community_id=comm_number,
            out_pth=fr"08_network_analysis\figures\{name}_main_community_{comm_number}.png"
        )

    ## Get sub-community number for company name and plot
    comm_numbers = get_community_id_of_company(
        df=df,
        name=name,
        comm_column="community_50"
    )
    for comm_number in comm_numbers:
        plot_network_graph_of_community_from_df(
            df=df,
            community_column="community_50",
            community_id=comm_number,
            out_pth=fr"08_network_analysis\figures\{name}_sub_community_{comm_number}.png"
        )

    ## BayWa Aktiengesellschaft
    comm_numbers = ["40_2", "40_1", "40_0"]
    for comm_number in comm_numbers:
        plot_network_graph_of_community_from_df(
            df=df,
            community_column="community_50",
            community_id=comm_number,
            out_pth=fr"08_network_analysis\figures\baywa aktiengesellschaft_sub_community_{comm_number}.png"
        )

    #################################### Some further but unused analysis ##############################################
    ## Get investors that invest in multiple communities but are below the threshold
    # # # get_common_owners(
    # # #     net_comm_pth=NETW_COMM_PTH,
    # # #     out_pth=""
    # # # )

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)



if __name__ == '__main__':
    main()
