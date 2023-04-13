# Clemens Jänicke
# github Repo: https://github.com/clejae

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#

import time
import math
import pandas as pd
import os
import json
import glob

# ------------------------------------------ USER INPUT ------------------------------------------------#
import helper_functions

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
UNINAMES_PTH = r"08_network_analysis\uni_names.csv"
POSSIBLE_MATCHES_FOLDER = r"08_network_analysis"
DAFNE_SEARCH_PTH = r"08_network_analysis\dafne_search_second_round.csv"


# AGRI_COMP_EXT_PTH = r"07_owner_name_cleaning\matching_results_alkis_farmsubs_futtermittel_extension.csv"

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

    #################################### Network analysis ####################################

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)



if __name__ == '__main__':
    main()
