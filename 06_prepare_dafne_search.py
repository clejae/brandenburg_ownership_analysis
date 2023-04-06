# Clemens Jänicke
# github Repo: https://github.com/clejae

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import pandas as pd
import math
import glob
import os
import time

## Project library
import helper_functions
# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input path
COMPANY_OWNERS_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified_private-companies.csv"
NONPROF_OWNERS_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified_non-profit-etc.csv"
# MATCHES_FOLDER = r"04_owner_class_reclassification"
NAME_COL = 'owner_clean'

## Output (the output folder will be created automatically!)
OUT_FOLDER = r"06_prepare_dafne_search"

# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#

def prepare_lists(owners_stretched_sub_pth, name_column, out_folder, sub_descr):
    """
    Get all names of owners, clean them to increase the likelihood of DAFNE search htis and save them in lists of 1000s.
    DAFNE only allows batch searches of 1000 names.
    :param owners_stretched_sub_pth: Path to dataframe with owner names.
    :param name_column: Column name of owner names
    :param out_folder: Folder in which the lists should be saved
    :param sub_descr: Description of the lists for output files.
    :return:
    """

    ## Read list of companies
    df1 = pd.read_csv(owners_stretched_sub_pth, sep=';')
    df1 = df1[[name_column]]
    df1.drop_duplicates(subset=name_column, inplace=True)

    ## Read already matched companies
    lst = glob.glob(rf"{WD}\{OUT_FOLDER}\*matches_v1.xlsx")
    df_lst = [pd.read_excel(pth) for pth in lst]
    df_done = pd.concat(df_lst, axis=0)

    ## Clean names of companies
    clean_name_col = 'clean_name'

    df1[clean_name_col] = df1[name_column].apply(helper_functions.remove_search_term_and_followup, substring='mit dem sitz in')
    df1[clean_name_col] = df1[clean_name_col].apply(helper_functions.remove_search_term_and_followup, substring='mit sitz in')
    df1[clean_name_col] = df1[clean_name_col].apply(helper_functions.remove_search_term_and_followup, substring=', sitz')

    address_code_words = ['straße', 'strasse', 'weg', ' zum ', 'dorf', 'ausbau', 'chausee', ' ot ', ' am ', ' an ']
    df1[clean_name_col] = df1[clean_name_col].apply(helper_functions.remove_address_part, delimiter=',', address_code_words=address_code_words)

    df1.drop_duplicates(subset=clean_name_col, inplace=True)
    df1 = df1.loc[~df1[clean_name_col].isin(df_done['Unternehmensname'])].copy()

    x = int(lst[-1].split('_')[-3])
    x = 9

    num_rows = len(df1)
    num_lists = math.ceil(num_rows / 1000)

    s = 0
    for i in range(1, num_lists + 1):
        e = i * 1000
        df_out = df1[s:e]
        s += 1000
        # out_pth = rf"{WD}\{OUT_FOLDER}\batch_search_companies_{i+x:02d}.txt"
        out_pth = rf"{out_folder}\batch_search_{sub_descr}_{i + x:02d}.txt"
        df_out[[clean_name_col]].to_csv(out_pth, sep=';', header=None, index=None)


def prepare_delivery(out_folder):

    """
    I can't remember for what reason I did this function.
    :return:
    """

    ## Read already matched companies
    lst = glob.glob(rf"{out_folder}\*matches_v1.xlsx")
    df_lst = [pd.read_excel(pth) for pth in lst]
    for i, df in enumerate(df_lst):
        df["ID"] = i
    df_done = pd.concat(df_lst, axis=0)
    df_done.loc[df_done["ID"] < 8, "type"] = "company"
    df_done.loc[df_done["ID"] >= 8, "type"] = "vereine, stiftungen, non-profit, etc"
    df_done.drop(columns=["ID"], inplace=True)

    df_done.drop_duplicates(subset="Unternehmensname", inplace=True)

    df_matches = df_done.loc[df_done['Gematchte BvD ID'].notna()].copy()
    df_misses = df_done.loc[df_done['Gematchte BvD ID'].isna()].copy()

    df_misses.reset_index(inplace=True)
    df_matches.reset_index(inplace=True)
    df_misses.drop(columns="index", inplace=True)
    df_matches.drop(columns="index", inplace=True)

    len(df_done.loc[df_done['type'] == 'company'])

    out_pth = rf"{out_folder}\Unternehmen.xlsx"
    with pd.ExcelWriter(out_pth) as writer:
        df_matches.to_excel(writer, sheet_name="Matches")
        df_misses.to_excel(writer, sheet_name="Misses")


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    helper_functions.create_folder(OUT_FOLDER)

    prepare_lists(
        owners_stretched_sub_pth=NONPROF_OWNERS_PTH,
        name_column=NAME_COL,
        out_folder=OUT_FOLDER,
        sub_descr="non_profit")
    prepare_lists(
        owners_stretched_sub_pth=COMPANY_OWNERS_PTH,
        name_column=NAME_COL,
        out_folder=OUT_FOLDER,
        sub_descr="companies")

    prepare_delivery()

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)


if __name__ == '__main__':
    main()