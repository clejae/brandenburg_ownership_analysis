# Author:
# github repository:

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import os
import time
import pandas as pd
import json

## Project library
import helper_functions

# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input
OWNERS_STRETCHED_PTH = r"02_identify_unique_owner_address_combinations\02_owners_and_addresses_stretched.csv"
CLASS_IDS_PTH = r'00_data\classification\class_ids_l3_preliminary.json'

## Output (the output folder will be created automatically!)
OWNERS_PRELIM_CLASSIFIED_PTH = r"03_owner_name_classification\03_owners_stretched_preliminary_classication.csv"
OUT_STATS = r"03_owner_name_classification\03_owners_stretched_preliminary_classication_stats.csv"
# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#

def classify_owner_names_into_classes(owner_stretched_pth, class_id_json_pth, owners_prelim_classified_pth):
    """
    Classifies owners into several preliminary classes that are specified in a json-file with help of code word lists.
    Classes are private people, different public institutions, company forms, religious institutions, non-profit
    organizations and other associations and foundations.
    :param owner_stretched_pth: Path to stretched owner dataframe.
    :param class_id_json_pth: Path to json, holding a dictionary with preliminary classes and their IDs.
    :param owners_prelim_classified_pth: Output path for dataframe with classified owners.
    :return:
    """
    print("Classify owner names into preliminary classes.")
    ## Read input
    print("\tRead data.")
    df = pd.read_csv(owner_stretched_pth, sep=';')
    with open(class_id_json_pth) as json_file:
        class_ids = json.load(json_file)

    ## Clean owner string
    df['owner_clean'] = df['owner_names'].str.lower()

    ## Remove strings
    # Don't remove '&' and '\+', because they are needed to identify mixed legal forms (& co ..)
    df['owner_clean'] = df['owner_clean'].apply(
        helper_functions.replace_characters_in_string,
        char_lst= ['.', ';', r'\\', '/', ':', '"', "'", '(', ')', '\+', '?'],
        replace_string='')

    df['owner_clean'] = df['owner_clean'].str.replace('ä', 'ae')
    df['owner_clean'] = df['owner_clean'].str.replace('ö', 'oe')
    df['owner_clean'] = df['owner_clean'].str.replace('ü', 'ue')
    df['owner_clean'] = df['owner_clean'].str.replace('ß', 'ss')

    ## Replace common word chains with typical abbreviations
    df['owner_clean'] = df['owner_clean'].str.replace('gesellschaft buergerlichen rechts', 'gbr')
    df['owner_clean'] = df['owner_clean'].str.replace('gesellschaft mit beschraenkter haftung', 'gmbh')
    df['owner_clean'] = df['owner_clean'].str.replace('mit beschraenkter haftung', 'mbh')
    df['owner_clean'] = df['owner_clean'].str.replace('kommanditgesellschaft', 'kg')
    df['owner_clean'] = df['owner_clean'].str.replace('kommandit-gesellschaft', 'kg')
    df['owner_clean'] = df['owner_clean'].str.replace('eingetragene genossenschaft', 'eg')
    df['owner_clean'] = df['owner_clean'].str.replace('aktiengesellschaft', 'ag')
    df['owner_clean'] = df['owner_clean'].str.replace('deutsche bahn', 'db')
    df['owner_clean'] = df['owner_clean'].str.replace('bundesrepublik', 'brd')
    df['owner_clean'] = df['owner_clean'].str.replace('eigentum des volkes', 'edv')
    df['owner_clean'] = df['owner_clean'].str.replace('gemeinnuetzige gmbh', 'ggmbh')
    df['owner_clean'] = df['owner_clean'].apply(helper_functions.remove_search_term_and_followup, search_term='mit dem sitz in')
    df['owner_clean'] = df['owner_clean'].apply(helper_functions.remove_search_term_and_followup, search_term='mit sitz in')
    df['owner_clean'] = df['owner_clean'].apply(helper_functions.remove_search_term_and_followup, search_term=', sitz')

    address_code_words = ['straße', 'strasse', 'weg', ' zum ', 'dorf', 'ausbau', 'chausee', ' ot ', ' am ', ' an ', 'str ']
    df['owner_clean'] = df['owner_clean'].apply(helper_functions.remove_address_part, delimiter=',',
                                                address_code_words=address_code_words)

    df.loc[df['owner_clean'].isna(), 'owner_clean'] = 'unbekannt'

    ## Classify all private persons based on asterisk
    df['category'] = 0
    df.loc[df['owner_clean'].str.count('\*') == 1, 'category'] = 1  # private person
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())

    ## Classify different "Gesellschaftsformen" with help of code words
    # All mixed legal forms
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['& co', ' &co', '& c o', '6 co', '+ co', '& go', 'und co', ' u co'],
        return_code=class_ids['mixed'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    ## gemeinnützige
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['ggmbh', 'gemeinnuetzige'],
        return_code=class_ids['gemeinn'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # BVVG
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['bvvg', 'treuhandanstalt', 'bodenverwert'],
        return_code=class_ids['bvvg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # GbR
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['gbr'],
        return_code=class_ids['gbr'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # OHG
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['ohg'],
        return_code=class_ids['ohg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # KG
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['kg', 'kommanditgesellschaft'],
        return_code=class_ids['kg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # EWIV
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['ewiv'],
        return_code=class_ids['ewiv'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # AG
    df.loc[df['category'] == 0, 'category']=df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['ag'],
        return_code=class_ids['ag'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['aktiengesellschaft', 'agraraktiengesellschaft'],
        return_code=class_ids['ag'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # GmbH
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['gmbh', 'mbh', 'mit beschraenkter haf'],
        return_code=class_ids['gmbh'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # ug
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['ug'],
        return_code=class_ids['ug'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['ug haftungsbeschraenkt', 'unternehmergesellschaft'],
        return_code=class_ids['ug'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # SE
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['se'],
        return_code=class_ids['se'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # limited
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['limited', 'sa', 'sarl', 'sárl', 'sàrl', 'holding', 'ltd'],
        return_code=class_ids['lim'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # eG
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['eg'],
        return_code=class_ids['eg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=[' e g'],
        return_code=class_ids['eg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Kirche
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2,
        search_terms=['christlich', 'evangelisch', 'katholisch', 'kirche', 'diakonie', 'pfarr', 'pfarrgemeinde', 'abtei',
                        'pfarrstelle', 'diakonisch',  'diaconat', 'diakonat', 'domstift', 'kantorat', 'predigerstelle',
                        'stift zum', 'juedisch', 'kirchgemeinde', 'hospital', 'jewish'],
        return_code=class_ids['kirch'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # e.V.
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['ev', 'verein'],
        return_code=class_ids['ev'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    df.loc[df['category'] == 0 , 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=[' e v', 'nabu', 'naturschutzbund'],
        return_code=class_ids['ev'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # w.V.
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['wv'],
        return_code=class_ids['wv'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Stiftungen
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['stiftung', 'naturschutzfonds', 'wwf'],
        return_code=class_ids['stift'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Nicht mehr existierende Institutionen
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['edv', 'rt', 'lpg'],
        return_code=class_ids['verg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['eigentum des volkes', 'des volk', 'separation', 'volkseigentum', 'rezess'],
        return_code=class_ids['verg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Unbekannt
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['unbekannt', 'nicht ermittelt', 'separation', 'aufgegeben', 'verstorben',
                                             'herrenlos', 'ermittel', 'nicht erfasst', 'seperation'],
        return_code=class_ids['unbek'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Gemeinden
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['gemeinde'],
        return_code=class_ids['gem'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2,
        search_terms=['stadtgemeinde', 'dorfgemeinde', 'landgemeinde', 'gemeindemitglieder', 'geimeindeverwaltung',
                        'anlieger', 'anliegenden', 'angrenzenden', 'oeffentlich'],
        return_code=class_ids['gem'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Land
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['land brandenburg', 'land brandenbung', 'landkreis', 'land berlin',
                                             'freistaat bayern', 'land hessen', 'land niedersachsen', 'freistaat thueringen',
                                             'landesstrassenverwaltung', 'landesbetrieb', 'landesregierung',
                                             'landesvermessung', 'nordrhein-westfalen', 'land baden',
                                             'land sachsen-anhalt'],
        return_code=class_ids['land'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Bund
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['bundesrepublik', 'bundesanstalt', 'bundesstrassenverwaltung', 'brd',
                                             'bundesfinanzverwaltung', 'bundesministerium', 'bunderepublik', ],
        return_code=class_ids['bund'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Zweckverband
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['zweckverband', 'wasserverband', 'verband'],
        return_code=class_ids['zweck'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['verband'],
        return_code=class_ids['zweck'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Erbengemeinschaften
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['erbengemein'],
        return_code=class_ids['erben'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())

    ## Second round, now confusion is less likely and check_occ_of_words_v2 can be applied for all
    # GbR
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['bgb', 'gbr', 'gesellschaft buergerlichen r', 'landwirtschaftsbetrieb'],
        return_code=class_ids['gbr'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Gemeinde
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['stadt', 'amt'],
        return_code=class_ids['gem'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # EG
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['agrargenossenschaft', 'produktions', 'bauerngenossenschaft',
                                             'weidegenossenschaft', 'waldgenossenschaft', 'fischergenossenschaft',
                                             'gaertnergenossenschaft', 'huefnergenossenschaft', 'ackerbuergergenossenschaft'],
        return_code=class_ids['eg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Vereine
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['verein', 'club'],
        return_code=class_ids['ev'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # limited
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['bv', 'holding'],
        return_code=class_ids['lim'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Zweckverband
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['verband', 'wohnungseigentuemergemeinschaft', 'interessengemeinschaft',
                                             'teilnehmergemeinschaft', 'guetergemeinschaft'],
        return_code=class_ids['zweck'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # OHG
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['ohg'],
        return_code=class_ids['ohg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Nicht mehr existierende Institutionen
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v1, search_terms=['ddr'],
        return_code=class_ids['verg'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())

    # Privatpersonen
    # All remaining owner names with a comma should be private owners, which don't have a birthdate
    df.loc[(df['category'] == 0) & (df['owner_clean'].str.count(',') > 0), 'category'] = class_ids['priv']
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # Unbekannt
    # All remaining owner names that have less than 10 characters are likely not identifiable
    df.loc[(df['category'] == 0) & (df['owner_clean'].str.len() < 10), 'category'] = class_ids['unbek']
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # All the rest will go into unbekannt
    df.loc[(df['category'] == 0), 'category'] = class_ids['unbek']
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())

    # Mixed forms of 'Kapitalgesellschaften'
    # gmbh & co kg
    df.loc[(df['category'] == class_ids['mixed']) & (df['owner_clean'].str.count('mbh') > 0) &
           (df['owner_clean'].str.count('co') > 0) & (df['owner_clean'].str.count('kg') > 0),
            'category'] = class_ids['gmco']
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # The following command sets all mixed cases that do not match to zero (important to keep in mind for next steps)
    df.loc[df['category'] == class_ids['mixed'], 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['mbh & co kg', 'mbh & cokg', 'mbh & coko', 'mbh & co ko', 'gmbh u cokg',
                                             'gmbh + co kg', 'mbh & co', 'haftung & co', 'gmbh und co', 'gmbh &co'],
        return_code=class_ids['gmco'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # ag & co kg
    df.loc[(df['category'] == 0) & (df['owner_clean'].str.count('ag') > 0) &
           (df['owner_clean'].str.count('co') > 0) & (df['owner_clean'].str.count('kg') > 0),
            'category'] = class_ids['agco']
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['ag & co kg', 'ag & go kg'],
        return_code=class_ids['agco'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # ug & co kg
    df.loc[(df['category'] == 0) & (df['owner_clean'].str.count('ug') > 0) &
            (df['owner_clean'].str.count('co') > 0) & (df['owner_clean'].str.count('kg') > 0),
            'category'] = class_ids['ugco']
    df.loc[df['category'] == 0, 'category'] = df['owner_clean'].apply(
        helper_functions.check_occ_of_words_v2, search_terms=['ug haftungsbeschraenkt & co kg', 'ughaftungsbeschraenkt & co kg', 'ug & co kg',
                                             'ug haftungebeschraenkt & cokg', 'ug & cokg'],
        return_code=class_ids['ugco'])
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())
    # other mixes
    df.loc[(df['category'] == 0), 'category'] = class_ids['andco']
    print("\tNo. classified owners:", len(df.loc[df['category'] != 0].copy()), df['category'].unique())

    ## Manual assignments
    df.loc[df['owner_clean'].str.count('berliner stadtgueter') > 0, 'category'] = class_ids['land']

    ## Write output
    print(f"\tWrite out to {owners_prelim_classified_pth}")
    df.to_csv(owners_prelim_classified_pth, sep=';', index=False)


def calculate_statistics(owners_prelim_classified_pth, class_id_json_pth, stats_out_pth):
    """
    Calculates some basic statistics on the classes.
    :param owners_prelim_classified_pth: Path to stretched owner dataframe.
    :param class_id_json_pth: Path to json, holding a dictionary with preliminary classes and their IDs.
    :param stats_out_pth: Output path to statistics.
    :return:
    """
    print("Calculate statistics.")
    print("\tRead data")
    df = pd.read_csv(owners_prelim_classified_pth, sep=';')

    with open(class_id_json_pth) as json_file:
        class_ids = json.load(json_file)

    ## Descriptive statistics
    df_owners = df.drop_duplicates(subset='owner_clean')
    df_count = df_owners[['owner_clean', 'category']].groupby('category').count().reset_index()
    df['p_count'] = 1
    df_pcount = df[['p_count', 'category']].groupby('category').sum().reset_index()
    df_area = df[['area', 'category']].groupby('category').sum().reset_index()
    df_area['area'] = df_area['area']/10000

    df_stats = pd.DataFrame([[k, v] for k, v in class_ids.items()], columns=['class_name', 'class_id'])
    df_stats = pd.merge(df_stats, df_count, left_on='class_id', right_on='category', how='left').drop(columns='category')
    df_stats = pd.merge(df_stats, df_pcount, left_on='class_id', right_on='category', how='left').drop(columns='category')
    df_stats = pd.merge(df_stats, df_area, left_on='class_id', right_on='category', how='left').drop(columns='category')
    df_stats.columns = ['class_name', 'class_id', 'number of owners', 'number of parcels', 'area [ha]']

    print("\tWrite out.")
    df_stats.to_csv(stats_out_pth, sep=';', index=False)


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    classify_owner_names_into_classes(
        owner_stretched_pth=OWNERS_STRETCHED_PTH,
        class_id_json_pth=CLASS_IDS_PTH,
        owners_prelim_classified_pth=OWNERS_PRELIM_CLASSIFIED_PTH
    )
    calculate_statistics(
        owners_prelim_classified_pth=OWNERS_PRELIM_CLASSIFIED_PTH,
        class_id_json_pth=CLASS_IDS_PTH,
        stats_out_pth=OUT_STATS
    )

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)

if __name__ == '__main__':
    main()
