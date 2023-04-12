# Clemens Jänicke
# github Repo: https://github.com/clejae

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import os
import time
import pandas as pd
import json
import datetime

## Project library
import helper_functions
# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input
OWNERS_PRELIM_CLASSIFIED_PTH = r"03_owner_name_classification\03_owners_stretched_preliminary_classication.csv"
CLASSIFIER_PTH = r'00_data\classification\class_ids_classifier.json'

## Output (the output folder will be created automatically!)
OWNERS_CLASSIFIED_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified.csv"

PRIVATE_OWNERS_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified_private-persons.csv"
PRIVATE_OWNERS_GROUPS_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified_private-persons-groups.csv"
PRIVATE_OWNERS_TEMP = r"04_owner_class_reclassification\04_owners_stretched_classified_private-persons-with-multiple-addresses.csv"

COMPANY_OWNERS_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified_private-companies.csv"
COMPANY_OWNERS_TEMP = r"04_owner_class_reclassification\04_owners_stretched_classified_private-companies-with-multiple-addresses.csv"

NONPROF_OWNERS_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified_non-profit-etc.csv"
NONPROF_OWNERS_TEMP = r"04_owner_class_reclassification\04_owners_stretched_classified_non-profit-etc-with-multiple-addresses.csv"

RELIGIOUS_OWNERS_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified_religious.csv"
RELIGIOUS_OWNERS_TEMP = r"04_owner_class_reclassification\04_owners_stretched_classified_religious-with-multiple-addresses.csv"

PUBLIC_OWNERS_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified_public.csv"
PUBLIC_OWNERS_TEMP = r"04_owner_class_reclassification\04_owners_stretched_classified_public-with-multiple-addresses.csv"

STATS_OUT_PTH = r'04_owner_class_reclassification\04_owners_stretched_classified_stats.xlsx'

## ToDo: for companies, non-profit, public and church: decide on an address if there are two or more addresses
## like its done for private people
# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#


def create_overall_stats(df, level_col, area_col, fid_col, owners_col):
    out_df = []
    uni_levels = df[level_col].unique()
    for lev in uni_levels:
        sub = df[df[level_col] == lev].copy()
        num_owners = len(sub[owners_col].unique())
        mean_area = round(sub[[area_col, owners_col]].groupby(owners_col).sum().mean().iloc[0] / 1000, 1)
        mean_num = sub[[fid_col, owners_col]].groupby(owners_col).count().mean().iloc[0]
        tot_area = round(sub[area_col].sum() / 1000, 1)
        out_df.append([lev, num_owners, mean_area, mean_num, tot_area])

    out_df = pd.DataFrame(out_df)
    out_df.columns = ['Category', 'Number of owners', 'Mean area per owner [ha]',
                      'Mean number of land parcels per owner', 'Total [ha]']
    out_df = out_df.sort_values(by=['Category'])
    return out_df


def overall_stats_wrapper(df, stats_out_pth, area_col, fid_col, owners_col):
    df_stats_l1 = create_overall_stats(df, 'level1', area_col, fid_col, owners_col)
    df_stats_l2 = create_overall_stats(df, 'level2', area_col, fid_col, owners_col)
    df_stats_l3 = create_overall_stats(df, 'level3', area_col, fid_col, owners_col)

    num_owners = len(df[owners_col].unique())
    mean_area = round(df[[area_col, owners_col]].groupby(owners_col).sum().mean().iloc[0] / 1000, 1)
    mean_num = df[[fid_col, owners_col]].groupby(owners_col).count().mean().iloc[0]
    tot_area = round(df[area_col].sum() / 1000, 1)
    df_stats_gen = pd.DataFrame([['Category', 'Number of owners', 'Mean area per owner [ha]',
                                  'Mean number of land parcels per owner', 'Total [ha]'],
                                 ['Overall', num_owners, mean_area, mean_num, tot_area]])
    new_header = df_stats_gen.iloc[0]
    df_stats_gen = df_stats_gen[1:]
    df_stats_gen.columns = new_header

    ## Export all category combinations and their respective owner strings
    writer = pd.ExcelWriter(stats_out_pth)
    for c, df in enumerate([df_stats_gen, df_stats_l1, df_stats_l2, df_stats_l3]):
        name = ['Overall', 'level1', 'level2', 'level3'][c]
        df.to_excel(writer, sheet_name=name, index=False)
    writer.save()


def reclassify_owner_classes(owners_prelim_classified_pth, classifier_json_pth, owners_classified_pth):
    """
    Classifies preliminary classification into a better classification scheme as provided in the classifier json.
    :param owners_prelim_classified_pth: Path to owner dataframe with preliminary classification.
    :param classifier_json_pth: Path to json with dictionary of classification scheme
    :param owners_classified_pth: Output path for dataframe with better classification.
    :return:
    """
    print("Reclassify owner classes.")
    print("\tRead data.")
    df = pd.read_csv(owners_prelim_classified_pth, sep=';')
    print('\tReclassify')
    with open(classifier_json_pth) as json_file:
        classifier = json.load(json_file)

    ## CLASSIFY LEVEL 1 TO 3
    ## Get unique class IDs from df
    uni_class_ids = list(df['category'].astype('str').unique())
    for class_id in uni_class_ids:
        df.loc[df['category'].astype('str') == class_id, 'level3'] = classifier[class_id][0]
        df.loc[df['category'].astype('str') == class_id, 'level2'] = classifier[class_id][1]
        df.loc[df['category'].astype('str') == class_id, 'level1'] = classifier[class_id][2]

    ## Get statistics for original df
    # AREA = 'AMTLFLSFL'
    # FID = 'OGC_FID'
    # OWNERS = 'owner_clean'
    # STATS_OUT_PTH = r'tables\ALKIS\results\01_overall_stats_all_parcels_stretched_owners_before_aggregation2_AMTLFLSFL.xlsx'
    # overall_stats_wrapper(df, STATS_OUT_PTH, AREA, FID, OWNERS)

    ## remove fill words and characters
    df['owner_merge'] = df['owner_clean'].str.replace('mit sitz in ', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace('sitz in ', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace(' mit sitz ', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace(' sitz ', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace(' in ', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace(' ', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace('-', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace(',', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace('&', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace('+', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace('`', '', regex=False)
    df['owner_merge'] = df['owner_merge'].str.replace('\n', '', regex=False)

    print(f"\tWrite out to {owners_classified_pth}")
    df.to_csv(owners_classified_pth, sep=";", index=False)


def choose_address_based_on_occurences(df, owner_names, out_pth):
    """
    For owners with multiple addresses, choose the address with the most occurences. In case of same occurences choose
    simply the first in the list.
    :param df:
    :param owner_names:
    :param out_pth:
    :return:
    """

    file = open(out_pth, "w+", encoding='ISO-8859-1')
    file.write("owner_merge;addresses;clean_address\n")

    for owner in owner_names:
        sub = df.loc[df['owner_merge'] == owner].copy()
        addresses = sub['clean_address'].tolist()
        uni_addresses = []
        for address in addresses:
            lst = address.split('_')
            for address in lst:
                if address != 'unbekannt':
                    uni_addresses.append(address)
        if uni_addresses:
            mf_address = helper_functions.most_frequent(uni_addresses)
        else:
            mf_address = 'unbekannt'
            uni_addresses = addresses
        file.write(f"{owner};{'_'.join(list(set(uni_addresses)))};{mf_address}\n")

    file.close()


def clean_company_identifiers(owners_classified_pth, company_owners_address_pth, company_owners_cleaned):
    """
    Clean the identifier of companies that occur multiple times in the dataframe, but have different spellings. With fuzzy
    matching one name is chosen. For this identifier, we assign the address with the most occurences.
    :param owners_classified_pth: Path to dataframe with classified owners.
    :param company_owners_address_pth: Output path to list of unique companies and the chosen address.
    :param company_owners_cleaned: Output path to dataframe with cleaned companies.
    :return:
    """
    print('Private companies')
    print("\tRead data.")
    df = pd.read_csv(owners_classified_pth, sep=';')

    df_comp = df[df['level1'] == 2].copy()

    ## Get df with unique company names and their frequencies (comes in ascending order)
    df_comp_uni = df_comp.groupby(['owner_merge']).size().reset_index(name='Freq').copy()

    ## Fuzzy match company names with forward looking moving window
    helper_functions.fuzzy_match_token_set_ratio(df_comp_uni, 'owner_merge', 'owner_new', 'Freq', 95, 20)

    ## Create dictionary with old names and new names
    own_dict = {}
    # for i in range(len(df_comp_uni)):
    for row in df_comp_uni.itertuples():
        owner_merge = row.owner_merge
        owner_new = row.owner_new
        own_dict[owner_merge] = owner_new

    ## Assign new names based on dictionary
    uni_owners = list(df_comp_uni['owner_merge'].astype('str').unique())
    for uni_own in uni_owners:
        df_comp.loc[df_comp['owner_merge'].astype('str') == uni_own, 'owner_merge'] = own_dict[uni_own]

    ## Decide on one address based on occurences
    print("\tNo. entries:", len(df_comp))
    owner_names = df_comp['owner_merge'].unique().tolist()

    choose_address_based_on_occurences(df=df_comp, owner_names=owner_names, out_pth=company_owners_address_pth)
    df_assigned = pd.read_csv(company_owners_address_pth, sep=';', encoding='ISO-8859-1')
    cols_left = ['OGC_FID', 'EIGENTUEME', 'AMTLFLSFL', 'area', 'owner_names', 'own_num', 'owner_clean', 'category',
                 'owner_merge', 'level3', 'level2', 'level1']

    # t = df_assigned.loc[df_assigned['owner_merge'].isin(owner_names)].copy()
    df_out = pd.merge(df_comp[cols_left], df_assigned, how='left', on='owner_merge')
    print("\tNo. entries:", len(df_out))

    print(f"\tWrite out to {company_owners_cleaned}")
    df_out.to_csv(company_owners_cleaned, sep=';', index=False)


def clean_non_profit_identifiers(owners_classified_pth, nonprof_owners_address_pth, nonprof_owners_cleaned):
    """
       Clean the identifier of non profit organization etc. that occur multiple times in the dataframe,
       but have different spellings. With fuzzy matching one name is chosen.
       For this name, we assign the address with the most occurences.
       :param owners_classified_pth: Path to dataframe with classified owners.
       :param nonprof_owners_address_pth: Output path to list of unique companies and the chosen address.
       :param nonprof_owners_cleaned: Output path to dataframe with cleaned companies.
       :return:
       """
    df = pd.read_csv(owners_classified_pth, sep=';')
    print('Non-profit')
    print("\tRead data.")
    df_nprof = df[df['level1'] == 3].copy()

    ## Aggregate manually import entities
    df_nprof.loc[df_nprof['owner_merge'].str.count('arbeitersamariter') > 0, 'owner_merge'] = 'arbeitersamariterbund'
    df_nprof.loc[df_nprof['owner_merge'].str.count('arbeiterwohlfahrt') > 0, 'owner_merge'] = 'arbeiterwohlfahrt'
    df_nprof.loc[df_nprof['owner_merge'].str.count('deutschesroteskreuz') > 0, 'owner_merge'] = 'deutschesroteskreuz'
    df_nprof.loc[df_nprof['owner_merge'].str.count('drk') > 0, 'owner_merge'] = 'deutschesroteskreuz'
    df_nprof.loc[
        df_nprof['owner_merge'].str.count('grosstrappenschutz') > 0, 'owner_merge'] = 'fordervereingrosstrappenschutz'
    df_nprof.loc[
        df_nprof['owner_merge'].str.count('heinzsielmann') > 0, 'owner_merge'] = 'heinzsielmannstiftunggutherbigshagen'
    df_nprof.loc[
        df_nprof['owner_merge'].str.count('landesanglerverbandbr') > 0, 'owner_merge'] = 'landesanglerverbandbrandenburgev'
    df_nprof.loc[df_nprof['owner_merge'].str.count(
        'vereinnuth') > 0, 'owner_merge'] = 'landschaftsfoerdervereinnuthenieplitzniederungev'
    df_nprof.loc[df_nprof['owner_merge'].str.count(
        'vereinoberesrhin') > 0, 'owner_merge'] = 'landschaftsfoerdervereinoberesrhinluchev'
    df_nprof.loc[df_nprof['owner_merge'].str.count('nabu') > 0, 'owner_merge'] = 'nabustiftungnationalesnaturerbe'
    df_nprof.loc[
        df_nprof['owner_merge'].str.count('naturschutzbund') > 0, 'owner_merge'] = 'nabustiftungnationalesnaturerbe'
    df_nprof.loc[df_nprof['owner_merge'].str.count('naturschutzfondsbr') > 0, 'owner_merge'] = 'naturschutzfondsbrandenburg'
    df_nprof.loc[df_nprof['owner_merge'].str.count(
        'naturparkschlaubetal') > 0, 'owner_merge'] = 'landschaftspflegeverbadnnaturparkschlaubetal'
    df_nprof.loc[df_nprof['owner_merge'].str.count(
        'edithmaryon') > 0, 'owner_merge'] = 'stiftungedithmaryonzurfoerderungsozialerwohnundarbeitsstaetten'
    df_nprof.loc[df_nprof['owner_merge'].str.count('stiftneuzelle') > 0, 'owner_merge'] = 'stiftungstiftneuzelle'
    df_nprof.loc[df_nprof['owner_merge'].str.count('wwf') > 0, 'owner_merge'] = 'umweltstiftungwwfdeutschland'
    df_nprof.loc[df_nprof['owner_merge'].str.count(
        'freundedesdeutschpolnisch') > 0, 'owner_merge'] = 'vereinfreundedesdeutschpolnischeneuropanationalparksunteresodertalev'
    df_nprof.loc[df_nprof['owner_merge'].str.count(
        'vogelschutzkomiteeev') > 0, 'owner_merge'] = 'vskvogelschutzkomiteeevgesellschaftzur'
    df_nprof.loc[df_nprof['owner_merge'].str.count(
        'zoologischegesellschaftfra') > 0, 'owner_merge'] = 'zoologischegesellschaftfrankfurtvon1858ev'

    ## Get df with unique nonprof names and their frequencies (comes in ascending order)
    df_nprof_uni = df_nprof.groupby(['owner_merge']).size().reset_index(name='Freq').copy()

    ## Fuzzy matching
    helper_functions.fuzzy_match_token_set_ratio(df_nprof_uni, 'owner_merge', 'owner_new', 'Freq', 97, 20)

    ## Create dictionary with old names and new names
    own_dict = {}
    for row in df_nprof_uni.itertuples():
        owner_merge = row.owner_merge
        owner_new = row.owner_new
        own_dict[owner_merge] = owner_new

    ## Assign new names based on dictionary
    uni_owners = list(df_nprof_uni['owner_merge'].astype('str').unique())
    for uni_own in uni_owners:
        df_nprof.loc[df_nprof['owner_merge'].astype('str') == uni_own, 'owner_merge'] = own_dict[uni_own]

    ## Decide on one address based on occurences
    print("\tNo. entries:", len(df_nprof))
    owner_names = df_nprof['owner_merge'].unique().tolist()

    choose_address_based_on_occurences(df=df_nprof, owner_names=owner_names, out_pth=nonprof_owners_address_pth)
    df_assigned = pd.read_csv(nonprof_owners_address_pth, sep=';', encoding='ISO-8859-1')
    cols_left = ['OGC_FID', 'EIGENTUEME', 'AMTLFLSFL', 'area', 'owner_names', 'own_num', 'owner_clean', 'category',
                 'owner_merge', 'level3', 'level2', 'level1']

    # t = df_assigned.loc[df_assigned['owner_merge'].isin(owner_names)].copy()
    df_out = pd.merge(df_nprof[cols_left], df_assigned, how='left', on='owner_merge')
    print("\tNo. entries:", len(df_out))
    print(f"\tWrite out to {nonprof_owners_cleaned}")
    df_out.to_csv(nonprof_owners_cleaned, sep=';', index=False)


def clean_religious_identifiers(owners_classified_pth, relig_owners_address_pth, relig_owners_cleaned):
    """
    Clean the identifier of non religious entities that occur multiple times in the dataframe,
    but have different spellings. With fuzzy  matching one name is chosen.
    For this name, we assign the address with the most occurences.
    :param owners_classified_pth:
    :param relig_owners_address_pth:
    :param relig_owners_cleaned:
    :return:
    """
    print('Religious')
    print("\tRead data.")
    df = pd.read_csv(owners_classified_pth, sep=';')

    df_rel = df[df['level1'] == 4].copy()

    ## Aggregate manually import entities
    df_rel.loc[
        df_rel['owner_merge'].str.count('conferenceonjewish') > 0, 'owner_merge'] = 'conferenceonjewishmaterialclaimsinc'

    ## Decide on one address based on occurences
    print("\tNo. entries:", len(df_rel))
    owner_names = df_rel['owner_merge'].unique().tolist()

    choose_address_based_on_occurences(df=df_rel, owner_names=owner_names, out_pth=relig_owners_address_pth)
    df_assigned = pd.read_csv(relig_owners_address_pth, sep=';', encoding='ISO-8859-1')
    cols_left = ['OGC_FID', 'EIGENTUEME', 'AMTLFLSFL', 'area', 'owner_names', 'own_num', 'owner_clean', 'category',
                 'owner_merge', 'level3', 'level2', 'level1']

    # t = df_assigned.loc[df_assigned['owner_merge'].isin(owner_names)].copy()
    df_out = pd.merge(df_rel[cols_left], df_assigned, how='left', on='owner_merge')

    print("\tNo. entries:", len(df_out))
    print(f"\tWrite out to {relig_owners_cleaned}")
    df_out.to_csv(relig_owners_cleaned, sep=';', index=False)

def clean_public_identifiers(owners_classified_pth, public_owners_address_pth, public_owners_cleaned):
    """
    Clean the identifier of non public entities that occur multiple times in the dataframe,
    but have different spellings.
    :param owners_classified_pth:
    :param public_owners_address_pth:
    :param public_owners_cleaned:
    :return:
    """
    print('Public')
    print("\tRead data.")
    df = pd.read_csv(owners_classified_pth, sep=';')

    df_pub = df[df['level1'] == 5].copy()
    df_pub_uni = df_pub.groupby(['owner_merge']).size().reset_index(name='Freq').copy()

    ## Aggregate manually import entities
    df_pub.loc[df_pub['owner_merge'].str.count('bvvg') > 0, 'owner_merge'] = 'bvvg'
    df_pub.loc[df_pub['owner_merge'].str.count('berlinerstadtgueter') > 0, 'owner_merge'] = 'berlinerstadtguetergmbh'
    df_pub.loc[df_pub['owner_merge'].str.count('bundesstrass') > 0, 'owner_merge'] = 'bundesstrassenverwaltung'
    df_pub.loc[df_pub['owner_merge'].str.count('bundeswasserstrass') > 0, 'owner_merge'] = 'bundeswasserstrassenverwaltung'
    df_pub.loc[df_pub['owner_merge'].str.count(
        'bundesanstaltfuerimmobilien') > 0, 'owner_merge'] = 'bundesanstaltfuerimmobilienaufgaben'
    df_pub.loc[df_pub['owner_merge'].str.count(
        'bundesanstaltfuervereinigungs') > 0, 'owner_merge'] = 'bundesanstaltfuervereinigungsbedingtesonderaufgaben'
    df_pub.loc[df_pub['owner_merge'].str.count('bundesfinanz') > 0, 'owner_merge'] = 'brdbundesfinanzverwaltung'
    df_pub.loc[df_pub['owner_merge'].str.count('landberli') > 0, 'owner_merge'] = 'landberlin'
    df_pub.loc[df_pub['owner_merge'].str.count('brandenburggrund') > 0, 'owner_merge'] = 'landbrandenburggrundstuecksfond'
    # df_pub.loc[df_pub['owner_merge'].str.count('landbra') > 0, 'owner_merge'] = 'landbrandenburg'

    df_pub_uni = df_pub.groupby(['owner_merge']).size().reset_index(name='Freq').copy()
    ## Fuzzy matching
    helper_functions.fuzzy_match_token_set_ratio(df_pub_uni, 'owner_merge', 'owner_new', 'Freq', 97, 20)

    ## Create dictionary with old names and new names
    own_dict = {}
    for row in df_pub_uni.itertuples():
        owner_merge = row.owner_merge
        owner_new = row.owner_new
        own_dict[owner_merge] = owner_new

    ## Assign new names based on dictionary
    uni_owners = list(df_pub_uni['owner_merge'].astype('str').unique())
    for uni_own in uni_owners:
        # print(uni_own)
        df_pub.loc[df_pub['owner_merge'].astype('str') == uni_own, 'owner_merge'] = own_dict[uni_own]

    ## Decide on one address based on occurences
    print("\tNo. entries:", len(df_pub))
    owner_names = df_pub['owner_merge'].unique().tolist()

    choose_address_based_on_occurences(df=df_pub, owner_names=owner_names, out_pth=public_owners_address_pth)
    df_assigned = pd.read_csv(public_owners_address_pth, sep=';', encoding='ISO-8859-1')
    cols_left = ['OGC_FID', 'EIGENTUEME', 'AMTLFLSFL', 'area', 'owner_names', 'own_num', 'owner_clean', 'category',
                 'owner_merge', 'level3', 'level2', 'level1']

    # t = df_assigned.loc[df_assigned['owner_merge'].isin(owner_names)].copy()
    df_out = pd.merge(df_pub[cols_left], df_assigned, how='left', on='owner_merge')
    print("\tNo. entries:", len(df_out))

    # df_out.loc[df_out['clean_address'].isnan(), 'clean_address'] = 'seeburger chaussee 2, 14476 potsdam'
    # df_out.loc[df_out['clean_address'].isnan(), 'owner_clean'] = 'land brandenburg landesnaturschutzflaechenverwaltung'
    # df_out.loc[df_out['clean_address'].isnan(), 'owner_merge'] = 'landbrandenburglandesnaturschutzflaechenverwaltung'

    print(f"\tWrite out to {public_owners_cleaned}")
    df_out.to_csv(public_owners_cleaned, sep=';', index=False)


def clean_private_owners_identifiers(owners_classified_pth, private_owners_address_pth, private_owner_groups_cleaned, private_owners_cleaned):
    """
    Clean the identifier of private people and groups of private people (e.g. partnerships etc.).
    :param owners_classified_pth:
    :param private_owners_address_pth:
    :param private_owner_groups_cleaned:
    :param private_owners_cleaned:
    :return:
    """
    print('Private')
    print("\tRead data.")
    df = pd.read_csv(owners_classified_pth, sep=';')
    ## CLEAN ALL PRIVATE OWNERS
    df_priv = df[df['level1'] == 1].copy()
    del df

    ## Split between single owners and groups of owners
    df_done = df_priv.loc[df_priv['category'] != 1].copy()
    df_work = df_priv.loc[df_priv['category'] == 1].copy()

    df_work_uni = df_work.drop_duplicates(subset=['owner_merge', 'clean_address'])
    df_name_count = df_work_uni[['owner_merge', 'clean_address']].groupby(by='owner_merge').count().reset_index()
    df_name_count.rename(columns={'clean_address': 'address_count'}, inplace=True)
    df_mult = df_name_count.loc[df_name_count['address_count'] > 1].copy()
    df_single = df_name_count.loc[df_name_count['address_count'] <= 1].copy()

    ## in df single, there are owners with multiple addresses that always occur with these addresses
    ## separate those from the others and assign first address to them
    df_temp = df_work_uni.loc[df_work_uni['owner_merge'].isin(df_single['owner_merge'])].copy()

    df_single_ambiguos = df_temp.loc[df_temp["clean_address"].str.count('_') > 0].copy()
    df_single = df_temp.loc[df_temp["clean_address"].str.count('_') == 0].copy()

    df_assigned3 = df_work_uni.loc[df_work_uni['owner_merge'].isin(df_single_ambiguos['owner_merge'])].copy()
    df_assigned3 = df_assigned3[["owner_merge", "addresses"]]
    df_assigned3["clean_address"] = df_assigned3.apply(lambda row: row.addresses.split('_')[0], axis=1)

    ##TODO remove clean addresses with '_'
    def most_frequent(in_list):
        return max(set(in_list), key=in_list.count)

    file = open(private_owners_address_pth, "w+")
    file.write("owner_merge;addresses;clean_address\n")

    for owner in df_mult['owner_merge']:
        sub = df_work.loc[df_work['owner_merge'] == owner].copy()
        addresses = sub['clean_address'].tolist()
        uni_addresses = []
        for address in addresses:
            lst = address.split('_')
            for address in lst:
                if address != 'unbekannt':
                    uni_addresses.append(address)
        mf_address = most_frequent(uni_addresses)
        file.write(f"{owner};{'_'.join(list(set(uni_addresses)))};{mf_address}\n")
        # out_dict["owner_merge"].append(owner)
        # out_dict["addresses"].append('_'.join(list(set(uni_addresses))))
        # out_dict["clean_address"] = mf_address
    file.close()

    df_assigned1 = df_work_uni.loc[df_work_uni['owner_merge'].isin(df_single['owner_merge'])].copy()
    df_assigned1 = df_assigned1[["owner_merge", "addresses", "clean_address"]]
    df_assigned2 = pd.read_csv(private_owners_address_pth, sep=';', encoding='ISO-8859-1')
    del df_work_uni, df_mult, df_single, df_name_count, df_temp

    df_assigned = pd.concat([df_assigned1, df_assigned2, df_assigned3], axis=0)
    del df_assigned1, df_assigned2, df_assigned3

    cols_left = ['OGC_FID', 'EIGENTUEME', 'AMTLFLSFL', 'area', 'owner_names', 'own_num', 'owner_clean', 'category',
                 'owner_merge', 'level3', 'level2', 'level1']
    df_work = pd.merge(df_work[cols_left], df_assigned, how='left', on='owner_merge')

    ## Reorder columns to original order
    original_order = list(df_priv.columns)
    if "Unnamed: 0" in original_order:
        original_order.remove("Unnamed: 0")
    df_work = df_work[original_order]

    df_work.rename(columns={'owner_merge': 'owner_merge0'}, inplace=True)
    df_work['owner_clean'] = df_work['owner_clean'].str.strip(',')
    df_work['owner_clean'] = df_work['owner_clean'].str.strip(' ')

    df_work['fam_name'] = df_work.apply(lambda row: row.owner_clean.split(',')[0], axis=1)
    df_work['owner_merge'] = df_work['fam_name'] + df_work['clean_address']
    df_work['owner_merge'] = df_work['owner_merge'].str.replace(' ', '', regex=False)
    df_work['owner_merge'] = df_work['owner_merge'].str.replace('-', '', regex=False)
    df_work['owner_merge'] = df_work['owner_merge'].str.replace(',', '', regex=False)
    df_work['owner_merge'] = df_work['owner_merge'].str.replace('&', '', regex=False)
    df_work['owner_merge'] = df_work['owner_merge'].str.replace('+', '', regex=False)
    df_work['owner_merge'] = df_work['owner_merge'].str.replace('`', '', regex=False)

    ## Filter by birthdate
    ############ uncomment if you want to start here ############
    ### df = pd.read_csv(OUT_PTH, sep=';')
    # df_priv = df[df['level1'] == 1].copy()
    #############################################################

    df_work['birthdate'] = df_work['owner_clean'].apply(helper_functions.get_birthdate)
    df_work['birthdate'] = df_work['birthdate'].str.replace(' ', '')
    df_work['birthdate'] = df_work['birthdate'].replace({'\*', ''}, regex=True)
    df_work['birthdate'] = df_work['birthdate'].str.replace('-', '')

    df_work['birthdate'] = pd.to_datetime(df_work['birthdate'], format='%Y%m%d', errors='coerce')
    date_accq = datetime.datetime(year=2020, month=11, day=15)
    # df['age'] = (now - df['dob']).astype('<m8[Y]')    # 3
    df_work['age_accq'] = (date_accq - df_work['birthdate']).astype('<m8[Y]')

    df_work.loc[df_work['age_accq'] > 101, 'owner_merge'] = 'unbekannt'
    df_work.loc[df_work['age_accq'] > 101, 'owner_clean'] = 'unbekannt'
    df_work.loc[df_work['age_accq'] > 101, 'clean_address'] = 'unbekannt'
    df_work.loc[df_work['age_accq'] > 101, 'category'] = 28
    df_work.loc[df_work['age_accq'] > 101, 'level3'] = '5_2_5'
    df_work.loc[df_work['age_accq'] > 101, 'level2'] = '5_2'
    df_work.loc[df_work['age_accq'] > 101, 'level1'] = 5

    df_work.loc[df_work['age_accq'].isna(), 'age_accq'] = -9999
    print(f"\tWrite out to {private_owner_groups_cleaned} and {private_owners_cleaned}")
    df_work.to_csv(private_owners_cleaned, sep=';', index=False)
    df_done.to_csv(private_owner_groups_cleaned, sep=';', index=False)


def merge_all_class_dfs(owners_prelim_classified_pth, owner_classified_pth, company_owners_cleaned, private_owners_cleaned,
                        private_owner_groups_cleaned, nonprof_owners_cleaned, religious_owners_cleaned,
                        public_owners_cleaned):

    """
    Combine all dataframes where the cleaned identifiers are cleaned and overwrite dataframe with classification.
    :param owners_prelim_classified_pth:
    :param company_owners_cleaned:
    :param private_owners_cleaned:
    :param private_owner_groups_cleaned:
    :param nonprof_owners_cleaned:
    :param religious_owners_cleaned:
    :param public_owners_cleaned:
    :return:
    """

    print("\tRead data.")
    df = pd.read_csv(owners_prelim_classified_pth, sep=';')
    df_priv = pd.read_csv(private_owners_cleaned, sep=';')
    df_priv_groups = pd.read_csv(private_owner_groups_cleaned, sep=';')
    df_comp = pd.read_csv(company_owners_cleaned, sep=';')
    df_nprof = pd.read_csv(nonprof_owners_cleaned, sep=';')
    df_rel = pd.read_csv(religious_owners_cleaned, sep=';')
    df_pub = pd.read_csv(public_owners_cleaned, sep=';')

    print("\tNo. original entries:", len(df))
    ## Bring all columns of all data frames into the same order
    columns = df_comp.columns
    df_priv = df_priv[columns]
    df_priv_groups = df_priv_groups[columns]
    df_nprof = df_nprof[columns]
    df_rel = df_rel[columns]
    df_pub = df_pub[columns]

    ## Merge all the separate dfs
    df_out = pd.concat([df_priv, df_priv_groups, df_comp, df_nprof, df_rel, df_pub])
    df_out.sort_values(by='OGC_FID', inplace=True)
    print("\tNo. final entries:", len(df_out))

    print(f"\tOverwrite {owner_classified_pth}")
    df_out.to_csv(owner_classified_pth, sep=";", index=False)


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    reclassify_owner_classes(
        owners_prelim_classified_pth=OWNERS_PRELIM_CLASSIFIED_PTH,
        classifier_json_pth=CLASSIFIER_PTH,
        owners_classified_pth=OWNERS_CLASSIFIED_PTH
    )
    clean_company_identifiers(
        owners_classified_pth=OWNERS_CLASSIFIED_PTH,
        company_owners_address_pth=COMPANY_OWNERS_TEMP,
        company_owners_cleaned=COMPANY_OWNERS_PTH
    )
    clean_non_profit_identifiers(
        owners_classified_pth=OWNERS_CLASSIFIED_PTH,
        nonprof_owners_address_pth=NONPROF_OWNERS_TEMP,
        nonprof_owners_cleaned=NONPROF_OWNERS_PTH
    )
    clean_religious_identifiers(
        owners_classified_pth=OWNERS_CLASSIFIED_PTH,
        relig_owners_address_pth=RELIGIOUS_OWNERS_TEMP,
        relig_owners_cleaned=RELIGIOUS_OWNERS_PTH
    )
    clean_public_identifiers(
        owners_classified_pth=OWNERS_CLASSIFIED_PTH,
        public_owners_address_pth=PUBLIC_OWNERS_TEMP,
        public_owners_cleaned=PUBLIC_OWNERS_PTH
    )
    clean_private_owners_identifiers(
        owners_classified_pth=OWNERS_CLASSIFIED_PTH,
        private_owners_address_pth=PRIVATE_OWNERS_TEMP,
        private_owner_groups_cleaned=PRIVATE_OWNERS_GROUPS_PTH,
        private_owners_cleaned=PRIVATE_OWNERS_PTH
    )
    merge_all_class_dfs(
        owners_prelim_classified_pth=OWNERS_PRELIM_CLASSIFIED_PTH,
        owner_classified_pth=OWNERS_CLASSIFIED_PTH,
        company_owners_cleaned=COMPANY_OWNERS_PTH,
        private_owners_cleaned=PRIVATE_OWNERS_PTH,
        private_owner_groups_cleaned=PRIVATE_OWNERS_GROUPS_PTH,
        nonprof_owners_cleaned=NONPROF_OWNERS_PTH,
        religious_owners_cleaned=RELIGIOUS_OWNERS_PTH,
        public_owners_cleaned=OWNERS_CLASSIFIED_PTH
    )

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)

if __name__ == '__main__':
    main()
