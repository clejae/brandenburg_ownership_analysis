# Clemens Jänicke
# github Repo: https://github.com/clejae

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import os
import time
import geopandas as gpd
import pandas as pd
import jaro

## Project library
import helper_functions
# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input
ALKIS_PTH = r"04_owner_class_reclassification\04_owners_stretched_classified.csv"
PARCELS_PTH = r"00_data\vector\ALKIS\v_eigentuemer_bb_reduced.shp"

ADRESSES_PATH = r"05_georeference_addresses\05_owners_stretched_addresses_not_geolocated.csv"  #r"05_calculate_geolocation\owner_addresses_missing.csv"
LOC_OSM_PATH = r"05_georeference_addresses\GER_places_post_code_lev4.shp"
ALREADY_GEOCODED_PATH = r"05_georeference_addresses\05_geocoded_addresses_preliminary.csv"

## Output (the output folder will be created automatically!)
OUTPATH_NOMINATIM_BACKUP = r"05_georeference_addresses\05_geocoded_addresses_after_nominatim_geocoding_backup.csv"
OUTPATH_GEOCODED_ALL = r"05_georeference_addresses\05_geocoded_addresses_complete.csv"

OUTPATH_GEOCODED_MISSES_L1 = r"05_georeference_addresses\05_missing_addresses_after_osm_geocoding_l1.csv"
OUTPATH_GEOCODED_MATCHES_L1 = r"05_georeference_addresses\05_geocoded_addresses_after_osm_geocoding_l1.csv"
OUTPATH_GEOCODED_MISSES_NOMINATIM = r"05_georeference_addresses\05_missing_addresses_after_nominatim_geocoding.csv"
OUTPATH_GEOCODED_MATCHES_NOMINATIM = r"05_georeference_addresses\05_geocoded_addresses_after_nominatim_geocoding.csv"
OUTPATH_GEOCODED_MISSES_FUZZY = r"05_georeference_addresses\05_missing_addresses_after_fuzzy_geocoding.csv"

ADDRESSES_SHP = r"05_georeference_addresses\05_geocoded_addresses_complete_4326.shp"
ADMINISTRATIVE_SHP = r"05_georeference_addresses\GER_adm_bound_lev4_4326.shp"
GEOCODED_ADDRESSES = r"05_georeference_addresses\05_geocoded_addresses_complete_with_administrative_names_4326.csv"

OUTPATH_MISSING = r"05_calculate_geolocation\05_owners_stretched_addresses_not_geolocated.csv"
OUTPATH_FINAL = r"05_calculate_geolocation\05_owners_stretched_addresses_geolocated.csv"
OUTPATH_DISTANCES = r"05_calculate_geolocation\05_owners_stretched_addresses_geolocated_distances.csv"

# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#

def transfer_geometry_from_geocoded_address_via_fuzzy_matching(df_addr, df_miss, plzcol_l, plzcol_r, addrcol_l, addrcol_r, geomcol_l, thresh):
    df_succ = pd.DataFrame(columns=list(df_addr.columns))

    if plzcol_l not in df_addr.columns:
        print(f'{plzcol_l} is not in columns of df_addr!')
        return df_succ

    if plzcol_r not in df_miss.columns:
        print(f'{plzcol_r} is not in columns of df_miss!')
        return df_succ

    if addrcol_l not in df_addr.columns:
        print(f'{addrcol_l} is not in columns of df_addr!')
        return df_succ

    if addrcol_r not in df_miss.columns:
        print(f'{addrcol_r} is not in columns of df_miss!')
        return df_succ

    df_addr[plzcol_l] = df_addr[plzcol_l].astype(str)
    df_miss[plzcol_r] = df_miss[plzcol_r].astype(str)

    df_addr['merge_address'] = df_addr[addrcol_l].str.replace(' ', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('.', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace(',', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('-', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('/', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('ä', 'ae', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('ü', 'ue', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('ö', 'oe', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('ß', 'ss', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('asse', '', regex=False)

    df_miss['merge_address'] = df_miss[addrcol_r].str.replace(' ', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('.', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace(',', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('-', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('/', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('ä', 'ae', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('ü', 'ue', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('ö', 'oe', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('ß', 'ss', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('asse', '', regex=False)

    df_miss.drop_duplicates(subset=['merge_address'], inplace=True)

    out_lst = []
    count = 0
    for i, address in enumerate(df_miss['merge_address']):
        plz = df_miss[plzcol_r].iloc[i]
        df_sub = df_addr.loc[df_addr[plzcol_l] == plz].copy()
        df_sub['jaro_winkler'] = df_sub['merge_address'].apply(jaro.jaro_winkler_metric, string2=address)
        df_sub = df_sub.loc[df_sub['jaro_winkler'] >= thresh].copy()
        max_val = df_sub['jaro_winkler'].max()
        df_sub = df_sub.loc[df_sub['jaro_winkler'] == max_val].copy()

        if not df_sub.empty:
            count += 1
            geom = df_sub[geomcol_l].iloc[0]
            point_addr = df_sub[addrcol_l].iloc[0]
            print(count, i, round(max_val, 4), address, '==', df_sub[addrcol_l].iloc[0])

            out_lst.append([address, geom, "fuzzy_matching", point_addr])
            # uni_adr.loc[uni_adr['address'] == address, 'geometry'] = geom
            # uni_adr.loc[uni_adr['address'] == address, 'geocoding'] = 'fuzzy_matching'
            # uni_adr.loc[uni_adr['address'] == address, 'point_address'] = point_addr

    df_succ = pd.DataFrame(out_lst)
    df_succ.columns = [addrcol_r, "geometry", "geocoding", "point_address"]

    return df_succ


def addresses_to_coordinates_with_osm_l1_data():

    ## read clean addresses
    df = pd.read_csv(ADRESSES_PATH, sep=';') #, dtype={'post_code': str})

    ## read localities of Germany
    gdf_loc = gpd.read_file(LOC_OSM_PATH)
    df_loc = pd.DataFrame(gdf_loc)[['geometry', 'fclass', 'plz', 'name']]
    df_loc['plz'] = df_loc['plz'].fillna('')
    df_loc['plz'] = df_loc['plz'].astype('str')

    ## clean localities string
    df_loc['name'] = df_loc['name'].str.lower()
    df_loc['name'] = df_loc['name'].apply(
        helper_functions.replace_characters_in_string, char_lst=['\n'], replace_string=',')
    df_loc['name'] = df_loc['name'].apply(
        helper_functions.replace_characters_in_string, char_lst=[' - '], replace_string=' | ')
    df_loc['name'] = df_loc['name'].apply(
        helper_functions.replace_characters_in_string, char_lst=['\x9a', '\x8a', '(', ')', '/'], replace_string='')
    df_loc['name'] = df_loc['name'].apply(
        helper_functions.replace_characters_in_string, char_lst=['-'], replace_string=' ')
    df_loc['name'] = df_loc['name'].apply(helper_functions.clean_city_text, search_terms=['ot', 'bei', '|'])

    ## create unique identifiers for level 1
    df_loc['pc_city'] = df_loc['plz'] + df_loc['name']
    df_loc['pc_city'] = df_loc['pc_city'].str.replace(' ', '')

    df_loc = df_loc.drop_duplicates(['pc_city'])

    ## identify street, post_code, city in df
    df = df[df['clean_address'].notna()]
    df['street'] = df['clean_address'].apply(helper_functions.identify_street)
    df['post_code'] = df['clean_address'].apply(helper_functions.identify_plz)
    df['city'] = df['clean_address'].apply(helper_functions.identify_city)
    df['address'] = df['clean_address']

    ## prepare df of unique addresses
    uni_adr = df[['street', 'post_code', 'city', 'address']].drop_duplicates(['address'])
    uni_adr.index = range(len(uni_adr['street']))

    ## check if address is complete
    ## 0: no address, 1: only street, 2: only city, 3: street + city,
    ## 4: only pc, 5: pc+ street, 6: pc + city, 7: full address
    uni_adr['sfull'] = 1
    uni_adr.loc[uni_adr['street'] == '', 'sfull'] = 0 ## uni_adr['sfull'][uni_adr['street'] == ''] = 0

    uni_adr['cfull'] = 2
    uni_adr.loc[uni_adr['city'] == '', 'cfull'] = 0 ## uni_adr['cfull'][uni_adr['city'] == ''] = 0

    uni_adr['pfull'] = 4
    uni_adr.loc[uni_adr['post_code'] == '', 'pfull'] = 0 ## uni_adr['pfull'][uni_adr['post_code'] == ''] = 0

    uni_adr['full_address'] = uni_adr['sfull'] + uni_adr['cfull'] + uni_adr['pfull']
    uni_adr = uni_adr.drop(columns=['sfull', 'cfull', 'pfull'])
    uni_adr['plz_len'] = uni_adr['post_code'].str.len()

    ## create unique identifier in df of addresses
    uni_adr['pc_city'] = uni_adr['post_code'] + uni_adr['city']
    uni_adr['pc_city'] = uni_adr['pc_city'].str.replace(' ', '')

    ## combine addresses with unique identifier of level 1
    ## divide between matches and misses
    df_comb = pd.merge(uni_adr, df_loc[['pc_city', 'fclass', 'geometry']], how='left', on='pc_city')
    df_clear = df_comb[df_comb['geometry'].notnull()].copy()
    df_clear['geocoding'] = 'osm_level1'
    df_miss = df_comb[df_comb['geometry'].isna()].copy()

    df_clear['point_address'] = df_clear['post_code'] + ' ' + df_clear['city']
    df_clear = df_clear[['street', 'post_code', 'city', 'address', 'full_address', 'plz_len',
                         'pc_city', 'fclass', 'geometry', 'geocoding', 'point_address']]

    df_miss.to_csv(OUTPATH_GEOCODED_MISSES_L1, sep=';', index=False)
    df_clear.to_csv(OUTPATH_GEOCODED_MATCHES_L1, sep=';', index=False)


def addresses_to_coordinates_with_nominatim():
    df_miss = pd.read_csv(OUTPATH_GEOCODED_MISSES_L1, sep=';')
    ## geocode addresses that could not be matched with level 1
    ## dived between matches and misses
    df_miss['geometry'] = df_miss['address'].apply(helper_functions.address_to_coordinates_nominatim, out_pth=OUTPATH_NOMINATIM_BACKUP)

    df_clear_nomi = df_miss[df_miss['geometry'].notnull()].copy()
    df_clear_nomi['geocoding'] = 'nominatim'
    df_clear_nomi['geometry'] = df_clear_nomi['geometry'].apply(helper_functions.transform_point, in_sr=4326, out_sr=25832)

    # point_address_df = pd.read_csv(OUTPATH_NOMINATIM_BACKUP, sep=';', encoding='ISO 8859-1')
    # point_address_df.drop_duplicates(subset='in_address', inplace=True)
    # df_clear_nomi = pd.merge(df_clear_nomi, point_address_df, how='left', left_on='address', right_on='in_address')
    # df_clear_nomi.drop(columns=[' point', 'in_address'], inplace=True)
    # df_clear_nomi.rename(columns={' point_address': 'point_address'}, inplace=True)
    # df_clear_nomi = df_clear_nomi[['street', 'post_code', 'city', 'address', 'full_address', 'plz_len',
    #                                'pc_city', 'fclass', 'geometry', 'geocoding', 'point_address']]

    df_miss_nomi = df_miss[df_miss['geometry'].isnull()].copy()

    # df_clear_nomi.to_csv(OUTPATH_GEOCODED_MATCHES_NOMINATIM, sep=';', index=False)

    df_clear_nomi.drop(columns='fclass', inplace=True)
    df_clear_nomi['point_address'] = df_clear_nomi['address']

    ## open already geocoded addresses and append newly geocoded addresses to them
    df_pre = pd.read_csv(ALREADY_GEOCODED_PATH, sep=';')
    df_pre = pd.concat([df_pre, df_clear_nomi], axis=0)

    df_pre.to_csv(ALREADY_GEOCODED_PATH, sep=';', index=False)
    df_miss_nomi.to_csv(OUTPATH_GEOCODED_MISSES_NOMINATIM, sep=';', index=False)


def address_to_coordinates_fuzzy_matching():

    df_addr = pd.read_csv(ALREADY_GEOCODED_PATH, sep=';')
    df_miss = pd.read_csv(OUTPATH_GEOCODED_MISSES_NOMINATIM, sep=';')

    df_clear = transfer_geometry_from_geocoded_address_via_fuzzy_matching(df_addr=df_addr,
                                                                          df_miss=df_miss,
                                                                          plzcol_l='post_code',
                                                                          plzcol_r='post_code',
                                                                          addrcol_l='point_address',
                                                                          addrcol_r='address',
                                                                          geomcol_l='geometry',
                                                                          thresh=0.93)

    addrcol_r = 'address'
    df_miss['merge_address'] = df_miss[addrcol_r].str.replace(' ', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('.', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace(',', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('-', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('/', '', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('ä', 'ae', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('ü', 'ue', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('ö', 'oe', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('ß', 'ss', regex=False)
    df_miss['merge_address'] = df_miss['merge_address'].str.replace('asse', '', regex=False)

    df_miss.drop_duplicates(subset=['merge_address'], inplace=True)

    df_miss.drop(columns=['fclass', 'geometry'], inplace=True)
    df_miss = pd.merge(df_miss, df_clear, how='left', left_on='merge_address', right_on='address')
    df_clear = df_miss.loc[df_miss['geometry'].notna()].copy()
    df_clear.drop(columns=['address_y', 'merge_address'], inplace=True)
    df_clear.rename(columns={"address_x":"address"}, inplace=True)
    df_addr.drop(columns=['merge_address'], inplace=True)

    df_miss = df_miss.loc[df_miss['geometry'].isna()].copy()
    df_miss.drop(columns=['merge_address', 'address_y', 'geometry', 'geocoding', 'point_address'], inplace=True)
    df_miss.rename(columns={"address_x": "address"}, inplace=True)

    df_addr = pd.concat([df_addr, df_clear], axis=0)
    df_addr.to_csv(ALREADY_GEOCODED_PATH, sep=';', index=False)
    df_miss.to_csv(OUTPATH_GEOCODED_MISSES_FUZZY, sep=';', index=False)


def addresses_to_coordinates_with_osm_l2_data():
    df_miss_nomi = pd.read_csv(OUTPATH_GEOCODED_MISSES_FUZZY, sep=';')

    ## read localities of Germany
    gdf_loc = gpd.read_file(LOC_OSM_PATH)
    df_loc = pd.DataFrame(gdf_loc)[['geometry', 'fclass', 'plz', 'NAME_4']]
    df_loc['plz'] = df_loc['plz'].fillna('')
    df_loc['plz'] = df_loc['plz'].astype('str')

    ## clean level 2 localities string (column name: NAME_4)
    df_loc['NAME_4'] = df_loc['NAME_4'].str.lower()
    df_loc['NAME_4'] = df_loc['NAME_4'].apply(
        helper_functions.replace_characters_in_string, char_lst=['\n'], replace_string=',')
    df_loc['NAME_4'] = df_loc['NAME_4'].apply(
        helper_functions.replace_characters_in_string, char_lst=[' - '], replace_string=' | ')
    df_loc['NAME_4'] = df_loc['NAME_4'].apply(
        helper_functions.replace_characters_in_string, char_lst=['\x9a', '\x8a', '(', ')', '/'], replace_string='')
    df_loc['NAME_4'] = df_loc['NAME_4'].apply(
        helper_functions.replace_characters_in_string, char_lst=['-'], replace_string=' ')
    df_loc['NAME_4'] = df_loc['NAME_4'].apply(helper_functions.clean_city_text, search_terms=['ot', 'bei', '|'])

    ## create unique identifiers for level 2
    df_loc['pc_city_l2'] = df_loc['plz'] + df_loc['NAME_4']
    df_loc['pc_city_l2'] = df_loc['pc_city_l2'].str.replace(' ', '')

    df_loc_l2 = df_loc.drop_duplicates(['pc_city_l2'])

    ## combine misses with unique identifier of level 2
    ## divide between matches and misses
    df_comb_l2 = pd.merge(df_miss_nomi[['street', 'post_code', 'city', 'address', 'full_address', 'plz_len', 'pc_city']],
                          df_loc_l2[['pc_city_l2', 'fclass', 'geometry']],
                          how='left', left_on='pc_city', right_on='pc_city_l2')
    df_clear_l2 = df_comb_l2[df_comb_l2['geometry'].notnull()].copy()
    df_clear_l2['geocoding'] = 'osm_level2'
    df_clear_l2 = df_clear_l2.drop(columns='pc_city_l2')
    df_miss_l2 = df_comb_l2[df_comb_l2['geometry'].isna()].copy()

    df_clear_l2['point_address'] = df_clear_l2['post_code'].astype(int).astype(str) + ' ' + df_clear_l2['city']
    df_clear_l2['geocoding'] = 'osm_level2'

    df_clear_l2 = df_clear_l2[['street', 'post_code', 'city', 'address', 'full_address', 'plz_len',
                               'pc_city', 'geometry', 'geocoding', 'point_address']]

    df_miss_l2 = df_comb_l2[df_comb_l2['geometry'].isna()].copy()
    df_miss_l2['geocoding'] = 'not possible'
    df_miss_l2['point_address'] = None
    df_miss_l2 = df_miss_l2[['street', 'post_code', 'city', 'address', 'full_address', 'plz_len',
                             'pc_city', 'geometry', 'geocoding', 'point_address']]

    ## open already geocoded addresses and append newly geocoded addresses to them
    df_pre = pd.read_csv(ALREADY_GEOCODED_PATH, sep=';')
    df_pre = pd.concat([df_pre, df_clear_l2, df_miss_l2], axis=0)

    df_pre.to_csv(OUTPATH_GEOCODED_ALL, sep=';', index=False)


def combine_addresses_with_administrative_levels():
    helper_functions.print_red("!! DO THIS ONLY AFTER YOU CONVERTED THE ADDRESSES TO WGS 84 COORDINATES IN QGIS !!")

    shp_addr = gpd.read_file(rf"{ADDRESSES_SHP}")
    shp_adm = gpd.read_file(rf"{ADMINISTRATIVE_SHP}")

    addr_with_adm = gpd.sjoin(shp_addr, shp_adm[['NAME_1', 'NAME_4', 'geometry']], how="inner", op='intersects')

    df_miss = shp_addr.loc[~shp_addr['address'].isin(addr_with_adm['address'])].copy()
    df_miss['NAME_1'] = 'Ausland'
    df_miss['NAME_4'] = 'Ausland'
    df_miss.loc[df_miss['geometry'].isna(), 'NAME_1'] = 'Unbekannt'
    df_miss.loc[df_miss['geometry'].isna(), 'NAME_4'] = 'Unbekannt'

    df_out = pd.DataFrame(addr_with_adm)
    df_out.drop(columns='index_right', inplace=True)

    df_out = pd.concat([df_out, df_miss], axis=0)

    df_out.rename(columns={"NAME_1": "fstateofowner", "NAME_4": "parishofowner"}, inplace=True)

    df_out.to_csv(GEOCODED_ADDRESSES, sep=';', index=False)

def combine_alkis_addresses_and_geometries(df_alkis, df_addr, addr_col_left, addr_col_right, geocoding=None,
                                           subset_addr_df=None):
    print(f"\nGeocoding: {geocoding}")
    print(f"Alkis merge column: {addr_col_left}; Address merge column: {addr_col_right}")

    if geocoding:
        df_addr = df_addr.loc[df_addr['geocoding'] == geocoding]

    if subset_addr_df:
        df_addr = df_addr[subset_addr_df]

    df_succ = pd.DataFrame(columns=list(df_alkis.columns))

    if df_addr.empty:
        print(f'There are no addresses that were geocoded with {geocoding}!')
        return df_succ, df_alkis

    if addr_col_right not in df_addr.columns:
        print(f'{addr_col_right} is not in columns of df_addr!')
        return df_succ, df_alkis

    if addr_col_left not in df_alkis.columns:
        print(f'{addr_col_left} is not in columns of df_alkis!')
        return df_succ, df_alkis

    df_addr.drop_duplicates(subset=[addr_col_right], inplace=True)

    df_addr['merge_address'] = df_addr[addr_col_right].str.replace(' ', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('.', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace(',', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('-', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('/', '', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('ä', 'ae', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('ü', 'ue', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('ö', 'oe', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('ß', 'ss', regex=False)
    df_addr['merge_address'] = df_addr['merge_address'].str.replace('asse', '', regex=False)

    df_alkis['merge_address'] = df_alkis[addr_col_left].str.replace(' ', '', regex=False)
    df_alkis['merge_address'] = df_alkis['merge_address'].str.replace('.', '', regex=False)
    df_alkis['merge_address'] = df_alkis['merge_address'].str.replace(',', '', regex=False)
    df_alkis['merge_address'] = df_alkis['merge_address'].str.replace('-', '', regex=False)
    df_alkis['merge_address'] = df_alkis['merge_address'].str.replace('/', '', regex=False)
    df_alkis['merge_address'] = df_alkis['merge_address'].str.replace('ä', 'ae', regex=False)
    df_alkis['merge_address'] = df_alkis['merge_address'].str.replace('ü', 'ue', regex=False)
    df_alkis['merge_address'] = df_alkis['merge_address'].str.replace('ö', 'oe', regex=False)
    df_alkis['merge_address'] = df_alkis['merge_address'].str.replace('ß', 'ss', regex=False)
    df_alkis['merge_address'] = df_alkis['merge_address'].str.replace('asse', '', regex=False)

    df_addr.drop_duplicates(subset=['merge_address'], inplace=True)

    df_merge = pd.merge(df_alkis, df_addr, how='left', on='merge_address')
    df_succ = df_merge.loc[df_merge.geometry.notna()].copy()
    df_fail = df_merge.loc[~df_merge.geometry.notna()].copy()

    df_fail.drop(columns=subset_addr_df + ['merge_address'], inplace=True)
    df_succ.drop(columns=['merge_address'], inplace=True)

    print(f"Number entries original df: {len(df_alkis)}\n"
          f"Number entries success: {len(df_succ)}\n"
          f"Number entries missed: {len(df_fail)}\n"
          f"Original - Success = {len(df_alkis) - len(df_succ)} Difference to missed: {(len(df_alkis) - len(df_succ)) - len(df_fail)}")

    if df_succ.empty:
        print('No entries could be matched!\n')
        df_alkis.drop(columns=['merge_address'], inplace=True)
        return df_succ, df_alkis

    return df_succ, df_fail


def combine_addresses_with_geolocations():
    os.chdir(WD)

    df_alkis = pd.read_csv(ALKIS_PTH, sep=';')

    length_input = len(df_alkis)

    df_addr = pd.read_csv(GEOCODED_ADDRESSES, sep=';', low_memory=False)
    df_addr_na = df_addr.loc[df_addr['geometry'].isna()].copy()
    df_addr.dropna(subset=["geometry"], inplace=True)
    df_addr.rename(columns={"point_addr": "point_address", "full_addre": "full_address"}, inplace=True)

    addr_cols = ['address', 'full_address', 'geometry', 'geocoding', 'point_address', 'fstateofowner', 'parishofowner']


    df_lst = []
    for geocoding in ['google_api', 'nominatim', 'osm_level1']:
        for left_col in ['clean_address']:  #, 'addresses'
            for right_col in ['address', 'point_address']:
                df_succ, df_alkis = combine_alkis_addresses_and_geometries(df_alkis=df_alkis,
                                                                          df_addr=df_addr,
                                                                          addr_col_left=left_col,
                                                                          addr_col_right=right_col,
                                                                          geocoding=geocoding,
                                                                          subset_addr_df=addr_cols)

                if not df_succ.empty:
                    df_out = df_succ.copy()
                    df_lst.append(df_out)

    ## This is done separately and after the others, because for fuzzy matching you should not use point address
    for geocoding in ['fuzzy_matching']:
        for left_col in ['clean_address']:  #, 'addresses'
            for right_col in ['address']:
                df_succ, df_alkis = combine_alkis_addresses_and_geometries(df_alkis=df_alkis,
                                                                          df_addr=df_addr,
                                                                          addr_col_left=left_col,
                                                                          addr_col_right=right_col,
                                                                          geocoding=geocoding,
                                                                          subset_addr_df=addr_cols)

                if not df_succ.empty:
                    df_out = df_succ.copy()
                    df_lst.append(df_out)

    for geocoding in ['osm_level2']:
        for left_col in ['clean_address']:  #, 'addresses'
            for right_col in ['address', 'point_address']:
                df_succ, df_alkis = combine_alkis_addresses_and_geometries(df_alkis=df_alkis,
                                                                          df_addr=df_addr,
                                                                          addr_col_left=left_col,
                                                                          addr_col_right=right_col,
                                                                          geocoding=geocoding,
                                                                          subset_addr_df=addr_cols)

                if not df_succ.empty:
                    df_out = df_succ.copy()
                    df_lst.append(df_out)

    df_miss = df_alkis.dropna(subset=['clean_address']).copy()
    df_miss = df_miss.loc[~df_miss['clean_address'].str.contains('00000')].copy()
    df_miss = df_miss.loc[~df_miss['clean_address'].str.contains('unknown')].copy()
    df_miss = df_miss.loc[~df_miss['clean_address'].str.contains('unbekannt')].copy()
    df_miss.drop_duplicates(subset=['clean_address'], inplace=True)

    df_miss.to_csv(OUTPATH_MISSING, sep=';', index=False)

    print("Addresses with no geolocation:", len(df_miss))

    for col in addr_cols:
        df_alkis[col] = None

    df_lst.append(df_alkis)
    df_out = pd.concat(df_lst, axis=0)

    print(f"Number of original entries: {length_input}\n"
          f"Number of output entries: {len(df_out)}")

    df_out.loc[df_out['geocoding'].isna(), 'geocoding'] = 'not_possible'
    df_out.to_csv(OUTPATH_FINAL, sep=';', index=False)


def calculate_owner_distance_to_parcel():
    print("!! Transform owner_df_pth to shapefile before in QGIS. And transform geometry from EPSG 25832 to 4326 !!")

    owner_loc = pd.read_csv(OUTPATH_FINAL, sep=';')
    parcels = gpd.read_file(PARCELS_PTH)

    print("Calculate centroids of parcels")
    parcels['parcel_loc'] = parcels['geometry'].centroid.to_crs(epsg=4326)
    parcels.rename(columns={'geometry': 'owner_loc'}, inplace=True)

    owner_df = pd.merge(owner_loc, parcels[['OGC_FID', 'parcel_loc']], how='left')
    ## This step doesn't make sense
    distance_df = owner_df.drop_duplicates(subset=['owner_merge', 'clean_address']).copy()
    distance_df['distance'] = distance_df.apply(lambda row: helper_functions.wkt_point_distance(row.parcel_loc, row.owner_loc), axis=1)

    owner_df = pd.merge(owner_df, distance_df[['owner_merge', 'clean_address', 'distance']], how='left',
                        on=['owner_merge', 'clean_address'])

    owner_df.to_csv(OUTPATH_DISTANCES, sep=';', index=False)


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    helper_functions.create_folder("05_georeference_addresses")

    addresses_to_coordinates_with_osm_l1_data()
    addresses_to_coordinates_with_nominatim()
    address_to_coordinates_fuzzy_matching()
    addresses_to_coordinates_with_osm_l2_data()

    ## DO THIS ONLY AFTER YOU CONVERTED THE ADDRESSES TO WGS 84 COORDINATES IN QGIS
    combine_addresses_with_administrative_levels()

    combine_addresses_with_geolocations()
    calculate_owner_distance_to_parcel()

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)

if __name__ == '__main__':
    main()


