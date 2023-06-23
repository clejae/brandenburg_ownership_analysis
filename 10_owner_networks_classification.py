# Author:
# github repository:

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import json
import time
import os
import pandas as pd
import numpy as np
import geopandas as gpd
import wordcloud

## Project library
import helper_functions
# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input pathes
## Miscellaneous
BB_SHP = r"00_data\vector\administrative\BB_municipalities.shp"
STOPWORD_LIST = r"00_data\miscellaneous\stopwords_german.txt"
CURR_FOLDER = "10_owner_network_classification"

## Parcels
ALKIS_IACS_PTH = r"09_alkis_intersection_with_other_layers\alkis_iacs_inters.shp"
ALKIS_IACS_MUNICIP_PTH = r"09_alkis_intersection_with_other_layers\alkis_munic_iacs_inter.shp"
ALK_MUNICIP_IACS_UNION_PTH = r"09_alkis_intersection_with_other_layers\alkis_munic_iacs_union.shp"

## Owners
OWNERS_W_THRESH_PTH = r"08_network_analysis\owners+communities_thr{0}\08_owners_stretched+comm_thr{0}.csv"
COMMUNITY_DAFNE_CSV_W_THRESH = r"08_network_analysis\network_connections_with_community.csv"
COMMUNITY_INFO_DICT_W_THRESH = r"08_network_analysis\owners+communities_thr{0}\08_alkis_owners_comm_thr{0}_info_dict.json"
NETW_CONN_W_TRESH_PTH = r"08_network_analysis\owners+communities_thr{0}\08_comp_comm_mcomp_thr{0}_dict.json"

NETW_BRANCH_PTH = r"08_network_analysis\all_companies_branches_and_locations.csv"
NETW_COMM_PTH = r"08_network_analysis\network_connections_with_community.csv"

## Output
## Owners paths
OWNERS_W_THRESH_DIST_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-dist.csv"
OWNERS_W_THRESH_CLEANED_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-dist+cleaned.csv"
COMMUNITY_INFO_FROM_DAFNE_W_THRESH = r"10_owner_network_classification\community_infos_from_dafne_thr{0}.csv"
OWNERS_W_THRESH_REDUCED_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-reduced.csv"
OWNERS_W_THRESH_AND_LOC_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-dist+cleaned+loc.csv"
OWNERS_W_THRESH_AND_LOC_CLASSIFIED_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-dist+cleaned+loc+class.csv"

## Classification
COMMUNITY_CLASSIFICATION_PTH = r"10_owner_network_classification\10_comm_classific_and_info_thr{0}.csv"

## Miscellaneous
MUNICIP_NEIGHBOR_DICT_PTH = r"10_owner_network_classification\municipality_neighbors_dict.json"
CHARACTERISTICS_OF_COMMS_PTH = r"10_owner_network_classification\10_characteristics_comms_from_alkis_thr{0}.csv"

# ------------------------------------------ LOAD DATA & PROCESSING ------------------------------------------#

def calculate_distances(owner_df_pth, out_pth):
    """
    Calculates distances of parcels to owner locations.
    Args:
        owner_df_pth: Path to owner data frame.
        out_pth: Output path to owner data frame with distances.

    Returns:

    """
    print("Calculate distances of parcels to owner locations!")
    print("\tRead owner data")
    df_owners = pd.read_csv(owner_df_pth, sep=";")

    print("\tCalculate distances")
    df_owners['distance'] = df_owners.apply(
        lambda row: helper_functions.wkt_point_distance(row.parcel_loc, row.geometry), axis=1)

    print("\tWrite out")
    df_owners.to_csv(out_pth, sep=";", index=False)


def unify_attribute_in_alkis(owners_pth, owners_pth_cleaned):
    """
    Unifies all attributes across all entries of a unique owner ID (owner merge) for further processing.
    Args:
        owners_pth: Path to owner dataframe.
        owners_pth_cleaned: Output path to cleaned owner dataframe.

    Returns:

    """
    print("Unify all attributes of all owner_merge for further processing.")
    print("\tRead owner information")
    df_alk = pd.read_csv(owners_pth, sep=";")

    cols = ["level1", "level2", "level3", "clean_address", "full_address", "geocoding", "point_address",
            "fstateofowner", "parishofowner", "agric", "city", "community_50", "mother_company"]

    ## Unify attributes per owner_merge
    print("\tGet all attributes per owner_merge")
    df_agg = df_alk.groupby("owner_merge").agg(
        level1=pd.NamedAgg("level1", list),
        level2=pd.NamedAgg("level2", list),
        level3=pd.NamedAgg("level2", list),
        clean_address=pd.NamedAgg("clean_address", list),
        full_address=pd.NamedAgg("full_address", list),
        geocoding=pd.NamedAgg("geocoding", list),
        point_address=pd.NamedAgg("point_address", list),
        fstateofowner=pd.NamedAgg("fstateofowner", list),
        parishofowner=pd.NamedAgg("parishofowner", list),
        agric=pd.NamedAgg("agric", list),
        city=pd.NamedAgg("city", list),
        community_50=pd.NamedAgg("community_50", list),
        mother_company=pd.NamedAgg("mother_company", list),
    ).reset_index()

    print("\tAssign most frequent value as final value.")

    for col in cols:
        df_agg[col] = df_agg[col].apply(helper_functions.get_most_frequent_item)

        clean_dict = dict(zip(df_agg.owner_merge, df_agg[col]))
        df_alk[col] = df_alk["owner_merge"].map(clean_dict)

    print("\tWrite results out.")
    df_alk.to_csv(owners_pth_cleaned, sep=";", index=False)


def reduce_alkis_data_to_necessary_information(owners_pth, comm_col, out_pth, alkis_iacs_inters_pth=None):
    """
    Reduces ALKIS data frame to unique entries, i.e. all owners should only occur once.
    Args:
        owners_pth: Path to owner data frame.
        comm_col: Column name of community IDs
        out_pth: Output path to reduce data frame.
        alkis_iacs_inters_pth: Optional: provide the parcels intersected with the IACS data to only look at parcels that are actually used for agriculture.

    Returns:

    """

    print("\tRead owner information")
    df_alk = pd.read_csv(owners_pth, sep=";")

    if alkis_iacs_inters_pth:
        print("\tRead parcels")
        gdf = gpd.read_file(alkis_iacs_inters_pth)
        gdf['area'] = gdf['geometry'].area
        gdf['area'] = gdf['area'] / 10000

        print("\tCombine both")
        df_alk = helper_functions.combine_parcels_with_owners(gdf[["OGC_FID", "EIGENTUEME"]], df_alk)

    df_alk.drop_duplicates(subset=[comm_col, "owner_merge"], inplace=True)
    df_alk.drop(columns=["OGC_FID"], inplace=True)

    df_alk.to_csv(out_pth, sep=";", index=False)


def create_neighbor_dictionary_for_municipalities(municipality_shp_pth, out_pth):
    """
    Creates a dictionary providing all neighbour municipalities of a municipality
    Args:
        municipality_shp_pth: Path to shapefile with muncipalities.
        out_pth: Output path to dictionary (json) with neigbourhood information.

    Returns:

    """

    gdf = gpd.read_file(municipality_shp_pth)

    out_dict = {}
    for index, row in gdf.iterrows():
        neighbors = gdf[gdf.geometry.touches(row['geometry'])].RS.tolist()
        if row.RS in neighbors:
            neighbors = neighbors.remove(row.RS)
        out_dict[row.RS] = neighbors

    with open(out_pth, "w") as file:
        json.dump(out_dict, file, indent=4)


def plot_correlation_between_measures(df_measures, cols, out_pth):

    import seaborn as sns

    df_measures = df_measures[cols]

    corr = df_measures.corr()
    for col in corr:
        corr.loc[corr[col].abs() < 0.2, col] = None
    ax = sns.heatmap(
        corr,
        vmin=-1, vmax=1, center=0,
        cmap=sns.diverging_palette(20, 220, n=200),
        square=True
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45,
                       horizontalalignment='right')

    fig = ax.get_figure()
    fig.savefig(out_pth)

    print()


def calculate_area_per_community(alkis_iacs_inters_pth, owners_pth, comm_col):
    print("Calculate areas.")
    print("Combine parcels with owner data")
    gdf = gpd.read_file(alkis_iacs_inters_pth)
    gdf["area"] = gdf["geometry"].area
    gdf['area'] = gdf['area'] / 10000
    df = pd.read_csv(owners_pth, sep=';', dtype={comm_col: str})

    df_alk = helper_functions.combine_parcels_with_owners(gdf[["OGC_FID", "EIGENTUEME", "BTNR", "area"]],
                                          df[["OGC_FID", "mother_company", "owner_merge", "level_c_category", "level3", comm_col]])

    ## Get number of original entries for "Unbekannt"
    df_unknown = df_alk.loc[df_alk["mother_company"] == 'unbekannt'].copy()
    num_unkown_owners = len(df_unknown["EIGENTUEME"].unique())
    print("number of original entries that are now classified into 'unbekannt':", num_unkown_owners)
    # print("Examples:")
    # for i in t["EIGENTUEME"].unique()[:10]:
    #     print(i)

    ## Get the area per community, add number and names of all owners from ALKIS
    df_area = df_alk[["area", comm_col]].groupby(comm_col).sum().reset_index()
    df_area.sort_values(by="area", ascending=False, inplace=True)
    total_area = df_alk['area'].sum()
    df_area['share'] = df_area['area'] / total_area * 100
    df_area['cum_share'] = df_area['share'].cumsum()
    df_area['area'] = round(df_area['area'], 5)

    return df_area


def derive_community_information_from_dafne_data(dafne_pth, info_dict_pth, comm_col, net_branch_pth, stopword_lst_pth, out_pth):
    """
    Derives some information on the networks only from the DAFNE data.
    Args:
        dafne_pth: Path to DAFNE data with community IDs.
        info_dict_pth: Path to info-dictionary about networks.
        comm_col: Column name of community/network IDs.
        net_branch_pth: Path to data frame with company names and economic branches.
        stopword_lst_pth: Path to text file with stopwords.
        out_pth: Output path to data frame with network information derived from DAFNE data.

    Returns:

    """
    print("Derive network information from DAFNE data.")

    print("\tRead data.")
    df = pd.read_csv(dafne_pth, sep=',', dtype={comm_col: str})
    df_branch = pd.read_csv(net_branch_pth, dtype={"postcode": str})
    with open(info_dict_pth) as json_file:
        comm_info_dict = json.load(json_file)

    ## Get no. subcompanies and mother company name dictionaries
    num_subcomps_dict = comm_info_dict["num_subcomps_dict"]
    # num_subcomps_dict = {float(key): int(float(num_subcomps_dict[key])) for key in num_subcomps_dict}
    name_mcomp_dict = comm_info_dict["comm_dict"]

    ## Concatenate left and right connections "vertically" to get unique network members
    df_l = df[["conn_left", comm_col]].copy()
    df_r = df[["conn_right", comm_col]].copy()
    df_r.rename(columns={"conn_right": "conn_left"}, inplace=True)
    df_lo = pd.concat([df_l, df_r], axis=0)
    df_lo.drop_duplicates(inplace=True)

    ## Add DAFNE information to network members
    df_lo = pd.merge(df_lo, df_branch[["company_name", "main_branch", "agric", "fstate", "self_descr", "international"]],
                      how="left", left_on="conn_left", right_on="company_name")
    df_lo.drop(columns="company_name", inplace=True)

    t = df_lo.loc[df_lo["main_branch"].isna()].copy()
    df_lo.loc[df_lo["main_branch"].isna(), "main_branch"] = "unkown"
    df_lo.loc[df_lo["fstate"].isna(), "fstate"] = "unkown"
    df_lo.loc[df_lo["self_descr"].isna(), "self_descr"] = "unkown"

    ## Derive network information by aggregating per community ID
    print("\tAggregate by community.")
    def count_people(group):
        group = pd.DataFrame(group)
        col = group.columns[0]
        sub_group = group.loc[group[col].str.count('_') > 0].copy()
        return len(sub_group)

    def count_companies(group):
        group = pd.DataFrame(group)
        col = group.columns[0]
        sub_group = group.loc[group[col].str.count('_') == 0].copy()
        return len(sub_group)

    def get_main_attribute(group):
        group = pd.DataFrame(group)
        col = group.columns[0]
        sub_group = group.loc[group[col] != "unkown"].copy()
        out = sub_group[col].mode()
        out = '_'.join(out.tolist())
        return out

    def get_all_attributes(group):
        group = pd.DataFrame(group)
        col = group.columns[0]
        sub_group = group.loc[group[col] != "unkown"].copy()
        sub_group.drop_duplicates(subset=col, inplace=True)
        sub_group.sort_values(by=col, inplace=True)
        out = '_'.join(sub_group[col].tolist())
        return out

    df_stat = df_lo.groupby(by=[comm_col]).agg(
        number_members=pd.NamedAgg(column="conn_left", aggfunc="count"),
        number_people=pd.NamedAgg(column="conn_left", aggfunc=count_people),
        number_companies=pd.NamedAgg(column="conn_left", aggfunc=count_companies),
        number_agric_comps=pd.NamedAgg(column="agric", aggfunc="sum"),
        number_intern_comps=pd.NamedAgg(column="international", aggfunc="sum"),
        main_branch=pd.NamedAgg(column="main_branch", aggfunc=get_main_attribute),
        all_branches=pd.NamedAgg(column="main_branch", aggfunc=get_all_attributes),
        main_fstate=pd.NamedAgg(column="fstate", aggfunc=get_main_attribute),
        all_fstates=pd.NamedAgg(column="fstate", aggfunc=get_all_attributes),
        all_descr=pd.NamedAgg(column="self_descr", aggfunc=lambda x: ' '.join(x))
    ).reset_index()

    df_stat["share_agric_comps"] = round(df_stat["number_agric_comps"] / df_stat["number_companies"], 1)
    df_stat["num_owners_in_alkis"] = df_stat[comm_col].map(num_subcomps_dict)
    df_stat.loc[df_stat["num_owners_in_alkis"].isna(), "num_owners_in_alkis"] = 0
    df_stat["num_owners_in_alkis"] = df_stat["num_owners_in_alkis"].astype(int)
    df_stat["name_mcomp"] = df_stat[comm_col].map(name_mcomp_dict)

    # t = df_stat.groupby(["bb_located", "only_bb", "single_company", "agric_related", "international"]).count().reset_index()

    ## Derive information from descriptions
    with open(stopword_lst_pth, "r") as file:
        stopwords = file.readlines()

    stopwords = [word.replace('\n', '') for word in stopwords if len(word) > 2]
    stopwords += ["insbesondere", "für", "unkown", "sonstige", "stehend", "stehenden", "genannt", "genannten", "bzw",
                  "bzw.", "sontige", "sonstigen"]

    def get_main_descriptors(description):
        description = description.replace("unkown", "")
        if description.replace(" ", "") == "":
            return "unkown"
        wc = wordcloud.WordCloud(stopwords=stopwords).generate(description)
        wc_words = [item[0][0] for item in wc.layout_]
        wc_words = '|'.join(wc_words)
        return wc_words

    print("\tGet main descriptors.")
    df_stat["main_words"] = df_stat["all_descr"].apply(get_main_descriptors)

    print("\tWrite out")
    df_stat.to_csv(out_pth, index=False)

    # ## First classify manually
    # number_communities = len(df_stat)
    # print("Number communities:", number_communities)
    #
    # df_stat["class"] = "0_unclassified"
    # df_stat.loc[df_stat["number_companies"] == 1, "class"] = "1_single_companies"
    # t = df_stat.loc[df_stat["class"] == "1_single_companies"].copy()
    # number_single_companies = len(t)
    # print("Number single companies:", number_single_companies)



    # print(df_stat.describe(percentiles=[.05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95]))
    # plot_correlation_between_measures(df_measures=df_stat,
    #                                   cols=["number_members", "number_people", "number_companies"],
    #                                   out_pth=r"13_community_classification\community_characteristics_correlation_w_thresh.png")
    ## Select data for analysis
    # features = ["number_members", "number_people", "number_companies", "num_owners_in_alkis"]
    #
    # x = df_stat.loc[:, features].values  # Separating out the target
    #
    # ## Standardize data
    # # x = StandardScaler().fit_transform(x)
    #
    # ## Perform K means clustering on standardized data
    # kmeans = KMeans(n_clusters=3)
    # kmeans.fit(x)
    # y_kmeans = kmeans.predict(x)
    #
    # ## Plot Scatterplot with K-Means groups
    # df_plt = pd.DataFrame(x, columns=features)
    # df_plt['kmeans_groups'] = y_kmeans
    # g = sns.PairGrid(df_plt, hue='kmeans_groups')
    # g.map_diag(sns.histplot)
    # g.map_offdiag(sns.scatterplot, s=3)
    # out_pth_fig = fr"11_community_classification\kmeans_clustering_scatterplots_{comm_col}.png"
    # g.savefig(out_pth_fig)
    #
    # fig, axs = plt.subplots(nrows=4)
    # sns.boxplot(y='number_members', data=df_plt, x='kmeans_groups', showfliers=False, ax=axs[0])
    # sns.boxplot(y='number_people', data=df_plt, x='kmeans_groups', showfliers=False, ax=axs[1])
    # sns.boxplot(y='number_companies', data=df_plt, x='kmeans_groups', showfliers=False, ax=axs[2])
    # sns.boxplot(y='num_owners_in_alkis', data=df_plt, x='kmeans_groups', showfliers=False, ax=axs[3])
    # out_pth_fig = rf"11_community_classification\kmeans_boxplots_scatterplots_{comm_col}.png"
    # fig.savefig(out_pth_fig)
    #
    # df_stat['kmeans_groups'] = y_kmeans
    # df_stat.to_csv(out_pth, index=False)


def classify_communities_from_info_dict(info_dict_pth, alkis_reduced_pth, netw_conn_pth, comm_col, owners_pth, out_pth,
                                        curr_folder):
    """
    Classifies the networks/communities with help of the information dictionary (creates "level_c_category").
    In between some manual input might be necessary.
    Args:
        info_dict_pth: Path to information dictionary.
        alkis_reduced_pth: Path to reduced owner data frame.
        netw_conn_pth: Path to csv with network connections and community ID.
        comm_col: Column name of community ID
        owners_pth: Path to cleaned owner dataframe.
        out_pth: Output path to owner data frame with classified networks and locations of mother companies.
        curr_folder: Path to current folder for intermediate outputs.

    Returns:

    """
    print("Classify communities into level-c-category.")

    df_netw = pd.read_csv(netw_conn_pth, dtype={comm_col: str})

    print("\tRead owner data")
    df_alk = pd.read_csv(alkis_reduced_pth, sep=';', dtype={comm_col: str})

    print("\tRead community info dictionary")
    with open(info_dict_pth) as json_file:
        comm_info_dict = json.load(json_file)

    comm_dict = comm_info_dict["comm_dict"]
    num_subcomps_dict = comm_info_dict["num_subcomps_dict"]
    name_subcomps_dict = comm_info_dict["name_subcomps_dict"]

    print("\tNumber of communities in total:", len(num_subcomps_dict))

    single_owner_comms = [key for key in num_subcomps_dict if num_subcomps_dict[key] == 1]
    multi_owner_comms = [key for key in num_subcomps_dict if (num_subcomps_dict[key] > 1) and (float(key.split("_")[0]) < 10000)]
    family_owner_comms = [key for key in num_subcomps_dict if (num_subcomps_dict[key] > 1) and (float(key.split("_")[0]) >= 10000)]

    print("\tNumber of single owners:", len(single_owner_comms))
    print("\tNumber of family groups:", len(family_owner_comms))
    print("\tNumber of multi-owner groups:", len(multi_owner_comms))
    print("\tCheck if numbers add up:", len(num_subcomps_dict) == len(single_owner_comms) + len(family_owner_comms) + len(multi_owner_comms))
    print("\n")

    single_private_people = [key for key in single_owner_comms if '_' in name_subcomps_dict[key]]
    single_non_people = [key for key in single_owner_comms if '_' not in name_subcomps_dict[key]]

    print("\tNumber of single private people:", len(single_private_people))
    print("\tNumber of single entities that are not private people:", len(single_non_people))
    print("\n")

    ## Create column for categorisation of communities and assign no data value
    df_alk["level_c_category"] = '99'

    ## First identify all owners that occur in communities with multiple owners
    df_alk.loc[df_alk[comm_col].isin(multi_owner_comms), "level_c_category"] = '90'

    ## Get all possible combinations of level1 that can occur in the communities
    df_mul = df_alk.loc[df_alk["level_c_category"] == '90'].copy()
    print("\tNumber of communities with multiple owner that are not just private people:", len(df_mul))
    df_mul["level1"] = df_mul["level1"].astype(str)
    def create_level1_combinations(group):
        uni_levels = list(set(group["level1"].tolist()))
        uni_levels.sort()
        out_str = '_'.join(uni_levels)
        return out_str

    df_comb = df_mul.groupby(comm_col).apply(create_level1_combinations).reset_index()
    df_comb.columns = [comm_col, "level3_combs"]
    df_comb_stats = df_comb.groupby("level3_combs").count().reset_index()
    print(df_comb_stats)

    ## Based on the combination of level3, classify the owner groups
    df_mul = pd.merge(df_mul, df_comb, how="left", on=comm_col)

    ## First separate groups out, that consist of only companies or companies with private owners
    df_mul_12 = df_mul.loc[df_mul["level3_combs"] == '1.0_2.0'].copy()
    df_mul_rest = df_mul.loc[df_mul["level3_combs"] != '1.0_2.0'].copy()

    df_mul_rest.loc[df_mul_rest["level3_combs"] == '2.0', "level_c_category"] = '2_9_1'
    df_mul_rest.loc[df_mul_rest["level3_combs"] == '2.0_4.0', "level_c_category"] = '4_9_1'
    df_mul_rest.loc[df_mul_rest["level3_combs"] == '2.0_5.0', "level_c_category"] = '5_9_1'
    df_mul_rest.loc[df_mul_rest["level3_combs"] == '3.0', "level_c_category"] = '3_9_1'

    df_mul_rest.loc[(df_mul_rest["level3_combs"] == '2.0_3.0') & (
                df_mul_rest["mother_company"].str.count("stiftung") > 0), "level_c_category"] = '3_9_1'
    df_mul_rest.loc[(df_mul_rest["level3_combs"] == '2.0_3.0') & (
            df_mul_rest["mother_company"].str.count("verein") > 0), "level_c_category"] = '3_9_1'
    df_mul_rest.loc[(df_mul_rest["level3_combs"] == '2.0_3.0') & (
            df_mul_rest["mother_company"].str.count("landkreis") > 0), "level_c_category"] = '5_9_1'
    df_mul_rest.loc[(df_mul_rest["level3_combs"] == '2.0_3.0') & (
            df_mul_rest["mother_company"].str.count("stadt") > 0), "level_c_category"] = '5_9_1'

    ## All cases where companies occur with private people have to be treated differently
    ## I want to separate communities with more than 1 company from communities with 1 company and 1 or more people
    ## This can only be done in DAFNE, as the ALKIS data don't hold all members of a community
    ## First count companies in DAFNE networks and classify them
    def count_occ_of_value_in_column(group, value, column):
        out = group.loc[group[column] == value]
        return len(out)

    df_netw_comps = df_netw.groupby(comm_col).apply(count_occ_of_value_in_column, "Unternehmen", "isperson").reset_index()
    df_netw_comps.columns = [comm_col, "num_comps"]
    df_netw_comps["num_comps"] += 1

    df_netw_comps.loc[df_netw_comps["num_comps"] > 1, "level_c_category"] = '2_9_1'
    df_netw_comps.loc[df_netw_comps["num_comps"] == 1, "level_c_category"] = '2_9_2'

    ## Add this classification to the groups of only companies or companies with private people
    df_mul_12.drop(columns="level_c_category", inplace=True)
    df_mul_12 = pd.merge(df_mul_12, df_netw_comps[[comm_col, "level_c_category"]], how="left", on=comm_col)

    community_to_levelc = {}

    l = df_mul_12[comm_col].tolist()
    t = df_mul_rest.loc[df_mul_rest[comm_col].isin(l)].copy()

    l = df_mul_rest[comm_col].tolist()
    t = df_mul_12.loc[df_mul_12[comm_col].isin(l)].copy()

    ind = list(df_mul_rest.columns).index(comm_col) + 1
    for row in df_mul_rest.itertuples():
        community = row[ind]
        levelc = row.level_c_category
        if community not in community_to_levelc:
            community_to_levelc[community] = levelc

    ind = list(df_mul_12.columns).index(comm_col) + 1
    for row in df_mul_12.itertuples():
        community = row[ind]
        levelc = row.level_c_category
        community_to_levelc[community] = levelc

    level_c_cat_dict = {
        '1_9_1': 'family network',
        '2_9_1': 'groups of companies',
        '2_9_2': 'company-owner groups',
        '3_9_1': 'groups of non-profit entities',
        '4_9_1': 'church related group of owners',
        '5_9_1': 'companies belonging to public entities'
    }

    ## Assign IDs to single private owners and families as they can be easily identified
    df_alk.loc[df_alk[comm_col].isin(single_private_people), "level_c_category"] = '1_1_1'
    df_alk.loc[df_alk[comm_col].isin(family_owner_comms), "level_c_category"] = '1_9_1'

    ## Now assign to the remaining owners the class IDs of level 3
    # ## For exploration separate them into separate dfs
    # df_rest = df_alk.loc[df_alk["level_c_category"] == '99'].copy()
    # df_rest_comp = df_rest.loc[df_rest["level1"] == 2].copy()
    # df_rest_priv = df_rest.loc[df_rest["level1"] == 1].copy()
    # df_rest_nonp = df_rest.loc[df_rest["level1"] == 3].copy()
    # df_rest_chur = df_rest.loc[df_rest["level1"] == 4].copy()
    # df_rest_publ = df_rest.loc[df_rest["level1"] == 5].copy()
    df_alk.loc[df_alk["level_c_category"] == '99', "level_c_category"] = \
        df_alk.loc[df_alk["level_c_category"] == '99', "level3"]

    ind = list(df_alk.columns).index(comm_col) + 1
    for row in df_alk.loc[df_alk["level_c_category"] != '90'].itertuples():
        community = row[ind]
        levelc = row.level_c_category
        if community not in community_to_levelc:
            community_to_levelc[community] = levelc

    df_alk["level_c_category"] = df_alk[comm_col].map(community_to_levelc)

    ####################################################################################################################
    ## Classify communities into agricultural related (at least one member is agricultural) and non-agricultural related
    level_c_categories = list(df_alk["level_c_category"].unique())

    def get_sum_of_attribute_per_subgroup_in_groups(groups, subgroup_col, attribute_col):

        df_subgroup = groups[[subgroup_col, attribute_col]].groupby(subgroup_col).sum().reset_index()
        df_subgroup.columns = [subgroup_col, f"sum_{attribute_col}"]

        return df_subgroup

    agri = df_alk.groupby("level_c_category").apply(get_sum_of_attribute_per_subgroup_in_groups, comm_col, "agric").reset_index()
    agri["comm_agri_related"] = np.where(agri['sum_agric'] > 0, 1, 0)
    df_alk = pd.merge(df_alk, agri[[comm_col, "comm_agri_related"]], "left", comm_col)


    ####################################################################################################################
    ## Classify mother companies into agricultural/non agricultural and foreign/domestic

    fstates = ['Baden-Württemberg', 'Brandenburg', 'Schleswig-Holstein', 'Sachsen-Anhalt', 'Mecklenburg-Vorpommern',
               'Berlin', 'Niedersachsen', 'Hessen', 'Rheinland-Pfalz', 'Nordrhein-Westfalen', 'Bayern', 'Thüringen',
               'Sachsen', 'Hamburg', 'Saarland', 'Bremen']

    df_alk1 = df_alk.loc[~df_alk["level_c_category"].isin(["2_9_1", "2_9_2", "3_9_1", "4_9_1", "5_9_1"])].copy()
    df_alk2 = df_alk.loc[df_alk["level_c_category"].isin(["2_9_1", "2_9_2", "3_9_1", "4_9_1", "5_9_1"])].copy()
    print("\tSplit df_alk correct:", len(df_alk) == len(df_alk1) + len(df_alk2))

    df_alk1["city_mcomp"] = df_alk1["city"]
    df_alk1["country_mcomp"] = df_alk1["fstateofowner"]
    df_alk1["fstate_mcomp"] = df_alk1["fstateofowner"]
    df_alk1.loc[df_alk1["fstate_mcomp"].isin(fstates), "country_mcomp"] = "Deutschland"
    df_alk1.loc[df_alk1["fstate_mcomp"] == "Ausland", "country_mcomp"] = df_alk1.loc[
        df_alk1["fstate_mcomp"] == "Ausland", "fstate_mcomp"]

    mcomps_all = list(df_alk2["mother_company"].unique())
    mcomps_done = df_alk2.loc[df_alk2["conn_left"].isin(mcomps_all), "mother_company"].tolist()
    mcomps_work = list(set(mcomps_all) - set(mcomps_done))

    df_alk2_done = df_alk2.loc[df_alk2["mother_company"].isin(mcomps_done)].copy()
    df_alk2_work = df_alk2.loc[df_alk2["mother_company"].isin(mcomps_work)].copy()

    df_alk2_done["city_mcomp"] = df_alk2_done["city"]
    df_alk2_done["country_mcomp"] = df_alk2_done["fstateofowner"]
    df_alk2_done["fstate_mcomp"] = df_alk2_done["fstateofowner"]
    df_alk2_done.loc[df_alk2_done["fstate_mcomp"].isin(fstates), "country_mcomp"] = "Deutschland"
    df_alk2_done.loc[df_alk2_done["fstate_mcomp"] == "Ausland", "country_mcomp"] = df_alk2_done.loc[
        df_alk2_done["fstate_mcomp"] == "Ausland", "fstate_mcomp"]

    lst = list(df_alk2_work["mother_company"].unique())
    ands = [item for item in lst if 'AND' in item]
    lst = list(set(lst) - set(ands))
    companies = [item for item in lst if '_' not in item]
    people = [item for item in lst if '_' in item]
    for item in ands:
        items = item.split("AND")
        for subitem in items:
            if '_' in subitem:
                people.append(subitem)
            else:
                companies.append(subitem)

    for comp in companies:
        with open(rf"{curr_folder}\mcomps_without_location_companies_{comm_col}.txt", "a") as file:
            file.write(f"{comp}\n")

    helper_functions.print_red("Assign locations to companies!")

    for person in people:
        with open(rf"{curr_folder}\mcomps_without_location_people_{comm_col}.txt", "a") as file:
            file.write(f"{person};\n")

    for item in ands:
        with open(rf"{curr_folder}\mcomps_without_location_ands_{comm_col}.txt", "a") as file:
            file.write(f"{item};\n")

    ## read classification table
    df_loc = pd.read_excel(rf"{curr_folder}\mother_companies_with_locations.xlsx")
    df_alk2_work["mother_company"] = df_alk2_work["mother_company"].str.strip()
    df_alk2_work = pd.merge(df_alk2_work, df_loc, "left", "mother_company")

    df_alk_out = pd.concat([df_alk1, df_alk2_done, df_alk2_work])
    df_alk_out = df_alk_out[['mother_company', 'level_c_category',
                             'comm_agri_related', 'city_mcomp', 'country_mcomp', 'fstate_mcomp']]

    ids = df_alk_out["mother_company"]
    t = df_alk_out[ids.isin(ids[ids.duplicated()])]
    df_alk_out.drop_duplicates(subset=['mother_company'], inplace=True)

    df_alk_out["fstate_mcomp"] = df_alk_out["fstate_mcomp"].str.strip()
    df_alk_out["country_mcomp"] = df_alk_out["country_mcomp"].str.strip()

    for col in ["city_mcomp", "country_mcomp", "fstate_mcomp"]:
        df_alk_out.loc[df_alk_out[col].isna(), col] = "Unbekannt"

    df_alk_complete = pd.read_csv(owners_pth, sep=";")

    for col in ["fstateofowner"]:
        df_alk_complete.loc[df_alk_complete[col].isna(), col] = "Unbekannt"

    df_alk_complete = pd.merge(df_alk_complete, df_alk_out, "left", on=['mother_company'])

    # print("\tLevel C category of community 5_2:", df_alk_complete.loc[df_alk_complete[comm_col].isin(["5_2"]), "level_c_category"].unique())
    # t = df_alk_complete.loc[df_alk_complete[comm_col] == "5_5"].copy()

    print("\tCalculate distance of mother companies to parcels.")
    df_cn = df_alk_complete.loc[df_alk_complete["level_c_category"].isin(['2_9_1', '2_9_2', '3_9_1', '4_9_1', '5_9_1'])].copy()
    df_cn.drop_duplicates(subset=["mother_company", "city_mcomp"], inplace=True)
    df_cn = df_cn.loc[(df_cn["city"] != df_cn["city_mcomp"]) & (df_cn["city_mcomp"] != "Unbekannt")].copy()
    print("\tNo. of entries for which the distances will be calculated:", len(df_cn))

    ## Create a list of addresses for all global ultimate owners with address from : "city, country"
    mcomp_addrs = []
    for row in df_cn.itertuples():
        mcomp_addrs.append(f"{row.city_mcomp}, {row.country_mcomp}")
    df_cn["mcomp_addr"] = mcomp_addrs

    ## Get geolocation of theses addresses with Nominatim
    df_cn['mcomp_loc'] = df_cn['mcomp_addr'].apply(helper_functions.address_to_coordinates_nominatim)

    ## For all cases where not geolocation could be derived, try again with city name
    df_cn_misses = df_cn.loc[df_cn["mcomp_loc"].isna()].copy()
    df_cn_misses["mcomp_loc"] = df_cn_misses['city_mcomp'].apply(helper_functions.address_to_coordinates_nominatim)

    ## Combine information
    df_cn = df_cn.loc[df_cn["mcomp_loc"].notna()].copy()
    df_cn = pd.concat([df_cn, df_cn_misses])
    df_cn.loc[df_cn["mcomp_loc"].isna(), "mcomp_loc"] = np.nan
    df_cn.loc[df_cn["city_mcomp"] == "Ledge", "mcomp_loc"] = "POINT (11.8383514 52.9103928)"

    ## Add information to ALKIS data frame
    geom_dict = {row.mother_company: row.mcomp_loc for row in df_cn.itertuples()}
    df_alk_complete["mcomp_loc"] = df_alk_complete["mother_company"].map(geom_dict)

    ## Calculate the distances of the global ultimate owners to the parcels
    df_alk_complete.loc[df_alk_complete['mcomp_loc'].notna(), "mcomp_dist"] = df_alk_complete.loc[df_alk_complete['mcomp_loc'].notna()].apply(
        lambda row: helper_functions.wkt_point_distance(row.parcel_loc, row.mcomp_loc), axis=1)
    df_alk_complete.loc[df_alk_complete["mcomp_loc"].isna(), "mcomp_loc"] = df_alk_complete.loc[df_alk_complete["mcomp_loc"].isna(), "geometry"]
    df_alk_complete.loc[df_alk_complete["mcomp_dist"].isna(), "mcomp_dist"] = df_alk_complete.loc[
        df_alk_complete["mcomp_dist"].isna(), "distance"]

    print("\tNumber of entries before export:", len(df_alk_complete))
    df_alk_complete.to_csv(out_pth, sep=';', index=False)


def get_characteristics_of_communities_from_alkis(alkis_iacs_municip_inters_pth, owners_pth, comm_col,
                                                  municips_covered_by_comm_pth, out_pth):
    """
    Derives characteristics of communities/networks from the ALKIS data.
    Args:
        alkis_iacs_municip_inters_pth: Path to ALKIS-municipality-IACS intesection/union shapefile.
        owners_pth: Path to owner data frame with classified networks and locations of mother companies
        comm_col: Column name of community/network IDs.
        municips_covered_by_comm_pth: Path to municipality-neighbourhood dictionary (json)
        out_pth: Output path to data frame with network information derived from ALKIS data.

    Returns:

    """
    print("Get characteristics of communities from ALKIS data (municips covered, farms covered etc.)")
    print("\tCalculate areas.")
    print("\tCombine parcels with owner data")
    gdf = gpd.read_file(alkis_iacs_municip_inters_pth)
    gdf["area"] = gdf["geometry"].area
    gdf['area'] = gdf['area'] / 10000
    df = pd.read_csv(owners_pth, sep=';', dtype={comm_col: str})

    df_alk = helper_functions.combine_parcels_with_owners(gdf[["OGC_FID", "RS", "BTNR", "area"]],
                                         df[["OGC_FID", "mother_company", "agric", "fstateofowner", "level_c_category", comm_col]])

    # def get_municips(community_group):
    #     community_group = pd.DataFrame(community_group)
    #     col = community_group.columns[0]
    #     unique_municips = list(community_group[col].unique())
    #     unique_municips.sort()
    #     return unique_municips

    print("\tGet characteristics")
    ## Get geographic information on communities
    def get_all_unique_attributes_in_list(group):
        out = list(set(list(group)))
        return out

    df_area = df_alk.groupby(by=[comm_col]).agg(
        area=pd.NamedAgg(column="area", aggfunc="sum"),
        municips_covered=pd.NamedAgg(column="RS", aggfunc=get_all_unique_attributes_in_list),
        num_farms_covered=pd.NamedAgg(column="BTNR", aggfunc=get_all_unique_attributes_in_list),
        # number_agric_parcels=pd.NamedAgg(column="agric", aggfunc="sum"),
        all_fstates_owners=pd.NamedAgg(column="fstateofowner", aggfunc=get_all_unique_attributes_in_list),
        level_c_category=pd.NamedAgg(column="level_c_category", aggfunc=get_all_unique_attributes_in_list)
    ).reset_index()

    ## Get information on agricultural owners, do this on entire dataset as some owners
    ## are not included in the combined owner-parcel dataset due to the intersection with IACS data
    df_agri = df.groupby(by=[comm_col]).agg(
        number_agric_parcels=pd.NamedAgg(column="agric", aggfunc="sum")
    ).reset_index()

    df_area = pd.merge(df_area, df_agri, "left", on=comm_col)

    df_area["agric_owners"] = np.where(df_area["number_agric_parcels"] > 0, 1, 0)
    df_area["num_farms_covered"] = [len(lst) for lst in df_area["num_farms_covered"]]
    df_area["level_c_category"] = [lst[0] for lst in df_area["level_c_category"]]
    df_area["all_fstates_owners"] = [[i for i in lst if type(i) == str] for lst in df_area["all_fstates_owners"]]
    df_area["all_fstates_owners"] = ['_'.join(lst) for lst in df_area["all_fstates_owners"]]

    with open(municips_covered_by_comm_pth, "r") as file:
        municip_dict = json.load(file)

    compare_lst = [municip_dict[key] for key in municip_dict]
    for i, key in enumerate(municip_dict):
        compare_lst[i].append(key)

    def sublist_of_any_lst(lst, list_of_lists):
        out = any(set(lst) <= set(lst2) for lst2 in list_of_lists)
        return out

    df_area["regional"] = df_area["municips_covered"].apply(sublist_of_any_lst, list_of_lists=compare_lst)
    df_area["num_municips_covered"] = df_area["municips_covered"].apply(lambda x: len(x))
    df_area.sort_values(by="area", ascending=False, inplace=True)

    print("\tWrite out.")
    df_area.to_csv(out_pth, index=False)


def combine_information_of_communities_and_classify(df_info_alkis_pth, df_info_dafne_pth, owners_pth, comm_col,
                                                    out_pth_classsification, alkis_out_pth):

    """
    Combines the the owner entries in ALKIS with network information derived from DAFNE and ALKIS data.
    Args:
        df_info_alkis_pth: Path to data frame with network information derived from ALKIS data.
        df_info_dafne_pth: Path to data frame with network information derived from DAFNE data.
        owners_pth: Path to owner data frame.
        comm_col: Column name of community IDs
        out_pth_classsification: Output path to data frame with community IDs and classification results.
        alkis_out_pth: Output path to owner data frame with community classification.

    Returns:

    """

    print("Combine all information of communities and derive new classification.")

    ## Read data
    print("\tRead data.")
    df_area = pd.read_csv(df_info_alkis_pth)
    df_stat = pd.read_csv(df_info_dafne_pth)
    df_alk = pd.read_csv(owners_pth, sep=";")

    print("\tCreate new classification.")
    ## Differentiate between communities of single companies and multiple companies, excluding private people in this regard
    df_stat["single_company"] = 1
    df_stat.loc[(df_stat["number_companies"] > 1), "single_company"] = 0

    ## Are there international owners in ALKIS data?
    df_area["intern_owners"] = 0
    df_area.loc[df_area["all_fstates_owners"].str.count("Ausland") > 0, "international"] = 1

    ## Combine all dfs
    df_comb = pd.merge(df_area, df_stat, how="left", on=comm_col)
    print("\tNumber of owners before classification:", len(df_comb))

    ## Combine information of different dfs
    df_comb.loc[df_comb["all_fstates"].isna(), "all_fstates"] = ""
    df_comb.loc[df_comb["all_fstates_owners"].isna(), "all_fstates_owners"] = ""
    df_comb["all_fstates_comb"] = df_comb["all_fstates"] + '_' + df_comb["all_fstates_owners"]
    df_comb.loc[df_comb["all_fstates_comb"].notna(), "all_fstates_comb"] = df_comb.loc[
        df_comb["all_fstates_comb"].notna(), "all_fstates_comb"].apply(lambda x: list(set(x.split('_'))))
    df_comb.loc[df_comb["all_fstates_comb"].notna(), "all_fstates_comb"] = df_comb.loc[
        df_comb["all_fstates_comb"].notna(), "all_fstates_comb"].apply(lambda x: '_'.join(x))
    df_comb["all_fstates_comb"] = df_comb["all_fstates_comb"].str.replace("__", "_")
    df_comb["all_fstates_comb"] = df_comb["all_fstates_comb"].str.strip("_")
    df_comb["fstates_count"] = df_comb["all_fstates_comb"].str.count('_') + 1

    ## Is community located in BB (i.e. is any of the network members located in BB)
    df_comb["bb_located"] = df_comb["all_fstates_comb"].str.count("Brandenburg")

    ## Is community only located in BB
    df_comb["only_bb"] = 0
    df_comb.loc[(df_comb["bb_located"] == 1) & (df_comb["fstates_count"] == 1), "only_bb"] = 1

    ## Is any of the members agricultural active?
    df_comb["agric_related"] = 0
    df_comb.loc[(df_comb["agric_owners"] > 0) | (df_comb["number_agric_comps"] > 0), "agric_related"] = 1

    ## Is the community also internationally active (i.e. companies registered intern., any member located in there?)
    df_comb["international"] = 0
    df_comb.loc[(df_comb["intern_owners"] > 0) | (df_comb["number_intern_comps"] > 0), "international"] = 1

    ## Is the community only regional
    df_comb["regional"] = df_comb["regional"].map({True: 1, False: 0})
    df_comb.loc[(df_comb["regional"] == 1) & (df_comb["only_bb"] == 0), "regional"] = 0

    ## Is community single company
    df_comb.loc[df_comb["single_company"].isna(), "single_company"] = 0

    ## Just to account for all communities
    df_comb["new_category"] = df_comb["level_c_category"]
    df_comb.loc[df_comb["new_category"].isna(), "new_category"] = '99'

    ## Reclassify wrong combinations
    ## 1. All private persons that are agricultural likely have a company
    ## Alternatively: differentiate between agricultural and non-agricultural private people
    df_comb.loc[
        (df_comb["new_category"].isin(
            ['1_1', '1_1_1', '1_2', '1_9_1'])) &
        (df_comb["agric_related"] == 1), "new_category"] = '1a'
    ## 2. Correct all companies that were wrongly classified into single companies
    df_comb.loc[
        (df_comb["new_category"].isin(
            ['2_1', '2_2', '2_9_1', '2_9_2'])) &
        (df_comb["number_companies"] > 1), "new_category"] = '2_9_1'
    df_comb.loc[
        (df_comb["new_category"].isin(
            ['90'])) &
        (df_comb["number_companies"] > 1), "new_category"] = '2_9_1'
    df_comb.loc[(df_comb["new_category"] == '2_9_1') & (df_comb["agric_related"] == 1), "new_category"] = '2_9_1a'
    df_comb.loc[(df_comb["new_category"] == '2_9_1') & (df_comb["agric_related"] != 1), "new_category"] = '2_9_1n'
    ## 3. Merge all single companies that are not agricultural active
    df_comb.loc[
        (df_comb["new_category"].isin(
            ['2_1_1', '2_1_2', '2_1_3', '2_1_4', '2_1_5', '2_1_6', '2_1_7', '2_1_8', '2_1', '2_2', '2_2_2', '2_9_2'])) &
        (df_comb["number_companies"] < 2), "new_category"] = '2_1'
    df_comb.loc[
        (df_comb["new_category"].isin(
            ['90'])) &
        (df_comb["number_companies"] < 2), "new_category"] = '2_1'

    df_comb.loc[(df_comb["new_category"] == '2_2'), "new_category"] = '2_1'

    df_comb.loc[(df_comb["new_category"].isin(['2_1'])) & ((df_comb["agric_related"] == 1)), "new_category"] = '2_1a'
    df_comb.loc[(df_comb["new_category"].isin(['2_1'])) & (df_comb["agric_related"] != 1), "new_category"] = '2_1n'

    ## 4. Merge all remaining private entities of level1 == 1
    df_comb.loc[(df_comb["new_category"].str.slice(0, 1) == '1') & (df_comb["new_category"] != '1a'), "new_category"] = '1_0_0'
    ## 5. Merge all public entities of level1 == 1 that are not BVVG
    df_comb.loc[(df_comb["new_category"].str.slice(0, 1) == '5') & (df_comb["new_category"] != '5_1_1'), "new_category"] = '5_0_0'
    ## 6. Merge all church entities of level1 == 1
    df_comb.loc[(df_comb["new_category"].str.slice(0, 1) == '4'), "new_category"] = '4_0_0'
    ## 6. Merge all non-profit entities of level1 == 1
    df_comb.loc[(df_comb["new_category"].str.slice(0, 1) == '3'), "new_category"] = '3_0_0'

    mean_num_comp = df_comb.loc[df_comb["new_category"] == '2_9_1', "number_companies"].mean()

    ## Divide between smaller and bigger company networks
    # df_comb.loc[(df_comb["new_category"] == "2_9_1") & (df_comb["number_companies"] < mean_num_comp), "new_category"] = '2_9_1'
    # df_comb.loc[(df_comb["new_category"] == "2_9_1") & (df_comb["number_companies"] > mean_num_comp), "new_category"] = '2_9_2'

    df_comb["only_bb"] = df_comb["only_bb"].map({0: 'in more states active', 1: 'only in BB active'})
    df_comb["regional"] = df_comb["regional"].map({0: 'supraregional ownership', 1: 'regional ownership'})
    df_comb["agric_related"] = df_comb["agric_related"].map(
        {0: 'not agricultural', 1: 'agricultural'})
    df_comb["label"] = df_comb["only_bb"] + ' & ' + df_comb["regional"] + ' & ' + df_comb["agric_related"]

    print("\tNumber of owners after classification:", len(df_comb))
    print("\tUnique categories:", df_comb["new_category"].unique())

    ## Assess the categories
    print("\tNumber of owners per category")
    t = df_comb.groupby("new_category").agg(
        number_owners=pd.NamedAgg(comm_col, "count")
    ).reset_index()
    print(t)

    print("\tAssess whether all agricultural owners fall into agricultural classes")
    t = df_comb.loc[df_comb["agric_owners"] > 0, "new_category"].unique()
    print("\tCategories of owners, where agric_owners==1", t)
    t = df_comb.loc[df_comb["agric_related"] == "agricultural", "new_category"].unique()
    print("\tCategories of owners, where agric_related==agricultural", t)

    df_comb["new_category"] = df_comb["new_category"].map(
        {'5_0_0': "PUBLIC",
         # '2_9_2': "laCNET",
         '2_9_1a': "aCONETW",
         '2_9_1n': "nCONETW",
         '3_0_0': "NONPRO",
         '2_1n': "siCOMP",
         '2_1a': "a_siCOMP",
         '1_0_0': "noagPR",
         '1a': "agriPR",
         '4_0_0': "CHURCH",
         # '2_2': "COOPER",
         '5_1_1': "PUBLIC"
         })
    print("\tUnique categories after reclassification:", df_comb["new_category"].unique())

    df_comb["new_category_ext"] = df_comb["new_category"] + ' - ' + df_comb["label"]

    ## Get overview of combinations of characteristics
    print("\tWrite Classification out.")
    df_comb.to_csv(out_pth_classsification, index=False)

    # t = df_comb.loc[df_comb["new_category"] == 'siCOMP'].copy()

    print("\tAttach to ALKIS.")
    df_alk = pd.merge(df_alk, df_comb[[comm_col, "new_category", "new_category_ext"]], how="left", on=comm_col)
    print("\tUnique categories in new classification:", df_alk["new_category"].unique())
    # t = df_alk.loc[df_alk["new_category"].isna()].copy()

    print("\tAssess again whether all agricultural owners fall into agricultural classes")
    t = df_alk.loc[df_alk["agric"] > 0, "new_category"].unique()
    t2 = df_alk.loc[df_alk["agric"] > 0].copy()
    print("\tCategories of owners, where agric_owners==1", t)
    t = df_alk.loc[df_alk["comm_agri_related"] > 0, "new_category"].unique()
    print("\tCategories of owners, where comm_agri_related==1", t)
    t = df_alk.loc[(df_alk["agric"] > 0) & (df_alk["new_category"] == "siCOMP")].copy()
    t = df_alk.loc[(df_alk["agric"] > 0) & (df_alk["new_category"] == "nCONETW")].copy()

    # t = df_alk.loc[df_alk["new_category"] == 'siCOMP'].copy()
    # t = df_alk.loc[df_alk[f"community_{50}"].isin(
    #     ["2046_0", "1888_0", "1506_0", "1108_0", "1043_0", "1021_0", "1061_0"]), "new_category"].unique()
    # for key in ["2046_0", "1888_0", "1506_0", "1108_0", "1043_0", "1021_0", "1061_0"]:
    #     print(key, df_alk.loc[df_alk[f"community_{50}"] == key, "new_category"].unique())
    # print(t)

    print("\tWrite ALKIS out.")
    df_alk.to_csv(alkis_out_pth, sep=";", index=False)


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    # create_neighbor_dictionary_for_municipalities(
    #     municipality_shp_pth=BB_SHP,
    #     out_pth=MUNICIP_NEIGHBOR_DICT_PTH
    # )

    for threshold in [50]:
        print(f"Threshold {threshold}\n ")
        comm_col = f"community_{threshold}"

        # calculate_distances(
        #     owner_df_pth=OWNERS_W_THRESH_PTH.format(threshold),
        #     out_pth=OWNERS_W_THRESH_DIST_PTH.format(threshold)
        # )

        # unify_attribute_in_alkis(
        #     owners_pth=OWNERS_W_THRESH_DIST_PTH.format(threshold),
        #     owners_pth_cleaned=OWNERS_W_THRESH_CLEANED_PTH.format(threshold)
        # )

        # reduce_alkis_data_to_necessary_information(
        #     # alkis_iacs_inters_pth=ALKIS_IACS_PTH.format(threshold),
        #     owners_pth=OWNERS_W_THRESH_CLEANED_PTH.format(threshold),
        #     comm_col=comm_col,
        #     out_pth=OWNERS_W_THRESH_REDUCED_PTH.format(threshold)
        # )

        # derive_community_information_from_dafne_data(
        #     dafne_pth=COMMUNITY_DAFNE_CSV_W_THRESH.format(threshold),
        #     info_dict_pth=COMMUNITY_INFO_DICT_W_THRESH.format(threshold),
        #     comm_col=comm_col,
        #     net_branch_pth=NETW_BRANCH_PTH,
        #     stopword_lst_pth=STOPWORD_LIST,
        #     out_pth=COMMUNITY_INFO_FROM_DAFNE_W_THRESH.format(threshold)
        # )

        # classify_communities_from_info_dict(
        #     info_dict_pth=COMMUNITY_INFO_DICT_W_THRESH.format(threshold),
        #     alkis_reduced_pth=OWNERS_W_THRESH_REDUCED_PTH.format(threshold),
        #     netw_conn_pth=COMMUNITY_DAFNE_CSV_W_THRESH,
        #     comm_col=comm_col,
        #     owners_pth=OWNERS_W_THRESH_CLEANED_PTH.format(threshold),
        #     out_pth=OWNERS_W_THRESH_AND_LOC_PTH.format(threshold),
        #     curr_folder=CURR_FOLDER
        # )

        ## Ab hier wird es etwas messy. Da ich hier nur die intersection der IACS+ALKIS Daten nutze, charakterisiere ich nur
        ## die communities, die Land auf den IACS Flächen besitzen. Dadurch entstehen NAs in der "new_category"!
        get_characteristics_of_communities_from_alkis(
            alkis_iacs_municip_inters_pth=ALKIS_IACS_MUNICIP_PTH,
            owners_pth=OWNERS_W_THRESH_AND_LOC_PTH.format(threshold),
            comm_col=comm_col,
            municips_covered_by_comm_pth=MUNICIP_NEIGHBOR_DICT_PTH,
            out_pth=CHARACTERISTICS_OF_COMMS_PTH.format(threshold))

        combine_information_of_communities_and_classify(
            df_info_alkis_pth=CHARACTERISTICS_OF_COMMS_PTH.format(threshold),
            df_info_dafne_pth=COMMUNITY_INFO_FROM_DAFNE_W_THRESH.format(threshold),
            comm_col=comm_col,
            out_pth_classsification=COMMUNITY_CLASSIFICATION_PTH.format(threshold),
            owners_pth=OWNERS_W_THRESH_AND_LOC_PTH.format(threshold),
            alkis_out_pth=OWNERS_W_THRESH_AND_LOC_CLASSIFIED_PTH.format(threshold)
        )

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)


if __name__ == '__main__':
    main()
