# Author:
# github repository:

## ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import os
import time
import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib
import math
import seaborn as sns
import networkx as nx
import json
from community import community_louvain

from collections import defaultdict
import networkx as nx
import plotly.graph_objects as go
from typing import Any, List, Dict, Tuple, Union, Callable


import helper_functions
import plotting_lib
## ------------------------------------------ USER INPUT ------------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input:
NETW_COMM_PTH = r"08_network_analysis\network_connections_with_community.csv"
COMMUNITY_MAX_DIST_DICT_PTH = r"08_network_analysis\community{0}_max_dist_dict.json"
OWNERS_W_THRESH_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-dist+cleaned+loc+class.csv"
OWNER_DF_FOR_PLOTTING = r"14_paper_figures\14_owners_w_parcels+comm_w_thr{0}.csv"
ALK_IACS_INTERS_PTH = r"09_alkis_intersection_with_other_layers\alkis_iacs_inters.shp"
BB_SHP = r"00_data\vector\administrative\BB_municipalities.shp"
CONC_MEASURES_MW_GRID_VERSIONS_PTH = r"11_ownership_concentration\mw_grid_buffers\mw_conc_meas-{0}.csv"
#r"11_ownership_concentration\mw_grid\mw_mean_conc_meas-{0}.csv"
#,

## Output:
OUT_SHP_FOLDER_W_THRESH = r"15_additional_analysis\largest_owners\_vector"
OUT_FIG_FOLDER_W_THRESH = r"15_additional_analysis\largest_owners\_community_distribution"


## Other output paths are defined in the functions below in this script.


## ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#
def owner_origin_statistics():
    helper_functions.create_folder("15_additional_analysis\owner_origins")

    threshold = 50

    df = pd.read_csv(OWNER_DF_FOR_PLOTTING.format(threshold), sep=";")
    df_owner = df.groupby("mother_company").agg(
        fstate_mcomp=pd.NamedAgg("fstate_mcomp", "first"),
        area=pd.NamedAgg("area", "sum")
    ).reset_index()

    df_fstate = df_owner.groupby("fstate_mcomp").agg(
        num_owner=pd.NamedAgg("mother_company", "nunique"),
        area=pd.NamedAgg("area", "sum")
    ).reset_index()

    df_fstate["origin"] = df_fstate["fstate_mcomp"].map({
        "Ausland": "Foreign",
        "Baden-Württemberg": "Old Länder",
        "Bayern": "Old Länder",
        "Berlin": "Berlin-Brandenburg",
        "Brandenburg": "Berlin-Brandenburg",
        "Bremen": "Old Länder",
        "Hamburg": "Old Länder",
        "Hessen": "Old Länder",
        "Kägiswil": "Foreign",
        "Mecklenburg-Vorpommern": "New Länder",
        "Niedersachsen": "Old Länder",
        "Rheinland-Pfalz": "Old Länder",
        "Saarland": "Old Länder",
        "Sachsen": "New Länder",
        "Sachen-Anhalt": "New Länder",
        "Schleswig-Holstein": "Old Länder",
        "Thüringen": "New Länder",
        "Unbekannt": "Unbekannt",
        "Unbekant": "Unbekannt"
    })

    df_origin = df_fstate.groupby("origin").sum().reset_index()
    df_origin["share_owners"] = df_origin["num_owner"] / df_origin["num_owner"].sum() * 100
    df_origin["share_area"] = df_origin["area"] / df_origin["area"].sum() * 100

    df_origin.to_excel("15_additional_analysis\owner_origins\owner_origins.xlsx", index=False)


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
    elif df['conn_left'].str.contains(name).any():
        # Extract the values from comm_column where conn_left contains the name
        comms = df.loc[df['conn_left'].str.contains(name), comm_column].tolist()
        return comms
    elif df['conn_right'].str.contains(name).any():
        # Extract the values from comm_column where conn_left contains the name
        comms = df.loc[df['conn_right'].str.contains(name), comm_column].tolist()
        return comms
    else:
        print(f"{name} not found in df.")


def save_network_members_to_df(df, df_areas, community_column, community_id, out_pth):
    """
    Saves a df of all network members to out_pth.
    Args:
        df: Network dataframe with companies and community IDs.
        community_column: Column name for community IDs.
        community_id: ID of community for which the df should be generated.
        out_pth: Output path for df.

    Returns:

    """
    df = df.loc[df[community_column] == community_id].copy()
    df2 = df[["conn_left"]].copy()
    df2["isperson"] = "Unternehmen"
    df2.rename(columns={"conn_left": "member"}, inplace=True)
    df.rename(columns={"conn_right": "member"}, inplace=True)
    df2 = pd.concat([df2, df[["member", "isperson"]]])
    df2.drop_duplicates(subset="member", inplace=True)
    df2.sort_values(by="isperson", ascending=False, inplace=True)

    df3 = df_areas.loc[df_areas[community_column] == community_id].copy()
    df3 = df3.groupby("owner_clean").agg(area=pd.NamedAgg("area", "sum")).reset_index()

    with pd.ExcelWriter(out_pth, engine='openpyxl') as writer:
        # Write each DataFrame to a different sheet
        df.to_excel(writer, sheet_name='community_information', index=False)
        df2.to_excel(writer, sheet_name='unique_members', index=False)
        df3.to_excel(writer, sheet_name='members_with_land', index=False)


def create_graph_from_df(df):
    function_to_color = {"Aktionär":  "#ffeda0",  # gelb
                         "Aufsichtsrat":  "#0f0f0f",  # grau
                         "Geschäftsführender Direktor":  "#a1d99b",  # grün
                         "Geschäftsführer":  "#a1d99b",  # grün
                         "Gesellschafter":  "#feb24c",  # orange
                         "Kommanditaktionär":  "#ffeda0",  # gelb
                         "Kommanditist":  "#feb24c",  # orange
                         "Komplementär":  "#feb24c",  # orange
                         "Liquidator":  "#0f0f0f",  # grau
                         "mother company":  "#f03b20",  # darkorange
                         "persönlich haftender Gesellschafter":  "#feb24c",  # orange
                         "Prokurist":  "#0f0f0f",  # grau
                         "subsidary":  "#f03b20",  # darkorange
                         "Verwaltungsrat":  "#0f0f0f",  # grau
                         "Vorstand":  "#0f0f0f"}  # grau

    df["color"] = df["function"].map(function_to_color)
    df["conn_right"] = df["conn_right"].apply(lambda x: x.split('_')[0])
    df["conn_left"] = df["conn_left"].apply(lambda x: x.split('_')[0])

    Graph = nx.from_pandas_edgelist(df=df, source="conn_right", target="conn_left", edge_attr=["color", "function", "share"])

    return Graph


def draw_network(Graph, debug=False):

    ## All this was found at https://towardsdatascience.com/visualizing-networks-in-python-d70f4cbeb259
    import dash
    import visdcc
    import dash_core_components as dcc
    import dash_html_components as html
    from dash.dependencies import Input, Output


    #########################################################################################################

    ## create app
    app = dash.Dash()

    nodes = [{'id': node_name, 'label': node_name, 'shape': 'dot', 'size': 7} for i, node_name in enumerate(Graph.nodes)]
    ## create edges from df
    function_to_color = {"Aktionär": "yellow", # "#ffeda0",  # gelb
                         "Aufsichtsrat": "grey", # "#bdbdbd",  # grau
                         "Geschäftsführender Direktor": "green", #" "#a1d99b",  # grün
                         "Geschäftsführer": "green", # "#a1d99b",  # grün
                         "Gesellschafter": "orange", #" "#feb24c",  # orange
                         "Kommanditaktionär": "yellow", # "#ffeda0",  # gelb
                         "Kommanditist": "orange", #" "#feb24c",  # orange
                         "Komplementär": "orange", # "#feb24c",  # orange
                         "Liquidator": "grey", # "#bdbdbd",  # grau
                         "mother company": "orange", # "#f03b20",  # darkorange
                         "persönlich haftender Gesellschafter": "orange", #" "#feb24c",  # orange
                         "Prokurist": "grey", # "#bdbdbd",  # grau
                         "subsidary": "orange", #" "#f03b20",  # darkorange
                         "Verwaltungsrat": "grey", # "#bdbdbd",  # grau
                         "Vorstand": "grey"} #" "#bdbdbd"}  # grau
    edges = []
    for edge in Graph.edges:
        source, target = edge
        weight = Graph[source][target]["share"]
        function = Graph[source][target]["function"]
        color = function_to_color.get(function, "grey")
        edges.append({
            "id": source + '__' + target,
            "edge_label": weight,
            "from": source,
            "to": target,
            "width": weight / 10 + 1,
            "color": color
        })

    ## app.layout is what specifies how your Dash app looks
    app.layout = html.Div([
        visdcc.Network(id='net',
                       data={'nodes': nodes, 'edges': edges},
                       options=dict(height='800px', width='100%'))
    ])

    ## define callback
    ## allows your graph to be interactive and listen to user events such as click or selection
    @app.callback(
        Output('net', 'options'),
        [Input('color', 'value')])
    def myfun(selection):
        if selection:
            selected_node = selection.get('nodes', [None])[0]
            if selected_node:
                return {'nodes': {'borderWidthSelected': 5}}
        return {}

    app.run_server(debug=debug)
    # app.run_server(host='0.0.0.0', port=8050, debug=True)
    # app.run_server(host='0.0.0.0', debug=True)


def plot_concentration_histogramm_per_owner(alkis_iacs_inters_pth, df_conc_pth, owners_pth, bb_shp_pth, comm_col, community_lst, out_shp_folder, out_fig_folder):

    print("Combine parcels with owner data")
    gdf = gpd.read_file(alkis_iacs_inters_pth)
    gdf["area"] = gdf["geometry"].area
    gdf['area'] = gdf['area'] / 10000
    df = helper_functions.read_table_to_df(owners_pth)

    grid_res = 4
    version_4km = 1
    gdf_grid = gpd.read_file(rf"00_data\vector\grids\square_grid_{grid_res}km_v{version_4km:02d}_with_12km_POLYIDs.shp")
    df_conc = pd.read_csv(df_conc_pth)

    gdf_conc = pd.merge(gdf_grid, df_conc, how='inner', left_on='POLYID', right_on='id_sp_unit')

    df_alk = helper_functions.combine_parcels_with_owners(
        gdf[["OGC_FID", "EIGENTUEME", "BTNR", "area", "geometry"]],
        df[["OGC_FID", "mother_company", "owner_merge", "level_c_category", comm_col]])

    df_alk = gpd.GeoDataFrame(df_alk)

    out_lst = [["mother_company", "comm", "min", "mean", "max"]]

    for i, comm in enumerate(community_lst):
        sub = df_alk.loc[df_alk[comm_col] == comm].copy()

        gdf_conc_touch = gpd.sjoin(gdf_conc, sub, op='intersects')
        df_plt = gdf_conc.loc[gdf_conc["POLYID"].isin(gdf_conc_touch["POLYID"].unique())].copy()

        # helper_functions.create_folder(out_shp_folder)
        # out_shp_pth = f"{out_shp_folder}/community_{comm}_hexagons.shp"
        # if not os.path.exists(out_shp_pth):
        #     sub.to_file(out_shp_pth)

        mother_company = sub["mother_company"].iloc[0]

        shp_bb = gpd.read_file(bb_shp_pth)

        folder = f"{out_fig_folder}/{i:02}_community"

        helper_functions.create_folder(folder)
        out_fig_pth = f"{folder}/{i:02}_community_{comm}_hexagons.png"

        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=plotting_lib.cm2inch(12, 12))
        shp_bb.plot(edgecolor='none', facecolor='#bebebe', ax=axs, zorder=0)
        df_plt.plot(edgecolor='none', facecolor='#0000FF', ax=axs, zorder=1)
        fig.suptitle(f'{comm}_{mother_company}', fontsize=9)
        plt.savefig(out_fig_pth, dpi=300)
        plt.close()

        plotting_lib.histogramm(
            df=df_plt,
            col="cr1",
            out_pth=f"{folder}/{i:02}_community_{comm}_histogramm.png"
        )

        sub_lst = [mother_company, comm, df_plt["cr1"].min(), df_plt["cr1"].mean(), df_plt["cr1"].max()]
        out_lst.append(sub_lst)

    df_out = pd.DataFrame(out_lst)
    df_out.to_csv(f"{out_fig_folder}/cr1_statistics.csv")

def get_largest_owners_in_concentrated_areas(df_conc_pth, df_conc_pth_v, conc_col, thresh_min, thresh_max=-9999):

    ## Read grid shapefiles
    grid_res = 4
    version_4km = 1
    gdf_grid = gpd.read_file(rf"00_data\vector\grids\square_grid_{grid_res}km_v{version_4km:02d}_with_12km_POLYIDs.shp")

    ## Read mean concentration measures
    df_conc = pd.read_csv(df_conc_pth)

    ## Read concentration measures of grid version 5 to identify the largest owner per grid cell
    df_top_owner = pd.read_csv(df_conc_pth_v)

    ## read 4km to 12km polyid dictionary
    owner_names = []
    for version in range(1, 10):
        target_unit_id_col = f"v{version:02d}_POLYID"

        ## Combine dataframes
        gdf_conc = pd.merge(gdf_grid, df_conc, how='inner', left_on='POLYID', right_on='POLYID')
        gdf_conc = pd.merge(gdf_conc, df_top_owner[['id_sp_unit', 'owner1']], how='left', left_on=target_unit_id_col,
                            right_on="id_sp_unit")

        if thresh_max == -9999:
            thresh_max = gdf_conc[conc_col].max()

        gdf_conc2 = gdf_conc.loc[(gdf_conc[conc_col] > thresh_min) & (gdf_conc[conc_col] <= thresh_max)].copy()

        temp_lst = list(gdf_conc2["owner1"].unique())
        owner_names += temp_lst

    owner_names = [name for name in owner_names if type(name) == str]

    return set(list(owner_names))


def get_largest_share_of_specific_owner_from_mw_calculation(conc_of_grid_version_pth, threshold, owner_name, out_pth=None, descr=None):
    grid_res = 4

    res_lst = []
    for version in range(1, 10):
        id_col = f"v{version:02d}_POLYID"
        if not descr:
            pth = conc_of_grid_version_pth.format(grid_res, version, threshold)
        else:
            pth = conc_of_grid_version_pth.format(grid_res, version, descr)

        df = pd.read_csv(pth, sep=',')

        df_sub = df.loc[df['owner1'] == owner_name].copy()

        if not df_sub.empty:
            df_sub["version"] = f"v{version:02d}"
            res_lst.append(df_sub)

    if res_lst:
        df_res = pd.concat(res_lst)

        if out_pth:
            df_res.to_csv(out_pth)

        return df_res


def get_largest_share_of_specific_owner_from_mw_12km_buffers(conc_pth, threshold, owner_name, out_pth=None):
    grid_res = 4

    pth = conc_pth.format(grid_res, threshold)

    df = pd.read_csv(pth, sep=',')
    df_sub = df.loc[df['owner1'] == owner_name].copy()

    if out_pth:
        df_sub.to_csv(out_pth)

    if not df_sub.empty:
        return df_sub

## ------------------------------------------ RUN PROCESSES ---------------------------------------------------#
def main():
    s_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + s_time)
    os.chdir(WD)

    # owner_origin_statistics()

    #################################### Visualization of network areas and network compositions ########################
    df = pd.read_csv(NETW_COMM_PTH)
    thresh = 50
    # df_areas = pd.read_csv(OWNERS_W_THRESH_PTH.format(thresh), sep=";")
    helper_functions.create_folder(r"15_additional_analysis\largest_owners")

    with open(COMMUNITY_MAX_DIST_DICT_PTH.format(thresh)) as json_file:
        no_layers_dict = json.load(json_file)

    ## Plot histograms of grid cells in which the largest owners are active
    df_largest_owners = pd.read_csv(r"13_general_statistics\ALKIS_IACS_intersection\area_per_community_w_thr50_alkis_iacs.csv")
    comm_col = f"community_{thresh}"
    community_lst = list(df_largest_owners[comm_col][:30])
    owner_lst = list(df_largest_owners['mother_company'][:16])
    owner_lst.remove('unbekannt')
    descr = f"mother_companies-comm_w_thr50-iacs_areas"

    # plot_concentration_histogramm_per_owner(
    #     alkis_iacs_inters_pth=ALK_IACS_INTERS_PTH,
    #     df_conc_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH.format(descr),
    #     owners_pth=OWNERS_W_THRESH_PTH.format(thresh),
    #     bb_shp_pth=BB_SHP,
    #     comm_col=comm_col,
    #     community_lst=community_lst,
    #     out_shp_folder=OUT_SHP_FOLDER_W_THRESH,
    #     out_fig_folder=OUT_FIG_FOLDER_W_THRESH
    # )

    res_lst = []
    for i, owner in enumerate(owner_lst):
        res = get_largest_share_of_specific_owner_from_mw_12km_buffers(
            conc_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH.format(descr),
            threshold=thresh,
            owner_name=owner)
        if isinstance(res, pd.DataFrame):
            res["rank"] = i+1
            res_lst.append(res)

    res_df = pd.concat(res_lst)
    res_df.to_csv(fr"{OUT_FIG_FOLDER_W_THRESH}\polygons_largest_15_owners.csv", index=False)
    # summary = res_df[["owner1", "rank", "cr1", "cr4", "gini_coeff", "hhi"]].groupby(['owner1', 'rank']).describe().reset_index()
    # summary.sort_values(by='rank', inplace=True)
    # summary.to_csv(fr"{OUT_FIG_FOLDER_W_THRESH}\concentrations_largest_30_owners.csv", index=False)

    ## Get largest owners in most concentrated areas
    # descr = f"mother_companies-comm_w_thr50-iacs_areas"
    # owner_names = get_largest_owners_in_concentrated_areas(
    #     df_conc_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH.format(descr),
    #     df_conc_pth_v=rf"11_ownership_concentration\mw_grid\mw_conc_meas-grid_4km_v05-{descr}.csv",
    #     conc_col="cr1",
    #     thresh_min=30,
    #     thresh_max=70)
    #
    # print("Owners in areas with CR1 larger than 30:", owner_names)
    #
    # descr = f"mother_companies-comm_w_thr50-iacs_areas"
    # owner_names = get_largest_owners_in_concentrated_areas(
    #     df_conc_pth=CONC_MEASURES_MW_GRID_VERSIONS_PTH.format(descr),
    #     df_conc_pth_v=rf"11_ownership_concentration\mw_grid\mw_conc_meas-grid_4km_v05-{descr}.csv",
    #     conc_col="cr1",
    #     thresh_min=25,
    #     thresh_max=30)
    #
    # print("Owners in areas with CR1 between 25 amd 30:", owner_names)
    #
    # # ## Examples
    # # # 'lindhorst', 'gross neuendorf-letschin eg', 'bioboden genossenschaft eg', 'christian olearius_hamburg', 'frank gross'
    # # # 'karsten alexander maria laubrock', 'bernd schmidt-ankum_glebitzsch', 'janine stratmann_nordwestuckermark', 'antje winkelmann'
    # # # 'muenchener rueckversicherungs-gesellschaft aktiengesellschaft in muenchen'
    #
    # names = ['werner kotzte']
    #
    # for name in names:
    #     out_folder = fr"15_additional_analysis\largest_owners\{name}"
    #     helper_functions.create_folder(out_folder)
    #
    #     ## Get community numbers for company name and plot
    #     comm_numbers = get_community_id_of_company(
    #         df=df,
    #         name=name,
    #         comm_column="community_50"
    #     )
    #
    #     comm_numbers = list(set(comm_numbers))
    #     comm_numbers = [i for i in comm_numbers if not pd.isna(i)]
    #
    #     for comm_number in comm_numbers:
    #
    #         plot_network_graph_of_community_from_df(
    #             df=df,
    #             community_column="community_50",
    #             community_id=comm_number,
    #             out_pth=fr"{out_folder}\{name}_main_community_{comm_number}.png"
    #         )
    #
    #         save_network_members_to_df(
    #             df=df,
    #             df_areas=df_areas,
    #             community_column="community_50",
    #             community_id=comm_number,
    #             out_pth=fr"{out_folder}\{name}_main_community_{comm_number}.xlsx")
    #
    #         no_layers = no_layers_dict[comm_number]
    #         print(f"Number of layers in network: {no_layers}")
    #
    #         Graph = create_graph_from_df(df.loc[df["community_50"] == comm_number].copy())
    #         ## This does not work in debug mode of PyCharm
    #         draw_network(Graph=Graph)

    e_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + s_time)
    print("end: " + e_time)


if __name__ == '__main__':
    main()
