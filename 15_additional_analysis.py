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
from community import community_louvain

import helper_functions
import plotting_lib
## ------------------------------------------ USER INPUT ------------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input:
NETW_COMM_PTH = r"08_network_analysis\network_connections_with_community.csv"

## Output:
OWNER_DF_FOR_PLOTTING = r"14_paper_figures\14_owners_w_parcels+comm_w_thr{0}.csv"
OWNERS_W_THRESH_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-dist+cleaned+loc+class.csv"

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


## ------------------------------------------ RUN PROCESSES ---------------------------------------------------#
def main():
    s_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + s_time)
    os.chdir(WD)

    # owner_origin_statistics()

    #################################### Visualization of network areas and network compositions ########################
    df = pd.read_csv(NETW_COMM_PTH)
    thresh = 50
    df_areas = pd.read_csv(OWNERS_W_THRESH_PTH.format(thresh), sep=";")
    helper_functions.create_folder(r"15_additional_analysis\largest_owners")

    ## Examples
    names = ['lindhorst', 'gross neuendorf-letschin eg', 'bioboden genossenschaft eg', 'christian olearius_hamburg', 'frank gross']

    for name in names:
        out_folder = fr"15_additional_analysis\largest_owners\{name}"
        helper_functions.create_folder(out_folder)

        ## Get community numbers for company name and plot
        comm_numbers = get_community_id_of_company(
            df=df,
            name=name,
            comm_column="community_50"
        )
        comm_numbers = list(set(comm_numbers))
        comm_numbers = [i for i in comm_numbers if not pd.isna(i)]
        for comm_number in comm_numbers:
            plot_network_graph_of_community_from_df(
                df=df,
                community_column="community_50",
                community_id=comm_number,
                out_pth=fr"{out_folder}\{name}_main_community_{comm_number}.png"
            )
            save_network_members_to_df(
                df=df,
                df_areas=df_areas,
                community_column="community_50",
                community_id=comm_number,
                out_pth=fr"{out_folder}\{name}_main_community_{comm_number}.xlsx")



    e_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + s_time)
    print("end: " + e_time)


if __name__ == '__main__':
    main()
