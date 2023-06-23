def lac(x, w=None, interval=False):
    """
    Coefficient of asymmetry in the Lorenz curve.
    LAC above 1 means that curve is steep towards right end, which means that large entities have high shares
    LAC below 1 means that curve is rather flat towards right end, which means a more equal distribution

    Adapted from R function in package ineq
    https://search.r-project.org/CRAN/refmans/ineq/html/Lasym.html
    :param x:
    :param w:
    :param interval:
    :return:
    """

    import numpy as np

    x = np.asarray(x)
    o = x.argsort() #np.argsort(x)
    x = x[o]

    if w is None:
        w = np.repeat([1], len(x))
    w = w[o]

    mu = np.mean(x * w)
    xlow = x < mu
    m = np.sum(w[xlow])
    n = np.sum(w)
    Lm = np.sum(w[xlow] * x[xlow])
    Ln = np.sum(w * x)

    xeq = x == mu

    if np.any(xeq):
        a = np.sum(w[xeq])
        Lma = np.sum(w[xlow | xeq] * x[xlow | xeq])
        Lac = (m / n + Lm / Ln, (m + a) / n + Lma / Ln)
        if not interval:
            Lac = np.mean(Lac)
    else:
        xm = np.max(x[xlow])
        xm1 = np.min(x[xlow==False])
        delta = (mu - xm) / (xm1 - xm)
        Lac = (m + delta) / n + (Lm + delta * xm1) / Ln

    return Lac


def gini(x, w=None):
    """
    Found @
    https://stackoverflow.com/questions/48999542/more-efficient-weighted-gini-coefficient-in-python/48999797#48999797

    :param x:
    :param w:
    :return:
    """

    import numpy as np

    # The rest of the code requires numpy arrays.
    x = np.asarray(x)
    if w is not None:
        w = np.asarray(w)
        sorted_indices = np.argsort(x)
        sorted_x = x[sorted_indices]
        sorted_w = w[sorted_indices]
        # Force float dtype to avoid overflows
        cumw = np.cumsum(sorted_w, dtype=float)
        cumxw = np.cumsum(sorted_x * sorted_w, dtype=float)
        return (np.sum(cumxw[1:] * cumw[:-1] - cumxw[:-1] * cumw[1:]) /
                (cumxw[-1] * cumw[-1]))
    else:
        sorted_x = np.sort(x)
        n = len(x)
        cumx = np.cumsum(sorted_x, dtype=float)
        # The above formula, with all weights equal to 1 simplifies to:
        return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n


def palma_v1(x):
    import numpy as np

    x = np.asarray(x)
    percentiles = np.percentile(x, np.arange(0, 100, 1))
    palma = np.sum(x[x > percentiles[99]]) / np.sum(x[x < percentiles[40]])
    return palma


def palma_v2(x):
    import numpy as np

    x = np.asarray(x)
    ventiles = np.percentile(x, np.arange(0, 100, 1))
    palma = np.sum(x[x > ventiles[99]]) / np.sum(x[x < ventiles[90]])
    return palma


def palma_v3(x):
    import numpy as np

    x = np.asarray(x)
    ventiles = np.percentile(x, np.arange(0, 100, 5))
    palma = np.sum(x[x > ventiles[19]]) / np.sum(x[x < ventiles[8]])
    return palma


def hhi_index(firm_sizes, weights=None):
    r"""Herfindahl–Hirschman Index
    A common measure of market concentration, defined as
    $$
    H = \sum_{i=1}^{N} s_i^2
    $$
    where $s_i$ is firm $i$'s market share in the industry of $N$ firms.
    Args:
        firm_sizes (np.ndarray): (n_firms,) array of firm sizes suchs as sales, used to compute market shares
        weights (np.ndarray): (n_firms,) array of the weights given to each firm's market share. Defaults to equal weights.
    !!! note
        If `weights` are provided, the HHI index is computed as:
        $$
        H = \sum_{i=1}^{N} s_i^2 \times w_i
        $$
        where $w_i$ is the weight given to the market share of firm $i$.
    Returns:
        float: HHI-index for the industry
    References:
        - [Wikipedia](https://en.wikipedia.org/wiki/Herfindahl%E2%80%93Hirschman_Index)
    Found at: https://github.com/mgao6767/frds/blob/master/frds/measures/func_hhi_index.py
    """
    import numpy as np

    if weights is None:
        weights = np.ones(firm_sizes.shape)
    mkt_shares = firm_sizes / np.sum(firm_sizes)
    return np.sum(np.square(mkt_shares) * weights) * 10000


def rosenbluth_index(x):
    ## Source: https://welt-der-bwl.de/Rosenbluth-Index
    import numpy as np

    x = np.array(x)

    ## sort descending
    x[::-1].sort()
    ## get ranks
    ranks = np.arange(1, len(x)+1, 1)
    ## calculate shares
    shares = x/sum(x)
    ## get sum of ranks*shares
    y = np.sum(ranks*shares)
    ## get rosenbluth index
    r_index = 1 / (2*y - 1)

    return r_index


def theil(x):
    import inequality
    return inequality.theil.Theil(x).T


def inequality_gini(x):
    import inequality
    return inequality.gini.Gini(x).g


def calculate_concentration_measures_from_df(df, target_unit_id_col, area_col, owner_col, owner_class_col=None, exclude_col=None, out_pth=None):
    """
    :param df: data frame of land owners, with polygons and attributes of polygons and owners
    :param target_unit_id_col: Attribute that assigns each polygon/parcel to a specific regions (e.g. municipality)
    :param area_col: Column holding area information.
    :param owner_col: Column with unique owner names.
    :param exclude_col: Dictionary that provides a column and the
    :param out_pth: Path for output.
    :return: Output dataframe.
    """

    import pandas as pd
    import numpy as np
    import inequality

    ## Calculate concentration measures per spatial unit (municipality or grid polygon)

    ## get unique IDs, remove nan
    unit_ids = pd.Series(df[target_unit_id_col].unique())
    unit_ids = unit_ids.dropna()

    out_dict = {
        "id_sp_unit": [],
        "total_area": [],
        "gini_coeff": [],
        "lac": [],
        "num_owners": [],
        "palma_v1": [],
        "palma_v2": [],
        # "num_parcels": [],
        "cr1": [],
        "cr3": [],
        "cr5": [],
        "hhi": [],
        "rosenbluth_index": [],
        "theil_index": [],
        "theil_bg": [], # betweem group
        "theil_wg": [], # within group
        "theil_bs": [], # between share
        "share1": [],
        "share2": [],
        "share3": [],
        "owner1": [],
        "owner2": [],
        "owner3": [],
        "owner4": [],
        "owner5": [],
        "share_p100": [], # Share of top 1%
        "share_p95_99": [], # Share of following 4%
        "share_v19": [], # Share of 2nd top 5%
        "share_m50": [], # Share of middle 50%
        "share_b40": [] # Share of bottom 40%
    }


    for unit_id in unit_ids:
        sub = df[df[target_unit_id_col] == unit_id].copy()

        if exclude_col:
            sub = sub[sub[exclude_col] != unit_id].copy()

        if len(sub) > 0:
            # num_parcels = len(sub)

            if owner_class_col:
                agg = sub[['area', owner_col, owner_class_col]].groupby([owner_col, owner_class_col]).sum().reset_index()
            else:
                agg = sub[['area', owner_col]].groupby([owner_col]).sum().reset_index()

            agg = agg.sort_values(by=[area_col])
            x_lorenz1 = np.array(agg[area_col]).cumsum() / agg[area_col].sum()
            x_lorenz1 = np.concatenate([np.array([0.0]), x_lorenz1])
            x1 = np.array(agg[area_col])
            g1 = gini(x1)

            lac1 = lac(x1)
            palma1 = palma_v1(x1)
            palma2 = palma_v2(x1)
            hhi = hhi_index(x1)
            rbi = rosenbluth_index(x1)
            theil_i = theil(x1)
            if owner_class_col:
                theil_dr = inequality.theil.TheilD(x1, agg[owner_class_col])
                theil_bg = theil_dr.bg[0]
                theil_wg = theil_dr.wg[0]
                theil_bs = theil_bg / theil_i
            else:
                theil_bg = np.nan
                theil_wg = np.nan
                theil_bs = np.nan

            ventiles = np.percentile(x1, np.arange(0, 100, 5))
            percentiles = np.percentile(x1, np.arange(0, 100, 1))
            share_p100 = np.sum(x1[x1 >= percentiles[99]]) / np.sum(x1) * 100
            share_p95_99 = np.sum(x1[(x1 >= percentiles[95]) & (x1 < percentiles[99])]) / np.sum(x1) * 100
            share_v19 = np.sum(x1[(x1 >= ventiles[18]) & (x1 < ventiles[19])]) / np.sum(x1) * 100
            share_m50 = np.sum(x1[(x1 >= ventiles[8]) & (x1 < ventiles[18])]) / np.sum(x1) * 100
            share_b40 = np.sum(x1[(x1 < ventiles[8])]) / np.sum(x1) * 100

            owner_num = len(agg[owner_col])

            owner_lst = agg[-5:][owner_col].tolist()
            owner_lst += ['None'] * (5 - len(owner_lst))
            owner5 = owner_lst[0]
            owner4 = owner_lst[1]
            owner3 = owner_lst[2]
            owner2 = owner_lst[3]
            owner1 = owner_lst[4]
            total_area = agg[area_col].sum()

            cr1 = round((agg[-1:]['area'].sum() / total_area) * 100, 2)
            cr3 = round((agg[-3:]['area'].sum() / total_area) * 100, 2)
            cr5 = round((agg[-5:]['area'].sum() / total_area) * 100, 2)

            share1 = round((agg['area'].iloc[-1] / total_area) * 100, 2)
            if len(agg) > 1:
                share2 = round((agg['area'].iloc[-2] / total_area) * 100, 2)
            else:
                share2 = 0
            if len(agg) > 2:
                share3 = round((agg['area'].iloc[-3] / total_area) * 100, 2)
            else:
                share3 = 0

            out_dict["id_sp_unit"].append(unit_id)
            out_dict["total_area"].append(total_area)
            out_dict["gini_coeff"].append(g1)
            out_dict["lac"].append(lac1)
            out_dict["palma_v1"].append(palma1)
            out_dict["palma_v2"].append(palma2)
            out_dict["num_owners"].append(owner_num)
            # out_dict["num_parcels"].append(num_parcels)
            out_dict["cr1"].append(cr1)
            out_dict["cr3"].append(cr3)
            out_dict["cr5"].append(cr5)
            out_dict["hhi"].append(hhi)
            out_dict["rosenbluth_index"].append(rbi)
            out_dict["theil_index"].append(theil_i)
            out_dict["theil_bg"].append(theil_bg)
            out_dict["theil_wg"].append(theil_wg)
            out_dict["theil_bs"].append(theil_bs)
            out_dict["share1"].append(share1)
            out_dict["share2"].append(share2)
            out_dict["share3"].append(share3)
            out_dict["owner1"].append(owner1)
            out_dict["owner2"].append(owner2)
            out_dict["owner3"].append(owner3)
            out_dict["owner4"].append(owner4)
            out_dict["owner5"].append(owner5)
            out_dict["share_p100"].append(share_p100)
            out_dict["share_p95_99"].append(share_p95_99)
            out_dict["share_v19"].append(share_v19)
            out_dict["share_m50"].append(share_m50)
            out_dict["share_b40"].append(share_b40)

    df_out = pd.DataFrame(out_dict)

    bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    df_out["bins_cr1"] = pd.cut(df_out.cr1, bins)
    df_out["bins_cr3"] = pd.cut(df_out.cr3, bins)
    df_out["bins_cr5"] = pd.cut(df_out.cr5, bins)

    if out_pth:
        df_out.to_csv(out_pth, index=False)

    return df_out


def get_share_of_owner_categories_per_spatial_unit(df, target_unit_id_col, area_col, owner_col, category_col,
                                                   category_sub=None, shares_lst=None, out_pth=None):
    """
    Calculates for a each spatial unit how much share of land each unique category of an owner classification has.
    This is done for the top 1% and 10% owners if not specified otherwise. A category subset can also be
    provided.
    :param df:
    :param target_unit_id_col:
    :param area_col:
    :param category_col:
    :param category_sub:
    :param shares_lst:
    :param out_pth:
    :return:
    """

    import pandas as pd

    if not shares_lst:
        shares_lst = [1, 10, 100]

    ## create column labels for output df
    if not category_sub:
        category_sub = df[category_col].unique()
    f1 = lambda x, y: f"share_{x}_sh{y}"
    f2 = lambda x, y: f"count_{x}_sh{y}"
    labels = [f(cat, share) for cat in category_sub for share in shares_lst for f in (f1, f2)]

    ## create out dictionary
    out_dict = {"id_sp_unit": []}
    for l in labels:
        out_dict[l] = []

    ## get unique names, remove nan
    sp_unit_ids = pd.Series(df[target_unit_id_col].unique())
    sp_unit_ids = sp_unit_ids.dropna()

    ## for each spatial unit
    for sp_unit in sp_unit_ids:
        sub = df[df[target_unit_id_col] == sp_unit].copy()

        ## if no land in spatial unit skip the polygon
        if sub.empty:
            continue

        out_dict["id_sp_unit"].append(sp_unit)

        def get_all_unique_attributes_in_list(group):
            out = list(set(list(group)))
            return out

        ## First aggregate by owner, take categories with
        agg = sub.groupby(owner_col).agg(
            area=pd.NamedAgg(column=area_col, aggfunc="sum"),
            category=pd.NamedAgg(column=category_col, aggfunc=get_all_unique_attributes_in_list)
        ).reset_index()
        agg["category"] = agg["category"].apply(lambda x: x[0])
        agg["share"] = agg["area"]/agg["area"].sum()

        ## get the percentile that each owner in this spatial unit falls into
        agg.sort_values(by=area_col, ascending=False, inplace=True)
        agg["count"] = range(1, len(agg) + 1)
        agg["percentiles"] = pd.cut(agg["count"], bins=100, labels=range(1, 101))

        ## for each share threshold calculate the share of each category
        for share in shares_lst:
            ## only use owners that are below the threshold
            agg_sub = agg.loc[agg["percentiles"] <= share].copy()

            ## if categories are provided, only use these
            agg_sub = agg_sub.loc[agg_sub["category"].isin(category_sub)].copy()

            ## it is possible that after the last subsetting the df is empty, in this case return for all categories 0s
            if agg_sub.empty:
                for cat in category_sub:
                    key_name_share = f"share_{cat}_sh{share}"
                    key_name_count = f"count_{cat}_sh{share}"
                    sh = 0
                    ct = 0
                    out_dict[key_name_share].append(sh)
                    out_dict[key_name_count].append(ct)
                continue

            ## if not empty then now aggregate by category
            cat_agg = agg_sub[["share", "category"]].groupby("category").agg(
                share=pd.NamedAgg(column="share", aggfunc="sum"),
                count=pd.NamedAgg(column="category", aggfunc="count")
            )

            ## convert to dict to more easily access the values
            agg_dict = cat_agg.to_dict('index')

            ## for each category retrieve the values, if not in dictionary return 0
            for cat in category_sub:
                key_name_share = f"share_{cat}_sh{share}"
                key_name_count = f"count_{cat}_sh{share}"

                if cat in agg_dict:
                    sh = agg_dict[cat]["share"]
                    ct = agg_dict[cat]["count"]
                else:
                    sh = 0
                    ct = 0

                out_dict[key_name_share].append(sh)
                out_dict[key_name_count].append(ct)

    ## reformat to df
    df_out = pd.DataFrame(out_dict)

    ## write out/return
    if out_pth:
        df_out.to_csv(out_pth, index=False)
    return df_out


def get_count_of_owner_categories_per_spatial_unit(df, target_unit_id_col, area_col, owner_col, category_col,
                                                   category_sub=None, topx_lst=None, out_pth=None):
    """
    Calculates for a each spatial unit how often each unique category of an owner classification occurs in the
    largest x owners. This is done for the largest 1, 3, 5 owners if not specified otherwise.
    A category subset can also be provided.
    :param df:
    :param target_unit_id_col:
    :param area_col:
    :param category_col:
    :param category_sub:
    :param topx_lst:
    :param out_pth:
    :return:
    """

    import pandas as pd
    import numpy as np

    if not topx_lst:
        topx_lst = [1, 3, 5]

    ## create column labels for output df
    if not category_sub:
        category_sub = df[category_col].unique()

    def get_all_unique_attributes_in_list(group):
        out = list(set(list(group)))
        return out

    sp_units = df[target_unit_id_col].unique()

    df_area = df[[area_col, owner_col, category_col, target_unit_id_col]].groupby([target_unit_id_col, owner_col]).agg(
        area=pd.NamedAgg(column=area_col, aggfunc="sum"),
        category=pd.NamedAgg(column=category_col, aggfunc=get_all_unique_attributes_in_list)
    ).reset_index()
    df_area.sort_values(by=[target_unit_id_col, "area"], ascending=False, inplace=True)
    df_area["rank"] = df_area.groupby(target_unit_id_col)["area"].rank(method="first", ascending=False)

    df_area["category"] = df_area["category"].apply(lambda x: x[0])

    topx = 5
    df_lst = []
    for topx in topx_lst:
        df_sub = df_area.loc[df_area["rank"] <= topx].copy()
        df_agg = df_sub.groupby([target_unit_id_col, "category"]).agg(
            count=pd.NamedAgg(column="mother_company", aggfunc="count")
        ).reset_index()
        df_agg = df_agg.loc[df_agg["category"].isin(category_sub)].copy()
        df_agg = pd.pivot(df_agg, columns= "category", index=target_unit_id_col, values="count").reset_index()
        df_agg.fillna(0, inplace=True)
        df_agg.columns = ["id_sp_unit"] + [f"count_{col}_top{topx}" for col in df_agg.columns[1:]]
        df_lst.append(df_agg)

    from functools import reduce
    df_out = reduce(lambda df1, df2: pd.merge(df1, df2, how="outer", on="id_sp_unit"), df_lst)
    df_out.fillna(0, inplace=True)

    sp_units = [sp_unit for sp_unit in list(sp_units) if sp_unit not in df_out["id_sp_unit"].tolist()]
    app_dict = {"id_sp_unit":sp_units}
    for col in df_out.columns[1:]:
        app_dict[col] = np.repeat(0, len(sp_units))
    df_app = pd.DataFrame(app_dict)

    df_out = pd.concat([df_out, df_app], axis=0)
    df_out.sort_values(by="id_sp_unit", inplace=True)

    ## write out/return
    if out_pth:
        df_out.to_csv(out_pth, index=False)
    return df_out

def compare_two_df_with_same_columns(pth_df1, pth_df2, descr1, descr2, join_col1, join_col2, cols, out_pth_diff, out_pth_descr):
    import pandas as pd

    df1 = pd.read_csv(pth_df1)
    df2 = pd.read_csv(pth_df2)

    df1.columns = [f"{col}_{descr1}" for col in list(df1.columns)]
    df2.columns = [f"{col}_{descr2}" for col in list(df2.columns)]

    df_comb = pd.merge(df1, df2, how="left", left_on=f"{join_col1}_{descr1}", right_on=f"{join_col2}_{descr2}")
    df_comb.drop(columns=[f"{join_col2}_{descr2}"], inplace=True)
    df_comb.rename(columns={f"{join_col1}_{descr1}": "id_sp_unit"}, inplace=True)

    for col in cols:
        df_comb[f"{col}_diff"] = round(df_comb[f"{col}_{descr2}"] - df_comb[f"{col}_{descr1}"], 3)
        df_comb[f"{col}_incr"] = round((df_comb[f"{col}_{descr2}"] / df_comb[f"{col}_{descr1}"]) - 1, 3)

    df_comb.to_csv(out_pth_diff, index=False)

    # df_count1 = df_comb.loc[df_comb["cr1_diff"] > 0].copy()
    # counts = df_count1.groupby("owner1_comp").agg(
    #     Mean_diff=("cr1_diff", np.mean),
    #     Min_diff=("cr1_diff", np.min),
    #     Max_diff=("cr1_diff", np.max),
    #     Count=("cr1_diff", np.count_nonzero),
    #     Municips=("id_sp_unit", list)
    # ).reset_index()
    # counts.sort_values(by="Count", inplace=True, ascending=False)
    # out_pth = fr"{out_tables_folder}\count_owners-munic-diff_mcomp_w_thr{threshold}-owner_merge.csv"
    # counts.to_csv(out_pth, index=False)

    df_descr = df_comb.describe(percentiles=[.05, .25, .50, .75, .95])
    df_descr.to_csv(out_pth_descr)
