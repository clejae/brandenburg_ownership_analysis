# Author:
# github repository:

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import os
import time
import pandas as pd
import geopandas as gpd
import shapely
import math
import warnings
from osgeo import gdal
from osgeo import ogr
import numpy as np
from collections import Counter
import json

## Project library
import helper_functions
import conc_meas_lib
# ------------------------------------------ USER VARIABLES ------------------------------------------------#
WD = os.path.curdir

## Input
OWNERS_W_THRESH_PTH = r"10_owner_network_classification\10_owners_stretched+comm_w_thr{0}-dist+cleaned+loc+class.csv"
GRID_CENTROIDS_PTH = r"00_data\vector\grids\square_grid_4km_centroids_BB.gpkg"
PARCELS_PTH = r"09_alkis_intersection_with_other_layers\alkis_iacs_inters.shp"

## Output
PARCELS_RAS_PTH = r"00_data\raster\land_parcels_ids.tiff"
CONC_MEASURES_MW_GRID_BUFFERS_PTH = r"11_ownership_concentration\mw_grid_buffers\mw_conc_meas-{0}.csv"
COUNTS_MW_GRID_BUFFER_PTH = r"11_ownership_concentration\mw_grid_buffers\mw_counts_categories_in_topx-{0}.csv"

# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#


def rasterize_shape(in_shp_pth, out_ras_pth, attribute, extent, res, no_data_val, gdal_dtype):
    """
    This function rasterizes a shapefile based on a provided attribute of the shapefile.
    :param in_shp_pth: Path to input shapefile. String.
    :param out_ras_pth: Path to output raster, including file name and ".shp". String.
    :param attribute: Attribute (i.e. field) of shapefile that should be rasterized. If attr is an integer,
    then only the geometries of the shape will be rasterized with the provided integer as the burn value.
    :param extent: List of extent of raster. Will be checked if it fits to provided resolution.
    [xmin, xmax, ymin, ymax]
    :param res: Resolution of raster in units of projection of input shapefile.
    :param no_data_val: No data value of raster.
    :gdal_dtype = gdal data type of raster.
    :return: Output raster will be written to specified location.
    """

    import math
    from osgeo import gdal
    from osgeo import ogr
    from osgeo import os

    ## Determine raster extent
    ## Reassuring, that extent and resolution fit together
    ## Assuming that upper left corner is correct (x_min, y_max)
    x_min = extent[0]
    x_max = extent[1]
    y_min = extent[2]
    y_max = extent[3]
    cols = math.ceil((x_max - x_min) / res)
    rows = math.ceil((y_max - y_min) / res)
    x_max = x_min + cols * res
    y_min = y_max - rows * res

    ## If input shape exists, then start the rasterization
    if os.path.exists(in_shp_pth):
        shp = ogr.Open(in_shp_pth, 0)  # 0=read only, 1=writeabel
        lyr = shp.GetLayer()

        #### Transform spatial reference of input shapefiles into projection of raster
        sr = lyr.GetSpatialRef()
        pr = sr.ExportToWkt()

        #### Create output raster
        target_ds = gdal.GetDriverByName('GTiff').Create(out_ras_pth, cols, rows, 1, gdal_dtype,
                                                         options=['COMPRESS=DEFLATE'])  # gdal.GDT_Int16)#
        target_ds.SetGeoTransform((x_min, res, 0, y_max, 0, -res))
        target_ds.SetProjection(pr)
        band = target_ds.GetRasterBand(1)
        band.Fill(no_data_val)
        band.SetNoDataValue(no_data_val)
        band.FlushCache()

        if isinstance(attribute, str):
            option_str = "ATTRIBUTE=" + attribute
            gdal.RasterizeLayer(target_ds, [1], lyr, options=[option_str])
        elif isinstance(attribute, int):
            gdal.RasterizeLayer(target_ds, [1], lyr, burn_values = [attribute])
        else:
            print("Provided attribute is not of type str or int.")

        del target_ds
    else:
        print(in_shp_pth, "doesn't exist.")


def rasterize_parcels(parcels_pth, ras_col, out_pth):
    """
    Rasterizes parcels or fields from shapefile
    :param parcels_fields_pth: Path to shapefile with information on fields/parcels
    :param ras_col: Column name for rasterization.
    :param out_pth: Output path for rasterized version.
    :return:
    """

    ## Prepare shapefile
    print("Read shapefile")
    shp = ogr.Open(parcels_pth, 0)
    lyr = shp.GetLayer()
    x_min_ext, x_max_ext, y_min_ext, y_max_ext = lyr.GetExtent()

    ## enlarge extent, so that buffers that overlap the actual extent still fall into the raster
    print("Define characteristics")
    x_min_ext -= 12500
    x_max_ext += 12500
    y_min_ext -= 12500
    y_max_ext += 12500

    res = 2.5
    no_data_val = -9999

    print("Rasterize.")
    rasterize_shape(parcels_pth, out_pth, ras_col, [x_min_ext, x_max_ext, y_min_ext, y_max_ext],
                    res, no_data_val, gdal_dtype=gdal.GDT_Int32)

    print("Done.")


def write_array_to_raster(in_array, out_path, gt, pr, no_data_value, type_code=None, options=['COMPRESS=DEFLATE', 'PREDICTOR=1']):
    """
    Writes an array to a tiff-raster. If no type code of output is given, it will be extracted from the input array.
    As default a deflate compression is used, but can be specified by the user.
    :param in_array: Input array
    :param out_path: Path of output raster
    :param gt: GeoTransfrom of output raster
    :param pr: Projection of output raster
    :param no_data_value: Value that should be recognized as the no data value
    :return: Writes an array to a raster file on the disc.
    """

    from osgeo import gdal
    from osgeo import gdal_array

    if type_code == None:
        type_code = gdal_array.NumericTypeCodeToGDALTypeCode(in_array.dtype)

    if len(in_array.shape) == 3:
        nbands_out = in_array.shape[0]
        x_res = in_array.shape[2]
        y_res = in_array.shape[1]

        out_ras = gdal.GetDriverByName('GTiff').Create(out_path, x_res, y_res, nbands_out, type_code, options=options)
        out_ras.SetGeoTransform(gt)
        out_ras.SetProjection(pr)

        for b in range(0, nbands_out):
            band = out_ras.GetRasterBand(b + 1)
            arr_out = in_array[b, :, :]
            band.WriteArray(arr_out)
            band.SetNoDataValue(no_data_value)
            band.FlushCache()

        del (out_ras)

    if len(in_array.shape) == 2:
        nbands_out = 1
        x_res = in_array.shape[1]
        y_res = in_array.shape[0]

        out_ras = gdal.GetDriverByName('GTiff').Create(out_path, x_res, y_res, nbands_out, type_code, options=options)
        out_ras.SetGeoTransform(gt)
        out_ras.SetProjection(pr)

        band = out_ras.GetRasterBand(1)
        band.WriteArray(in_array)
        band.SetNoDataValue(no_data_value)
        band.FlushCache()

        del (out_ras)

        # Conversion dictionary:
        # NP2GDAL_CONVERSION = {
        #     "uint8": 1,
        #     "int8": 1,
        #     "uint16": 2,
        #     "int16": 3,
        #     "uint32": 4,
        #     "int32": 5,
        #     "float32": 6,
        #     "float64": 7,
        #     "complex64": 10,
        #     "complex128": 11,
        # }

def prepare_iacs_farm_centroids_v1(iacs_pth, out_pth):

    print("Read IACS.")
    iacs = gpd.read_file(iacs_pth)
    iacs = iacs.loc[iacs["BTNR"].str.slice(0, 2) == "12"].copy()
    iacs["centroid"] = iacs.geometry.centroid
    iacs["x"] = iacs["centroid"].x
    iacs["y"] = iacs["centroid"].y
    ## Get farm centroids
    farms = iacs.groupby("BTNR")[["x", "y"]].mean().reset_index()
    farms["geometry"] = gpd.points_from_xy(farms['x'], farms['y'])

    total_area_lst = []
    buffer_lst = []
    share_lst = []
    max_radius_lst = []

    f = 4210
    farm_id = farms["BTNR"].iloc[f]
    for f, farm_id in enumerate(list(farms["BTNR"])):
        print(f+1, len(farms))
        iacs_sub = iacs.loc[iacs["BTNR"] == farm_id].copy()
        minx, miny, maxx, maxy = iacs_sub.geometry.total_bounds
        max_radius = max([maxx-minx, maxy-miny]) / 2
        farm_centroid = farms.loc[farms["BTNR"] == farm_id].copy()
        farm_centroid.crs = 25833

        total_area = iacs_sub.area.sum()
        share = 0
        buffer_size = 100
        while share < 0.9:
            farm = farm_centroid.copy()
            farm.geometry = farm.buffer(buffer_size)
            intersection = gpd.overlay(iacs_sub, farm, how="intersection", keep_geom_type=False,
                                       make_valid=True)
            intersection_area = intersection.area.sum()
            share = intersection_area / total_area
            # print(buffer_size, round(share, 2))
            buffer_size += 100
        buffer_size -= 100

        # farm_centroid.to_file(
        #     fr"C:\Users\IAMO\Documents\work_data\ownership_paper\11_ownership_concentration\temp\farm_centroid_{farm_id}.gpkg",
        #     driver="GPKG")
        # farm.to_file(rf"C:\Users\IAMO\Documents\work_data\ownership_paper\11_ownership_concentration\temp\farm_buffer_{farm_id}.gpkg",
        #              driver="GPKG")
        # iacs_sub.drop(columns=["centroid", "x", "y"], inplace=True)
        # iacs_sub.to_file(rf"C:\Users\IAMO\Documents\work_data\ownership_paper\11_ownership_concentration\temp\iacs_sub_{farm_id}.gpkg",
        #              driver="GPKG")

        total_area_lst.append(total_area)
        buffer_lst.append(buffer_size)
        share_lst.append(share)
        max_radius_lst.append(max_radius)

    farms["buffer_radius"] = buffer_lst
    farms["total_area"] = total_area_lst
    farms["share"] = share_lst
    farms["max_radius"] = max_radius_lst

    farms.csr = 25833
    farms.to_file(out_pth, driver="GPKG")

# iacs_pth = r"00_data\vector\IACS\IACS_BB_2020.shp"
# out_pth = r"11_ownership_concentration\farm_centroids.gpkg"
def prepare_iacs_farm_centroids_v2(iacs_pth, out_pth):

    print("Read IACS.")
    iacs = gpd.read_file(iacs_pth)
    iacs = iacs.loc[iacs["BTNR"].str.slice(0, 2) == "12"].copy()
    iacs["centroid"] = iacs.geometry.centroid
    iacs["x"] = iacs["centroid"].x
    iacs["y"] = iacs["centroid"].y
    ## Get farm centroids
    farms = iacs.groupby("BTNR")[["x", "y"]].mean().reset_index()
    farms["geometry"] = gpd.points_from_xy(farms['x'], farms['y'])

    buffer_lst = []

    # f = 1
    # farm_id = farms["BTNR"].iloc[f]
    for f, farm_id in enumerate(list(farms["BTNR"])):
        print(f+1, len(farms))
        iacs_sub = iacs.loc[iacs["BTNR"] == farm_id].copy()
        minx, miny, maxx, maxy = iacs_sub.geometry.total_bounds
        farm_centroid = farms.loc[farms["BTNR"] == farm_id].copy()
        farm_centroid.crs = 25833

        def coord_lister(geom):
            if geom.geom_type == 'MultiPolygon':
                coords = []
                for x in geom.geoms:
                    coords += list(x.exterior.coords)
            else:
                coords = list(geom.exterior.coords)
            return (coords)

        coords = []
        for row in iacs_sub.itertuples():
           coords += coord_lister(row.geometry)

        coords_df = gpd.points_from_xy([i[0] for i in coords], [i[1] for i in coords])
        buffer_size = max([farm_centroid.distance(point).iloc[0] for point in coords_df])

        farm = farm_centroid.copy()
        farm.geometry = farm.buffer(buffer_size)

        # farm_centroid.to_file(
        #     fr"C:\Users\IAMO\Documents\work_data\ownership_paper\11_ownership_concentration\temp\farm_centroid_{farm_id}.gpkg",
        #     driver="GPKG")
        # farm.to_file(rf"C:\Users\IAMO\Documents\work_data\ownership_paper\11_ownership_concentration\temp\farm_buffer_{farm_id}.gpkg",
        #              driver="GPKG")
        # iacs_sub.drop(columns=["centroid", "x", "y"], inplace=True)
        # iacs_sub.to_file(rf"C:\Users\IAMO\Documents\work_data\ownership_paper\11_ownership_concentration\temp\iacs_sub_{farm_id}.gpkg",
        #              driver="GPKG")

        buffer_lst.append(buffer_size)

    farms["buffer_radius"] = buffer_lst

    farms.csr = 25833
    farms.to_file(out_pth, driver="GPKG")


def calculate_concentration_around_purchases_on_raster(rastered_parcels_pth, grid_centroids_pth, owner_pth, owner_col,
                                                       out_pth_concs, out_pth_counts, parcel_id_col, point_id_col,
                                                       feature_count=0, buffer_radius=12000, buffer_col=None):

    """
    Calculates the ownerhip or LU concentration in a buffer around the land purchases. It uses a rasterized version of
    the parcels (ALKIS) or fields (IACS) and rasterizes the purchase buffer as well for a faster calculation as gpd
    intersection is too slow.
    :param rastered_parcels_pth: Path to rasterized parcels.
    :param grid_centroids_pth: Path to shapefile with points for which the concentration should be calculated in a specific buffer
    :param owner_pth: Path to df with owner information. Needed for ownerhip concentration calculation.
    :param buffer_radius: Radius of buffer in m.
    :param out_pth_concs: Output path to results with concentration measures.
    :param parcel_id_col: Column name of parcels IDs (probably OGC_FID).
    :param feature_count: Can be set for debugging. Subsets the land purchase shapefile.
    :return:
    """

    ## Read input and get additional information
    parcels_ras = gdal.Open(rastered_parcels_pth)
    input_arr = parcels_ras.ReadAsArray()
    gt = parcels_ras.GetGeoTransform()
    x_ref = gt[0]
    y_ref = gt[3]

    owner_df = pd.read_csv(owner_pth, sep=";")
    point_ds = ogr.Open(grid_centroids_pth)
    point_lyr = point_ds.GetLayer()
    sr = point_lyr.GetSpatialRef()

    # point_lyr.SetAttributeFilter('"xxx" = 172')
    # point_lyr.SetAttributeFilter('"xxx" IN (3593.0, 3594.0, 3595.0)')

    ## If feature count is specified, then use only subset. If not then use all features.
    if feature_count == 0:
        feature_count = point_lyr.GetFeatureCount()

    print(f"Calculate concentration measures for {feature_count} features.")

    ## List for results
    df_lst = []
    df_count_lst = []

    ## Loop over features
    for f, feat in enumerate(point_lyr):
        ## Check if current feature is in subset
        if f >= feature_count:
            break

        stime = time.time()

        pointid = feat.GetField(point_id_col)
        if isinstance(buffer_col, str):
            buffer_radius = feat.GetField(buffer_col)
            if buffer_radius > 15000:
                buffer_radius = 15000

        # if pointid in [3593, 3594]:
        #     print(pointid)
        #     print("")
        ## Get Buffer
        geom = feat.GetGeometryRef()
        buffer = geom.Buffer(buffer_radius)
        geom_wkt = buffer.ExportToWkt()

        ## Get buffer extent
        extent = buffer.GetEnvelope()
        x_min_ext = extent[0]
        x_max_ext = extent[1]
        y_min_ext = extent[2]
        y_max_ext = extent[3]

        resolution = gt[1]

        ## align coordinates to parcel/field raster
        dist_x = x_min_ext - x_ref
        steps_x = (math.floor(dist_x / resolution))
        x_min_ali = x_ref + steps_x * resolution  # - 30

        dist_x = x_max_ext - x_ref
        steps_x = (math.floor(dist_x / resolution))
        x_max_ali = x_ref + steps_x * resolution  # + 30

        dist_y = y_ref - y_min_ext
        steps_y = (math.floor(dist_y / resolution))
        y_min_ali = y_ref - steps_y * resolution  # - 30

        dist_y = y_ref - y_max_ext
        steps_y = (math.floor(dist_y / resolution))
        y_max_ali = y_ref - steps_y * resolution  # + 30

        # slice input raster array to common dimensions
        px_min = int((x_min_ali - gt[0]) / gt[1])
        px_max = int((x_max_ali - gt[0]) / gt[1])
        # raster coordinates count from S to N, but array count from Top to Bottom, thus pymax = ymin
        py_max = int((y_min_ali - gt[3]) / gt[5])
        py_min = int((y_max_ali - gt[3]) / gt[5])
        geom_arr = input_arr[py_min: py_max, px_min: px_max].copy()

        # create memory layer for rasterization
        driver_mem = ogr.GetDriverByName('Memory')
        ogr_ds = driver_mem.CreateDataSource('wrk')
        ogr_lyr = ogr_ds.CreateLayer('poly', srs=sr)
        feat_mem = ogr.Feature(ogr_lyr.GetLayerDefn())
        feat_mem.SetGeometryDirectly(ogr.Geometry(wkt=geom_wkt))
        ogr_lyr.CreateFeature(feat_mem)

        # rasterize geom
        col_sub = px_max - px_min
        row_sub = py_max - py_min
        step_size_x = gt[1]
        step_size_y = gt[5]
        gt_mem = (x_min_ali, step_size_x, 0, y_max_ali, 0, step_size_y)
        target_ds = gdal.GetDriverByName('MEM').Create('', col_sub, row_sub, 1, gdal.GDT_Byte)
        # target_ds = gdal.GetDriverByName('GTiff').Create(rf"data\raster\temp\buffer{feat.GetField('field_1')}.tiff", col_sub, row_sub, 1, gdal.GDT_Byte, options=['COMPRESS=DEFLATE'])
        target_ds.SetProjection(parcels_ras.GetProjection())
        target_ds.SetGeoTransform(gt_mem)
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(-9999)
        gdal.RasterizeLayer(target_ds, [1], ogr_lyr, burn_values=[1])

        ## Read buffer array
        buffer_arr = band.ReadAsArray()

        target_ds.FlushCache()
        del target_ds
        del band

        # write_array_to_raster(
        #     in_array=geom_arr,
        #     out_path=rf"data\raster\temp\in_arr_{feat.GetField('field_1')}.tiff",
        #     gt=gt_mem,
        #     pr=parcels_ras.GetProjection(),
        #     no_data_value=-9999,
        #     type_code=None,
        #     options=['COMPRESS=DEFLATE', 'PREDICTOR=1'])

        ## Mask the sliced parcel array with buffer array
        geom_arr[buffer_arr == 0] = -9999

        mtime = time.time()

        ############################## FIRST CODE VERSION - FASTER? ########################################
        ## Get no pixels per id parcel, i.e. the area
        unique_vals = np.unique(geom_arr, return_counts=True)

        ## Transform into df for concentration calculation
        df_parcels = pd.DataFrame.from_dict(data={parcel_id_col: unique_vals[0], "pix_count": unique_vals[1]})
        df_parcels["area"] = df_parcels["pix_count"] * (resolution * resolution) / 10000
        df_parcels[point_id_col] = pointid
        # df_parcels.sort_values(by="area", inplace=True)

        ## Th df needs to be filled with information on landowners
        df_comb = helper_functions.combine_parcels_with_owners(
            df_parcels[[point_id_col, parcel_id_col, "area"]],
            owner_df[[parcel_id_col, "community_50", "mother_company", "new_category"]],
            id_col=parcel_id_col
        )

        df_comb = df_comb.loc[df_comb[parcel_id_col] != -9999].copy()
        # df_comb.to_csv(rf"00_data\raster\temp\df_comb{feat.GetField('field_1')}_raster.csv")

        ## Calculate concentration for buffer
        df_conc = conc_meas_lib.calculate_concentration_measures_from_df(
            df=df_comb,
            target_unit_id_col=point_id_col,
            area_col='area',
            owner_col=owner_col,
            out_pth=None)
        df_lst.append(df_conc)

        ## Calculate counts of owners for buffer
        if not df_comb.empty:
            df_count = conc_meas_lib.get_count_of_owner_categories_per_spatial_unit(
                df=df_comb,
                target_unit_id_col=point_id_col,
                area_col='area',
                owner_col=owner_col,
                category_col="new_category",
                out_pth=None
            )
            df_count_lst.append(df_count)

        etime = time.time()
        print(f"{f + 1}/{feature_count} done - field_1: {pointid}.", "Rasterization time share",
              round((mtime - stime) / (etime - stime), 2))
    point_lyr.ResetReading()

    print("Write out.")
    df_conc = pd.concat(df_lst)
    df_conc.index = range(len(df_conc))
    df_conc.to_csv(out_pth_concs, index=False)

    df_count = pd.concat(df_count_lst)
    df_count.index = range(len(df_count))
    df_count.to_csv(out_pth_counts, index=False)

def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    ########################### OWNERSHIP CONCENTRATION MEASURES ###########################
    # rasterize_parcels(
    #     parcels_pth=PARCELS_PTH,
    #     ras_col="OGC_FID",
    #     out_pth=PARCELS_RAS_PTH)

    # prepare_iacs_farm_centroids_v1(
    #     iacs_pth=r"00_data\vector\IACS\IACS_BB_2020.shp",
    #     out_pth=r"11_ownership_concentration\farm_centroids.gpkg")

    prepare_iacs_farm_centroids_v2(
        iacs_pth=r"00_data\vector\IACS\IACS_BB_2020.shp",
        out_pth=r"11_ownership_concentration\farm_centroids_max_radius.gpkg")

    # threshold = 50
    # descr = f"mother_companies-comm_w_thr{threshold}-iacs_areas"
    # for radius in [4]:
    #     print(f"\n###########################\nCALCULATE CONCENTRATION MEASURES IN {radius}KM RADIUS.\n###########################\n")
    #     ########################### OWNERSHIP CONCENTRATION MEASURES ###########################
    #     calculate_concentration_around_purchases_on_raster(
    #         rastered_parcels_pth=PARCELS_RAS_PTH,
    #         grid_centroids_pth=GRID_CENTROIDS_PTH,
    #         owner_pth=OWNERS_W_THRESH_PTH.format(threshold),
    #         buffer_radius=radius*1000,
    #         parcel_id_col="OGC_FID",
    #         owner_col='mother_company',
    #         # feature_count=300,
    #         point_id_col="POLYID",
    #         out_pth_concs=CONC_MEASURES_MW_GRID_BUFFERS_PTH.format(descr + f"_{radius}km"),
    #         out_pth_counts=COUNTS_MW_GRID_BUFFER_PTH.format(descr + f"_{radius}km")
    #     )
    #
    # print(
    #     f"\n###########################\nCALCULATE CONCENTRATION MEASURES IN FLEXIBLE RADIUS.\n###########################\n")
    # descr = f"mother_companies-comm_w_thr{threshold}-iacs_areas-flexible_farm_buffers"
    # calculate_concentration_around_purchases_on_raster(
    #     rastered_parcels_pth=PARCELS_RAS_PTH,
    #     grid_centroids_pth=r"11_ownership_concentration\farm_centroids.gpkg",
    #     owner_pth=OWNERS_W_THRESH_PTH.format(threshold),
    #     buffer_col="buffer_radius",
    #     parcel_id_col="OGC_FID",
    #     owner_col='mother_company',
    #     point_id_col="BTNR",
    #     # feature_count=10,
    #     out_pth_concs=CONC_MEASURES_MW_GRID_BUFFERS_PTH.format(descr),
    #     out_pth_counts=COUNTS_MW_GRID_BUFFER_PTH.format(descr)
    # )
    #
    # print(
    #     f"\n###########################\nCALCULATE CONCENTRATION MEASURES IN FLEXIBLE RADIUS.\n###########################\n")
    # descr = f"mother_companies-comm_w_thr{threshold}-iacs_areas-12km_farm_buffers"
    # calculate_concentration_around_purchases_on_raster(
    #     rastered_parcels_pth=PARCELS_RAS_PTH,
    #     grid_centroids_pth=r"11_ownership_concentration\farm_centroids.gpkg",
    #     owner_pth=OWNERS_W_THRESH_PTH.format(threshold),
    #     buffer_radius=12000,
    #     parcel_id_col="OGC_FID",
    #     owner_col='mother_company',
    #     point_id_col="BTNR",
    #     # feature_count=10,
    #     out_pth_concs=CONC_MEASURES_MW_GRID_BUFFERS_PTH.format(descr),
    #     out_pth_counts=COUNTS_MW_GRID_BUFFER_PTH.format(descr)
    # )


    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)


if __name__ == '__main__':
    main()