# Clemens Jänicke
# github Repo: https://github.com/clejae

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import os
import time
import geopandas as gpd

## Project library
import helper_functions

# ------------------------------------------ USER INPUT ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"
os.chdir(WD)

## Input paths
PTH_IACS = r"00_data\vector\IACS\IACS_BB_2020.shp"
PTH_ALKIS_REDUCED = r"00_data\vector\ALKIS\v_eigentuemer_bb_reduced.shp"
MUNICIPALITY_PTH = r"00_data\vector\administrative\BB_municipalities.shp"

## Output paths
GRIDS_12KM_FOLDER = r"00_data\vector\grids\12km_grids"
TEST_GRID_SIZES_FOLDER = r"00_data\vector\grids\test_grid_sizes"
FINAL_GRID_PTH = r"00_data\vector\grids\square_grid_{0}km_v{1:02d}_with_{2}km_POLYIDs.shp"
ALK_MUNIC_INTERS_PTH = r"09_alkis_intersection_with_other_layers\alkis_munic_inters.shp"
ALK_MUNICIP_IACS_UNION_PTH = r"09_alkis_intersection_with_other_layers\alkis_munic_iacs_union.shp"

ALK_GRID_INTERS_PTH = r"09_alkis_intersection_with_other_layers\alkis_grid_{0}km_v{1:02d}_inters.shp"
ALK_GRID_IACS_UNION_PTH = r"09_alkis_intersection_with_other_layers\alkis_grid_{0}km_v{1:02d}_iacs_union.shp"

# ALK_IACS_INTERS_PTH = r"09_alkis_intersection_with_other_layers\alkis_iacs_inters.shp"
# ALK_IACS_MUNIC_INTERS_PTH = r"09_alkis_intersection_with_other_layers\alkis_iacs_munic_inters.shp"
# IACS_MUNIC_INTERS_PTH = r"00_data\vector\IACS\IACS_municip_BB_2020_klassifiziert_25832.shp"

# ------------------------------------------ LOAD DATA & PROCESSING ------------------------------------------#

def create_polygon_grid(grid_res, ref_shp_pth, out_folder, offsets=[0]):
    """
    Creates a polygon grid over the extent of a reference shapefile.
    Args:
        grid_res: x-y-extent of polygon grid.
        ref_shp_pth: Path to reference shapefile from which extent and projection will be sourced.
        offsets: Shifts the extent in -x and -y by the given offset. If more than 1 offset is provided, the grid is
        generate no.offsets^2 times, as all combinations of x and y offsets are generated.

    Returns:

    """

    import os
    import math
    from osgeo import ogr

    print(f"Create polygon grid(s). Save to {out_folder}")

    v = 1

    shp = ogr.Open(ref_shp_pth)
    lyr = shp.GetLayer()
    ## xmin, xmax, ymin, ymax
    extent = lyr.GetExtent()
    sr = lyr.GetSpatialRef()

    helper_functions.create_folder(out_folder)

    for x_offset in offsets:
        for y_offset in offsets:

            x_ext = extent[1] + 12000 - extent[0] + 12000
            y_ext = extent[3] + 12000 - extent[2] + 12000

            grid_col = math.ceil(x_ext / grid_res)
            grid_row = math.ceil(y_ext / grid_res)

            out_shp_name = rf'{out_folder}/square_grid_{int(grid_res / 1000)}km_v{v:02d}.shp'
            drv = ogr.GetDriverByName('ESRI Shapefile')

            if os.path.exists(out_shp_name):
                drv.DeleteDataSource(out_shp_name)

            lyr_name = f'{int(grid_res/1000)}km_grid'
            ds = drv.CreateDataSource(out_shp_name)
            out_polygons = ds.CreateLayer(lyr_name, sr, ogr.wkbPolygon)

            # add fields
            out_polygons.CreateField(ogr.FieldDefn("POLYID", ogr.OFTString))

            curr_xmin = extent[0] - x_offset
            for col in range(0, grid_col):
                curr_ymin = extent[2] - y_offset

                for row in range(0, grid_row):
                    ring = ogr.Geometry(ogr.wkbLinearRing)
                    ring.AddPoint(curr_xmin, curr_ymin)
                    ring.AddPoint(curr_xmin + grid_res, curr_ymin)
                    ring.AddPoint(curr_xmin + grid_res, curr_ymin + grid_res)
                    ring.AddPoint(curr_xmin, curr_ymin + grid_res)
                    ring.AddPoint(curr_xmin, curr_ymin)

                    poly = ogr.Geometry(ogr.wkbPolygon)
                    poly.AddGeometry(ring)
                    poly.CloseRings()

                    poly_id = '{0:04d}_{1:04d}'.format(col + 1, row + 1)

                    out_defn = out_polygons.GetLayerDefn()
                    out_feat = ogr.Feature(out_defn)
                    out_feat.SetGeometry(poly)
                    out_polygons.CreateFeature(out_feat)

                    out_feat.SetField(0, poly_id)
                    out_polygons.SetFeature(out_feat)

                    curr_ymin += grid_res

                curr_xmin += grid_res

            ds = None
            v += 1


def cut_alkis_shp_with_grid_or_iacs(alkis_pth, grid_pth, out_pth, keep_cols_grid=None, keep_cols_alkis=None,
                                    shape_parameters=None, overlay_option="intersection"):
    """
    Overlays (i.e. intersects or cuts) ALKIS data with a grid or another polygon layer such as the IACS data.
    :param alkis_pth: Path to shapefile with parcel geometries from ALKIS.
    :param grid_pth: Path to grid or other polygon layer that should be used for the overlay
    :param keep_cols_grid: Columns to keep from the grid layer.
    :param keep_cols_alkis: Columns to keept from ALKIS.
    :param shape_parameters: needs to be a dictionary with "thin_ratio" and "area_thresh" as keys and float as values to clean the intersected layer.
    :param overlay_option: Option to use for geopandas overlay. Default is intersection.
    :param out_pth: Output path to overlayed ALKIS-grid shapefile.
    :return:
    """

    ## Read input
    print(f"Read {alkis_pth}")
    gdf_alk = gpd.read_file(alkis_pth)
    print(f"Read {grid_pth}")
    gdf_grid = gpd.read_file(grid_pth)

    ## Reproject IACS to ALKIS if necessary
    gdf_alk_epsg = int(gdf_alk.crs.srs.split(':')[1])
    gdf_iacs_epsg = int(gdf_grid.crs.srs.split(':')[1])
    if gdf_alk_epsg != gdf_iacs_epsg:
        print(f'\tInput shapefiles do not have the same projection. Reproject shp2 to epsg:{gdf_alk_epsg}')
        try:
            gdf_grid = gdf_grid.to_crs(gdf_alk_epsg)
        except:
            print(f'\tReprojection from {gdf_alk_epsg} to {gdf_iacs_epsg} failed!')
            return

    ## Subset ALKIS to necessary columns if provided
    # keep_cols_alkis = ['OGC_FID', 'geometry']
    if keep_cols_alkis:
        gdf_alk = gdf_alk[keep_cols_alkis].copy()

    ## Subset GRID shape to necessary information if provided
    # keep_cols_grid = ['BTNR', 'CODE', 'CODE_BEZ', 'ID', 'ID_KTYP', 'geometry']
    if keep_cols_grid:
        gdf_grid = gdf_grid[keep_cols_grid].copy()
        ## ToDo: Check if OBJECTID is in output
    gdf_grid['OBJECTID'] = range(len(gdf_grid))

    ## Intersect both shapefiles
    print("\tOverlay with opion:", overlay_option)
    gdf = gpd.overlay(gdf_grid, gdf_alk, how=overlay_option, keep_geom_type=False, make_valid=True)

    ## Recalculate area
    gdf['area'] = gdf['geometry'].area

    ## Calculate shape parameters
    if shape_parameters:
        thin_ratio = shape_parameters["thin_ratio"]#
        area_thresh = shape_parameters["area_thresh"]#

        gdf['perimeter'] = gdf['geometry'].length
        gdf['thin_ratio'] = 4 * 3.14 * (gdf['area'] / (gdf['perimeter'] * gdf['perimeter']))

        ## Clean sliver polygons
        # ## ToDo: verify somehow the threshold for the thinness ratio
        # gdf_slivers = gdf.loc[gdf['thin_ratio'] <= 0.01]
        # gdf_slivers.to_file(INTERSECTION_SLIVERS_PTH)

        # gdf = gdf.loc[gdf['thin_ratio'] > 0.01]
        # gdf = gdf.loc[gdf['area'] > 0.0]

        gdf = gdf.loc[(gdf['area'] >= area_thresh) & (gdf['thin_ratio'] > thin_ratio)].copy()

    no_entries = len(gdf)
    gdf = gdf.loc[(gdf.geometry.type == 'POLYGON')].copy()
    no_entries_clean = len(gdf)

    print("\tNo. of entries with no polygon:", no_entries - no_entries_clean)

    ## Write to disct
    print(f"\tWrite intersection to disc {out_pth}")
    gdf.to_file(out_pth)


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    ################################ INTERSECTION WITH GRIDS ################################
    ## We use a moving window approach of calculating the concentration measures. In our moving window approach, the
    ## concentration measures will be calculated 9 times, each time with a different 12x12km windows that shifts around
    ## a middle point. For this analysis we created a 4x4km grid, and nine 12x12km grids. For each grid cell in the
    ## 4x4km, we save the information in which 12x12km grid cell it lands. This way, we only need to intersect the ALKIS
    ## data once with the 4x4km grid and then can group the 4x4km grid cells with the IDs of the 12x12km grid cell IDs.
    ## 12km are used because of findings by Plogman et al. 2022 (https://doi.org/10.1016/j.landusepol.2022.106036)
    target_grid_res = 4

    ## Create 4km grid
    # create_polygon_grid(
    #     grid_res=4000,
    #     ref_shp_pth=PTH_IACS,
    #     offsets=[0],
    #     out_folder=GRIDS_12KM_FOLDER
    # )
    #
    # ## Create nine 12km grids
    # create_polygon_grid(
    #     grid_res=12000,
    #     ref_shp_pth=PTH_IACS,
    #     offsets=[12000, 8000, 4000],
    #     out_folder=GRIDS_12KM_FOLDER
    # )

    ## We first intersect the ALKIS data with the Grids to assign each parcel to a grid cell (or parts of parcels if
    ## they intersect with multiple grid cells. Then we add the information who uses the land (IACS) via the overlay-
    ## option "union" to the intersections. Not all grid-parcel intersections overlay with the land user information.
    ## They carry no information on the land users. It is important to do it in this order and to keep the parcels with
    ## no information on the land users, as otherwise in later steps, not the entire dataset could be used (e.g. in
    ## 10_owner_networks_classification --> get_characteristics_of_communities_from_alkis)

    # ## Add 12km IDs to 4km Shapefile.
    # target_grid = gpd.read_file(f"{GRIDS_12KM_FOLDER}\square_grid_{target_grid_res}km_v01.shp")
    # out_grid = target_grid.copy()
    # for v in range(1, 10):
    #     larger_grid = gpd.read_file(f"{GRIDS_12KM_FOLDER}\square_grid_{12}km_v{v:02d}.shp")
    #     larger_grid.rename(columns={"POLYID": f"v{v:02d}_POLYID"}, inplace=True)
    #     out_grid = gpd.overlay(larger_grid, out_grid, how="intersection", keep_geom_type=False, make_valid=True)
    #     out_grid = out_grid.loc[out_grid["geometry"].astype(str).str.contains("POLYGON") == True].copy()
    #     out_grid.index = range(1, len(out_grid)+1)
    # out_grid = out_grid[["POLYID", "v01_POLYID", "v02_POLYID", "v03_POLYID", "v04_POLYID", "v05_POLYID", "v06_POLYID",
    #                      "v07_POLYID", "v08_POLYID", "v09_POLYID", "geometry"]]
    # out_grid.to_file(FINAL_GRID_PTH.format(target_grid_res, 1, 12))
    #
    # ## Intersect ALKIS with 4x4km Grid
    # cut_alkis_shp_with_grid_or_iacs(
    #     alkis_pth=PTH_ALKIS_REDUCED,
    #     grid_pth=f"{GRIDS_12KM_FOLDER}\square_grid_{target_grid_res}km_v01.shp",
    #     keep_cols_grid=['POLYID', 'geometry'],
    #     keep_cols_alkis=['OGC_FID', 'EIGENTUEME', 'geometry'],
    #     out_pth=ALK_GRID_INTERS_PTH.format(target_grid_res, 1),
    #     shape_parameters=None,
    #     overlay_option='intersection')

    ## Overlay-union ALKIS-Grid intersection with IACS
    cut_alkis_shp_with_grid_or_iacs(
        alkis_pth=ALK_GRID_INTERS_PTH.format(target_grid_res, 1),
        grid_pth=PTH_IACS,
        keep_cols_grid=['BTNR', 'ID', 'geometry'],
        keep_cols_alkis=['OGC_FID', 'EIGENTUEME', 'POLYID', 'geometry'],
        out_pth=ALK_GRID_IACS_UNION_PTH.format(target_grid_res, 1),
        shape_parameters=None,
        overlay_option='union'
    )

    ## Remove areas that were added from IACS Shapefile, because the polygons don't overlap perfectly.
    print("Remove areas without an OGC_FID, i.e. that they are not part of the ALKIS data.")
    union = gpd.read_file(ALK_GRID_IACS_UNION_PTH.format(target_grid_res, 1))
    union = union.loc[union["OGC_FID"].notna()].copy()

    ## Check for missing parcels and missing owners
    orig_alkis = gpd.read_file(PTH_ALKIS_REDUCED)

    uni_ids_full = orig_alkis["OGC_FID"].unique()
    uni_ids = union["OGC_FID"].unique()
    miss = set(uni_ids_full).difference(uni_ids)
    miss_owners = list(orig_alkis.loc[orig_alkis["OGC_FID"].isin(miss), "EIGENTUEME"].unique())

    uni_owners = union["EIGENTUEME"].unique()
    owners_missing = set(miss_owners).difference(uni_owners)

    print("No. of missing owners:", len(owners_missing))

    union.to_file(ALK_GRID_IACS_UNION_PTH.format(target_grid_res, 1))

    ################################ INTERSECTION WITH MUNICIPALITIES ################################
    ## 1. Intersect ALKIS with Municipalities
    cut_alkis_shp_with_grid_or_iacs(
        alkis_pth=PTH_ALKIS_REDUCED,
        grid_pth=MUNICIPALITY_PTH,
        keep_cols_grid=['GEN', 'RS', 'geometry'],
        keep_cols_alkis=['OGC_FID', 'EIGENTUEME', 'geometry'],
        out_pth=ALK_MUNIC_INTERS_PTH,
        shape_parameters=None,
        overlay_option='intersection')

    ## 2. Overlay-union ALKIS-MUNICIP intersection with IACS
    cut_alkis_shp_with_grid_or_iacs(
        alkis_pth=ALK_MUNIC_INTERS_PTH,
        grid_pth=PTH_IACS,
        keep_cols_grid=['BTNR', 'CODE', 'ID', 'geometry'],
        keep_cols_alkis=['OGC_FID', 'EIGENTUEME', 'GEN', 'RS', 'geometry'],
        out_pth=ALK_MUNICIP_IACS_UNION_PTH,
        shape_parameters=None,
        overlay_option='union'
    )

    ## Remove areas that were added from IACS Shapefile, because the polygons don't overlap perfectly.
    print("Remove areas without an OGC_FID, i.e. that they are not part of the ALKIS data.")
    union = gpd.read_file(ALK_MUNICIP_IACS_UNION_PTH)
    union = union.loc[union["OGC_FID"].notna()].copy()

    ## Check for missing parcels and missing owners
    orig_alkis = gpd.read_file(PTH_ALKIS_REDUCED)

    uni_ids_full = orig_alkis["OGC_FID"].unique()
    uni_ids = union["OGC_FID"].unique()
    miss = set(uni_ids_full).difference(uni_ids)
    miss_owners = list(orig_alkis.loc[orig_alkis["OGC_FID"].isin(miss), "EIGENTUEME"].unique())

    uni_owners = union["EIGENTUEME"].unique()
    owners_missing = set(miss_owners).difference(uni_owners)

    print("No. of missing owners:", len(owners_missing))

    union.to_file(ALK_MUNICIP_IACS_UNION_PTH)


    # ################################ TEST DIFFERENT GRID SIZES ################################
    # target_grid_res = 1000
    # create_polygon_grid(
    #     grid_res=target_grid_res,
    #     ref_shp_pth=PTH_IACS,
    #     offsets=[0],
    #     out_folder=TEST_GRID_SIZES_FOLDER
    # )
    #
    # ## Generate grids for test of different grid sizes.
    # for i in range(2, 21):
    #     grid_res = i * 1000
    #     create_polygon_grid(
    #         grid_res=grid_res,
    #         ref_shp_pth=PTH_IACS,
    #         offsets=[0],
    #         out_folder=TEST_GRID_SIZES_FOLDER
    #         )
    #
    # ## Intersect ALKIS with 4x4km Grid
    # cut_alkis_shp_with_grid_or_iacs(
    #     alkis_pth=PTH_ALKIS_REDUCED,
    #     grid_pth=f"{GRID_FOLDER}\square_grid_{target_grid_res}km_v01.shp",
    #     keep_cols_grid=['POLYID', 'geometry'],
    #     keep_cols_alkis=['OGC_FID', 'EIGENTUEME', 'geometry'],
    #     out_pth=ALK_GRID_INTERS_PTH.format(target_grid_res, 1),
    #     shape_parameters=None,
    #     overlay_option='intersection')
    #
    # ## Overlay-union ALKIS-Grid intersection with IACS
    # cut_alkis_shp_with_grid_or_iacs(
    #     alkis_pth=ALK_GRID_INTERS_PTH.format(target_grid_res, 1),
    #     grid_pth=PTH_IACS,
    #     keep_cols_grid=['BTNR', 'ID', 'geometry'],
    #     keep_cols_alkis=['OGC_FID', 'EIGENTUEME', 'POLYID', 'geometry'],
    #     out_pth=ALK_GRID_IACS_UNION_PTH.format(target_grid_res, 1),
    #     shape_parameters=None,
    #     overlay_option='union'
    # )
    #
    # ## Remove areas that were added from IACS Shapefile, because the polygons don't overlap perfectly.
    # print("Remove areas without an OGC_FID, i.e. that they are not part of the ALKIS data.")
    # union = gpd.read_file(ALK_GRID_IACS_UNION_PTH.format(target_grid_res, 1))
    # union = union.loc[union["OGC_FID"].notna()].copy()
    #
    # ## Check for missing parcels and missing owners
    # orig_alkis = gpd.read_file(PTH_ALKIS_REDUCED)
    #
    # uni_ids_full = orig_alkis["OGC_FID"].unique()
    # uni_ids = union["OGC_FID"].unique()
    # miss = set(uni_ids_full).difference(uni_ids)
    # miss_owners = list(orig_alkis.loc[orig_alkis["OGC_FID"].isin(miss), "EIGENTUEME"].unique())
    #
    # uni_owners = union["EIGENTUEME"].unique()
    # owners_missing = set(miss_owners).difference(uni_owners)
    #
    # print("No. of missing owners:", len(owners_missing))
    #
    # union.to_file(ALK_GRID_IACS_UNION_PTH.format(target_grid_res, 1))

    ################################ OLD WORKFLOW ################################
    # ## 1. Intersect ALKIS with IACS
    # cut_alkis_shp_with_grid_or_iacs(alkis_pth=PTH_ALKIS_REDUCED,
    #                                 grid_pth=PTH_IACS,
    #                                 keep_cols_grid=['BTNR', 'CODE', 'ID', 'ID_KTYP', 'geometry'],
    #                                 keep_cols_alkis=['OGC_FID',  'EIGENTUEME', 'geometry'],
    #                                 out_pth=ALK_IACS_INTERS_PTH,
    #                                 shape_parameters=None)

    # ## 2. Intersect ALKIS_IACS with Municipalities
    # cut_alkis_shp_with_grid_or_iacs(alkis_pth=ALK_IACS_INTERS_PTH,
    #                                 grid_pth=MUNICIPALITY_PTH,
    #                                 keep_cols_grid=['RS', 'GEN', 'geometry'],
    #                                 keep_cols_alkis=['OGC_FID', 'EIGENTUEME', 'BTNR', 'geometry'],
    #                                 out_pth=ALK_IACS_MUNIC_INTERS_PTH,
    #                                 shape_parameters=None)

    ## 3. Intersect ALKIS_IACS with Grid
    # grid_res = 4
    # for i in range(1, 2):
    #     grid_pth = fr"{GRID_FOLDER}\square_grid_{grid_res}km_v{i:02d}.shp"
    #     out_pth = ALK_IACS_GRID_INTERS_PTH.format(grid_res, i)
    #     cut_alkis_shp_with_grid_or_iacs(alkis_pth=ALK_IACS_INTERS_PTH,
    #                                     grid_pth=grid_pth,
    #                                     keep_cols_grid=['POLYID', 'geometry'],
    #                                     keep_cols_alkis=['OGC_FID', 'EIGENTUEME', 'BTNR', 'geometry'],
    #                                     out_pth=out_pth,
    #                                     shape_parameters=None)

    ## 4. Intersect ALKIS with Municipalities
    # cut_alkis_shp_with_grid_or_iacs(alkis_pth=PTH_ALKIS_REDUCED,
    #                                 grid_pth=MUNICIPALITY_PTH,
    #                                 keep_cols_grid=['GEN', 'RS', 'geometry'],
    #                                 keep_cols_alkis=['OGC_FID', 'EIGENTUEME', 'geometry'],
    #                                 out_pth=ALK_MUNIC_INTERS_PTH,
    #                                 shape_parameters=None)

    ## 5. Intersect IACS with Municipalities
    # cut_alkis_shp_with_grid_or_iacs(alkis_pth=PTH_IACS,
    #                                 grid_pth=MUNICIPALITY_PTH,
    #                                 keep_cols_grid=['GEN', 'RS', 'geometry'],
    #                                 out_pth=IACS_MUNIC_INTERS_PTH,
    #                                 shape_parameters=None)

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)

if __name__ == '__main__':
    main()

