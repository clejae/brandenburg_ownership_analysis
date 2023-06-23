# Author:
# github repository:

# ------------------------------------------ LOAD PACKAGES ---------------------------------------------------#
import os
import time

## Project functions
import helper_functions

# ------------------------------------------ USER VARIABLES ------------------------------------------------#
WD = r"C:\Users\IAMO\Documents\work_data\ownership_paper"

# ------------------------------------------ DEFINE FUNCTIONS ------------------------------------------------#


def main():
    stime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    os.chdir(WD)

    ## Create folders for entire project
    helper_functions.create_folder("00_data")
    helper_functions.create_folder("01_clean_owner_strings")
    helper_functions.create_folder("02_identify_unique_owner_address_combinations")
    helper_functions.create_folder("03_owner_name_classification")
    helper_functions.create_folder("04_owner_class_reclassification")
    helper_functions.create_folder("05_georeference_addresses")
    helper_functions.create_folder("06_prepare_dafne_search")
    helper_functions.create_folder("07_owner_name_cleaning")
    helper_functions.create_folder("08_network_analysis")
    helper_functions.create_folder("09_alkis_intersection_with_other_layers")
    helper_functions.create_folder("10_owner_network_classification")
    helper_functions.create_folder("11_ownership_concentration_calculation")
    helper_functions.create_folder("12_ownership_concentration_plotting")
    helper_functions.create_folder("13_general_statistics")
    helper_functions.create_folder("14_paper_figures")

    etime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("start: " + stime)
    print("end: " + etime)

if __name__ == '__main__':
    main()

