## Data preparation 
The data were provided in the Norm-based Data Interface (NAS) format by the Landesvermessung und Geobasisinformation Brandenburg (LGB). They came as 44 separate zipped XML-files. These files were imported into a PostgreSQL database using PostgreSQL Version 11 and the norGIS-ALKIS tool (http://norbit.de/68/). Until this point, we followed the workflow of Müller et al. (2021), who also provide a detailed description of these steps (https://doi.org/10.22004/ag.econ.311013).

After creating the database, the parcel data with ownership information could be found under public --> v_eigentuemer. We exported the parcel data from the database using the algorithm "pgsql2shp" from BostonGIS, which can be downloaded from http://www.bostongis.com/pgsql2shp_shp2pgsql_quickguide_20.bqg. The following steps were taken:
	1. Open cmd
	2. Cd C:\Path\to\pgsql2shp.exe  
	3. Run:
	pgsql2shp -f "C:\Path\to\output_file\v_eigentuemer_bb_pgsql2shp.shp" -h localhost -u postgres -P password123 database_name public.v_eigentuemer 
The exported file was then used for further processing in Python.

## Description of original data
The shapefile had the following attributes:
OCG_FID = unique ID for parcels
EIGENTUEME = Owner information. For private people the names, birthdates and addresses were provided. For all others only names and addresses were provided. A single parcel could be owned by multiple different owners (i.e. shared ownership).
AMTLFLSFL = Area of parcel.
geomeotry = Geographic coordinates of parcel boundaries.
... and other columns that were not of interest to us.

## Data Preprocessing
The pre-processing and analysis of the data was done in several steps and mostly automatically. In some cases, however, manual input had to be made. Each major processing step was also divided into several sub-steps to allow repetition of the individual steps if necessary. We decided to save the results of all sub-steps as separate files. This made it possible to enter each processing step without going through the entire previous workflow again. However, this increased the storage space required.

### 1. Separate owner information from geometries and clean owner strings
There were several problems with the owner information. The same owner could appear with multiple addresses throughout the dataset. Names, dates of birth and addresses could be spelled differently due to spelling errors and different abbreviations, and an owner could own both individual parcels and parcels with common ownership. Therefore, we had to clean up the owner information. To reduce loading times, we separated the geometries and the owner information into a shapefile that contained only the geographical information and a csv file that contained only the owner information. In the csv file we separated the names of the owners (+birth dates in the case of private individuals) from the addresses.

### 2. Stretch dataframe and clean addresses
An owner could appear in the dataframe as a single owner of a parcel and as one owner of several owners of a parcel with common ownership. Furthermore, this owner could appear in different spellings. These problems made it necessary to standardise the different spellings (which was done later). To this end, we extended the data set by separating all owners of parcels with common ownership. Each of them received a new data entry in the extended dataset. After extending the dataset, we cleaned up all addresses.

### 3. Preliminary owner name classification
We used the list of all stretched owners to classify them into a set of 29 preliminary classes. These classes helped in the further processing steps to perform various operations on the different classes. The preliminary classification was based on a classification of business forms (https://www.praxis-agrar.de/betrieb/recht/rechtsformen-landwirtschaftlicher-unternehmen) and we extended it iteratively. An overview of the classes can be found in "00_data\classification\class_catalogue.xlsx". The classification was done with the help of code words. The order of classification was important because some code words could occur for different classes.

### 4. Improve the classification and create unique identifiers
The preliminary classification was then transformed into an improved hierarchical classification, which can also be found in "00_data\classification\class_catalogue.xlsx". In addition, we created a unique identifier for each owner by removing typical but superfluous words and characters (e.g. spaces, dots, "mit sitz in" etc.) from the name. For individuals, we used the family name and address to create a unique identifier. We used fuzzy matching to identify the same owner who appeared with slightly different spellings (e.g. due to misspellings or different abbreviations) to standardise their identifier. In some cases, the unique identifier was listed with the addresses due to different spellings, etc. We have chosen to assign the most frequent address to all occurrences of the identifier. 

### 5. Georeference the addresses
In order to be able to calculate the distances of the owners to their parcels later on, we georeferenced the addresses of all owners. The georeferencing was done in different ways (Google API, Nominatim, matching with OSM data, fuzzy matching of addresses) and the scripts for this are quite messy and will likely not work as they are provided in this repository. The reason for choosing different methods was that we didn't have enough money to georeference all addresses with the Google API, although it is probably the most accurate. Also, we didn't want to put too much load on the OSM server as they don't want that. Besides, it also took a long time with their API. So we downloaded all German OSM data to create a shapefile with all cities and villages in Germany and their postcodes. With the help of this shapefile, we also georeferenced a number of addresses at least to the city centre. All remaining addresses that were not previously georeferenced were also matched using a fuzzy matching approach of address names. We therefore recommend carrying out an analysis of distances only on the basis of distance classes.

### 6. Prepare search of companies and other organizations in DAFNE
To identify all possible relationships between private individuals, companies and other organisations, we used the DAFNE company database. The database contains all information about the corporate structures of companies registered in Germany. We used the names of private companies and non-profit organisations, associations, foundations and others for a search in DAFNE. In DAFNE we have researched all possible information on boards of directors, supervisory boards, managing directors, subsidiaries and shareholders of companies. 