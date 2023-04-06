############################################ String cleaning ###########################################################

def replace_characters_in_string(text, char_lst, replace_string):
    """
    Replaces multiple characters in a string with another string.
    :param text:
    :param char_lst:
    :param replace_string:
    :return:
    """
    ## check if string is text
    if type(text) != str:
        text = str(text)

    for char in char_lst:
        text = text.replace(char, replace_string)

    return text


def remove_search_term_and_followup_in_substrings(text, search_term):
    """
    Removes all words after a search term from a sub part of a text. Sub parts are indicated by commas.
    E.g. search_term = "geb.":
     Maria Mueller geb. Schmidt, Hauptstraße 1, 11111 Stadt --> Maria Mueller, Hauptstraße 1, 11111 Stadt
    :param text: Input string.
    :param search_term: Term after which all words will be deleted.
    :return:
    """
    lst = text.split(',')

    out_lst = []
    for s, sub, in enumerate(lst):
        sub_lst = sub.split(' ')

        ## if search term in sub part, then don't keep
        if search_term not in sub_lst:
            sub = sub.strip()
            out_lst.append(sub)
        else:
            pass

    out_text = ', '.join(out_lst)

    return out_text


def remove_search_term_and_followup(text, search_term):
    """
    Removes all words after search term including the search term.
    :param text: Input string.
    :param search_term: Term after which all words will be deleted.
    :return:
    """

    ## e.g. remove all OT extensions (e.g. 'OT Manchow')

    if search_term in text:
        start_index = text.find(search_term)
        out_text = text[:start_index]
    else:
        out_text = text
    return out_text


def remove_address(text):
    """
    Removes addresses from the input text with regex expressions.
    :param text: Input string.
    :return:
    """

    import re
    lst = text.split(',')
    out_lst = []

    ## Define regex expressions
    ## find 4- or 5-digit numbers and a subsequent word
    rgx1 = r'\d{4,5} \b\w+\b'
    ## find just 5-digit numbers
    rgx2 = r'\d{5}'
    ## find street names with regex expression found at (for Straße 1, Straße 1a, Straße 1 a)
    ## https://stackoverflow.com/questions/5014537/regex-parse-streetname-number/5014624
    ## similar, but not quite working rgx3 = '(?:[^ ]+ ){0,5}\d{1,4} [a-f]?$'
    rgx3 = '^(.+)\s(\d+(\s*[^\d\s]+)*)$'
    ## regex expression found at (but adapted, for Straß1)
    ## https://stackoverflow.com/questions/24132763/regex-how-many-strings-start-with-letters-and-end-with-numbers/24132850
    rgx4 = '[A-Za-z]+.*\d{1,4}$'

    for sub in lst:
        if (not re.search(rgx1, sub)) and (
                not re.search(rgx2, sub) and (not re.search(rgx3, sub)) and (not re.search(rgx4, sub))):
            sub = sub.strip()
            # print(sub)
            out_lst.append(sub)
        else:
            pass
    out_text = ', '.join(out_lst)
    return out_text


def identify_owners(text, return_count=False):
    """
    Identify owner names in a text. Return count of owners, optionally.
    :param text: Input string.
    :param return_count: True or False.
    :return:
    """
    ## check if text is not empty
    if (text == None) or (text == ''):
        text = 'unbekannt'

    ## clean text a little bit
    text = text.lower()
    text = text.strip(',')
    text = text.strip()

    ## split at breaks
    lst = text.split('\n')

    ## remove addresses and 'born as names'
    owner_lst = []
    for item in lst:
        item = remove_search_term_and_followup_in_substrings(item, 'geb.')
        item = remove_address(item)

        ## if there is an asterisk in the item, then the part after the asterisk is an adress that is cut off
        ## (e.g. 'lack', 'wolfgang', '* 1999-99-99', 'bah'] --> bah = bahnhofstraße
        ## this part should be removed
        cut_lst = item.split(',')
        if any("*" in item for item in cut_lst):
            ind = [i for i, item in enumerate(cut_lst) if '*' in item][0]
            cut_lst = cut_lst[:ind + 1]
        item = ', '.join(cut_lst)

        owner_lst.append(item)

    ## check if last owner, which is sometimes a string that is cut off,
    ## can be found in any of the previous items
    last_item = owner_lst[-1]
    if any(last_item in prev_item for prev_item in owner_lst[:-1]):
        owner_lst.remove(last_item)
    else:
        pass

    owner_lst = [item.strip() for item in owner_lst]

    ## derive unique owners or number of owners
    unique_owners = list(set(owner_lst))
    out = ' | '.join(unique_owners)
    out = out.replace('  ', ' ')
    if return_count == True:
        out = len(unique_owners)

    return out


def get_birthdate(str):
    """
    Get the birthdate in a string for cases where it is indicated by a asterisk,
    :param str: Input string.
    :return: String with birthdate.
    """

    str_lst = str.split('*')

    if len(str_lst) > 1:
        birthdate = str_lst[1]
        if ',' in birthdate:
            birthdate = birthdate.split(',')[0]
    else:
        birthdate = ''

    return birthdate

############################################ Address and location identication #########################################

def find_addresses_in_string(owners, owner_string):
    """
    Find all addresses in an owner string, by identifying first the owner names and then determine their addresses.
    :param owners: Names of owners.
    :param owner_string: String with owner names and addresses.
    :return:
    """
    import re
    # print(owners)
    # print(owner_string)

    owner_lst = owners.split('|')
    owner_lst = [item.strip() for item in owner_lst]

    ## Define regex expressions
    ## find 4- or 5-digit numbers and a subsequent word
    rgx1 = r'\d{4,5} \b\w+\b'
    ## find just 5-digit numbers
    rgx2 = r'\d{5}'
    ## find street names with regex expression found at (for Straße 1, Straße 1a, Straße 1 a)
    ## https://stackoverflow.com/questions/5014537/regex-parse-streetname-number/5014624
    ## similar, but not quite working rgx3 = '(?:[^ ]+ ){0,5}\d{1,4} [a-f]?$'
    rgx3 = '^(.+)\s(\d+(\s*[^\d\s]+)*)$'
    ## regex expression found at (but adapted, for Straß1)
    ## https://stackoverflow.com/questions/24132763/regex-how-many-strings-start-with-letters-and-end-with-numbers/24132850
    rgx4 = '[A-Za-z]+.*\d{1,4}$'

    ## get length of  the complete owner string so that it can be searched till the end
    length = len(owner_string)

    addresses_all = []
    # owner= owner_lst[0]
    for owner in owner_lst:
        # print(owner)
        addresses_owner = []
        indices = []

        ## find start and end indices of current owner name + birthdate in complete owner string
        ## replace special characters, so that re is not confused
        owner_search = owner.replace('*', '_')
        owner_search = owner_search.replace('(', '_')
        owner_search = owner_search.replace(')', '_')
        owner_search = owner_search.replace('/', '_')
        owner_search = owner_search.replace('+', '_')
        # owner_search = owner_search.replace('.', '_')
        # owner_search = owner_search.replace('-', '_')

        owner_string_search = owner_string.replace('*', '_')
        owner_string_search = owner_string_search.replace('(', '_')
        owner_string_search = owner_string_search.replace(')', '_')
        owner_string_search = owner_string_search.replace('/', '_')
        owner_string_search = owner_string_search.replace('+', '_')
        # owner_string_search = owner_string_search.replace('.', '_')
        # owner_string_search = owner_string_search.replace('-', '_')

        for match in re.finditer(owner_search, owner_string_search):
            indices.append((match.start(), match.end()))
        # print(indices)

        ## loop through indices and check if the next two sub parts of the owner string (defined by ',')
        ## are a street+housenumber and a postal code + city name
        # i=0
        for i in range(len(indices)):
            owner_string_sub = owner_string[indices[i][1]:length]
            owner_string_sub = owner_string_sub.strip(',')
            owner_string_sub = owner_string_sub.split(',')
            owner_string_sub = owner_string_sub[:2]

            street_num = ''
            # sub = owner_string_sub[0]
            for sub in owner_string_sub:
                if street_num == '':
                    if re.search(rgx3, sub):
                        street_num = sub
                if street_num == '':
                    if re.search(rgx4, sub):
                        street_num = sub

            city_plz = ''
            for sub in owner_string_sub:
                if city_plz == '':
                    if re.search(rgx1, sub):
                        city_plz = sub
                if city_plz == '':
                    if re.search(rgx2, sub):
                        city_plz = sub

            address = street_num + ',' + city_plz
            address = address.strip(',')
            address = address.strip(' ')
            if len(address) > 1:
                addresses_owner.append(address)
            else:
                ## ToDo:
                addresses_owner.append('unkown 1, 00000 unkown')

        addresses_owner = [item.strip() for item in addresses_owner]
        addresses_owner = list(set(addresses_owner))
        addresses_owner = '_'.join(addresses_owner)
        addresses_all.append(addresses_owner)

    addresses_all = [item.strip() for item in addresses_all]
    addresses_all = ' | '.join(addresses_all)
    # print(addresses_all)

    return addresses_all


def identify_plz(text):
    import re

    text = text.lower()

    ## identify street names and postal codes + city name
    lst = text.split(',')
    pc = ''
    for l, sub in enumerate(lst):
        ## look for a postcal code + city name, if there hasn't been found any yet
        if pc == '':
            ## look  for 5-digit number and a subsequent word
            if re.search(r'\d{5} \b\w+\b', sub):
                pc = re.findall(r'\d{5}', sub)[0]
        if pc == '':
            ## if there was no 5-digit number, try a 4-digit number + word (e.g. Switzerland has 4-digit postal codes)
            if re.search(r'\d{4} \b\w+\b', sub):
                pc = re.findall(r'\d{4}', sub)[0]
        if pc == '':
            ## if there was no 5-digit number + city name, try only a 5-digit number
            if re.search(r'\d{5}', sub):
                pc = re.findall(r'\d{5}', sub)[0]

    pc = pc.strip()

    return pc


def identify_street(text):
    import re

    text = text.lower()
    text = text.replace('\n', ',')

    ## identify street names and postal codes + city name
    lst = text.split(',')
    street = ''
    for l, sub in enumerate(lst):
        ## only look for street name if there hasn't been found any yet
        ## look for words with subseqent 1-, 2-, 3- or 4-digit numbers
        if street == '':
            if re.search('(?:[^ ]+ ){0,5}\d{1,4}$', sub):
                street = sub
        if street == '':
            if re.search('(?:[^ ]+ ){0,5}\d{1,4} [a-f]$', sub):
                street = sub
        if street == '':
            if re.search('(?:[^ ]+ ){0,5}\d{1,4}[a-f]$', sub):
                street = sub

    street = street.strip()

    return street


def identify_city(text):
    import re
    text = text.lower()

    ## identify street names and postal codes + city name
    lst = text.split(',')
    city = ''
    for l, sub in enumerate(lst):
        ## look for a postcal code + city name, if there hasn't been found any yet
        if city == '':
            ## look  for 5-digit number and a subsequent word
            if re.search(r'\d{5} \b\w+\b', sub):
                city = re.findall(r'\D{1,100}', sub)
                city = ' '.join(city)
        if city == '':
            ## if there was no 5-digit number, try a 4-digit number + word (e.g. Switzerland has 4-digit postal codes)
            if re.search(r'\d{4} \b\w+\b', sub):
                city = re.findall(r'\D{1,100}', sub)
                city = ' '.join(city)

    city = city.strip()

    return city


def clean_city_text(text, search_terms, num_words = 2):
    ## e.g. remove all OT extensions (e.g. 'OT Manchow')
    if text != None:

        str_lst = text.split(" ")
        for search_term in search_terms:
            if search_term in str_lst:
                i = str_lst.index(search_term)
                sub_words = str_lst[i:i + num_words]
                sub_words = ' '.join(sub_words)
            else:
                sub_words = ''

            text = text.replace(sub_words, '')
            text = text.strip()

    else:
        pass

    return text


def clean_address(address_str):
    import re
    ## clean
    address_str = address_str.replace('ß', 'ss')
    address_str = address_str.replace('ä', 'ae')
    address_str = address_str.replace('ö', 'oe')
    address_str = address_str.replace('ü', 'ue')
    address_str = address_str.replace('.', '')
    address_str = address_str.replace('-', ' ')
    address_str = address_str.replace(':', '')

    address_str = address_str.replace('strasse', 'str')

    addresses = address_str.split('_')
    addresses = [item.strip() for item in addresses]
    addresses = [address for address in addresses if '00000' not in address]
    addresses = [address for address in addresses if 'unbekannt' not in address]
    lst = []

    # address = addresses[1]
    for address in addresses:

        street = identify_street(address)
        if street != '':
            num = re.search('\d{1,4}\s?[a-f]?$', street)
            number = street[int(num.span()[0]): int(num.span()[1])]
            number = number.replace(' ', '')
            if number[0] == '0':
                number = number[1:]
            else:
                pass
            street = street[:int(num.span()[0])]
            street = street.strip()

            street = street + ' ' + number + ', '

        plz = identify_plz(address)

        city = identify_city(address)
        city = city.split('/')[0]
        city = clean_city_text(city, search_terms=['bei', 'ot'])

        if (plz + city) != '':
            clean_add = street + plz + ' ' + city
            lst.append(clean_add)

    lst = set(lst)
    if len(lst) > 1:
        sorted_lst = sorted(lst, key=len)
        lst = []

        for i, item in enumerate(sorted_lst[:-1]):
            if not any(item.replace(' ', '') in s_item.replace(' ', '') for s_item in sorted_lst[i+1:]):
                lst.append(item)
            else:
                pass
        lst.append(sorted_lst[-1])

    lst = '_'.join(lst)

    return lst


def remove_street(owner_str, lookup_lst):
    """
    Remove street names from a string with help of a look-up list of key words.
    :param owner_str: Input string
    :param lookup_lst: Look up list of key words.
    :return:
    """
    if owner_str != None:

        str_lst = owner_str.split(" ")

        if len(str_lst) > 1:
            last = str_lst[-1]
            for search_term in lookup_lst:
                if search_term in last:
                    owner_str = owner_str.replace(last, '')
                    break

    return owner_str


def remove_address_part(text, delimiter, address_code_words):
    """
    Remove an address part in a text with a code word list.
    :param text:  Input text
    :param delimiter: Delimiter separating addresses from other parts.
    :param address_code_words: Code word list.
    :return:
    """

    if delimiter in text:
        ## search for delimiter and get possible address part (i.e the last part of the text)
        text_lst = text.split(delimiter)
        address_part = text_lst[-1]

        ## check if any of the code words is in address part
        ## if so, remove the part from the text
        for address_code in address_code_words:
            if address_code in address_part:
                out_text = ','.join(text_lst[:-1])
                break
            else:
                out_text = text
    else:
        out_text = text

    return out_text

############################################ Fuzzy matching of names and addresses #####################################


def jaro_address(address_str, thresh):
    import jaro
    """
    Check if two addresses, separated by "_" are basically the same with the Jaro Winkler 
    metric and a provided threshold.
    :param address_str: String of addresses.
    :param thresh: Threshold above which the addresses are the same.
    :return: 
    """
    addresses = address_str.split('_')
    addresses = [item.strip() for item in addresses]
    addresses = sorted(addresses, key=len)

    address1 = addresses[0]
    address2 = addresses[1]

    # address1 = 'bahnhofstr 26, 16352 basdorf'#'borkumstr 2, 13189 berlin'#'schoenerlinder dorfstr 19, 16348 wandlitz'
    # address2 = 'bahnhofstr 26, 16348 wandlitz'#'hauptstr 120, 16352 berlin'#'schoenerlinder dorfstr 19, 16352 wandlitz'

    jw_value = jaro.jaro_winkler_metric(address1, address2)
    # j_value = SequenceMatcher(None, word, station).ratio()
    if jw_value > thresh:
        out = address2
    else:
        # out = '-'
        out = address_str

    return out


def fuzzy_match_token_set_ratio(df, look_col, new_col, freq_col, thresh, search_range):
    from fuzzywuzzy import fuzz
    """
    Look for the same owner with similar spellings with a forward looking search algorithm on an ordered list of 
    owner names. Match them with fuzzy match token set ratio. Changes the input dataframe.
    :param df: Input dataframe
    :param look_col: Column on which search is based on
    :param new_col: Column where the matched owner name will be written to
    :param freq_col: Column in which the frequency will be saved in.
    :param thresh: fuzzy match token set ratio treshhold.
    :param search_range: Length of forward looking window.
    :return:
    """
    for i in range(len(df)):
        num_rows = len(df)
        add = search_range
        if i > (num_rows - search_range):
            add = num_rows - i
        df_sub = df[i:i + add].copy()
        df_sub['index_old'] = df_sub.index
        curr_owner = df_sub[look_col][df_sub['index_old'] == i].iloc[0]
        other_owners = list(df_sub[look_col][df_sub['index_old'] != i])

        lst = []
        ind_lst = [i]
        for j, oth_owner in enumerate(other_owners):
            match_value = fuzz.token_set_ratio(curr_owner, oth_owner)
            if (match_value > thresh):  # & (match_value <= 98): #.98,.99
                # print(oth_owner)
                lst.append(oth_owner)
                ind_lst.append(i + (j + 1))

        df_sub2 = df.iloc[ind_lst].copy()
        corr_owner = df_sub2[look_col][df_sub2[freq_col] == df_sub2[freq_col].max()].iloc[0]
        df.loc[ind_lst, new_col] = corr_owner


############################################ General utility ###########################################################

def create_folder(directory):
    """
    Tries to create a folder at the specified location. Path should already exist (excluding the new folder).
    If folder already exists, nothing will happen.
    :param directory: Path including new folder.
    :return: Creates a new folder at the specified location.
    """

    import os
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def get_indices_of_item_in_list(in_list, search_item):
    """
    Searches an item in a list and returns all indeces where it occurs.
    :param in_list: A list.
    :param search_item: The item to search.
    :return: Returns a list of indeces, if item not in list, then list empty.
    """

    if search_item in in_list:
        indices = [i for i, item in enumerate(in_list) if item == search_item]
    else:
        indices = []

    return indices


def print_red(skk):
    print("\033[91m {}\033[00m".format(skk))


def get_most_frequent_item(lst):
    from collections import Counter
    data = Counter(lst)
    return data.most_common(1)[0][0]


def most_frequent(in_list):
    """
    Similar to get_most_frequent_item, just implemented in another script. Thus, I don't dare to replace one of them,
    because I'm afraid something unexpected might happen.
    :param in_list:
    :return:
    """
    return max(set(in_list), key=in_list.count)


def merge_list_of_dfs(df_lst, merge_col):
    from functools import reduce
    import pandas as pd
    df = reduce(lambda df1, df2: pd.merge(df1, df2, on=merge_col, how="left"), df_lst)
    return df


def read_table_to_df(in_pth):
    import os
    import pandas as pd

    ext = in_pth.split('.')[-1]

    if ext == 'xlsx' or ext == 'xls':
        df = pd.read_excel(in_pth)
        return df

    if ext == 'csv':
        try:
            df = pd.read_csv(in_pth, sep=',')
        except:
            df = pd.read_csv(in_pth, sep=';')
        return df

    else:
        print("Please provide path with extension xlsx, xls or csv.")


class MyDict(dict):
    """
    Very usefull dictionary if you want to map stuff on a column but the dictionary doesn't have all values of
    that column as keys.
    Does not work when values of one column are used to update values of another column.
    """
    def __missing__(self, key):
        return key


############################################ Georeferencing ############################################################
def address_to_coordinates_google(address, api_key):
    import requests

    params = {
        'key': api_key,
        'address': address
    }

    base_url = 'https://maps.googleapis.com/maps/api/geocode/json?'

    response = requests.get(base_url, params=params)
    response = response.json()

    if response['status'] == 'OK':
        geometry = response['results'][0]['geometry']
        lat = geometry['location']['lat']
        lon = geometry['location']['lng']
        point = 'POINT ({0} {1})'.format(lon, lat)
    else:
        point = None

    return point


def address_to_coordinates_nominatim(address, out_pth=None):
    from geopy.geocoders import Nominatim
    from geopy.extra.rate_limiter import RateLimiter

    addresses = address.split('_')

    point_lst = []
    addr_lst = []
    for addr in addresses:
        addr = addr.replace('oe', 'ö')
        addr = addr.replace('ae', 'ä')
        addr = addr.replace('ue', 'ü')

        geolocator = Nominatim(user_agent="pdlen583ngkdlrz")
        geocode = RateLimiter(geolocator.geocode, min_delay_seconds=0.1)
        location = geocode(addr)

        if location:
            lat = location.latitude
            lon = location.longitude
            point = 'POINT ({0} {1})'.format(lon, lat)
            point_lst.append(point)
            addr_lst.append(addr)
        else:
            point = None

    if point_lst:
        point = point_lst[0]
        addr = addr_lst[0]
    else:
        point = None

    if out_pth:
        with open(out_pth, "a", encoding='ISO-8859-1') as file:
            file.write(f"{address};{point};{addr}\n")

    return point


def transform_point(in_point, in_sr, out_sr):
    from osgeo import osr
    from osgeo import ogr

    ## Define Coordinate transformation
    source = osr.SpatialReference()
    source.ImportFromEPSG(in_sr)  ## The geocoding provides coordinates in WGS 84
    target = osr.SpatialReference()
    target.ImportFromEPSG(out_sr)  ## ALKIS and IACS are stored in ETRS89 / UTM zone 32
    transform = osr.CoordinateTransformation(source, target)

    point = in_point.replace('POINT (', '')
    point = point.replace(')', '')
    point = point.split(' ')
    lon_orig, lat_orig = float(point[0]), float(point[1])

    ## transform to EPSG 25832. AddPoint order lat lon:
    out_point = ogr.Geometry(ogr.wkbPoint)
    out_point.AddPoint(lat_orig, lon_orig)
    out_point.Transform(transform)
    lat_out = out_point.GetY()
    lon_out = out_point.GetX()
    out_point = 'POINT ({0} {1})'.format(lon_out, lat_out)

    return out_point


def wkt_point_distance(wkt1, wkt2):
    """
    WKTs must be in EPSG:4326
    :param wkt1:
    :param wkt2:
    :return:
    """
    import shapely
    from geopy.distance import distance

    if wkt1 == None or wkt2 == None or type(wkt2) == float or type(wkt1) == float:
        dist = None
    else:

        ## extract points
        if type(wkt1) is shapely.geometry.point.Point:
            point1 = wkt1.x, wkt1.y
        elif type(wkt1) is str:
            point1 = wkt1.replace('POINT (', '')
            point1 = point1.replace(')', '')
            point1 = point1.split(' ')
            point1 = float(point1[0]), float(point1[1])

        if type(wkt2) is shapely.geometry.point.Point:
            point2 = wkt2.x, wkt2.y
        elif type(wkt2) is str:
            point2 = wkt2.replace('POINT (', '')
            point2 = point2.replace(')', '')
            point2 = point2.split(' ')
            point2 = (float(point2[0]), float(point2[1]))

        for val in point1:
            if not -90 < val < 90:
                point1 = (0.0, 0.0)
            else:
                pass

        for val in point2:
            if not -90 < val < 90:
                point2 = (0.0, 0.0)
            else:
                pass

        if point1 == (0.0, 0.0) or point2 == (0.0, 0.0):
            dist = None
        else:
            dist = distance(point1, point2).m

    return dist


############################################ Special ALKIS functions ###################################################

def combine_parcels_with_owners(parcels, owner_df, id_col="OGC_FID"):
    """
    Combines parcels with owner information. Corrects area of parcels with number of owners, i.e. if a parcel
    has multiple owners, it will be duplicated in the merging process as all owners are listed separately. To
    avoid that areas are counted multiple times, the area will be divided by the number of original owners.
    The resulting dataset should not be used as spatial explicit dataset, only for summing up areas etc.
    :param parcels: Shapefile with polygons and area column.
    :param owner_df:
    :return:
    """
    from collections import Counter
    import pandas as pd

    ogc_counts = Counter(owner_df[id_col].tolist())

    # if "area" not in parcels.columns:
    #     parcels["area"] = parcels["geometry"].area

    df_comb = pd.merge(parcels, owner_df, on=id_col, how="left")
    df_comb["d"] = df_comb[id_col].map(ogc_counts)
    df_comb["area"] = df_comb["area"] / df_comb["d"]

    return df_comb

############################################ String manipulation #######################################################

def get_sub(x):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    sub_s = "ₐ₈CDₑբGₕᵢⱼₖₗₘₙₒₚQᵣₛₜᵤᵥwₓᵧZₐ♭꜀ᑯₑբ₉ₕᵢⱼₖₗₘₙₒₚ૧ᵣₛₜᵤᵥwₓᵧ₂₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎"
    res = x.maketrans(''.join(normal), ''.join(sub_s))
    return x.translate(res)