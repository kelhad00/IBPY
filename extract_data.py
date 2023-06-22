import json
import os
import pympi
import numpy
from os import listdir
from os.path import join
from .processing import tuple_to_sequence
import pympi.Elan as elan

def get_all_filepaths(root, ext, to_fill=None):
    """Return list of eaf files.
    
    Args:
        root (str): directory path to read from.
        ext (str): extension of files to read.
        to_fill ([type], optional): list to add the paths to. Defaults to None.
                                    If none, considered an empty list.
    Returns:
        [type]: [description]
    """
    if to_fill is None:
        to_fill = []
    if not isinstance(to_fill, list):
        raise AttributeError("to_fill parameter must be a list")
    paths = os.walk(root)
    for root, _, files in paths:
        for f in files:
            if f.endswith(ext):
                to_fill.append(os.path.join(root, f))
    return to_fill


def read_eaf_to_dict(filepath, mark=True, tiers=None):
    """Read elan EAF file in a dictionary.
    
    Args:
        filepath ([type]): path to single eaf file
        mark (bool, optional):  append same ID to parent/children tiers.
                                Defaults to True.
        tiers (list): list of tier names to keep and discard the rest.
    Returns:
        [dict]: {tier name: [(begin, end, value)]}
    """
    eaf = pympi.Elan.Eaf(filepath)
    dct = {}
    if tiers is None:
        tiers = list(eaf.get_tier_names())
    if mark:  # mark parents and children tiers with same ID
        # get parents and children
        par_child_dct = {}
        for tier in tiers:
            param = eaf.get_parameters_for_tier(tier)
            try:  # parents exist
                if param["PARENT_REF"] not in par_child_dct:
                    par_child_dct[param["PARENT_REF"]] = [tier]
                else:
                    par_child_dct[param["PARENT_REF"]].append(tier)
            except:  # no parents
                continue
        par_child_id = 0
        for parent in par_child_dct:
            for t in range(len(tiers)):
                if (tiers[t] in parent) or (tiers[t] in par_child_dct[parent]):
                    tiers[t] = "_".join([tiers[t], str(par_child_id)])
            par_child_id += 1
    # create final dct
    for tier in tiers:
        dct[tier] = eaf.get_annotation_data_for_tier(tier.split("_")[0])
    return dct


# def get_tier_from_file(filepath, tier, values=None):
#     """Return a dict of {tier:[(strt, stp, val),...]} if val in values
    
#     Args:
#         filepath (str): path to eaf file
#         tier (str): tier name
#         values (list, optional): list of values to keep. Defaults to None.
#     Returns:
#         [dict]: {tier: [(strt, stp, val),...]}
#     """
#     dct = read_eaf_to_dict(filepath)
#     if values is not None:  # if none keep all
#         if not isinstance(values, (list, tuple, numpy.ndarray, str)):
#             raise AttributeError("Values must be a list-like object or a string")
#         if isinstance(values, str):
#             values = [values]
    
#     if tier not in dct:
#         return {tier: []}
    
#     filtered_data = []
#     if values is None:
#         values = set([lab for _, _, lab in dct[tier]])
#     for annot in dct[tier]:
#         if annot[2] in values:
#             filtered_data.append(annot)

#     return {tier: filtered_data}

def get_tier_from_file(filepath, tier, values=None, filename=None):
    """Return a dict of {tier:[(strt, stp, val),...]} if val in values
    Args:
        filepath (str): path to eaf file
        tier (str): tier name
        values (list, optional): list of values to keep. Defaults to None.
        filename (str): name of the json file
    Returns:
        [dict]: {tier: [(strt, stp, val),...]}
    """
    """
    This function is used to get the tier from the eaf file and replace 
    the value of the tier with the value in the json file;
    It has been modified to include the filename parameter which is the name 
    of the json file as part of the CBA-toolkit development: https://github.com/kelhad00/CBA-toolkit.git 
    """
    try :
        with open(f'{filename}.json') as json_file:
                data = json.load(json_file)
        replace_value = data['TIER_LISTS'][tier]['Replace_Value']
        dct = read_eaf_to_dict(filepath)
        if values is not None:  # if none keep all
            if not isinstance(values, (list, tuple, numpy.ndarray, str)):
                raise AttributeError("Values must be a list-like object or a string")
            if isinstance(values, str):
                values = [values]
        if tier not in dct:
            return {tier: []}
        filtered_data = []
        if replace_value != "" :
            for i, item in enumerate(dct[tier]):
                if item[2] != "":
                    dct[tier][i] = (item[0], item[1], replace_value)
                else:
                    dct[tier][i] = (item[0], item[1], "No_" + replace_value)
        if values is None:
            values = set([lab for _, _, lab in dct[tier]])
        for annot in dct[tier]:
            if annot[2] in values:
                filtered_data.append(annot)
        return {tier: filtered_data}
    except :
        dct = read_eaf_to_dict(filepath)
        if values is not None:  # if none keep all
            if not isinstance(values, (list, tuple, numpy.ndarray, str)):
                raise AttributeError("Values must be a list-like object or a string")
            if isinstance(values, str):
                values = [values]
        if tier not in dct:
            return {tier: []}
        filtered_data = []
        if values is None:
            values = set([lab for _, _, lab in dct[tier]])
        for annot in dct[tier]:
            if annot[2] in values:
                filtered_data.append(annot)
        return {tier: filtered_data}


def replace_label(lst, to_replace, value=None, inplace=False, append=None):
    """Replace to_replace label with value in list lst.
    Note: Elements of lst should be in the (start time, stop time, label) format.

    Args:
        lst (list): list of (start time, stop time, label) tuples.
        to_replace (str): label to replace.
        value (str, optional): value to replace with. Defaults to None.
        inplace (bool, optional): if True, replace in place. Defaults to False. 
        append (str, optional): append to the label. Defaults to None.
    Raises:
        AttributeError: both value and append parameters cannot be None.
    Returns:
        [list]: list of (start time, stop time, label) tuples.
    """
    if not (value or append):
        raise AttributeError("both value and append parameters cannot be None.")
    if isinstance(to_replace, str):
        to_replace = [to_replace]
    if not inplace:
        newlst = lst[:]
    for l in range(len(newlst)):
        if newlst[l][2] in to_replace:
            if append:
                value = "_".join([newlst[l][2], append])
            newlst[l] = newlst[l][0], newlst[l][1], value
    return newlst

def get_time_eaf(folder, tiers=None):
    """Return a list of max end time of each eaf file in folder.    

    Args:
        folder (list): list of eaf file paths.
        tiers (list, optional): list of tiers to consider. Defaults to None.
    Returns:
        [list]: list of max end time of each eaf file in folder.
    """
    lst_time = []

    for file in folder:
        dct = read_eaf_to_dict(file, tiers=tiers)
        max_end_time = 0

        for tier in dct:
            annotations = dct[tier]
            for annotation in annotations:
                end_time = annotation[1]
                if end_time > max_end_time:
                    max_end_time = end_time

        lst_time.append(max_end_time / 1000)

    return lst_time


def get_tier_count(folder, tier_name):
    """Return a list of number of annotations in each eaf file in folder.

    Args:
        folder (list): list of eaf file paths.
        tier_name (list): list of tiers to consider.
    Returns:
        [list]: list of number of annotations in each eaf file in folder.
    """
    lst = []

    for file in folder:
        dct = read_eaf_to_dict(file)
        lst_count = []

        for tier in tier_name:
            if tier in dct:
                annotations = dct[tier]
                lst_count.append(len(annotations))
            else:
                lst_count.append(0)

        lst.append(lst_count)

    return lst

def get_max_min_time_tier(folder, tier) :
    """Return a list of max and min duration of each annotation in each eaf file in folder.

    Args:
        folder (list): list of eaf file paths.
        tier (str): tier to consider.
    Returns:
        [list]: list of max and min duration of each annotation in each eaf file in folder.
    """
    lst_max = []
    lst_min = []

    for file in folder :
        dct = read_eaf_to_dict(file)
        if tier not in dct:
            lst_min.append(0.0)
            lst_max.append(0.0)
            continue
        
        annotations = dct[tier]
        min_duration = float("inf")
        max_duration = 0.0
        
        for annotation in annotations:
            start_time = annotation[0]
            end_time = annotation[1]
            duration = end_time - start_time
            
            if duration < min_duration:
                min_duration = duration
            
            if duration > max_duration:
                max_duration = duration
        
        lst_min.append(min_duration/1000)
        lst_max.append(max_duration/1000)

    return lst_min, lst_max


def get_tier_intensities(folder, tier, intensities) :
    """Return a list of dict of {intensity: count} for each file in folder.""

    Args:
        folder (list): list of file paths.
        tier (string): tier name.
        intensities (list): list of intensities to count.
    Returns:
        list: list of dict of {intensity: count} for each file in folder.
    """
    lst = []

    for file in folder:
        dct = read_eaf_to_dict(file)
        if tier not in dct:
            lst.append(dict.fromkeys(intensities, 0))
            continue
        
        annotations = dct[tier]
        intensity_count = dict.fromkeys(intensities, 0)

        for annotation in annotations:
            intensity = annotation[2].strip()
            
            if intensity in intensity_count:
                intensity_count[intensity] += 1
    
        lst.append(intensity_count)
        
    return lst

def remove_label(lst, to_remove):
    """Return list not containing the labels in to_remove.
    
    Args:
        lst (list of tuples): list of the type [ ( _, _, label) ] from which data should be removed.
        to_remove (string or list of strings): label(s) that should be removed from lst.
    
    Returns:
        list: list of the same type as lst.
    """
    if isinstance(to_remove, str):
        to_remove = [to_remove]
    newlst = []
    for tup in lst:
        if tup[2] not in to_remove:
            newlst.append(tup)
    return newlst


def keep_only(lst, tier_to_keep, inplace=False):
    """Keep filepaths for files with annotation data in the tier tier_to_keep.

    Note: keep it whether file is entirely annotated or not.

    Args:
        lst (list): list of filepaths.
        tier_to_keep (list): list of tiers to keep.
        inplace (bool): if True, lst is modified in place.

    Returns:
        list: list of filepaths.
    """
    if inplace:
        filter_lst = lst
    else:
        filter_lst = lst[:]
    i = 0
    limit = len(filter_lst)
    while i < limit:
        eaf = pympi.Elan.Eaf(filter_lst[i])
        dct = read_eaf_to_dict(filter_lst[i])
        if tier_to_keep not in dct:
            filter_lst.pop(i)
            limit -= 1
            continue
        i += 1
    return filter_lst


#ADDED
def tuple_access(filepath, expression, where, index):
    """This function accesses the information of a tuple.
    
    Args: 
        filepath (str): path of the file.
        expression (str): name of the tier.
        where (int): position of the tuple to which we want to access.
        index (int): index of the element to be accessed in the tuple (0, 1 ou 2).
    Returns: The wanted information.
    """
    to_dict = read_eaf_to_dict(filepath, mark=True, tiers=None)
    m = to_dict[expression][where-1][index]

    return m

def find_indice_tuple(lst, value, idx):
    """Finds the index of a tuple in a list of multiple tuples.
    
    Args:
        lst (list): the list containing the tuples.
        value (int/str): value contained in the tuple.
        idx (int): position in the tuple corresponding to "value".
    Returns: 
        int: the index of the tuple where the value at position idx is.
    """
    indice = 0
    for _ in lst:
        indice += 1
        if value == _[idx]:
            return (indice)

def get_all_filenames(root, ext, to_fill=None):
    """Return list of eaf files name.

    Args:
        root (str): directory path to read from.
        ext (str): extension of files to read.
        to_fill ([type], optional): list to add the names to. Defaults to None.
                                    If none, considered an empty list.
    Returns:
        [filename1, filname2,....].
    """
    if to_fill is None:
        to_fill = []
    if not isinstance(to_fill, list):
        raise AttributeError("to_fill parameter must be a list")
    paths = os.walk(root)
    for root, _, files in paths:
        #print(root)
        for f in files:
            if f.endswith(ext):
                to_fill.append(os.path.join(f))
    return to_fill
 
def replace_intensity(lst):
    """This function replace intensity by numbers.

    Args: 
        lst (list): list of tuples (stt, stp, label).
    Returns: 
        list: list of tuples (stt, stp, label).
    """
    labels = set(i[2] for i in lst)
    # num_labels = len(labels)

    label_mapping = {}
    for index, label in enumerate(labels):
        label_mapping[label] = index + 1

    return [(stt, stp, label_mapping[label]) for stt, stp, label in lst]

def tuple_to_int_sequence(lst, width, shift):
    """This function convert tuple (stt, stp, label) into a sequence of int corresponding to the label.

    Args:
        lst (list): list of tuple (stt, stp, label).
        width (numeric): window width in ms.
        shift (numeric): window shift in ms.
    Returns:
        list: A sequence of int.
    """
    if len(lst) == 0:
        lst_int = []
    else:
        L = replace_intensity(lst)
        lst_int = tuple_to_sequence(L, width, shift)
        for i in range(len(lst_int)):
            if lst_int[i] is None:
                lst_int[i] = 0
    return lst_int

