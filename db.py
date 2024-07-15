import os
import re
from .utils import list_to_df
from .extract_data import get_all_filenames, get_all_filepaths

def form_pairs_ndcme(lst):
    """Return filename pairs [(), (), ...].

    Args:
        lst (list): list of filenames without the path.
    Returns:
        list: [(), (), ...].
    """
    lst_sorted = sorted(lst)
    i = 0
    final = []
    while i <= len(lst_sorted) - 1:
        curr = '_'.join(lst_sorted[i].split('_')[:2])
        nxt = '_'.join(lst_sorted[i + 1].split('_')[:2])
        if curr == nxt:
            final.append((lst_sorted[i], lst_sorted[i + 1]))
            i += 2
        else:
            i += 1
    return final

def form_pairs_ndc_me(lst):
    """Return filename pairs [(), (), ...].

    Args:
        lst (list): list of filenames without the path.
    Returns:
        list: [(), (), ...].
    """
    lst_sorted = sorted(lst)
    i = 0
    final = []
    while i <= len(lst_sorted) - 1:
        curr = '_'.join(lst_sorted[i].split('_')[:2])
        nxt = '_'.join(lst_sorted[i + 1].split('_')[:2])
        if curr == nxt:
            final.append((lst_sorted[i], lst_sorted[i + 1]))
            i += 2
        else:
            i += 1
    return final

def form_pairs_ccdb(lst):
    """Return filename pairs [(), (), ...].

    Args:
        lst (list): list of filenames without the path.
    Returns:
        list: [(),(),...].
    """
    final = []
    replace_with = {'dizzy': 'monk', 'monk': 'dizzy'}
    i = 0
    while i < len(lst):
        key = lst[i].split('_')[-1].split('.')[0]
        pair = lst[i].replace(key, replace_with[key])
        if pair in lst:
            final.append((lst[i], pair))
            lst.remove(lst[i])
            lst.remove(pair)
        else:
            i = +1
    for l in lst:
        continue
    return final

def form_pairs_ifadv(lst):
    """Return list of filename pairs.

    Args:
        lst (list): list of filenames.
    Returns:
        list: [(), (), ...].
    """
    lst1 = []
    lst2 = []
    for l in lst:
        ind = int(re.findall('\\d+', l)[0])
        if l[2] == 'A':
            lst1.append(l)
        else:
            lst2.append(l)
    new_lst2 = []
    for i in range(len(lst1)):
        pair_name = lst1[i].replace('A', 'B', 1)[:4] + '{}'
        pair_exists = any((name.startswith(pair_name[:4]) for name in lst2))
        if not pair_exists:
            lst1[i] = None
        else:
            index = [j for j, name in enumerate(lst2) if name.startswith(pair_name[:4])][0]
            new_lst2.append(lst2[index])
    lst1 = [name for name in lst1 if name is not None]
    lst1.sort()
    new_lst2.sort()
    return list(zip(lst1, new_lst2))

def form_pairs_ab(lst):
    """Return filename pairs [(), (), ...].

    Args:
        lst (list): list of filenames without the path.
    Returns:
        list: [(), (), ...].
    """
    pairs = []
    for i in range(len(lst)):
        filename = lst[i]
        name_parts = filename.split('_')
        prefix = name_parts[0]
        number = name_parts[1]
        for j in range(i + 1, len(lst)):
            filename2 = lst[j]
            name_parts2 = filename2.split('_')
            prefix2 = name_parts2[0]
            number2 = name_parts2[1]
            if prefix == 'A' and prefix2 == 'B' and (number == number2):
                pair2 = 'A_' + number2
                if pair2 == '_'.join(name_parts[:2]) and (filename, filename2) not in pairs:
                    pairs.append((filename, filename2))
        pairs = sorted(pairs)
    return pairs

def form_pairs(ROOT1, ROOT2, ROOT3):
    """Gives filespath of pairs for each dataset.

    Args:
        ROOT1 (str): path of ccdb directory.
        ROOT2 (str): path of ifadv directory.
        ROOT3 (str): path of ndc directory.         
    Returns:
        list: [pair1_A, pair1_B, pair2_A, pair2_B,....].
    """
    c = form_pairs_ccdb(get_all_filenames(ROOT1, 'eaf'))
    i = form_pairs_ifadv(get_all_filenames(ROOT2, 'eaf'))
    n = form_pairs_ndcme(get_all_filenames(ROOT3, 'eaf'))
    liste_ccdb = list(sum(c, ()))
    liste_ifadv = list(sum(i, ()))
    liste_ndc = list(sum(n, ()))
    pair_ccdb, pair_ifadv, pair_ndc = ([] for _ in range(3))
    for _ in liste_ccdb:
        pair_ccdb.append(os.path.join(ROOT1, _))
    for _ in liste_ifadv:
        pair_ifadv.append(os.path.join(ROOT2, _))
    for _ in liste_ndc:
        pair_ndc.append(os.path.join(ROOT3, _))
    return (pair_ccdb, pair_ifadv, pair_ndc)

def form_list_pairs(ROOT, foldername):
    """Gives filespath of pairs for your dataset.
    
    Args:
        ROOT(str): path of your directory.
        foldername(str): name of the folder.

    Returns:
        list: [pair1_A, pair1_B, pair2_A, pair2_B,....]
    """
    print('foldername', foldername)
    function_name = f'form_pairs_{foldername}'
    n = eval(function_name)(get_all_filenames(ROOT, 'eaf'))
    liste = list(sum(n, ()))
    pair = []
    for item in liste:
        pair.append(os.path.join(ROOT, item))
    return pair

def get_db_from_func_pair(dir, func, database, expression_choice, tier_lists):
    """This function takes a path as an argument and creates a dataset 
        based on the number of items in the folder using a chosen function.
    It takes into account pairs in our datasets.

    Args:
        dir (str) : Path of the folder containing all datasets.
        func (function): Function we want to use.
        database (str): Dataset we want to use.
        expression_choice (str): Expression we want to use.
        tier_lists (list): List of tiers we want to use.
    Returns:
        dataframe: A dataframe corresponding to the function chosen.
    """
    n = 0
    L = []
    dct = {db.lower(): globals()['form_list_pairs'] for db in database}
    for path in os.listdir(dir):
        if os.path.isdir(os.path.join(dir, path)):
            n += 1
            for i, j in zip(list(dct.keys()), list(dct.values())):
                if path == i:
                    L.append((j(os.path.join(dir, path), i), path))
    dg = []
    for i in range(len(L)):
        dg += func(L[i][0], L[i][1], expression_choice, tier_lists)[0]
    dg = list_to_df(dg, func(L[0][0], L[0][1], expression_choice, tier_lists)[1])
    return dg

def get_db_from_func_pair_tier(dir, func, database, tier1, tier2, entity1, entity2):
    """This function takes a path as an argument and creates a dataset 
        based on the number of items in the folder using a chosen function.
    It takes into account pairs in our datasets.

    Args:
        dir (str) : Path of the folder containing all datasets.
        func (function): Function we want to use.
        database (str): Dataset we want to use.
        tier1 (str): First tier we want to use.
        tier2 (str): Second tier we want to use.
        entity1 (str): Entity of tier1 we want to use.
        entity2 (str): Entity of tier2 we want to use.
    Returns:
        dataframe: A dataframe corresponding to the function chosen.
    """
    n = 0
    L = []
    dct = {database.lower(): globals()['form_list_pairs'] for _ in database}
    for path in os.listdir(dir):
        if os.path.isdir(os.path.join(dir, path)):
            n += 1
            for i, j in zip(list(dct.keys()), list(dct.values())):
                if path == i:
                    L.append((j(os.path.join(dir, path), i), path))
    dg = []
    for i in range(len(L)):
        dg += func(L[i][0], L[i][1], tier1, tier2, entity1, entity2)[0]
    dg = list_to_df(dg, func(L[0][0], L[0][1], tier1, tier2, entity1, entity2)[1])
    return dg

def get_db_from_func_no_pair(dir, func, database_names, tier):
    """This function takes a path as an argument and creates a dataset 
        based on the number of items in the folder using a chosen function.
    It doesn't take into account pairs in our datasets.

    Args:
        dir (str) : Path of the folder containing all datasets.
        func (function): Function we want to use.
        database_names (list): List of datasets names.
        tier (str): Tier name.
    Returns:
        dataframe: A dataframe corresponding to the function chosen.
    """
    n = 0
    L = []
    for path in os.listdir(dir):
        if os.path.isdir(os.path.join(dir, path)):
            n += 1
            for i in database_names:
                if path == i.lower():
                    L.append((get_all_filepaths(os.path.join(dir, path), 'eaf', None), path))
    dg = []
    for i in range(len(L)):
        dg += func(L[i][0], L[i][1], tier)[0]
    dg = list_to_df(dg, func(L[0][0], L[0][1], tier)[1])
    return dg

def get_db_from_func_no_pair_tier(dir, func, database_names, tier1, tier2, entity):
    """This function takes a path as an argument and creates a dataset 
        based on the number of items in the folder using a chosen function.
    It doesn't take into account pairs in our datasets.

    Args:
        dir (str) : Path of the folder containing all databases.
        func (function): Function we want to use.
        database_names (list): List of database names.
        tier1 (str): First tier name.
        tier2 (str): Second tier name.
        entity (str): Entity name of tier1.
    Returns:
        dataframe: A dataframe corresponding to the function chosen.
    """
    n = 0
    L = []
    for path in os.listdir(dir):
        if os.path.isdir(os.path.join(dir, path)):
            n += 1
            for i in database_names:
                if path == i.lower():
                    L.append((get_all_filepaths(os.path.join(dir, path), 'eaf', None), path))
    dg = []
    for i in range(len(L)):
        dg += func(L[i][0], L[i][1], tier1, tier2, entity)[0]
    dg = list_to_df(dg, func(L[0][0], L[0][1], tier1, tier2, entity)[1])
    return dg

def form_pairs_IB(lst):
    """Return filename pairs [(), (), ...].

    Args:
        lst (list): list of filenames without the path.
    Returns:
        list: [(), (), ...].
    """
    pairs = form_pairs_ab(lst)
    return pairs
