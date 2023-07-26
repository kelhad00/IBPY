import os 
import pandas as pd
from .utils import overlapping_dct_from_indices_to_vals
from .extract_data import get_tier_from_file

def get_overlapping_segments_ind(lstA, lstB):  
    """Get segments in A and B that overlap.
    
    Args:
        lstA (list of tuples): [(start time, stop time, lab),..].
        lstB (list of tuples): [(start time, stop time, lab),..].
    Returns
        dict: {index of segment in lstA: [indices of segments in lstB]}
    """
    indA = 0
    indB = 0
    dct = {}
    while indA < len(lstA) and indB < len(lstB):
        while lstA[indA][0] >= lstB[indB][1]:       
            indB += 1
            if indB >= len(lstB):
                return dct
        while lstA[indA][1] <= lstB[indB][0]:      
            indA += 1
            if indA >= len(lstA):
                return dct
        if (lstA[indA][1] > lstB[indB][1] > lstA[indA][0]) or (
            lstA[indA][1] > lstB[indB][0] > lstA[indA][0]
        ):
            while indB < len(lstB) and lstB[indB][1] < lstA[indA][1]:
                if indA in dct:
                    if indB in dct[indA]:           #added
                        pass                        #added
                    else:                           #added
                      dct[indA].append(indB)
                else:
                    dct[indA] = [indB]
                dct[indA].append(indB+1)            #added
                indB += 1
            indA += 1

        elif (lstB[indB][0] <= lstA[indA][0]) and (lstB[indB][1] >= lstA[indA][1]): 
            while indA < len(lstA) and (lstA[indA][1] < lstB[indB][1]):
                if indA in dct:
                    dct[indA].append(indB)
                else:
                    dct[indA] = [indB]
                indA += 1
            indB += 1

    return dct


def gos_ind(A, B): #get_overlapping_segments_ind - correct function
    """Get segments in A and B that overlap.
    
    Args:
        lstA (list of tuples): [(start time, stop time, lab),..].
        lstB (list of tuples): [(start time, stop time, lab),..].
    Returns:
        dict: {index of segment in lstA: [indices of segments in lstB]}.
    """
    indA = 0
    indB = 0
    dct = {}

    while indA < len(A) and indB < len(B):
        # Left bound for intersecting segment
        l = max(A[indA][0], B[indB][0])
         
        # Right bound for intersecting segment
        r = min(A[indA][1], B[indB][1])
         
        # If segment is valid, we add indices corresponding to A and B to the dictionnary
        if l <= r:
            if indA in dct:            #if A indice is already in the dict,
                dct[indA].append(indB)      #we add B indice to those who are aloready ther
            else :
                dct[indA] = [indB]            
 
        #If the endtime of the i-th interval of list A is smaller, we increment indice A.
        #Else, we increment indice B
        if A[indA][1] < B[indB][1]:
            indA += 1
        else:
            indB += 1

    return dct


def get_overlapping_segments(lstA, lstB, values_only=False):
    """Get segments in lstB overlapping with segments of lstA.
    
    Args:
        lstA ([type]): [(start, stop, label), etc.].
        lstB ([type]): [(start, stop, label), etc.].
        values_only (bool, optional): [description]. Defaults to False.
    Returns:
        dict: {Segments in A: [Segments in B]}.
    """
    dct_inds = gos_ind(lstA, lstB)
    if values_only:
        lstA_tempo = [val for b, e, val in lstA]
        lstB_tempo = [val for b, e, val in lstB]
    else:
        lstA_tempo = lstA[:]
        lstB_tempo = lstB[:]
    dct = overlapping_dct_from_indices_to_vals(dct_inds, lstA_tempo, lstB_tempo)
    return dct

# high levels
def count_mimicry(lstA, lstB, delta_t=0):
    """Count the occurences of B mimicking A by delta_t.
   
    This implementation counts mimicry based on method in [1]
    and also returns the instances of mimickry.

    The times in each of lstA and lstB cannot overlap internally.
    They have to be successive segments in each of lstA and lstB.

    [1] Feese, Sebastian, et al. "Quantifying behavioral mimicry by automatic 
    detection of nonverbal cues from body motion." 2012 International Conference 
    on Privacy, Security, Risk and Trust and 2012 International Conference on 
    Social Computing. IEEE, 2012.
    
    Args:
        lstA (list): list of tuples (start, stop, label) of expressions mimicked.
        lstB (list): list of tuples (start, stop, label) of expressions mimicking.
        delta_t (int, optional): Defaults to 0.
                                Time after which expression occuring still counts as mimicry.
                                Should be in the same unit as the times in lstA and lstB.    
    Returns:
        int: number of times B mimicked A (=len(the list described below)).
        list: [(indA, indB),...]
              where the indB element of B mimick the indA element of A
              following the definition of mimickry described in the reference above.
    """
    indA = 0
    indB = 0
    count = 0  # number of mimicry events
    mimic_ind = []  # indices of mimicry events in lstB
    if len(lstA) == 0 or len(lstB) == 0:  # if at least one of them has no expression
        return 0, []
    while indA < len(lstA) and indB < len(lstB):
        if lstB[indB][0] <= lstA[indA][0]:
            indB += 1
        elif (
            lstB[indB][0] > lstA[indA][0] and (lstB[indB][0] - delta_t) <= lstA[indA][1]
        ):
            # avoid double counting incase delta_t is > (lstA[indA+1][0] - lstA[indA][1])
            if (indA + 1) < len(lstA):
                if lstB[indB][0] > lstA[indA + 1][0]:
                    indA += 1  # skip to next lstA expression
                    continue
            count += 1
            mimic_ind.append((indA, indB))
            # if no double counting
            # check if several expressions from B overlap with A's
            while lstB[indB][1] <= lstA[indA][1]:
                indB += 1  # skip to the following expression untill no more overlapping
                if indB == len(lstB):
                    break
            indA += 1
        elif (lstB[indB][0] - delta_t) > lstA[indA][1]:
            indA += 1
    return count, mimic_ind


def count_mimicry_per_value_in_tier(ref, target, delta_t):
    """Return the number of times mimicry occured.
    
    Considers that all expresssions in ref and target are the same.
    So all are potential mimicry events.

    Args:
        ref (dict): dictionary of values in tier mimicked.
        target (dict): dictionary of values in tier potentially mimicking.
        delta_t (float): time after which expression occuring still counts as mimicry.
                        Should be in the same unit as the times in ref and target.
    Raises:
        AttributeError: if there are duplicate values in ref.
    Returns:
        dict: {ref_value: {target_value: (count, [(indA, indB),...])}}
    """
    final = {}
    if len(set(ref)) != len(ref):
        raise AttributeError("No parameter is allowed in the parameter ref")
    for r in ref:
        final[r] = {}
        for tar in target:
            final[r][tar] = count_mimicry(ref[r], target[tar], delta_t=delta_t)
    return final


def calculate_mimicking_ratio(total_mimicker_expressions, total_mimicked_expressions):
    """Return the ratio of the total number of expression that are mimicking to 
    the total number of a certain expression.
    
    Args:
        total_mimicker_expression (float): total number of expression that are mimicking.
        total_mimicked_expressions (float): total number of expression that are mimicked.
    Returns:
        float: ratio of the total number of expression that are mimicking to
    """

    return total_mimicked_expressions / total_mimicker_expressions


def following_expressions(lst, delta_t=0):
    """succession of expressions in tier
    
    Args:
        lst (list): list of tuples (start, stop, label) of expressions.
        delta_t (int, optional): Defaults to 0.
                                Time after which expression occuring still counts as mimicry.
                                Should be in the same unit as the times in lstA and lstB.
    Returns:
        dict: {label: [label, label, ...]}.
    """
    dct = {}
    for l in range(len(lst) - 1):
        if (lst[l + 1][0] - lst[l][1]) <= delta_t:
            if lst[l][2] in dct:
                dct[lst[l][2]].append(lst[l + 1])
            else:
                dct[lst[l][2]] = [lst[l + 1]]
        else:
            if lst[l][2] in dct:
                dct[lst[l][2]].append(None)
            else:
                dct[lst[l][2]] = [None]
    return dct


def count_vals_in_tier(lst, vals_to_count=None):
    """Count the number of times each value in vals_to_count appears in lst.

    Args:
        lst (list): list of tuples (start, stop, label) of expressions.
        vals_to_count (list, optional): Defaults to None.
                                        List of values to count.
                                        If None, all values are counted.
    Returns:
        dict: {label: count}.
    """
    dct = {}
    for lab in lst:
        if lab is not None:
            lab = lab[2]

        if lab in dct:
            dct[lab] += 1
        else:
            dct[lab] = 1
    return dct


def calculate_correlation(lstA, lstB):
    """Calculate the correlation between two tiers.

    Args:
        lstA (list): list of tuples (start, stop, label) of expressions.
        lstB (list): list of tuples (start, stop, label) of expressions.
    Returns:
        float: correlation between the two tiers.
    """

    pass #TODO


def count_following(lst, n, max_dist):
    """Count the number of times each label is followed by each other label in a list.

    Args:
        lst (list): list of tuples (start, stop, label) of expressions.
        n (int): number of following expressions to consider.
        max_dist (int): maximum distance between two expressions to consider them as following.
    Returns:
        dict: {label: {n: {label: count}}}.
    """
    labs = set([l for _, _, l in lst])
    dct = {}
    for l in labs:
        dct[l] = {}
    for i in range(len(lst) - n):
        if lst[i][2] not in dct:
            dct[lst[i][2]] = {}
        if (lst[i + n][0] - lst[i][1]) <= max_dist:
            for j in range(1, n + 1):
                dct[lst[i][2]][j] = {}
                if lst[i + j][2] not in dct[lst[i][2]][j]:
                    dct[lst[i][2]][j][lst[i + j][2]] = 0
                dct[lst[i][2]][j][lst[i + j][2]] += 1
    return dct


def get_next_n_exp(lst, n, max_dist, append_none=True):
    """return lists of n labels following each different label in a list.
    
    Args:
        lst (list of tuples): list of type [(start, stop, label)]
        n (int): number of elements.
        max_dist (int): maximum distance between elements, in number of elements.
                        After this distance, labels are not considered following the current one,
        append_none (bool, optional): fill with None if no more following label. Defaults to True.
    
    Returns:
        dict: {label: [followinglabels]}
    """
    dct = {}
    for l in range(len(lst) - n):  # skip the last n elements (cannot assume they are None)
        lab = lst[l][2]
        if lab not in dct:
            dct[lab] = []
        temp = []
        for ind_next in range(1, n + 1):
            next_close = lst[l + ind_next - 1][1]
            next_far = lst[l + ind_next][0]
            if (next_far - next_close) <= max_dist:
                temp.append(lst[l + ind_next][2])
            else:
                if append_none:
                    temp.extend([None] * (n - ind_next + 1))
                break
        if len(temp) == n:
            dct[lab].append(temp)
    return dct


def get_prev_n_exp(lst, n, max_dist, append_none=True):
    """return lists of n labels preceding each different label in a list.
    
    Args:
        lst (list of tuples): list of type [(start, stop, label)]
        n (int): number of elements.
        max_dist (int): maximum distance between elements, in number of elements.
                        After this distance, labels are not considered preceding the current one,
        append_none (bool, optional): fill with None if no more preceding label. Defaults to True.
    Returns:
        dict: {label: [preceding labels]}
    """
    dct = {}
    for l in range(
        n, len(lst)
    ):  # skip the first n elements (cannot assume they are None)
        lab = lst[l][2]
        if lab not in dct:
            dct[lab] = []
        temp = []
        for ind_next in range(1, n + 1):
            prev_close = lst[l - ind_next + 1][0]
            prev_far = lst[l - ind_next][1]
            if (prev_close - prev_far) <= max_dist:
                temp.append(lst[l - ind_next][2])
            else:
                if append_none:
                    temp.extend([None] * (n - ind_next + 1))
                break
        if len(temp) == n:
            dct[lab].append(temp)
    return dct


## ADDED ##

def get_overlapping_seg(A, B): 
    """Get segments in A and B that overlap.
    Same as get_overlapping_segments but here, the function makes directly intersection between the segments.
    
    Args:
        lstA (list of tuples): [(start time, stop time, lab),..].
        lstB (list of tuples): [(start time, stop time, lab),..]
    Returns
        list: [(startime overlap, endtime overlap, lab)]
    """
    indA = 0
    indB = 0
    lst = []

    while indA < len(A) and indB < len(B):
        # Left bound for intersecting segment
        l = max(A[indA][0], B[indB][0])
         
        # Right bound for intersecting segment
        r = min(A[indA][1], B[indB][1])
         
        # If segment is valid print it
        if l <= r:
            lst.append((l,r,B[indB][2]))
 
        # If i-th interval's right bound is
        # smaller increment i else increment j
        if A[indA][1] < B[indB][1]:
            indA += 1
        else:
            indB += 1

    return lst

def overlap_count_SL(databases_name, databases_pairs, databases_pair_paths, expression_pairs, expressions_track, choice):
    """
    This function returns the number and percentage of overlapping segments of smiles and laughs between two persons for each pair of files.
    When we talk about overlapping between person A to person B, this means that we are looking at all the 
    segments of person B which overlap a segment of A and this for all the segments of the tier concerned for 
    the "A" files. When we change direction (so person B to person A), we just switched the direction of the 
    files in the overlapping function. (A and B are the pair files)
        - person A to person B: {(segmentA: (segmentB n°1), (segmentB n°N), etc),...}
        - person B to person A: {(segmentB: (segmentA n°1), (segmentA n°N), etc),...}

    Args:
        databases_name (list): list of databases names
        databases_pairs (list): list of databases pairs names
        databases_pair_paths (dict): dictionary of databases pairs paths
        expression_pairs (list): list of expressions pairs names. For example : [("Smiles_0", "Smiles_0"), 
                                                                                 ("Smiles_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Smiles_0"),
                                                                                 ("Role", "Role")]
        expressions_track (list): list of expressions names. For example : [("Smiles_0", "Smiles_0"), 
                                                                            ("Smiles_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Smiles_0")]
        choice (str): choice between "A/B" for  or "B/A". For example : "A/B" for overlapping between person A to person B 
                      or "B/A" for the opposite.

    Returns:
        dataframes (dict): dictionary of dataframes
    """
    if choice == "A/B":
        person1 = "A"
        person2 = "B"
    elif choice == "B/A":
        person1 = "B"
        person2 = "A"

    lstA = {}
    lstB = {}
    overlapping_segments_dict = {}
    for i, database in enumerate(databases_name):
        if database == databases_pairs[i].replace('_pairs', '').upper():
            databases_list = databases_pair_paths[databases_pairs[i]]
            dataset_dict = {}
            for i in range(0, len(databases_list), 2):
                filepath_A = databases_list[i]
                filepath_B = databases_list[i+1] 
                pair_file_A = os.path.basename(filepath_A)
                pair_file_B = os.path.basename(filepath_B)

                if pair_file_A and pair_file_B:
                    if person1 == "A":
                        pair_name = f"{pair_file_A}_&_{pair_file_B}"
                    elif person1 == "B":
                        pair_name = f"{pair_file_B}_&_{pair_file_A}"

                pair_dict = {}
                overlapping_data = {}

                for tier_A, tier_B in expression_pairs:
                    lstA_tier = get_tier_from_file(filepath_A, tier_A)
                    lstB_tier = get_tier_from_file(filepath_B, tier_B)

                    if tier_A in lstA:
                        lstA[tier_A].extend(lstA_tier[tier_A])
                    else:
                        lstA[tier_A] = lstA_tier[tier_A]

                    if tier_B in lstB:
                        lstB[tier_B].extend(lstB_tier[tier_B])
                    else:
                        lstB[tier_B] = lstB_tier[tier_B]

                    if person1 == "A":
                        overlapping_segments = get_overlapping_segments(lstA_tier[tier_A], lstB_tier[tier_B])
                        overlapping_data[f"{tier_A} vs {tier_B}"] = {'Segments': overlapping_segments}

                        overlapping_data[f"{tier_B} count in lstB"] = 0

                        tiers_in_lstB = set(lstB_tier[tier_B])
                        overlapping_data[f"{tier_B} count in lstB"] = len(tiers_in_lstB)
                    elif person1 == "B":
                        overlapping_segments = get_overlapping_segments(lstB_tier[tier_B], lstA_tier[tier_A])
                        overlapping_data[f"{tier_B} vs {tier_A}"] = {'Segments': overlapping_segments}

                        overlapping_data[f"{tier_A} count in lstA"] = 0

                        tiers_in_lstA = set(lstA_tier[tier_A])
                        overlapping_data[f"{tier_A} count in lstA"] = len(tiers_in_lstA)

                dataset_dict[pair_name] = overlapping_data

            overlapping_segments_dict[database] = dataset_dict
    dataframes = {}
    overlap_segments_set = set()
    for database, dataset_dict in overlapping_segments_dict.items():
        overlap_count_list = []
        for pair_name, pair_dict in dataset_dict.items():
            overlap_count_spk_vs_lsn = 0
            overlap_count_lsn_vs_spk = 0
            percentage_spk_vs_lsn = 0
            percentage_lsn_vs_spk = 0
            if person1 == "A":
                count_smiles = pair_dict["Smiles_0 count in lstB"]
                count_laughs = pair_dict["Laughs_0 count in lstB"]
            elif person1 == "B":
                count_smiles = pair_dict["Smiles_0 count in lstA"]
                count_laughs = pair_dict["Laughs_0 count in lstA"]
            segments = pair_dict["Role vs Role"]["Segments"]
            expression = expressions_track
            if person1 == "A":
                for segmentA, segmentB in segments.items():
                    for segB in segmentB:
                        segment_key = f"{segB}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if A is "spk" and B is "lsn"
                            if (segmentA[2].replace(" ", "") == "spk" and segB[2].replace(" ", "") == "lsn"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        overlap_count_spk_vs_lsn += 1
                            # Check if A is "lsn" and B is "spk"
                            elif (segmentA[2].replace(" ", "") == "lsn" and segB[2].replace(" ", "") == "spk"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        overlap_count_lsn_vs_spk += 1
                        break
            elif person1 == "B":
                for segmentB, segmentA in segments.items():
                    for segA in segmentA:
                        segment_key = f"{segA}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if B is "spk" and A is "lsn"
                            if (segmentB[2].replace(" ", "") == "spk" and segA[2].replace(" ", "") == "lsn"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        overlap_count_spk_vs_lsn += 1
                            # Check if B is "lsn" and A is "spk"
                            elif (segmentB[2].replace(" ", "") == "lsn" and segA[2].replace(" ", "") == "spk"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        overlap_count_lsn_vs_spk += 1
                        break

            if count_smiles != 0 or count_laughs != 0:
                percentage_spk_vs_lsn = overlap_count_spk_vs_lsn / (count_smiles + count_laughs) * 100        
                percentage_lsn_vs_spk= overlap_count_lsn_vs_spk / (count_smiles + count_laughs) * 100

            overlap_count_list.append({
                'Database': database,
                'Pair': pair_name,
                f'Overlap Count for {person1} spk / {person2} lsn - S&L': overlap_count_spk_vs_lsn,
                f'Overlap Count for {person1} lsn / {person2} spk - S&L': overlap_count_lsn_vs_spk,
                f'Overlap Percentage for {person1} spk / {person2} lsn - S&L': percentage_spk_vs_lsn,
                f'Overlap Percentage for {person1} lsn / {person2} spk - S&L': percentage_lsn_vs_spk,
            }) 
        df_overlap_count = pd.DataFrame(overlap_count_list)
        dataframes[database] = df_overlap_count
    return dataframes

def overlap_count_SL_advanced(databases_name, databases_pairs, databases_pair_paths, expression_pairs, expressions_track, choice):
    """
    This function returns the number and percentage of overlapping segments of smiles vs smiles, laughs vs smiles, smiles vs laughs and laughs vs laughs between two persons for each pair of files.
    When we talk about overlapping between person A to person B, this means that we are looking at all the 
    segments of person B which overlap a segment of A and this for all the segments of the tier concerned for 
    the "A" files. When we change direction (so person B to person A), we just switched the direction of the 
    files in the overlapping function. (A and B are the pair files)
        - person A to person B: {(segmentA: (segmentB n°1), (segmentB n°N), etc),...}
        - person B to person A: {(segmentB: (segmentA n°1), (segmentA n°N), etc),...}

    Args:
        databases_name (list): list of databases names
        databases_pairs (list): list of databases pairs names
        databases_pair_paths (dict): dictionary of databases pairs paths
        expression_pairs (list): list of expressions pairs names. For example : [("Smiles_0", "Smiles_0"), 
                                                                                 ("Smiles_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Smiles_0"),
                                                                                 ("Role", "Role")]
        expressions_track (list): list of expressions names. For example : [("Smiles_0", "Smiles_0"), 
                                                                            ("Smiles_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Smiles_0")]
        choice (str): choice between "A/B" for  or "B/A". For example : "A/B" for overlapping between person A to person B 
                      or "B/A" for the opposite.

    Returns:
        dataframes (dict): dictionary of dataframes
    """
    if choice == "A/B":
        person1 = "A"
        person2 = "B"
    elif choice == "B/A":
        person1 = "B"
        person2 = "A"

    lstA = {}
    lstB = {}
    overlapping_segments_dict = {}
    for i, database in enumerate(databases_name):
        if database == databases_pairs[i].replace('_pairs', '').upper():
            databases_list = databases_pair_paths[databases_pairs[i]]
            dataset_dict = {}
            for i in range(0, len(databases_list), 2):
                filepath_A = databases_list[i]
                filepath_B = databases_list[i+1] 
                pair_file_A = os.path.basename(filepath_A)
                pair_file_B = os.path.basename(filepath_B)

                if pair_file_A and pair_file_B:
                    if person1 == "A":
                        pair_name = f"{pair_file_A}_&_{pair_file_B}"
                    elif person1 == "B":
                        pair_name = f"{pair_file_B}_&_{pair_file_A}"

                pair_dict = {}
                overlapping_data = {}

                for tier_A, tier_B in expression_pairs:
                    lstA_tier = get_tier_from_file(filepath_A, tier_A)
                    lstB_tier = get_tier_from_file(filepath_B, tier_B)

                    if tier_A in lstA:
                        lstA[tier_A].extend(lstA_tier[tier_A])
                    else:
                        lstA[tier_A] = lstA_tier[tier_A]

                    if tier_B in lstB:
                        lstB[tier_B].extend(lstB_tier[tier_B])
                    else:
                        lstB[tier_B] = lstB_tier[tier_B]

                    if person1 == "A":
                        overlapping_segments = get_overlapping_segments(lstA_tier[tier_A], lstB_tier[tier_B])
                        overlapping_data[f"{tier_A} vs {tier_B}"] = {'Segments': overlapping_segments}

                        overlapping_data[f"{tier_B} count in lstB"] = 0

                        tiers_in_lstB = set(lstB_tier[tier_B])
                        overlapping_data[f"{tier_B} count in lstB"] = len(tiers_in_lstB)
                    elif person1 == "B":
                        overlapping_segments = get_overlapping_segments(lstB_tier[tier_B], lstA_tier[tier_A])
                        overlapping_data[f"{tier_B} vs {tier_A}"] = {'Segments': overlapping_segments}

                        overlapping_data[f"{tier_A} count in lstA"] = 0

                        tiers_in_lstA = set(lstA_tier[tier_A])
                        overlapping_data[f"{tier_A} count in lstA"] = len(tiers_in_lstA)

                dataset_dict[pair_name] = overlapping_data

            overlapping_segments_dict[database] = dataset_dict
    dataframes = {}
    overlap_segments_set = set()
    for database, dataset_dict in overlapping_segments_dict.items():
        overlap_count_list = []
        for pair_name, pair_dict in dataset_dict.items():
            overlap_count_spk_vs_lsn_Smiles_0_vs_Smiles_0 = 0
            overlap_count_spk_vs_lsn_Laughs_0_vs_Laughs_0 = 0
            overlap_count_lsn_vs_spk_Smiles_0_vs_Smiles_0 = 0
            overlap_count_lsn_vs_spk_Laughs_0_vs_Laughs_0 = 0
            overlap_count_lsn_vs_spk_Smiles_0_vs_Laughs_0 = 0
            overlap_count_spk_vs_lsn_Smiles_0_vs_Laughs_0 = 0
            overlap_count_lsn_vs_spk_Laughs_0_vs_Smiles_0 = 0
            overlap_count_spk_vs_lsn_Laughs_0_vs_Smiles_0 = 0
            percentage_spk_vs_lsn_Smiles_0_vs_Smiles_0 = 0
            percentage_spk_vs_lsn_Laughs_0_vs_Laughs_0 = 0
            percentage_lsn_vs_spk_Smiles_0_vs_Smiles_0 = 0
            percentage_lsn_vs_spk_Laughs_0_vs_Laughs_0 = 0
            percentage_lsn_vs_spk_Smiles_0_vs_Laughs_0 = 0
            percentage_spk_vs_lsn_Smiles_0_vs_Laughs_0 = 0
            percentage_lsn_vs_spk_Laughs_0_vs_Smiles_0 = 0
            percentage_spk_vs_lsn_Laughs_0_vs_Smiles_0 = 0
            if person1 == "A":
                count_smiles = pair_dict["Smiles_0 count in lstB"]
                count_laughs = pair_dict["Laughs_0 count in lstB"]
            elif person1 == "B":
                count_smiles = pair_dict["Smiles_0 count in lstA"]
                count_laughs = pair_dict["Laughs_0 count in lstA"]
            segments = pair_dict["Role vs Role"]["Segments"]
            expression = expressions_track
            if person1 == "A":
                for segmentA, segmentB in segments.items():
                    for segB in segmentB:
                        segment_key = f"{segB}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if A is "spk" and B is "lsn"
                            if (segmentA[2].replace(" ", "") == "spk" and segB[2].replace(" ", "") == "lsn"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        globals()[f"overlap_count_spk_vs_lsn_{tierA}_vs_{tierB}"] += 1
                            # Check if A is "lsn" and B is "spk"
                            elif (segmentA[2].replace(" ", "") == "lsn" and segB[2].replace(" ", "") == "spk"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        globals()[f"overlap_count_lsn_vs_spk_{tierA}_vs_{tierB}"] += 1
                        break
            elif person1 == "B":
                for segmentB, segmentA in segments.items():
                    for segA in segmentA:
                        segment_key = f"{segA}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if B is "spk" and A is "lsn"
                            if (segmentB[2].replace(" ", "") == "spk" and segA[2].replace(" ", "") == "lsn"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        globals()[f"overlap_count_spk_vs_lsn_{tierB}_vs_{tierA}"] += 1
                            # Check if B is "lsn" and A is "spk"
                            elif (segmentB[2].replace(" ", "") == "lsn" and segA[2].replace(" ", "") == "spk"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        globals()[f"overlap_count_lsn_vs_spk_{tierB}_vs_{tierA}"] += 1
                        break

            if count_smiles != 0:
                percentage_spk_vs_lsn_Smiles_0_vs_Smiles_0 = overlap_count_spk_vs_lsn_Smiles_0_vs_Smiles_0 / count_smiles * 100        
                percentage_lsn_vs_spk_Smiles_0_vs_Smiles_0 = overlap_count_lsn_vs_spk_Smiles_0_vs_Smiles_0 / count_smiles * 100
                percentage_lsn_vs_spk_Laughs_0_vs_Smiles_0 = overlap_count_lsn_vs_spk_Laughs_0_vs_Smiles_0 / count_smiles * 100
                percentage_spk_vs_lsn_Laughs_0_vs_Smiles_0 = overlap_count_spk_vs_lsn_Laughs_0_vs_Smiles_0 / count_smiles * 100
            if count_laughs != 0:  
                percentage_spk_vs_lsn_Laughs_0_vs_Laughs_0 = overlap_count_spk_vs_lsn_Laughs_0_vs_Laughs_0 / count_laughs * 100
                percentage_lsn_vs_spk_Laughs_0_vs_Laughs_0 = overlap_count_lsn_vs_spk_Laughs_0_vs_Laughs_0 / count_laughs * 100
                percentage_lsn_vs_spk_Smiles_0_vs_Laughs_0 = overlap_count_lsn_vs_spk_Smiles_0_vs_Laughs_0 / count_laughs * 100
                percentage_spk_vs_lsn_Smiles_0_vs_Laughs_0 = overlap_count_spk_vs_lsn_Smiles_0_vs_Laughs_0 / count_laughs * 100

            overlap_count_list.append({
                'Database': database,
                'Pair': pair_name,
                f'Overlap Percentage for {person1} spk / {person2} lsn - Smiles_0 vs Smiles_0': percentage_spk_vs_lsn_Smiles_0_vs_Smiles_0,
                f'Overlap Percentage for {person1} spk / {person2} lsn - Laughs_0 vs Laughs_0': percentage_spk_vs_lsn_Laughs_0_vs_Laughs_0,
                f'Overlap Percentage for {person1} spk / {person2} lsn - Smiles_0 vs Laughs_0': percentage_spk_vs_lsn_Smiles_0_vs_Laughs_0,
                f'Overlap Percentage for {person1} spk / {person2} lsn - Laughs_0 vs Smiles_0': percentage_spk_vs_lsn_Laughs_0_vs_Smiles_0,
                f'Overlap Percentage for {person1} lsn / {person2} spk - Smiles_0 vs Smiles_0': percentage_lsn_vs_spk_Smiles_0_vs_Smiles_0,
                f'Overlap Percentage for {person1} lsn / {person2} spk - Laughs_0 vs Laughs_0': percentage_lsn_vs_spk_Laughs_0_vs_Laughs_0,
                f'Overlap Percentage for {person1} lsn / {person2} spk - Smiles_0 vs Laughs_0': percentage_lsn_vs_spk_Smiles_0_vs_Laughs_0,
                f'Overlap Percentage for {person1} lsn / {person2} spk - Laughs_0 vs Smiles_0': percentage_lsn_vs_spk_Laughs_0_vs_Smiles_0
            }) 
        df_overlap_count = pd.DataFrame(overlap_count_list)
        dataframes[database] = df_overlap_count
    return dataframes

def overlap_total_duration_B(databases_name, databases_pairs, databases_pair_paths, expression_pairs, expressions_track, choice):
    """
    This function returns the duration and percentage of overlapping segments of smiles and laughs between two persons for each pair of files.
    The total duration is calculated by adding all the segments of the tier concerned of the file B (or the second files if you prefered).
    When we talk about overlapping between person A to person B, this means that we are looking at all the 
    segments of person B which overlap a segment of A and this for all the segments of the tier concerned for 
    the "A" files. When we change direction (so person B to person A), we just switched the direction of the 
    files in the overlapping function. (A and B are the pair files)
        - person A to person B: {(segmentA: (segmentB n°1), (segmentB n°N), etc),...}
        - person B to person A: {(segmentB: (segmentA n°1), (segmentA n°N), etc),...}

    Args:
        databases_name (list): list of databases names
        databases_pairs (list): list of databases pairs names
        databases_pair_paths (dict): dictionary of databases pairs paths
        expression_pairs (list): list of expressions pairs names. For example : [("Smiles_0", "Smiles_0"), 
                                                                                 ("Smiles_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Smiles_0"),
                                                                                 ("Role", "Role")]
        expressions_track (list): list of expressions names. For example : [("Smiles_0", "Smiles_0"), 
                                                                            ("Smiles_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Smiles_0")]
        choice (str): choice between "A/B" for  or "B/A". For example : "A/B" for overlapping between person A to person B 
                      or "B/A" for the opposite.

    Returns:
        dataframes (dict): dictionary of dataframes
    """
    if choice == "A/B":
        person1 = "A"
        person2 = "B"
    elif choice == "B/A":
        person1 = "B"
        person2 = "A"

    lstA = {}
    lstB = {}
    overlapping_segments_dict = {}
    for i, database in enumerate(databases_name):
        if database == databases_pairs[i].replace('_pairs', '').upper():
            databases_list = databases_pair_paths[databases_pairs[i]]
            dataset_dict = {}
            for i in range(0, len(databases_list), 2):
                filepath_A = databases_list[i]
                filepath_B = databases_list[i+1] 
                pair_file_A = os.path.basename(filepath_A)
                pair_file_B = os.path.basename(filepath_B)

                if pair_file_A and pair_file_B:
                    if person1 == "A":
                        pair_name = f"{pair_file_A}_&_{pair_file_B}"
                    elif person1 == "B":
                        pair_name = f"{pair_file_B}_&_{pair_file_A}"

                pair_dict = {}
                overlapping_data = {}

                for tier_A, tier_B in expression_pairs:
                    lstA_tier = get_tier_from_file(filepath_A, tier_A)
                    lstB_tier = get_tier_from_file(filepath_B, tier_B)

                    if tier_A in lstA:
                        lstA[tier_A].extend(lstA_tier[tier_A])
                    else:
                        lstA[tier_A] = lstA_tier[tier_A]

                    if tier_B in lstB:
                        lstB[tier_B].extend(lstB_tier[tier_B])
                    else:
                        lstB[tier_B] = lstB_tier[tier_B]

                    if person1 == "A":
                        overlapping_segments = get_overlapping_segments(lstA_tier[tier_A], lstB_tier[tier_B])
                        overlapping_data[f"{tier_A} vs {tier_B}"] = {'Segments': overlapping_segments}

                        overlapping_data[f"{tier_B} duration in lstB"] = 0
                        tiers_in_lstB = set(lstB_tier[tier_B])
                        for seg in tiers_in_lstB:
                            overlapping_data[f"{tier_B} duration in lstB"] += seg[1] - seg[0]
                    elif person1 == "B":
                        overlapping_segments = get_overlapping_segments(lstB_tier[tier_B], lstA_tier[tier_A])
                        overlapping_data[f"{tier_B} vs {tier_A}"] = {'Segments': overlapping_segments}

                        overlapping_data[f"{tier_A} duration in lstA"] = 0
                        tiers_in_lstA = set(lstA_tier[tier_A])
                        for seg in tiers_in_lstA:
                            overlapping_data[f"{tier_A} duration in lstA"] += seg[1] - seg[0]

                dataset_dict[pair_name] = overlapping_data

            overlapping_segments_dict[database] = dataset_dict
    dataframes = {}
    overlap_segments_set = set()
    for database, dataset_dict in overlapping_segments_dict.items():
        overlap_percentage_list = []
        for pair_name, pair_dict in dataset_dict.items():
            overlap_duration_spk_vs_lsn = 0
            overlap_duration_lsn_vs_spk = 0
            percentage_spk_vs_lsn = 0
            percentage_lsn_vs_spk = 0
            if person1 == "A":
                duration_smiles = pair_dict["Smiles_0 duration in lstB"]
                duration_laughs = pair_dict["Laughs_0 duration in lstB"]
            elif person1 == "B":
                duration_smiles = pair_dict["Smiles_0 duration in lstA"]
                duration_laughs = pair_dict["Laughs_0 duration in lstA"]
            segments = pair_dict["Role vs Role"]["Segments"]
            expression = expressions_track
            if person1 == "A":
                for segmentA, segmentB in segments.items():
                    for segB in segmentB:
                        segment_key = f"{segB}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if A is "spk" and B is "lsn"
                            if (segmentA[2].replace(" ", "") == "spk" and segB[2].replace(" ", "") == "lsn"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - b[0]
                            # Check if A is "lsn" and B is "spk"
                            elif (segmentA[2].replace(" ", "") == "lsn" and segB[2].replace(" ", "") == "spk"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - b[0]
                        break
            elif person1 == "B":
                for segmentB, segmentA in segments.items():
                    for segA in segmentA:
                        segment_key = f"{segA}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if B is "spk" and A is "lsn"
                            if (segmentB[2].replace(" ", "") == "spk" and segA[2].replace(" ", "") == "lsn"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - a[0] 
                            # Check if B is "lsn" and A is "spk"
                            elif (segmentB[2].replace(" ", "") == "lsn" and segA[2].replace(" ", "") == "spk"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - a[0]
                        break

            if duration_smiles != 0 or duration_laughs !=0 :
                percentage_spk_vs_lsn = overlap_duration_spk_vs_lsn / (duration_smiles + duration_laughs) * 100
                percentage_lsn_vs_spk = overlap_duration_lsn_vs_spk / (duration_smiles + duration_laughs) * 100
            overlap_percentage_list.append({
                'Database': database,
                'Pair': pair_name,
                f'Overlap Duration for {person1} spk / {person2} lsn - S&L': overlap_duration_spk_vs_lsn,
                f'Overlap Percentage for {person1} spk / {person2} lsn - S&L': percentage_spk_vs_lsn,
                f'Overlap Duration for {person1} lsn / {person2} spk - S&L': overlap_duration_lsn_vs_spk,
                f'Overlap Percentage for {person1} lsn / {person2} spk - S&L': percentage_lsn_vs_spk,
            })
        df_overlap_count = pd.DataFrame(overlap_percentage_list)
        dataframes[database] = df_overlap_count
    return dataframes

def overlap_total_duration_A(databases_name, databases_pairs, databases_pair_paths, expression_pairs, expressions_track, choice):
    """
    This function returns the duration and percentage of overlapping segments of smiles and laughs between two persons for each pair of files.
    The total duration is calculated by adding all the segments of the tier concerned of the file A (or the first files if you prefered).
    When we talk about overlapping between person A to person B, this means that we are looking at all the 
    segments of person B which overlap a segment of A and this for all the segments of the tier concerned for 
    the "A" files. When we change direction (so person B to person A), we just switched the direction of the 
    files in the overlapping function. (A and B are the pair files)
        - person A to person B: {(segmentA: (segmentB n°1), (segmentB n°N), etc),...}
        - person B to person A: {(segmentB: (segmentA n°1), (segmentA n°N), etc),...}

    Args:
        databases_name (list): list of databases names
        databases_pairs (list): list of databases pairs names
        databases_pair_paths (dict): dictionary of databases pairs paths
        expression_pairs (list): list of expressions pairs names. For example : [("Smiles_0", "Smiles_0"), 
                                                                                 ("Smiles_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Smiles_0"),
                                                                                 ("Role", "Role")]
        expressions_track (list): list of expressions names. For example : [("Smiles_0", "Smiles_0"), 
                                                                            ("Smiles_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Smiles_0")]
        choice (str): choice between "A/B" for  or "B/A". For example : "A/B" for overlapping between person A to person B 
                      or "B/A" for the opposite.

    Returns:
        dataframes (dict): dictionary of dataframes
    """
    if choice == "A/B":
        person1 = "A"
        person2 = "B"
    elif choice == "B/A":
        person1 = "B"
        person2 = "A"

    lstA = {}
    lstB = {}
    overlapping_segments_dict = {}
    for i, database in enumerate(databases_name):
        if database == databases_pairs[i].replace('_pairs', '').upper():
            databases_list = databases_pair_paths[databases_pairs[i]]
            dataset_dict = {}
            for i in range(0, len(databases_list), 2):
                filepath_A = databases_list[i]
                filepath_B = databases_list[i+1] 
                pair_file_A = os.path.basename(filepath_A)
                pair_file_B = os.path.basename(filepath_B)

                if pair_file_A and pair_file_B:
                    if person1 == "A":
                        pair_name = f"{pair_file_A}_&_{pair_file_B}"
                    elif person1 == "B":
                        pair_name = f"{pair_file_B}_&_{pair_file_A}"

                pair_dict = {}
                overlapping_data = {}

                for tier_A, tier_B in expression_pairs:
                    lstA_tier = get_tier_from_file(filepath_A, tier_A)
                    lstB_tier = get_tier_from_file(filepath_B, tier_B)

                    if tier_A in lstA:
                        lstA[tier_A].extend(lstA_tier[tier_A])
                    else:
                        lstA[tier_A] = lstA_tier[tier_A]

                    if tier_B in lstB:
                        lstB[tier_B].extend(lstB_tier[tier_B])
                    else:
                        lstB[tier_B] = lstB_tier[tier_B]

                    if person1 == "A":
                        overlapping_segments = get_overlapping_segments(lstA_tier[tier_A], lstB_tier[tier_B])
                        overlapping_data[f"{tier_A} vs {tier_B}"] = {'Segments': overlapping_segments}

                        overlapping_data[f"{tier_A} duration in lstA"] = 0
                        tiers_in_lstA = set(lstA_tier[tier_A])
                        for seg in tiers_in_lstA:
                            overlapping_data[f"{tier_A} duration in lstA"] += seg[1] - seg[0]
                    elif person1 == "B":
                        overlapping_segments = get_overlapping_segments(lstB_tier[tier_B], lstA_tier[tier_A])
                        overlapping_data[f"{tier_B} vs {tier_A}"] = {'Segments': overlapping_segments}

                        overlapping_data[f"{tier_B} duration in lstB"] = 0
                        tiers_in_lstB = set(lstB_tier[tier_B])
                        for seg in tiers_in_lstB:
                            overlapping_data[f"{tier_B} duration in lstB"] += seg[1] - seg[0]

                dataset_dict[pair_name] = overlapping_data

            overlapping_segments_dict[database] = dataset_dict
    dataframes = {}
    overlap_segments_set = set()
    for database, dataset_dict in overlapping_segments_dict.items():
        overlap_percentage_list = []
        for pair_name, pair_dict in dataset_dict.items():
            overlap_duration_spk_vs_lsn = 0
            overlap_duration_lsn_vs_spk = 0
            percentage_spk_vs_lsn = 0
            percentage_lsn_vs_spk = 0
            if person1 == "A":
                duration_smiles = pair_dict["Smiles_0 duration in lstA"]
                duration_laughs = pair_dict["Laughs_0 duration in lstA"]
            elif person1 == "B":
                duration_smiles = pair_dict["Smiles_0 duration in lstB"]
                duration_laughs = pair_dict["Laughs_0 duration in lstB"]
            segments = pair_dict["Role vs Role"]["Segments"]
            expression = expressions_track
            if person1 == "A":
                for segmentA, segmentB in segments.items():
                    for segB in segmentB:
                        segment_key = f"{segB}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if A is "spk" and B is "lsn"
                            if (segmentA[2].replace(" ", "") == "spk" and segB[2].replace(" ", "") == "lsn"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - b[0]
                            # Check if A is "lsn" and B is "spk"
                            elif (segmentA[2].replace(" ", "") == "lsn" and segB[2].replace(" ", "") == "spk"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - b[0]
                        break
            elif person1 == "B":
                for segmentB, segmentA in segments.items():
                    for segA in segmentA:
                        segment_key = f"{segA}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if B is "spk" and A is "lsn"
                            if (segmentB[2].replace(" ", "") == "spk" and segA[2].replace(" ", "") == "lsn"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - a[0] 
                            # Check if B is "lsn" and A is "spk"
                            elif (segmentB[2].replace(" ", "") == "lsn" and segA[2].replace(" ", "") == "spk"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - a[0]
                        break

            if duration_smiles != 0 or duration_laughs !=0 :
                percentage_spk_vs_lsn = overlap_duration_spk_vs_lsn / (duration_smiles + duration_laughs) * 100
                percentage_lsn_vs_spk = overlap_duration_lsn_vs_spk / (duration_smiles + duration_laughs) * 100
            overlap_percentage_list.append({
                'Database': database,
                'Pair': pair_name,
                f'Overlap Duration for {person1} spk / {person2} lsn - S&L': overlap_duration_spk_vs_lsn,
                f'Overlap Percentage for {person1} spk / {person2} lsn - S&L': percentage_spk_vs_lsn,
                f'Overlap Duration for {person1} lsn / {person2} spk - S&L': overlap_duration_lsn_vs_spk,
                f'Overlap Percentage for {person1} lsn / {person2} spk - S&L': percentage_lsn_vs_spk,
            })
        df_overlap_count = pd.DataFrame(overlap_percentage_list)
        dataframes[database] = df_overlap_count
    return dataframes

def overlap_total_duration_union(databases_name, databases_pairs, databases_pair_paths, expression_pairs, expressions_track, choice):
    """
    This function returns the duration and percentage of overlapping segments of smiles and laughs between two persons for each pair of files.
    The total duration is calculated by adding all the segments of the tier concerned of the two pair files (we do the union of the segments).
    When we talk about overlapping between person A to person B, this means that we are looking at all the 
    segments of person B which overlap a segment of A and this for all the segments of the tier concerned for 
    the "A" files. When we change direction (so person B to person A), we just switched the direction of the 
    files in the overlapping function. (A and B are the pair files)
        - person A to person B: {(segmentA: (segmentB n°1), (segmentB n°N), etc),...}
        - person B to person A: {(segmentB: (segmentA n°1), (segmentA n°N), etc),...}

    Args:
        databases_name (list): list of databases names
        databases_pairs (list): list of databases pairs names
        databases_pair_paths (dict): dictionary of databases pairs paths
        expression_pairs (list): list of expressions pairs names. For example : [("Smiles_0", "Smiles_0"), 
                                                                                 ("Smiles_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Laughs_0"), 
                                                                                 ("Laughs_0", "Smiles_0"),
                                                                                 ("Role", "Role")]
        expressions_track (list): list of expressions names. For example : [("Smiles_0", "Smiles_0"), 
                                                                            ("Smiles_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Smiles_0")]
        choice (str): choice between "A/B" for  or "B/A". For example : "A/B" for overlapping between person A to person B 
                      or "B/A" for the opposite.

    Returns:
        dataframes (dict): dictionary of dataframes
    """
    if choice == "A/B":
        person1 = "A"
        person2 = "B"
    elif choice == "B/A":
        person1 = "B"
        person2 = "A"

    lstA = {}
    lstB = {}
    overlapping_segments_dict = {}
    for i, database in enumerate(databases_name):
        if database == databases_pairs[i].replace('_pairs', '').upper():
            databases_list = databases_pair_paths[databases_pairs[i]]
            dataset_dict = {}
            for i in range(0, len(databases_list), 2):
                filepath_A = databases_list[i]
                filepath_B = databases_list[i+1] 
                pair_file_A = os.path.basename(filepath_A)
                pair_file_B = os.path.basename(filepath_B)

                if pair_file_A and pair_file_B:
                    if person1 == "A":
                        pair_name = f"{pair_file_A}_&_{pair_file_B}"
                    elif person1 == "B":
                        pair_name = f"{pair_file_B}_&_{pair_file_A}"

                pair_dict = {}
                overlapping_data = {}

                for tier_A, tier_B in expression_pairs:
                    lstA_tier = get_tier_from_file(filepath_A, tier_A)
                    lstB_tier = get_tier_from_file(filepath_B, tier_B)

                    if tier_A in lstA:
                        lstA[tier_A].extend(lstA_tier[tier_A])
                    else:
                        lstA[tier_A] = lstA_tier[tier_A]

                    if tier_B in lstB:
                        lstB[tier_B].extend(lstB_tier[tier_B])
                    else:
                        lstB[tier_B] = lstB_tier[tier_B]

                    if person1 == "A":
                        overlapping_segments = get_overlapping_segments(lstA_tier[tier_A], lstB_tier[tier_B])
                        overlapping_data[f"{tier_A} vs {tier_B}"] = {'Segments': overlapping_segments}

                    elif person1 == "B":
                        overlapping_segments = get_overlapping_segments(lstB_tier[tier_B], lstA_tier[tier_A])
                        overlapping_data[f"{tier_B} vs {tier_A}"] = {'Segments': overlapping_segments}

                    overlapping_data[f"Total duration"] = 0
                    tiers_in_lstA = set(lstA_tier[tier_A])
                    tiers_in_lstB = set(lstB_tier[tier_B])
                    for segA, segB in zip(tiers_in_lstA, tiers_in_lstB):
                        if segA[0] < segB[0] and segA[1] > segB[1]:
                            overlapping_data[f"Total duration"] += segA[1] - segA[0]
                        elif segB[0] < segA[1] and segB[1] > segA[0]:
                            overlapping_data[f"Total duration"] += segB[1] - segB[0]
                        elif segA[0] > segB[1] or segA[1] < segB[0]:
                            overlapping_data[f"Total duration"] += segA[1] - segA[0]
                        elif segB[0] > segA[1] or segB[1] < segA[0]:
                            overlapping_data[f"Total duration"] += segB[1] - segB[0]
                        elif segA[0] < segB[0] and segA[1] < segB[1]:
                            overlapping_data[f"Total duration"] += (segB[1] - segA[0]) - (segA[1] - segB[0])
                        elif segB[0] < segA[0] and segB[1] < segA[1]:
                            overlapping_data[f"Total duration"] += (segA[1] - segB[0]) - (segB[1] - segA[0])

                dataset_dict[pair_name] = overlapping_data

            overlapping_segments_dict[database] = dataset_dict
    dataframes = {}
    overlap_segments_set = set()
    for database, dataset_dict in overlapping_segments_dict.items():
        overlap_percentage_list = []
        for pair_name, pair_dict in dataset_dict.items():
            overlap_duration_spk_vs_lsn = 0
            overlap_duration_lsn_vs_spk = 0
            percentage_spk_vs_lsn = 0
            percentage_lsn_vs_spk = 0
            total_duration = pair_dict["Total duration"]
            segments = pair_dict["Role vs Role"]["Segments"]
            expression = expressions_track
            if person1 == "A":
                for segmentA, segmentB in segments.items():
                    for segB in segmentB:
                        segment_key = f"{segB}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if A is "spk" and B is "lsn"
                            if (segmentA[2].replace(" ", "") == "spk" and segB[2].replace(" ", "") == "lsn"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - b[0]
                            # Check if A is "lsn" and B is "spk"
                            elif (segmentA[2].replace(" ", "") == "lsn" and segB[2].replace(" ", "") == "spk"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - b[0]
                        break
            elif person1 == "B":
                for segmentB, segmentA in segments.items():
                    for segA in segmentA:
                        segment_key = f"{segA}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if B is "spk" and A is "lsn"
                            if (segmentB[2].replace(" ", "") == "spk" and segA[2].replace(" ", "") == "lsn"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - a[0] 
                            # Check if B is "lsn" and A is "spk"
                            elif (segmentB[2].replace(" ", "") == "lsn" and segA[2].replace(" ", "") == "spk"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - a[0]
                        break

            if total_duration != 0 :
                percentage_spk_vs_lsn = overlap_duration_spk_vs_lsn / total_duration * 100
                percentage_lsn_vs_spk = overlap_duration_lsn_vs_spk / total_duration * 100
            overlap_percentage_list.append({
                'Database': database,
                'Pair': pair_name,
                f'Overlap Duration for {person1} spk / {person2} lsn - S&L': overlap_duration_spk_vs_lsn,
                f'Overlap Percentage for {person1} spk / {person2} lsn - S&L': percentage_spk_vs_lsn,
                f'Overlap Duration for {person1} lsn / {person2} spk - S&L': overlap_duration_lsn_vs_spk,
                f'Overlap Percentage for {person1} lsn / {person2} spk - S&L': percentage_lsn_vs_spk,
            })
        df_overlap_count = pd.DataFrame(overlap_percentage_list)
        dataframes[database] = df_overlap_count
    return dataframes

def overlap_percentage_B(databases_name, databases_pairs, databases_pair_paths, expression_pairs, expressions_track, choice):
    '''
    This function calculates the overlap percentage between two person for each pair of files.
    The total duration is calculated by adding the duration of the overlapping segments of the second person.
    When we talk about overlapping between person A to person B, this means that we are looking at all the 
    segments of person B which overlap a segment of A and this for all the segments of the tier concerned for 
    the "A" files. When we change direction (so person B to person A), we just switched the direction of the 
    files in the overlapping function. (A and B are the pair files)
        - person A to person B: {(segmentA: (segmentB n°1), (segmentB n°N), etc),...}
        - person B to person A: {(segmentB: (segmentA n°1), (segmentA n°N), etc),...}

    Args:
        databases_name: list of databases names.
        databases_pairs: list of databases pairs names.
        databases_pair_paths: dictionary of databases pairs paths.
        expression_pairs: list of expression pairs to study. For example : [("Smiles_0", "Smiles_0"), 
                                                                            ("Smiles_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Smiles_0")].
        expressions_track: list of expression pairs to study. For example : [("Smiles_0", "Smiles_0"),
                                                                            ("Smiles_0", "Laughs_0"),
                                                                            ("Laughs_0", "Laughs_0"),
                                                                            ("Laughs_0", "Smiles_0")].
        choice: chchoice between "A/B" or "B/A". For example : "A/B" for overlapping between person A to person B 
                or "B/A" for the opposite.
        
    Returns:
        dataframes: dictionary of dataframes containing the overlap percentage and duration for each pair of 
                    files.
    '''

    if choice == "A/B":
        person1 = "A"
        person2 = "B"
    elif choice == "B/A":
        person1 = "B"
        person2 = "A"

    lstA = {}
    lstB = {}

    overlapping_segments_dict = {}

    for i, database in enumerate(databases_name):
        if database == databases_pairs[i].replace('_pairs', '').upper():
            databases_list = databases_pair_paths[databases_pairs[i]]
            dataset_dict = {}

            for i in range(0, len(databases_list), 2):
                filepath_A = databases_list[i]
                filepath_B = databases_list[i+1]
                pair_file_A = os.path.basename(filepath_A)
                pair_file_B = os.path.basename(filepath_B)

                if pair_file_A and pair_file_B:
                    if person1 == "A":
                        pair_name = f"{pair_file_A}_&_{pair_file_B}"
                    elif person1 == "B":
                        pair_name = f"{pair_file_B}_&_{pair_file_A}"
                pair_dict = {}
                overlapping_data = {}
                for tier_A, tier_B in expression_pairs:
                    lstA_tier = get_tier_from_file(filepath_A, tier_A)
                    lstB_tier = get_tier_from_file(filepath_B, tier_B)

                    if tier_A in lstA:
                        lstA[tier_A].extend(lstA_tier[tier_A])
                    else:
                        lstA[tier_A] = lstA_tier[tier_A]

                    if tier_B in lstB:
                        lstB[tier_B].extend(lstB_tier[tier_B])
                    else:
                        lstB[tier_B] = lstB_tier[tier_B]
                    
                    if person1 == "A":
                        overlapping_segments = get_overlapping_segments(lstA_tier[tier_A], lstB_tier[tier_B])
                        overlapping_data[f"{tier_A} vs {tier_B}"] = {'Segments': overlapping_segments}
                    else:
                        overlapping_segments = get_overlapping_segments(lstB_tier[tier_B], lstA_tier[tier_A])
                        overlapping_data[f"{tier_B} vs {tier_A}"] = {'Segments': overlapping_segments}


                dataset_dict[pair_name] = overlapping_data

            overlapping_segments_dict[database] = dataset_dict

    dataframes = {}
    overlap_segments_set = set()
    for database, dataset_dict in overlapping_segments_dict.items():
        overlap_percentage_list = []
        for pair_name, pair_dict in dataset_dict.items():
            overlap_duration_spk_vs_lsn = 0
            overlap_duration_lsn_vs_spk = 0
            percentage_spk_vs_lsn = 0
            percentage_lsn_vs_spk = 0
            duration = 0
            segments = pair_dict["Role vs Role"]["Segments"]
            expression = expressions_track
            if person1 == "A":
                for segmentA, segmentB in segments.items():
                    for segB in segmentB:
                        segment_key = f"{segB}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if A is "spk" and B is "lsn"
                            if (segmentA[2].replace(" ", "") == "spk" and segB[2].replace(" ", "") == "lsn"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        duration += b[1] - b[0]
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - b[0]
                            # Check if A is "lsn" and B is "spk"
                            elif (segmentA[2].replace(" ", "") == "lsn" and segB[2].replace(" ", "") == "spk"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        duration += b[1] - b[0]
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - b[0]     
                        break
            else: 
                for segmentB, segmentA in segments.items():
                    for segA in segmentA:
                        segment_key = f"{segA}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if B is "spk" and A is "lsn"
                            if (segmentB[2].replace(" ", "") == "spk" and segA[2].replace(" ", "") == "lsn"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        duration += a[1] - a[0]
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - a[0]
                            # Check if B is "lsn" and A is "spk"
                            elif (segmentB[2].replace(" ", "") == "lsn" and segA[2].replace(" ", "") == "spk"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        duration += a[1] - a[0]
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - a[0]       
                        break
            if duration != 0 :
                percentage_spk_vs_lsn = overlap_duration_spk_vs_lsn / (duration) * 100
                percentage_lsn_vs_spk = overlap_duration_lsn_vs_spk / (duration) * 100
            overlap_percentage_list.append({
                'Database': database,
                'Pair': pair_name,
                f'Overlap Duration for {person1} spk / {person2} lsn - S&L': overlap_duration_spk_vs_lsn,
                f'Overlap Percentage for {person1} spk / {person2} lsn - S&L': percentage_spk_vs_lsn,
                f'Overlap Duration for {person1} lsn / {person2} spk - S&L': overlap_duration_lsn_vs_spk,
                f'Overlap Percentage for {person1} lsn / {person2} spk - S&L': percentage_lsn_vs_spk,
            })
        df_overlap_percentage = pd.DataFrame(overlap_percentage_list)
        dataframes[database] = df_overlap_percentage
    return dataframes

def overlap_percentage_A(databases_name, databases_pairs, databases_pair_paths, expression_pairs, expressions_track, choice):
    '''
    This function calculates the overlap percentage between two person for each pair of files.
    The total duration is calculated by adding the duration of the overlapping segments of the first person.
    When we talk about overlapping between person A to person B, this means that we are looking at all the 
    segments of person B which overlap a segment of A and this for all the segments of the tier concerned for 
    the "A" files. When we change direction (so person B to person A), we just switched the direction of the 
    files in the overlapping function. (A and B are the pair files)
        - person A to person B: {(segmentA: (segmentB n°1), (segmentB n°N), etc),...}
        - person B to person A: {(segmentB: (segmentA n°1), (segmentA n°N), etc),...}

    Args:
        databases_name: list of databases names.
        databases_pairs: list of databases pairs names.
        databases_pair_paths: dictionary of databases pairs paths.
        expression_pairs: list of expression pairs to study. For example : [("Smiles_0", "Smiles_0"), 
                                                                            ("Smiles_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Smiles_0")].
        expressions_track: list of expression pairs to study. For example : [("Smiles_0", "Smiles_0"),
                                                                            ("Smiles_0", "Laughs_0"),
                                                                            ("Laughs_0", "Laughs_0"),
                                                                            ("Laughs_0", "Smiles_0")].
        choice: chchoice between "A/B" or "B/A". For example : "A/B" for overlapping between person A to person B 
                or "B/A" for the opposite.
        
    Returns:
        dataframes: dictionary of dataframes containing the overlap percentage and duration for each pair of 
                    files.
    '''

    if choice == "A/B":
        person1 = "A"
        person2 = "B"
    elif choice == "B/A":
        person1 = "B"
        person2 = "A"

    lstA = {}
    lstB = {}

    overlapping_segments_dict = {}

    for i, database in enumerate(databases_name):
        if database == databases_pairs[i].replace('_pairs', '').upper():
            databases_list = databases_pair_paths[databases_pairs[i]]
            dataset_dict = {}

            for i in range(0, len(databases_list), 2):
                filepath_A = databases_list[i]
                filepath_B = databases_list[i+1]
                pair_file_A = os.path.basename(filepath_A)
                pair_file_B = os.path.basename(filepath_B)

                if pair_file_A and pair_file_B:
                    if person1 == "A":
                        pair_name = f"{pair_file_A}_&_{pair_file_B}"
                    elif person1 == "B":
                        pair_name = f"{pair_file_B}_&_{pair_file_A}"
                pair_dict = {}
                overlapping_data = {}
                for tier_A, tier_B in expression_pairs:
                    lstA_tier = get_tier_from_file(filepath_A, tier_A)
                    lstB_tier = get_tier_from_file(filepath_B, tier_B)

                    if tier_A in lstA:
                        lstA[tier_A].extend(lstA_tier[tier_A])
                    else:
                        lstA[tier_A] = lstA_tier[tier_A]

                    if tier_B in lstB:
                        lstB[tier_B].extend(lstB_tier[tier_B])
                    else:
                        lstB[tier_B] = lstB_tier[tier_B]
                    
                    if person1 == "A":
                        overlapping_segments = get_overlapping_segments(lstA_tier[tier_A], lstB_tier[tier_B])
                        overlapping_data[f"{tier_A} vs {tier_B}"] = {'Segments': overlapping_segments}
                    elif person1 == "B":
                        overlapping_segments = get_overlapping_segments(lstB_tier[tier_B], lstA_tier[tier_A])
                        overlapping_data[f"{tier_B} vs {tier_A}"] = {'Segments': overlapping_segments}


                dataset_dict[pair_name] = overlapping_data

            overlapping_segments_dict[database] = dataset_dict

    dataframes = {}
    overlap_segments_set = set()
    for database, dataset_dict in overlapping_segments_dict.items():
        overlap_percentage_list = []
        for pair_name, pair_dict in dataset_dict.items():
            overlap_duration_spk_vs_lsn = 0
            overlap_duration_lsn_vs_spk = 0
            percentage_spk_vs_lsn = 0
            percentage_lsn_vs_spk = 0
            duration = 0
            segments = pair_dict["Role vs Role"]["Segments"]
            expression = expressions_track
            if person1 == "A":
                for segmentA, segmentB in segments.items():
                    for segB in segmentB:
                        segment_key = f"{segB}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if A is "spk" and B is "lsn"
                            if (segmentA[2].replace(" ", "") == "spk" and segB[2].replace(" ", "") == "lsn"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        duration += A[1] - A[0]
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - b[0]
                            # Check if A is "lsn" and B is "spk"
                            elif (segmentA[2].replace(" ", "") == "lsn" and segB[2].replace(" ", "") == "spk"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        duration += A[1] - A[0]
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - A[0]
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - b[0]     
                        break
            elif person1 == "B": 
                for segmentB, segmentA in segments.items():
                    for segA in segmentA:
                        segment_key = f"{segA}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if B is "spk" and A is "lsn"
                            if (segmentB[2].replace(" ", "") == "spk" and segA[2].replace(" ", "") == "lsn"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        duration += B[1] - B[0]
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - a[0]
                            # Check if B is "lsn" and A is "spk"
                            elif (segmentB[2].replace(" ", "") == "lsn" and segA[2].replace(" ", "") == "spk"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        duration += B[1] - B[0]
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - B[0]
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - a[0]       
                        break
            if duration != 0 :
                percentage_spk_vs_lsn = overlap_duration_spk_vs_lsn / (duration) * 100
                percentage_lsn_vs_spk = overlap_duration_lsn_vs_spk / (duration) * 100
            overlap_percentage_list.append({
                'Database': database,
                'Pair': pair_name,
                f'Overlap Duration for {person1} spk / {person2} lsn - S&L': overlap_duration_spk_vs_lsn,
                f'Overlap Percentage for {person1} spk / {person2} lsn - S&L': percentage_spk_vs_lsn,
                f'Overlap Duration for {person1} lsn / {person2} spk - S&L': overlap_duration_lsn_vs_spk,
                f'Overlap Percentage for {person1} lsn / {person2} spk - S&L': percentage_lsn_vs_spk,
            })
        df_overlap_percentage = pd.DataFrame(overlap_percentage_list)
        dataframes[database] = df_overlap_percentage
    return dataframes

def overlap_percentage_union(databases_name, databases_pairs, databases_pair_paths, expression_pairs, expressions_track, choice):
    '''
    This function calculates the overlap percentage between two person for each pair of files.
    The total duration is calculated by adding the duration of the two overlapping segments of the two person.
    When we talk about overlapping between person A to person B, this means that we are looking at all the 
    segments of person B which overlap a segment of A and this for all the segments of the tier concerned for 
    the "A" files. When we change direction (so person B to person A), we just switched the direction of the 
    files in the overlapping function. (A and B are the pair files)
        - person A to person B: {(segmentA: (segmentB n°1), (segmentB n°N), etc),...}
        - person B to person A: {(segmentB: (segmentA n°1), (segmentA n°N), etc),...}

    Args:
        databases_name: list of databases names.
        databases_pairs: list of databases pairs names.
        databases_pair_paths: dictionary of databases pairs paths.
        expression_pairs: list of expression pairs to study. For example : [("Smiles_0", "Smiles_0"), 
                                                                            ("Smiles_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Laughs_0"), 
                                                                            ("Laughs_0", "Smiles_0")].
        expressions_track: list of expressions names. For example : [("Smiles_0", "Smiles_0"),
                                                                    ("Smiles_0", "Laughs_0"),               
                                                                    ("Laughs_0", "Laughs_0"),
                                                                    ("Laughs_0", "Smiles_0")].
        choice: chchoice between "A/B" or "B/A". For example : "A/B" for overlapping between person A to person B 
                or "B/A" for the opposite.
        
    Returns:
        dataframes: dictionary of dataframes containing the overlap percentage and duration for each pair of 
                    files.
    '''

    if choice == "A/B":
        person1 = "A"
        person2 = "B"
    elif choice == "B/A":
        person1 = "B"
        person2 = "A"

    lstA = {}
    lstB = {}

    overlapping_segments_dict = {}

    for i, database in enumerate(databases_name):
        if database == databases_pairs[i].replace('_pairs', '').upper():
            databases_list = databases_pair_paths[databases_pairs[i]]
            dataset_dict = {}

            for i in range(0, len(databases_list), 2):
                filepath_A = databases_list[i]
                filepath_B = databases_list[i+1]
                pair_file_A = os.path.basename(filepath_A)
                pair_file_B = os.path.basename(filepath_B)

                if pair_file_A and pair_file_B:
                    if person1 == "A":
                        pair_name = f"{pair_file_A}_&_{pair_file_B}"
                    elif person1 == "B":
                        pair_name = f"{pair_file_B}_&_{pair_file_A}"
                pair_dict = {}
                overlapping_data = {}
                for tier_A, tier_B in expression_pairs:
                    lstA_tier = get_tier_from_file(filepath_A, tier_A)
                    lstB_tier = get_tier_from_file(filepath_B, tier_B)

                    if tier_A in lstA:
                        lstA[tier_A].extend(lstA_tier[tier_A])
                    else:
                        lstA[tier_A] = lstA_tier[tier_A]

                    if tier_B in lstB:
                        lstB[tier_B].extend(lstB_tier[tier_B])
                    else:
                        lstB[tier_B] = lstB_tier[tier_B]
                    
                    if person1 == "A":
                        overlapping_segments = get_overlapping_segments(lstA_tier[tier_A], lstB_tier[tier_B])
                        overlapping_data[f"{tier_A} vs {tier_B}"] = {'Segments': overlapping_segments}
                    else:
                        overlapping_segments = get_overlapping_segments(lstB_tier[tier_B], lstA_tier[tier_A])
                        overlapping_data[f"{tier_B} vs {tier_A}"] = {'Segments': overlapping_segments}


                dataset_dict[pair_name] = overlapping_data

            overlapping_segments_dict[database] = dataset_dict

    dataframes = {}
    overlap_segments_set = set()
    for database, dataset_dict in overlapping_segments_dict.items():
        overlap_percentage_list = []
        for pair_name, pair_dict in dataset_dict.items():
            overlap_duration_spk_vs_lsn = 0
            overlap_duration_lsn_vs_spk = 0
            percentage_spk_vs_lsn = 0
            percentage_lsn_vs_spk = 0
            duration = 0
            segments = pair_dict["Role vs Role"]["Segments"]
            expression = expressions_track
            if person1 == "A":
                for segmentA, segmentB in segments.items():
                    for segB in segmentB:
                        segment_key = f"{segB}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if A is "spk" and B is "lsn"
                            if (segmentA[2].replace(" ", "") == "spk" and segB[2].replace(" ", "") == "lsn"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - b[0]
                                                            duration += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - A[0]
                                                            duration += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_spk_vs_lsn += b[1] - A[0]
                                                            duration += (A[1] - b[0]) - (b[1] - A[0])
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_spk_vs_lsn += A[1] - b[0]
                                                            duration += (b[1] - A[0]) - (A[1] - b[0])
                            # Check if A is "lsn" and B is "spk"
                            elif (segmentA[2].replace(" ", "") == "lsn" and segB[2].replace(" ", "") == "spk"):
                                for tierA, tierB in expression:
                                    segments_tier = pair_dict[f"{tierA} vs {tierB}"]["Segments"]   
                                    for A, B in segments_tier.items():
                                        if A[0] < segB[1] and A[1] > segB[0]:
                                            for b in B:
                                                if b[0] < segB[1] and b[1] > segB[0]:
                                                    tier_key = f"{b}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if b[0] > A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - b[0]
                                                            duration += A[1] - A[0]
                                                        elif b[0] < A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - A[0]
                                                            duration += b[1] - b[0]
                                                        elif b[0] < A[0] and b[1] < A[1]:
                                                            overlap_duration_lsn_vs_spk += b[1] - A[0]
                                                            duration += (A[1] - b[0]) - (b[1] - A[0])
                                                        elif b[0] > A[0] and b[1] > A[1]:
                                                            overlap_duration_lsn_vs_spk += A[1] - b[0]    
                                                            duration += (b[1] - A[0]) - (A[1] - b[0])     
                        break
            else: 
                for segmentB, segmentA in segments.items():
                    for segA in segmentA:
                        segment_key = f"{segA}"
                        if segment_key not in overlap_segments_set:
                            overlap_segments_set.add(segment_key)
                            # Check if B is "spk" and A is "lsn"
                            if (segmentB[2].replace(" ", "") == "spk" and segA[2].replace(" ", "") == "lsn"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - a[0]
                                                            duration += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - B[0]
                                                            duration += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_spk_vs_lsn += a[1] - B[0]
                                                            duration += (B[1] - a[0]) - (a[1] - B[0])
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_spk_vs_lsn += B[1] - a[0]
                                                            duration += (a[1] - B[0]) - (B[1] - a[0])
                            # Check if B is "lsn" and A is "spk"
                            elif (segmentB[2].replace(" ", "") == "lsn" and segA[2].replace(" ", "") == "spk"):
                                for tierB, tierA in expression:
                                    segments_tier = pair_dict[f"{tierB} vs {tierA}"]["Segments"]   
                                    for B, A in segments_tier.items():
                                        if B[0] < segA[1] and B[1] > segA[0]:
                                            for a in A:
                                                if a[0] < segA[1] and a[1] > segA[0]:
                                                    tier_key = f"{a}"
                                                    if tier_key not in overlap_segments_set:
                                                        overlap_segments_set.add(tier_key)
                                                        if a[0] > B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - a[0]
                                                            duration += B[1] - B[0]
                                                        elif a[0] < B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - B[0]
                                                            duration += a[1] - a[0]
                                                        elif a[0] < B[0] and a[1] < B[1]:
                                                            overlap_duration_lsn_vs_spk += a[1] - B[0]
                                                            duration += (B[1] - a[0]) - (a[1] - B[0])
                                                        elif a[0] > B[0] and a[1] > B[1]:
                                                            overlap_duration_lsn_vs_spk += B[1] - a[0] 
                                                            duration += (a[1] - B[0]) - (B[1] - a[0])        
                        break
            if duration != 0 :
                percentage_spk_vs_lsn = overlap_duration_spk_vs_lsn / (duration) * 100
                percentage_lsn_vs_spk = overlap_duration_lsn_vs_spk / (duration) * 100
            overlap_percentage_list.append({
                'Database': database,
                'Pair': pair_name,
                f'Overlap Duration for {person1} spk / {person2} lsn - S&L': overlap_duration_spk_vs_lsn,
                f'Overlap Percentage for {person1} spk / {person2} lsn - S&L': percentage_spk_vs_lsn,
                f'Overlap Duration for {person1} lsn / {person2} spk - S&L': overlap_duration_lsn_vs_spk,
                f'Overlap Percentage for {person1} lsn / {person2} spk - S&L': percentage_lsn_vs_spk,
            })
        df_overlap_percentage = pd.DataFrame(overlap_percentage_list)
        dataframes[database] = df_overlap_percentage
    return dataframes

## TO IMPROVE ##
def overlap_count(databases_name, databases_pairs, databases_pair_paths, choice, tier="Role"):
    """
    This function returns the number of overlapping segments of a specific tier between two persons for each pair of files.
    When we talk about overlapping between person A to person B, this means that we are looking at all the 
    segments of person B which overlap a segment of A and this for all the segments of the tier concerned for 
    the "A" files. When we change direction (so person B to person A), we just switched the direction of the 
    files in the overlapping function. (A and B are the pair files)
        - person A to person B: {(segmentA: (segmentB n°1), (segmentB n°N), etc),...}
        - person B to person A: {(segmentB: (segmentA n°1), (segmentA n°N), etc),...}

    Args:
        databases_name (list): list of databases names
        databases_pairs (list): list of databases pairs names
        databases_pair_paths (dict): dictionary of databases pairs paths
        choice (str): choice between "A/B" for  or "B/A". For example : "A/B" for overlapping between person A to person B 
                      or "B/A" for the opposite.
        tier (str): tier name. Default: "Role"

    Returns:
        dataframes (dict): dictionary of dataframes
    """
    if choice == "A/B":
        person1 = "A"
        person2 = "B"
    elif choice == "B/A":
        person1 = "B"
        person2 = "A"

    lstA = {}
    lstB = {}
    overlapping_segments_dict = {}
    for i, database in enumerate(databases_name):
        if database == databases_pairs[i].replace('_pairs', '').upper():
            databases_list = databases_pair_paths[databases_pairs[i]]
            dataset_dict = {}
            for i in range(0, len(databases_list), 2):
                filepath_A = databases_list[i]
                filepath_B = databases_list[i+1] 
                pair_file_A = os.path.basename(filepath_A)
                pair_file_B = os.path.basename(filepath_B)

                if pair_file_A and pair_file_B:
                    if person1 == "A":
                        pair_name = f"{pair_file_A}_&_{pair_file_B}"
                    elif person1 == "B":
                        pair_name = f"{pair_file_B}_&_{pair_file_A}"

                pair_dict = {}
                lstA_tier = get_tier_from_file(filepath_A, "Role")
                lstB_tier = get_tier_from_file(filepath_B, "Role")

                if tier in lstA:
                    lstA[tier].extend(lstA_tier[tier])
                else:
                    lstA[tier] = lstA_tier[tier]

                if tier in lstB:
                    lstB[tier].extend(lstB_tier[tier])
                else:
                    lstB[tier] = lstB_tier[tier]

                if person1 == "A":
                    overlapping_segments = get_overlapping_segments(lstA_tier[tier], lstB_tier[tier])
                elif person1 == "B":
                    overlapping_segments = get_overlapping_segments(lstB_tier[tier], lstA_tier[tier])
                pair_dict[tier] = {'Segments': overlapping_segments}

                dataset_dict[pair_name] = pair_dict

            overlapping_segments_dict[database] = dataset_dict
    dataframes = {}
    overlap_segments_set = set()
    overlap_count_dict = {}
    for database, dataset_dict in overlapping_segments_dict.items():
        overlap_count_list = []
        for pair_name, pair_dict in dataset_dict.items():
            overlap_count = 0
            overlap_count_spk_lsn = 0
            overlap_count_lsn_spk = 0
            segments = pair_dict[tier]["Segments"]
            if person1 == "A":
                for segmentA, segmentB in segments.items():
                    for segB in segmentB:
                        segment_key = f"{segB}"
                        if segment_key not in overlap_segments_set:
                            overlap_count += 1
                            overlap_segments_set.add(segment_key)
                            if (segmentA[2].replace(" ", "") == "spk" and segB[2].replace(" ", "") == "lsn"):
                                overlap_count_spk_lsn += 1
                            elif (segmentA[2].replace(" ", "") == "lsn" and segB[2].replace(" ", "") == "spk"):
                                overlap_count_lsn_spk += 1
                        break 
            elif person1 == "B":
                for segmentB, segmentA in segments.items():
                    for segA in segmentA:
                        segment_key = f"{segA}"
                        if segment_key not in overlap_segments_set:
                            overlap_count += 1
                            overlap_segments_set.add(segment_key)
                            if (segmentB[2].replace(" ", "") == "spk" and segA[2].replace(" ", "") == "lsn"):
                                overlap_count_spk_lsn += 1
                            elif (segmentB[2].replace(" ", "") == "lsn" and segA[2].replace(" ", "") == "spk"):
                                overlap_count_lsn_spk += 1
                        break
            overlap_count_list.append({
                'Database': database,
                'Pair': pair_name,
                f'Overlap Count for {tier} in general': overlap_count,
                f'Overlap Count for {person1} spk / {person2} lsn': overlap_count_spk_lsn,
                f'Overlap Count for {person1} lsn / {person2} spk': overlap_count_lsn_spk
            })
        
        df_overlap_count = pd.DataFrame(overlap_count_list)
        dataframes[database] = df_overlap_count
    return dataframes