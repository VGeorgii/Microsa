import numpy as np
from skimage import io
from matplotlib import colors
from skimage.morphology import medial_axis
from skimage.measure import regionprops
from skimage.measure import label
from scipy import ndimage
from skimage.graph import route_through_array
from scipy.ndimage import binary_closing, binary_hit_or_miss
from skimage.filters import frangi
from scipy.spatial import distance
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import time
import concurrent.futures
import gc


def cells_cells_neigh(cells_coords_list, radius):
    ###
    # Function returns list of cell neighbors indexes

    # arguments:
    # - cells_coords_list: list of cells coords
    # - radius: radius of neighborhood outline

    # function returns:
    # 'cell_neigh_index': array that contains indexes of neighboring cells
    # 'len_cell_neigh_num': array that contains number of neighboring cells
    ###

    cells_cells_dist_bool = []

    for i in np.arange(len(cells_coords_list)):
        cells_cells_dist = np.logical_and(
            (np.sqrt(np.sum((np.array(cells_coords_list[i]) - np.array(cells_coords_list)) ** 2, axis=1))) > 0,
            (np.sqrt(np.sum((np.array(cells_coords_list[i]) - np.array(cells_coords_list)) ** 2, axis=1))) < radius)

        cells_cells_dist_bool.append(cells_cells_dist)

    cell_neigh_index = []

    for j in np.arange(len(cells_cells_dist_bool)):

        neigh_index = []

        cell_neigh_index.append(neigh_index)

        for i in np.arange(len(cells_cells_dist_bool[j])):

            if cells_cells_dist_bool[j][i] == True:
                neigh_index.append(i)

    cell_neigh_num = []

    for i in np.arange(len(cell_neigh_index)):
        neigh_num = len(cell_neigh_index[i])

        cell_neigh_num.append(neigh_num)

    return ({'cell_neigh_index': cell_neigh_index,
             'cell_neigh_num': cell_neigh_num})


def cell_cell_type_neigbor(cells_neigh_index, cell_type_list, cell_type):
    ###
    # Function returns list of cell neighbors indexes

    # arguments:
    # - cells_neigh_index: array that contains indexes of neighboring cells
    # - cell_type_list: list of cell types
    # - cell_type: cell type of interest

    # function returns:
    # 'cell_cell_type': array that contain indexes of cells of certain type,
    # 'len_cell_cell_type': array that contains total number of neighboring cells of certain type
    ###

    cell_cell_type = []

    for j in np.arange(len(cells_neigh_index)):

        c_type = []

        cell_cell_type.append(c_type)

        for i in cells_neigh_index[j]:

            if cell_type_list[i] == cell_type:
                c_type.append(i)

    len_cell_cell_type = []

    for i in np.arange(len(cell_cell_type)):
        len_cell_type = len(cell_cell_type[i])

        len_cell_cell_type.append(len_cell_type)

    return ({'cell_cell_type': cell_cell_type,
             'len_cell_cell_type': len_cell_cell_type})


def cell_cell_neigh_features(cell_neigh_index, cell_feature_list):
    ###
    # Function returns average feature of interest of neighboring cells

    # arguments:
    # - cell_neigh_index: array with indexis of neighboring cells
    # - cell_feature_list: array with features of segmented cells

    # function returns:
    # - mean_result: average feature of interest of neighboring cells
    ###

    feature_result = []

    for j in np.arange(len(cell_neigh_index)):

        feature_list = []

        feature_result.append(feature_list)

        for i in cell_neigh_index[j]:
            feature = cell_feature_list[i]

            feature_list.append(feature)

    mean_result = []

    for i in feature_result:
        result = np.mean(np.array(i))

        mean_result.append(result)

    return (mean_result)


def cell_cell_spatial(cells_coords_list, radius, cell_type='None', cell_type_list='None', cell_feature_list='None'):
    ###
    # Function returns dataframe which contains information about neighboring cells, their type, and features

    # arguments:
    # - cells_coords_list: list of cells coords
    # - radius: radius of neighborhood outline
    # - cell_type: default 'None', list of cell types for spatial analysis. 'All' make function perform calculation for all types of cell
    # - cell_type_list: list(column) with cell types
    # - cell_feature_list: list(colum) of feature of cells

    # function returns:
    # - pd.DataFrame with calculatedd spatial information of neighboring cells
    ###

    neighbors = cells_cells_neigh(cells_coords_list, radius)

    dicts_neigh_result = {}

    dicts_feature_result = {}

    if cell_type != 'None':

        if cell_type == 'All':

            types = []

            type_index_list = []

            type_neighbors_list = []

            for i in set(cell_type_list):
                type_neigbors = cell_cell_type_neigbor(neighbors['cell_neigh_index'], cell_type_list, i)

                type_index_list.append(type_neigbors['cell_cell_type'])

                type_neighbors_list.append(type_neigbors['len_cell_cell_type'])

                types.append(i)

            dicts_typs = {}

            for i, j in zip(types, type_neighbors_list):
                dicts_typs[i] = j

            dicts_neighbors = {'cells_number': neighbors['cell_neigh_num']}

            dicts_neigh_result = {**dicts_neighbors, **dicts_typs}


        else:

            types = []

            type_index_list = []

            type_neighbors_list = []

            for i in cell_type:
                type_neighbors = cell_cell_type_neigbor(neighbors['cell_neigh_index'], cell_type_list, i)

                type_index_list.append(type_neighbors['cell_cell_type'])

                type_neighbors_list.append(type_neighbors['len_cell_cell_type'])

                types.append(i)

            dicts_typs = {}

            for i, j in zip(types, type_neighbors_list):
                dicts_typs[i] = j

            dicts_neighbors = {'cell_neighbors': neighbors['cell_neigh_num']}

            dicts_neigh_result = {**dicts_neighbors, **dicts_typs}

    if cell_type == 'None':
        dicts_neigh_result = {'cell_neighbors': neighbors['cell_neigh_num']}

    if cell_feature_list != 'None':

        if cell_type != 'None':

            if cell_type == 'All':

                feature_result = []

                var_names = []

                for i in np.arange(len(cell_feature_list)):
                    neighbors_feature = cell_cell_neigh_features(neighbors['cell_neigh_index'], cell_feature_list[i])

                    feature_result.append(neighbors_feature)

                    var_name = 'common_feature_' + '%s' % i

                    var_names.append(var_name)

                dicts_feature_result = {}

                for i, j in zip(var_names, feature_result):
                    dicts_feature_result[i] = j


            else:

                feature_result = []

                var_names = []

                for i in np.arange(len(cell_feature_list)):

                    for j in np.arange(len(type_index_list)):
                        neighbors_feature = cell_cell_neigh_features(type_index_list[j], cell_feature_list[i])

                        feature_result.append(neighbors_feature)

                        var_name = '%s' % types[j] + '_feature_' + '%s' % i

                        var_names.append(var_name)

                dicts_feature_result = {}

                for i, j in zip(var_names, feature_result):
                    dicts_feature_result[i] = j

        if cell_type == 'None':

            feature_result = []

            var_names = []

            for i in np.arange(len(cell_feature_list)):
                neighbors_feature = cell_cell_neigh_features(neighbors['cell_neigh_index'], cell_feature_list[i])

                feature_result.append(neighbors_feature)

                var_name = 'common_feature_' + '%s' % i

                var_names.append(var_name)

            dicts_feature_result = {}

            for i, j in zip(var_names, feature_result):
                dicts_feature_result[i] = j

    if cell_feature_list == 'None':
        pass

    return (pd.DataFrame({**dicts_neigh_result, **dicts_feature_result}))


def cell_fibs_neigh(executed_fibs, cells_coords_list, radius):
    ###
    # Function returns number of fibs neighbors

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers
    # - cells_coords_list: list of cells coords
    # - radius: radius of neighborhood outline

    # function returns:
    # 'fibs_neigh_index': array with indexes of neighboring fibers,
    # 'fibs_neigh_num': array with number of neighboring fibers,
    ###

    fibers_coords = []

    for j in np.arange(len(executed_fibs['props_pruned'])):
        coords = np.array(executed_fibs['props_pruned'][j].coords)

        fibers_coords.append(coords)

    cells_fibcoords_dist_bool = []

    for i in np.arange(len(cells_coords_list)):

        cells_fibcoords_dist = []

        cells_fibcoords_dist_bool.append(cells_fibcoords_dist)

        for j in np.arange(len(fibers_coords)):
            cells_fibcoords = np.sqrt(
                np.sum((np.array(cells_coords_list[i]) - np.array(fibers_coords[j])) ** 2, axis=1)) < radius

            cells_fibcoords_dist.append(cells_fibcoords)

    cell_fibs_neigh_index = []

    for k in np.arange(len(cells_fibcoords_dist_bool)):

        cell_fibs_neigh = []

        cell_fibs_neigh_index.append(cell_fibs_neigh)

        for j in np.arange(len(cells_fibcoords_dist_bool[k])):

            if any(cells_fibcoords_dist_bool[k][j]):
                cell_fibs_neigh.append(j)

    cell_fibs_neigh_num = []

    for i in np.arange(len(cell_fibs_neigh_index)):
        cell_fibs_neigh = len(cell_fibs_neigh_index[i])

        cell_fibs_neigh_num.append(cell_fibs_neigh)

    return ({'fibs_neigh_index': cell_fibs_neigh_index,
             'fibs_neigh_num': cell_fibs_neigh_num})


def cell_fibs_neigh_length(executed_fibs, cell_fibs_neigh_index):
    ###
    # Function calculates average length of neighboring fibers

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers
    # - cell_fibs_neigh_index: array with indexes of neighboring fibers

    # function returns:
    # cell_fibs_neigh_length array with mean length of neighboring fibers
    ###

    fibs_length = []

    for i in np.arange(len(executed_fibs['props_pruned'])):
        cell_fibs = executed_fibs['props_pruned'][i].perimeter

        fibs_length.append(cell_fibs)

    length = []

    for j in np.arange(len(cell_fibs_neigh_index)):

        pre_length = []

        length.append(pre_length)

        if len(cell_fibs_neigh_index[j]) > 0:

            for i in cell_fibs_neigh_index[j]:
                pre_length.append(fibs_length[i])

    mean_result = []

    for i in length:
        result = np.mean(np.array(i))

        mean_result.append(result)

    return (mean_result)


def cell_fibs_neigh_angle(executed_fibs, cell_fibs_neigh_index):
    ###
    # Function returns average angle of neighboring fibers

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers
    # - cell_fibs_neigh_index: array with indexes of neighboring fibers

    # function returns:
    # cell_fibs_neigh_angle
    ###

    cell_fibs_neigh_angle = []

    for j in np.arange(len(cell_fibs_neigh_index)):

        cell_fibs_angle = []

        cell_fibs_neigh_angle.append(cell_fibs_angle)

        for i in cell_fibs_neigh_index[j]:
            cell_fibs = np.rad2deg(np.arctan2((executed_fibs['props_pruned'][i].coords[
                                                   np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[
                                                       -1]][0] - executed_fibs['props_pruned'][i].coords[
                                                   np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[
                                                       0]][0]),
                                              (executed_fibs['props_pruned'][i].coords[
                                                   np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[
                                                       -1]][1] - executed_fibs['props_pruned'][i].coords[
                                                   np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[
                                                       0]][1])))

            cell_fibs_angle.append(cell_fibs)

    mean_result = []

    for i in cell_fibs_neigh_angle:
        result = np.mean(np.array(i))

        mean_result.append(result)

    return (mean_result)


def cell_fibs_neigh_strightness(executed_fibs, cell_fibs_neigh_index):
    ###
    # Function returns average straightness of neighboring fibers

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers
    # - cell_fibs_neigh_index: array with indexes of neighboring fibers

    # function returns:
    # cell_fibs_neigh_strightness
    ###

    cell_fibs_neigh_strightness = []

    for j in np.arange(len(cell_fibs_neigh_index)):

        cell_fibs_strightness = []

        cell_fibs_neigh_strightness.append(cell_fibs_strightness)

        for i in cell_fibs_neigh_index[j]:
            fibs_strightness = distance.euclidean(executed_fibs['props_pruned'][i].coords[np.lexsort(
                np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]],
                                                  executed_fibs['props_pruned'][i].coords[np.lexsort(
                                                      np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]]) / \
                               executed_fibs['props_pruned'][i].perimeter

            cell_fibs_strightness.append(fibs_strightness)

    mean_result = []

    for i in cell_fibs_neigh_strightness:
        result = np.mean(np.array(i))

        mean_result.append(result)

    return (mean_result)


def cell_fibs_neigh_thikness(executed_fibs, cell_fibs_neigh_index):
    ###
    # Function returns average thickness of neighboring fibers

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers
    # - cell_fibs_neigh_index: array with indexes of neighboring fibers

    # function returns:
    # cell_fibs_neigh_thikness
    ###

    thikness_matrix = []

    for i in np.arange(len(executed_fibs['skeleton'])):
        thikness = executed_fibs['skeleton'][i] * executed_fibs['distance'][i]

        thikness_matrix.append(thikness)

    labels_thikness_row = []

    for j in np.arange(len(executed_fibs['props_pruned'])):

        labels = []

        labels_thikness_row.append(labels)

        for i in np.arange(len(executed_fibs['props_pruned'][j].coords)):
            labels.append(thikness_matrix[executed_fibs['props_pruned'][j].coords[i][0]][
                              executed_fibs['props_pruned'][j].coords[i][1]])

    labels_thikness = []

    for i in np.arange(len(executed_fibs['props_pruned'])):
        thikness = np.sum(labels_thikness_row[i]) / executed_fibs['props_pruned'][i].area

        labels_thikness.append(thikness)

    mean_thikness = []

    for j in np.arange(len(cell_fibs_neigh_index)):

        fibers = []

        mean_thikness.append(fibers)

        if len(cell_fibs_neigh_index[j]) > 0:

            for i in cell_fibs_neigh_index[j]:
                fibers.append(labels_thikness[i])

    mean_result = []

    for i in mean_thikness:
        result = np.mean(np.array(i))

        mean_result.append(result)

    return (mean_result)


def cell_fibs_neigh_linearity(executed_fibs, cell_fibs_neigh_index):
    ###
    # Function returns average linearity of neighboring fibers

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers
    # - cell_fibs_neigh_index: array with indexes of neighboring fibers

    # function returns:
    # cell_fibs_neigh_linearity
    ###

    fibs_cos = []

    for i in np.arange(len(cell_fibs_neigh_index)):

        fibs_linearity_mean = []

        fibs_cos.append(fibs_linearity_mean)

        for j in cell_fibs_neigh_index[i]:
            linearity_list = np.cos(abs(np.rad2deg(np.arctan2(
                executed_fibs['props_pruned'][i].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]][0] -
                executed_fibs['props_pruned'][i].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]][0],
                executed_fibs['props_pruned'][i].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]][1] -
                executed_fibs['props_pruned'][i].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]][1]))
                                        - np.rad2deg(np.arctan2(
                executed_fibs['props_pruned'][j].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[-1]][0] -
                executed_fibs['props_pruned'][j].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[0]][0],
                executed_fibs['props_pruned'][j].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[-1]][1] -
                executed_fibs['props_pruned'][j].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[0]][1]))))

            fibs_linearity_mean.append(linearity_list)

    mean_fibs_linearity = []

    for i in fibs_cos:
        mean_lin = np.mean(list(map(abs, i)))

        mean_fibs_linearity.append(mean_lin)

    return (mean_fibs_linearity)


def cell_fibs_spatial(executed_fibs, cells_coords_list, radius):
    ###
    # Function returns dataframe which contains information about neighboring fibers and their features

    # arguments:
    # - executed_fibs: dictionary executed with fiber_executer function that contains information about segmented fibers
    # - cells_coords_list: list with cells coordinates (y,x)
    # - radius: radius of neighborhood outline

    # function returns:
    # - pd.DataFrame with calculatedd spatial information of neighboring fibers
    ###

    cell_fibs = cell_fibs_neigh(executed_fibs, cells_coords_list, radius)

    cell_length = cell_fibs_neigh_length(executed_fibs, cell_fibs['fibs_neigh_index'])

    cell_fibs_angle = cell_fibs_neigh_angle(executed_fibs, cell_fibs['fibs_neigh_index'])

    cell_fibs_strightness = cell_fibs_neigh_strightness(executed_fibs, cell_fibs['fibs_neigh_index'])

    cell_fibs_thikness = cell_fibs_neigh_thikness(executed_fibs, cell_fibs['fibs_neigh_index'])

    cell_fibs_linearity = cell_fibs_neigh_linearity(executed_fibs, cell_fibs['fibs_neigh_index'])

    return (pd.DataFrame({'fibs_neigh_num': cell_fibs['fibs_neigh_num'],
                          'length': cell_length,
                          'angle': cell_fibs_angle,
                          'strightness': cell_fibs_strightness,
                          'thickness': cell_fibs_thikness,
                          'linearity': cell_fibs_linearity}))


def endPoints(skel):
    ##
    # Function returns boolen array that contains True values in the coordinates that fit as the endpoint of the fibers.

    # arguments:
    # - ‘skel’: array with skeletonized fibers

    # function returns:
    # - 'endp': boolen array that contains True values in the coordinates that fit as the endpoint of the fibers
    ##

    epoint1 = np.array([[0, 0, 0], [0, 1, 0], [0, 1, 0]])

    epoint2 = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1]])

    epoint3 = np.array([[0, 0, 0], [0, 1, 1], [0, 0, 0]])

    epoint4 = np.array([[0, 0, 1], [0, 1, 0], [0, 0, 0]])

    epoint5 = np.array([[0, 1, 0], [0, 1, 0], [0, 0, 0]])

    epoint6 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 0]])

    epoint7 = np.array([[0, 0, 0], [1, 1, 0], [0, 0, 0]])

    epoint8 = np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0]])

    epoint9 = np.array([[0, 0, 0], [0, 1, 0], [1, 1, 0]])

    epoint10 = np.array([[0, 0, 0], [0, 1, 0], [0, 1, 1]])

    epoint11 = np.array([[0, 0, 0], [0, 1, 1], [0, 0, 1]])

    epoint12 = np.array([[0, 0, 1], [0, 1, 1], [0, 0, 0]])

    epoint13 = np.array([[0, 1, 1], [0, 1, 0], [0, 0, 0]])

    epoint14 = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 0]])

    epoint15 = np.array([[1, 0, 0], [1, 1, 0], [0, 0, 0]])

    epoint16 = np.array([[0, 0, 0], [1, 1, 0], [1, 0, 0]])

    endp1 = binary_hit_or_miss(skel, structure1=epoint1)

    endp2 = binary_hit_or_miss(skel, structure1=epoint2)

    endp3 = binary_hit_or_miss(skel, structure1=epoint3)

    endp4 = binary_hit_or_miss(skel, structure1=epoint4)

    endp5 = binary_hit_or_miss(skel, structure1=epoint5)

    endp6 = binary_hit_or_miss(skel, structure1=epoint6)

    endp7 = binary_hit_or_miss(skel, structure1=epoint7)

    endp8 = binary_hit_or_miss(skel, structure1=epoint8)

    endp9 = binary_hit_or_miss(skel, structure1=epoint9)

    endp10 = binary_hit_or_miss(skel, structure1=epoint10)

    endp11 = binary_hit_or_miss(skel, structure1=epoint11)

    endp12 = binary_hit_or_miss(skel, structure1=epoint12)

    endp13 = binary_hit_or_miss(skel, structure1=epoint13)

    endp14 = binary_hit_or_miss(skel, structure1=epoint14)

    endp15 = binary_hit_or_miss(skel, structure1=epoint15)

    endp16 = binary_hit_or_miss(skel, structure1=epoint16)

    endp = endp1 + \
           endp2 + \
           endp3 + \
           endp4 + \
           endp5 + \
           endp6 + \
           endp7 + \
           endp8 + \
           endp9 + \
           endp10 + \
           endp11 + \
           endp12 + \
           endp13 + \
           endp14 + \
           endp15 + \
           endp16

    return endp


def endpoints_storage(skel_labels, nlabels):
    ###
    # Function returns array that include cooridnates of fiber endpoints.

    # arguments:
    # - ‘skel_labels’: array that include labels of skeletonized fibers
    # - ‘nlabels’: total number of labels

    # function returns:
    # - 'endp': array that contains coordinates of each endpoint of the fiber
    ###

    endpoints_storage = np.where(endPoints(np.array(skel_labels)) == True)

    endpoints_list = []

    for i in np.arange(nlabels + 1):
        endpoints_list.append([])

    for i in np.arange(len(endpoints_storage[0])):
        endpoints_list[np.array(skel_labels[endpoints_storage[0][i]][endpoints_storage[1][i]])].append(
            [endpoints_storage[0][i], endpoints_storage[1][i]])

    return (endpoints_list)


def array_transform(array):
    ###
    # Function transformed array where values that are equal to 0 transformed to big numbers to performe graph algorithm.

    # arguments:
    # - ‘array’: array that include labels of skeletonized fibers

    # function returns:
    # - 'endp': transformed array where values that are equal to 0 transformed to big numbers
    ###

    array[array == 0] = 2 * 400000

    return (array)


def path_list(array, endpoints_list, nlabels):
    ###
    # Function returns array that include labels of pruned (without additional branches) fibers

    # arguments:
    # - ‘array’: array that include labels of skeletonized fibers
    # - ‘endpoints_list’: array that contains coordinates of each endpoint of the fiber
    # - ‘nlabels’: total number of labels

    # function returns:
    # - 'endp': transformed array where values that are equal to 0 transformed to big numbers
    ###

    path_list = []

    for i in np.arange(nlabels):

        if len(endpoints_list[i]) > 0:
            path = route_through_array(np.array(array),
                                       endpoints_list[i][np.lexsort(np.array(endpoints_list[i]).T[:])[0]],
                                       endpoints_list[i][np.lexsort(np.array(endpoints_list[i]).T[:])[-1]],
                                       fully_connected=True, geometric=True)[0]

            path_list.append(path)

    pruned_array = np.zeros((len(array), len(array[0])))

    path_list = list(map(np.array, path_list))

    for i in np.arange(len(path_list)):
        pruned_array[path_list[i][:, 0], path_list[i][:, 1]] = i + 1

    return pruned_array


def fibers_executor(gray):
    ##
    # Function mainly performs segmentation of fibers and extract parameters that are required for other functions like number of labels, distance between centroids etc.

    # arguments:
    # image – gray scale image to process

    # function returns:
    # - 'skeleton': image with skeletonized objects,
    # - 'distance': distance from skeleton to the edge of original object (radius),
    # - 'skel_labels_pruned': pruned and labeled skeletonized objects,
    # - 'props_pruned': properties of labeled objects,
    # - 'nlabels_pruned': number of labeled objects
    ##

    skeleton, distance = medial_axis(gray, return_distance=True)

    skel_labels = label(skeleton, neighbors=8)

    props = regionprops(skel_labels)

    nlabels = len(props)

    endpoints = endpoints_storage(skel_labels, nlabels)

    transformed = array_transform(skel_labels)

    pruning = path_list(transformed, endpoints, nlabels)

    skel_labels_pruned = label(pruning, neighbors=8)

    props_pruned = regionprops(skel_labels_pruned)

    return {'skeleton': skeleton,
            'distance': distance,
            'skel_labels_pruned': skel_labels_pruned,
            'props_pruned': props_pruned}


def fibs_cell_neigh(cells_coords_list, executed_fibs, radius):
    ###
    # Function returns array that include fiber and indexes of cells’ (or other objects). The distance is calculated between cells centroid and closest pixel of fibers skeleton.

    # arguments:
    # - ‘cells_coords’: list of cells coords
    # - ‘pruned_skel’: pruned and labeled skeletonized objects
    # - ‘radius’: radius of neighborhood outline

    # function returns:
    # - 'fibs_cell_neigh': array that contains indexes of neighboring cells,
    # - 'cells_neigh_num': array that contains total number of neighboring cells
    ###

    fibers_coords = []

    for i in np.arange(len(executed_fibs['props_pruned'])):
        coords = np.array(executed_fibs['props_pruned'][i].coords)

        fibers_coords.append(coords)

    cells_fibcoords_dist_bool = []

    for i in np.arange(len(cells_coords_list)):

        cells_fibcoords_dist_list = []

        cells_fibcoords_dist_bool.append(cells_fibcoords_dist_list)

        for j in np.arange(len(fibers_coords)):
            cells_fibcoords_dist = np.sqrt(
                np.sum((np.array(cells_coords_list[i]) - np.array(fibers_coords[j])) ** 2, axis=1)) < 25

            cells_fibcoords_dist_list.append(list(cells_fibcoords_dist))

    fibs_neigh_index = []

    for k in np.arange(len(cells_fibcoords_dist_bool)):

        fibs_neigh = []

        fibs_neigh_index.append(fibs_neigh)

        for j in np.arange(len(cells_fibcoords_dist_bool[k])):

            if any(cells_fibcoords_dist_bool[k][j]):
                fibs = j

                fibs_neigh.append(fibs)

    fibs_cell_neigh = []

    for j in np.arange(len(fibers_coords)):

        fibs_cell = []

        fibs_cell_neigh.append(fibs_cell)

        for i in np.arange(len(fibs_neigh_index)):

            if j in fibs_neigh_index[i]:
                cell = i

                fibs_cell.append(cell)

    cells_neigh_num = []

    for i in np.arange(len(fibs_cell_neigh)):
        cells_neigh = len(fibs_cell_neigh[i])

        cells_neigh_num.append(cells_neigh)

    return ({'fibs_cell_neigh': fibs_cell_neigh,
             'cells_neigh_num': cells_neigh_num})


def fibs_cell_type_neighbor(fibs_cell_neigh, cell_type_list, cell_type):
    ###
    # Function use the variable that was returned by fibs_cell_neigh function to classify neighbor objects (for example different type of cells)

    # arguments:
    # - fibs_neigh_index: array that contains indexes of neighboring cells
    # - cell_type_list: list of cell types
    # - cell_type: cell type of interest

    # function returns:
    # - 'fibs_cell_type': array that contain indexes of cells of certain type,
    # - 'len_fibs_cell_type': array that contains total number of neighboring cells of certain type
    ###

    fibs_cell_type = []

    for j in np.arange(len(fibs_cell_neigh)):

        c_cell_type = []

        fibs_cell_type.append(c_cell_type)

        for i in fibs_cell_neigh[j]:

            if cell_type_list[i] == cell_type:
                c_type = i

                c_cell_type.append(c_type)

    cell_type_number = []

    for i in fibs_cell_type:
        result = len(i)

        cell_type_number.append(result)

    return ({'fibs_cell_type': fibs_cell_type,
             'cell_type_number': cell_type_number})


def fibs_spatial(cells_coords_list, executed_fibs, radius, cell_type='None', cell_type_list='None'):
    ###
    # Function returns dataframe which contains information about neighboring cells, their type, and features

    # arguments:
    # - cells_coords_list: list of cells coords
    # - executed_fibs: array with calculated features of labeled fibers
    # - radius: radius of neighborhood outline
    # - cell_type: default 'None', list of cell types for spatial analysis. 'All' make function perform calculation for all types of cell
    # - cell_type_list: default 'None', list(column) with cell types

    # function returns:
    # - pd.DataFrame with calculatedd spatial information of neighboring cells
    ###

    neighbors = fibs_cell_neigh(cells_coords_list, executed_fibs, radius)

    res_dict = {}

    if cell_type != 'None':

        if cell_type == 'All':

            types = []

            type_neighbors_list = []

            for i in set(cell_type_list):
                type_neighbors = fibs_cell_type_neighbor(neighbors['fibs_cell_neigh'], cell_type_list, i)

                type_neighbors_list.append(type_neighbors['cell_type_number'])

                types.append(i)

            dicts_typs = {}

            for i, j in zip(types, type_neighbors_list):
                dicts_typs[i] = j

            dicts_neighbors = {'cell_neighbors': neighbors['cells_neigh_num']}

            res_dict = {**dicts_neighbors, **dicts_typs}


        else:

            types = []

            type_neighbors_list = []

            for i in cell_type:
                type_neighbors = fibs_cell_type_neighbor(neighbors['fibs_cell_neigh'], cell_type_list, i)

                type_neighbors_list.append(type_neighbors['cell_type_number'])

                types.append(i)

            dicts_typs = {}

            for i, j in zip(types, type_neighbors_list):
                dicts_typs[i] = j

            dicts_neighbors = {'cell_neighbors': neighbors['cells_neigh_num']}

            res_dict = {**dicts_neighbors, **dicts_typs}


    else:

        res_dict = {'cell_neighbors': neighbors['cells_neigh_num']}

    return (pd.DataFrame(res_dict))


def fibs_total_number(executed_fibs):
    number = len(executed_fibs['props_pruned'])

    return (number)


def fibs_fibs_neigh(executed_fibs, radius):
    ###
    # Function calculates number of other fibers neighbor in the given radius. The distance is calculated between fibers centroids.

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers
    # - radius: radius of neighborhood outline

    # function returns:
    # 'fibs_fibs_dist_index': indexes of neighboring fibers
    # 'number_of_fibs_neighbors': number of neighboring fibers
    ###

    fibers_coords = []

    for i in np.arange(len(executed_fibs['props_pruned'])):
        coords = np.array(executed_fibs['props_pruned'][i].centroid)

        fibers_coords.append(coords)

    fibs_fibs_dist_bool = []

    for i in np.arange(len(fibers_coords)):
        fibs_fibs_dist = np.logical_and(
            (np.sqrt(np.sum((np.array(fibers_coords[i]) - np.array(fibers_coords)) ** 2, axis=1))) > 0,
            (np.sqrt(np.sum((np.array(fibers_coords[i]) - np.array(fibers_coords)) ** 2, axis=1))) < radius)

        fibs_fibs_dist_bool.append(fibs_fibs_dist)

    fibs_fibs_dist_index = [[i for i in np.arange(len(fibs_fibs_dist_bool[j])) if fibs_fibs_dist_bool[j][i] == True]
                            for j in np.arange(len(fibs_fibs_dist_bool))]

    fibs_neigh_index = []

    for j in np.arange(len(fibs_fibs_dist_bool)):

        neigh_index = []

        fibs_neigh_index.append(neigh_index)

        for i in np.arange(len(fibs_fibs_dist_bool[j])):

            if fibs_fibs_dist_bool[j][i] == True:
                neigh_index.append(i)

    len_fibs_neigh_num = []

    for i in np.arange(len(fibs_neigh_index)):
        len_fibs_neigh = len(fibs_neigh_index[i])

        len_fibs_neigh_num.append(len_fibs_neigh)

    return ({'fibs_neigh_index': fibs_neigh_index,
             'len_fibs_neigh_num': len_fibs_neigh_num})


def fibs_length(executed_fibs):
    ###
    # Function returns list of fibers length

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers

    # function returns:
    # fibs_length mean length of labeled fibers
    ###

    fibs_length = []

    for i in np.arange(len(executed_fibs['props_pruned'])):
        length = executed_fibs['props_pruned'][i].perimeter

        fibs_length.append(length)

    return (fibs_length)


def fibs_angle(executed_fibs):
    ###
    # Function returns list of fibers angle to against X axis

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers

    # function returns:
    # fibs_angle mean angle of labeled fibers against X axis
    ###

    fibs_angle = []

    for fibs in np.arange(len(executed_fibs['props_pruned'])):
        angle = np.rad2deg(np.arctan2(
            (executed_fibs['props_pruned'][fibs].coords[-1][0] - executed_fibs['props_pruned'][fibs].coords[0][0]),
            (executed_fibs['props_pruned'][fibs].coords[-1][1] - executed_fibs['props_pruned'][fibs].coords[0][1])))

        fibs_angle.append(angle)

    return (fibs_angle)


def fibs_strightness(executed_fibs):
    ###
    # Functions performs calculation of straightness of fibers

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers

    # function returns:
    # fibs_straightness mean straightness of labeled fibers
    ###

    fibs_strightness = []

    for i in np.arange(len(executed_fibs['props_pruned'])):
        strightness = distance.euclidean(executed_fibs['props_pruned'][i].coords[
                                             np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]],
                                         executed_fibs['props_pruned'][i].coords[
                                             np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]]) / \
                      executed_fibs['props_pruned'][i].perimeter

        fibs_strightness.append(strightness)

    return (fibs_strightness)


def fibs_thikness(executed_fibs):
    ###
    # Function returns list of fibers thickness

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers

    # function returns:
    # fibs_thikness mean thickness of labeled fibers
    ###

    thikness_matrix = []

    for i in np.arange(len(executed_fibs['skeleton'])):
        thikness = np.array(executed_fibs['skeleton'][i] * executed_fibs['distance'][i])

        thikness_matrix.append(thikness)

    labels_thikness = [np.sum(
        [thikness_matrix[executed_fibs['props_pruned'][j].coords[i][0]][executed_fibs['props_pruned'][j].coords[i][1]]
         for i in np.arange(len(executed_fibs['props_pruned'][j].coords))]) / executed_fibs['props_pruned'][j].area
                       for j in np.arange(len(executed_fibs['props_pruned']))]

    labels_thikness_row = []

    for j in np.arange(len(executed_fibs['props_pruned'])):

        labels = []

        labels_thikness_row.append(labels)

        for i in np.arange(len(executed_fibs['props_pruned'][j].coords)):
            labels.append(thikness_matrix[executed_fibs['props_pruned'][j].coords[i][0]][
                              executed_fibs['props_pruned'][j].coords[i][1]])

    labels_thikness = []

    for i in np.arange(len(executed_fibs['props_pruned'])):
        thikness = np.sum(labels_thikness_row[i]) / executed_fibs['props_pruned'][i].area

        labels_thikness.append(thikness)

    return (labels_thikness)


def fibs_linearity(executed_fibs, fibs_fibs_dist_index):
    ###
    # Function returns linearity between fiber of interest and other closest fiber neighbors.

    # arguments:
    # - props_pruned: array with calculated features of labeled fibers

    # function returns:
    # fibs_linearity mean linearity of labeled fibers
    ###

    fibs_cos = []

    for i in np.arange(len(fibs_fibs_dist_index)):

        fibs_linearity_mean = []

        fibs_cos.append(fibs_linearity_mean)

        for j in fibs_fibs_dist_index[i]:
            linearity_list = np.cos(abs(np.rad2deg(np.arctan2(
                executed_fibs['props_pruned'][i].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]][0] -
                executed_fibs['props_pruned'][i].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]][0],
                executed_fibs['props_pruned'][i].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]][1] -
                executed_fibs['props_pruned'][i].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]][1]))
                                        - np.rad2deg(np.arctan2(
                executed_fibs['props_pruned'][j].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[-1]][0] -
                executed_fibs['props_pruned'][j].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[0]][0],
                executed_fibs['props_pruned'][j].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[-1]][1] -
                executed_fibs['props_pruned'][j].coords[
                    np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[0]][1]))))

            fibs_linearity_mean.append(linearity_list)

    mean_fibs_linearity = []

    for i in fibs_cos:
        mean_lin = np.mean(list(map(abs, i)))

        mean_fibs_linearity.append(mean_lin)

    return (mean_fibs_linearity)


def fibs_geom(executed_fibs, radius):
    ###
    # Function returns dataframe which contains information about fibers

    # arguments:
    # - executed_fibs: array with calculated features of labeled fibers
    # - radius: radius of neighborhood outline

    # function returns:
    # - pd.DataFrame with calculated fibers features (number, length, angle, strightness, thickness, linearity)
    ###

    fibs_fibs = fibs_fibs_neigh(executed_fibs, radius)

    number = fibs_total_number(executed_fibs)

    length = fibs_length(executed_fibs)

    angle = fibs_angle(executed_fibs)

    strightness = fibs_strightness(executed_fibs)

    thickness = fibs_thikness(executed_fibs)

    linearity = fibs_linearity(executed_fibs, fibs_fibs['fibs_neigh_index'])

    return (pd.DataFrame({'number': number,
                          'length': length,
                          'angle': angle,
                          'strightness': strightness,
                          'thickness': thickness,
                          'linearity': linearity}))


def cell_map(list_of_coords, list_of_types, gray, dilation):
    ###
    # Function returns array with marked cells and their type

    # arguments:
    # - cells_coords_list: list with cells coordinates (y,x)
    # - list_of_types: list with cells types
    # - gray: original image in gray scale
    # - dilation: size of dot which indicates cell localization

    # function returns:
    # - array with marked cells and their type
    ###

    array = np.zeros((len(gray), len(gray[0])))

    score = 0

    indexis = []

    for i in np.arange(len(list_of_types)):
        index = []

        indexis.append(index)

    for i in set(list_of_types):

        for j in np.arange(len(list_of_types)):

            if list_of_types[j] == i:
                indexis[j] = score

        score += 1

    list_of_coords_types = np.concatenate((list(list_of_coords), [[i] for i in indexis]), axis=1)

    for i in list_of_coords_types[:, 2]:
        array[list_of_coords_types[list_of_coords_types[:, 2] == i][:, 0], list_of_coords_types[
                                                                               list_of_coords_types[:, 2] == i][:,
                                                                           1]] = i + 2

    bigger_points = ndimage.grey_dilation(array, size=(5, 5), structure=np.ones((dilation, dilation)))

    return (bigger_points)


def fibs_map(array):
    ###
    # Function returns array that contains fibers skeleton

    # arguments:
    # - array: array with labeled fibers

    # function returns:
    # - array that contains fibers skeleton
    ###

    fib_map = np.multiply(array > 0, 1)

    return (fib_map)


def pict(array, size, colorlist=['black', 'red', 'green', 'yellow', 'blue', 'purple']):
    ###
    # Function visualize combined cellular and fibrilar maps

    # arguments:
    # - array: array to be visualized
    # - size: desired size of image in pixels
    # - colorlist: list of colors that will be used to indicate different types of cells

    # function returns:
    # - visualized array with marked cells and fibers
    ###

    color = matplotlib.colors.ListedColormap(colorlist)

    fig, ax = plt.subplots(figsize=(size, size))

    ax = plt.imshow(array, cmap=color)

    plt.axis('off')

    return (plt.show())


def visualization(image, list_of_coords, cell_type_list, executed_fibs, dilation, size=100):
    ###
    # Function visualize combined cellular and fibrillar maps

    # arguments:
    # - image: original image in gray scale
    # - list_of_coords: list with cells coordinates (y,x)
    # - cell_type_list: list with cells types
    # - executed_fibs: dictionary executed with fiber_executer function that contains information about segmented fibers
    # - dilation: size of dot which indicates cell localization
    # - size: default 100, desired size of image in pixels

    # function returns:
    # - visualized array with marked cells and fibers
    ###

    cell_m = cell_map(list_of_coords, cell_type_list, image, dilation)

    fib_mape = fibs_map(executed_fibs['skel_labels_pruned'])

    visual_combin = fib_mape + cell_m

    visual = pict(visual_combin, size)

    return (visual)

import glob


def radius_pass(x):
    return fibs_geom(x, 25)


if __name__ == '__main__':

    path_5 = 'C:/Users/microsa_update/Pict/Test/5'
    path_10 = 'C:/Users/microsa_update/Pict/Test/10'
    path_15 = 'C:/Users/microsa_update/Pict/Test/15'
    path_20 = 'C:/Users/microsa_update/Pict/Test/20'
    path_25 = 'C:/Users/microsa_update/Pict/Test/25'
    path_30 = 'C:/Users/microsa_update/Pict/Test/30'
    path_35 = 'C:/Users/microsa_update/Pict/Test/35'
    path_40 = 'C:/Users/microsa_update/Pict/Test/40'
    path_45 = 'C:/Users/microsa_update/Pict/Test/45'
    path_50 = 'C:/Users/microsa_update/Pict/Test/50'

    files_list_5 = glob.glob(path_5 + '/*.jpg')
    files_list_10 = glob.glob(path_10 + '/*.jpg')
    files_list_15 = glob.glob(path_15 + '/*.jpg')
    files_list_20 = glob.glob(path_20 + '/*.jpg')
    files_list_25 = glob.glob(path_25 + '/*.jpg')
    files_list_30 = glob.glob(path_30 + '/*.jpg')
    files_list_35 = glob.glob(path_35 + '/*.jpg')
    files_list_40 = glob.glob(path_40 + '/*.jpg')
    files_list_45 = glob.glob(path_45 + '/*.jpg')
    files_list_50 = glob.glob(path_50 + '/*.jpg')


    gray_frangi_bit_10 = []

    for i in files_list_10:
        img = io.imread(i)
        gray_frangi = np.array(frangi(img, sigmas=range(4, 6, 10), gamma=25, black_ridges=False))
        sub_gray_frangi_bit = (gray_frangi / 255) > 0.000000001
        gray_frangi_bit_10.append(sub_gray_frangi_bit)

    print('Calculation is in process')



    start = time.perf_counter()


    r_array = gray_frangi_bit_10

    chunk = 0

    step = 13

    fibere_exe_res = []

    while chunk < len(r_array):

        if chunk + step <= len(r_array):

            print(chunk)

            sub_res = []

            for i in np.arange(chunk, chunk + step):

                sub_res.append(r_array[i])

            with concurrent.futures.ProcessPoolExecutor() as executor:

                fibere_exe = list(executor.map(fibers_executor, sub_res))

                fibere_geom = list(executor.map(radius_pass, fibere_exe))

            for i in np.arange(chunk + step):

                fibere_geom[i].to_csv('C:/Users/microsa_update/Time test/' + files_list_10[i][files_list_10[i].index('10') + 3 :files_list_10[i].index('.jpg')] + '.csv')

            chunk += step
            del (fibere_exe)
            del (fibere_geom)
            del (sub_res)
            gc.collect()


        if chunk + step > len(r_array):

            sub_res = []

            print(chunk)

            for i in np.arange(chunk, len(r_array)):

                sub_res.append(r_array[i])

            with concurrent.futures.ProcessPoolExecutor() as executor:

                fibere_exe = list(executor.map(fibers_executor, sub_res))

                fibere_geom = list(executor.map(radius_pass, fibere_exe))

            for i in np.arange(len(r_array) - chunk):

                fibere_geom[i].to_csv('C:/Users/microsa_update/Time test/' + files_list_10[i][files_list_10[i].index('10') + 3 :files_list_10[i].index('.jpg')] + '.csv')

            chunk += len(r_array) - chunk
            del (fibere_exe)
            del (fibere_geom)
            del (sub_res)
            gc.collect()


    finish = time.perf_counter()
    print('Finished in {} second(s)'.format(round(finish - start, 2)))