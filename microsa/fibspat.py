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
        length = executed_fibs['props_pruned'][i].area

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
        strightness = dist_cheb.chebyshev(executed_fibs['props_pruned'][i].coords[
                                              np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]],
                                          executed_fibs['props_pruned'][i].coords[
                                              np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]]) / \
                      executed_fibs['props_pruned'][i].area

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


def executed_fibs(executed_fibs, radius):
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