def cells_cells_neigh (cells_coords_list, radius):


    ###
    #Function returns list of cell neighbors indexes

    #arguments:
    #- cells_coords_list: list of cells coords
    #- radius: radius of neighborhood outline

    #function returns:
    #'cell_neigh_index': array that contains indexes of neighboring cells
    #'len_cell_neigh_num': array that contains number of neighboring cells
    ###

    cells_coords_list = list(cells_coords_list)

    cells_cells_dist_bool = []
    
    for i in np.arange(len(cells_coords_list)):

        cells_cells_dist  = np.logical_and((np.sqrt(np.sum((np.array(cells_coords_list[i]) - np.array(cells_coords_list)) ** 2, axis=1))) > 0,
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
        
        cell_neigh_num .append(neigh_num)

    
    return ({'cell_neigh_index': cell_neigh_index, 
             'cell_neigh_num': cell_neigh_num})




def cell_cell_type_neigbor (cells_neigh_index, cell_type_list, cell_type):

    
    ###
    #Function returns list of cell neighbors indexes

    #arguments:
    #- cells_neigh_index: array that contains indexes of neighboring cells
    #- cell_type_list: list of cell types
    #- cell_type: cell type of interest

    #function returns:
    #'cell_cell_type': array that contain indexes of cells of certain type, 
    #'len_cell_cell_type': array that contains total number of neighboring cells of certain type
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




def cell_cell_neigh_features (cell_neigh_index, cell_feature_list):
    
    
    ###
    #Function returns average feature of interest of neighboring cells

    #arguments:
    #- cell_neigh_index: array with indexis of neighboring cells
    #- cell_feature_list: array with features of segmented cells

    #function returns:
    #- mean_result: average feature of interest of neighboring cells
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


def cell_cell_spatial (cells_coords_list, radius, cell_type = 'None', cell_type_list = 'None', cell_feature_list = 'None'):
    
    
    ###
    #Function returns dataframe which contains information about neighboring cells, their type, and features

    #arguments:
    #- cells_coords_list: list of cells coords
    #- radius: radius of neighborhood outline
    #- cell_type: default 'None', list of cell types for spatial analysis. 'All' make function perform calculation for all types of cell
    #- cell_type_list: list(column) with cell types
    #- cell_feature_list: list(colum) of feature of cells    

    #function returns:
    # - pd.DataFrame with calculatedd spatial information of neighboring cells
    ###
    
    
    
    neighbors = cells_cells_neigh (cells_coords_list, radius)
    
    
    dicts_neigh_result = {}
    
    
    dicts_feature_result = {}
    
    
    if cell_type != 'None':
        
         
        if cell_type == 'All':
            
        
            types = []
            
            
            type_index_list = []
            
            
            type_neighbors_list = []
            
        
            for i in set(cell_type_list):
                
                
                type_neigbors = cell_cell_type_neigbor (neighbors['cell_neigh_index'], cell_type_list, i)
                
                
                type_index_list.append(type_neigbors['cell_cell_type'])
                
                
                type_neighbors_list.append(type_neigbors['len_cell_cell_type'])
                
                
                types.append(i)
                
            
            dicts_typs = {}
            
            
            for i, j in zip(types, type_neighbors_list):
                
                
                dicts_typs[i] = j
                    
                    
            dicts_neighbors = {'cells_number' : neighbors['cell_neigh_num']}
            
            
            dicts_neigh_result = {**dicts_neighbors, **dicts_typs}
                
        
        else:
            
            
            types = []
            
            
            type_index_list = []
            
            
            type_neighbors_list = []
                        
                        
            for i in cell_type:
                
            
                type_neighbors = cell_cell_type_neigbor (neighbors['cell_neigh_index'], cell_type_list, i)


                type_index_list.append(type_neighbors['cell_cell_type'])
                
                
                type_neighbors_list.append(type_neighbors['len_cell_cell_type'])
                

                types.append(i)
                
                
            dicts_typs = {}
            
            
            for i, j in zip(types, type_neighbors_list):
                
                
                dicts_typs[i] = j
                
                
            dicts_neighbors = {'cell_neighbors' : neighbors['cell_neigh_num']}
            
            
            dicts_neigh_result = {**dicts_neighbors, **dicts_typs}
        
    
    
    if cell_type == 'None':
        
        
        dicts_neigh_result = {'cell_neighbors' : neighbors['cell_neigh_num']}
        
    
    
    
    if cell_feature_list != 'None':
        
        
        if cell_type != 'None':
        
         
            if cell_type == 'All': 
        
        
                feature_result = []


                var_names = []


                for i in np.arange (len(cell_feature_list)):


                    neighbors_feature = cell_cell_neigh_features (neighbors['cell_neigh_index'], cell_feature_list[i])


                    feature_result.append(neighbors_feature)


                    var_name = 'common_feature_' + '%s' % i


                    var_names.append(var_name)
                    
                
                dicts_feature_result = {}
                
                
                for i, j in zip(var_names, feature_result):

                    
                    dicts_feature_result[i] = j
                    
                    
            else:
                
                
                feature_result = []


                var_names = []
                

                for i in np.arange (len(cell_feature_list)):

                    
                    for j in np.arange (len(type_index_list)):
                    
                    
                        neighbors_feature = cell_cell_neigh_features (type_index_list[j], cell_feature_list[i])


                        feature_result.append(neighbors_feature)


                        var_name = '%s' %types[j] + '_feature_' + '%s' % i


                        var_names.append(var_name)


                dicts_feature_result = {}
                
                
                for i, j in zip(var_names, feature_result):

                    
                    dicts_feature_result[i] = j
                    
    
        if cell_type == 'None':


            feature_result = []


            var_names = []


            for i in np.arange (len(cell_feature_list)):


                neighbors_feature = cell_cell_neigh_features (neighbors['cell_neigh_index'], cell_feature_list[i])


                feature_result.append(neighbors_feature)


                var_name = 'common_feature_' + '%s' % i


                var_names.append(var_name)


            dicts_feature_result = {}


            for i, j in zip(var_names, feature_result):


                dicts_feature_result[i] = j
                
    
    if cell_feature_list == 'None':
        
        
        pass
                    
    
    return (pd.DataFrame({**dicts_neigh_result, **dicts_feature_result}))




def cell_fibs_neigh (executed_fibs, cells_coords_list, radius):


    ###
    #Function returns number of fibs neighbors

    #arguments:
    #- props_pruned: array with calculated features of labeled fibers
    #- cells_coords_list: list of cells coords
    #- radius: radius of neighborhood outline

    #function returns:
    #'fibs_neigh_index': array with indexes of neighboring fibers,
    #'fibs_neigh_num': array with number of neighboring fibers,
    ###


    fibers_coords = []

    for j in np.arange(len(executed_fibs['props_pruned'])):

        coords =  np.array(executed_fibs['props_pruned'][j].coords)
        
        fibers_coords.append(coords)


    cells_fibcoords_dist_bool = []
    
    for i in np.arange(len(cells_coords_list)):

        cells_fibcoords_dist = []

        cells_fibcoords_dist_bool.append(cells_fibcoords_dist)

        for j in np.arange(len(fibers_coords)):

            cells_fibcoords = np.sqrt(np.sum((np.array(cells_coords_list[i]) - np.array(fibers_coords[j])) ** 2, axis=1)) < radius

            cells_fibcoords_dist.append(cells_fibcoords)
 
        
    cell_fibs_neigh_index = []

    for k in np.arange(len(cells_fibcoords_dist_bool)):

        cell_fibs_neigh = []
        
        cell_fibs_neigh_index.append(cell_fibs_neigh)

        for j in np.arange(len(cells_fibcoords_dist_bool[k])):
            
            if any (cells_fibcoords_dist_bool[k][j]):

                    cell_fibs_neigh.append(j)


    cell_fibs_neigh_num = []

    for i in np.arange(len(cell_fibs_neigh_index)):

        cell_fibs_neigh = len(cell_fibs_neigh_index[i])

        cell_fibs_neigh_num.append(cell_fibs_neigh) 
        
  
    return ({'fibs_neigh_index': cell_fibs_neigh_index,
             'fibs_neigh_num': cell_fibs_neigh_num})



def cell_fibs_neigh_length (executed_fibs, cell_fibs_neigh_index):


    ###
    #Function calculates average length of neighboring fibers

    #arguments:
    #- props_pruned: array with calculated features of labeled fibers
    #- cell_fibs_neigh_index: array with indexes of neighboring fibers

    #function returns:
    #cell_fibs_neigh_length array with mean length of neighboring fibers
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



def cell_fibs_neigh_angle (executed_fibs, cell_fibs_neigh_index):


    ###
    #Function returns average angle of neighboring fibers


    #arguments:
    #- props_pruned: array with calculated features of labeled fibers
    #- cell_fibs_neigh_index: array with indexes of neighboring fibers

    #function returns:
    #cell_fibs_neigh_angle
    ###

    
    cell_fibs_neigh_angle = []

    for j in np.arange(len(cell_fibs_neigh_index)):

        cell_fibs_angle = []

        cell_fibs_neigh_angle.append(cell_fibs_angle)
        
        for i in cell_fibs_neigh_index[j]:

            cell_fibs = np.rad2deg(np.arctan2((executed_fibs['props_pruned'][i].coords[np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]][0] - executed_fibs['props_pruned'][i].coords[np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]][0]),
                                              (executed_fibs['props_pruned'][i].coords[np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]][1] - executed_fibs['props_pruned'][i].coords[np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]][1])))

            cell_fibs_angle.append(cell_fibs)
            
            
    mean_result = []
    
    for i in cell_fibs_neigh_angle:
        
        result = np.mean(np.array(i))
        
        mean_result.append(result)


    return (mean_result)



def cell_fibs_neigh_strightness (executed_fibs, cell_fibs_neigh_index):


    ###
    #Function returns average straightness of neighboring fibers

    #arguments:
    #- props_pruned: array with calculated features of labeled fibers
    #- cell_fibs_neigh_index: array with indexes of neighboring fibers

    #function returns:
    #cell_fibs_neigh_strightness
    ###


    cell_fibs_neigh_strightness = []

    for j in np.arange(len(cell_fibs_neigh_index)):

        cell_fibs_strightness = []

        cell_fibs_neigh_strightness.append(cell_fibs_strightness)

        for i in cell_fibs_neigh_index[j]:

            fibs_strightness = distance.euclidean(executed_fibs['props_pruned'][i].coords[np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[0]],
                                                   executed_fibs['props_pruned'][i].coords[np.lexsort(np.array(executed_fibs['props_pruned'][i].coords).T[:])[-1]]) / executed_fibs['props_pruned'][i].perimeter

            cell_fibs_strightness.append(fibs_strightness)

    
    
    mean_result = []
    
    for i in cell_fibs_neigh_strightness:
        
        result = np.mean(np.array(i))
        
        mean_result.append(result)


    return (mean_result)




def cell_fibs_neigh_thikness (executed_fibs, cell_fibs_neigh_index):


    ###
    #Function returns average thickness of neighboring fibers

    #arguments:
    #- props_pruned: array with calculated features of labeled fibers
    #- cell_fibs_neigh_index: array with indexes of neighboring fibers

    #function returns:
    #cell_fibs_neigh_thikness
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
        
            labels.append(thikness_matrix[executed_fibs['props_pruned'][j].coords[i][0]][executed_fibs['props_pruned'][j].coords[i][1]])
            
    
    
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




def cell_fibs_neigh_alignment(executed_fibs, cell_fibs_neigh_index):


    ###
    #Function returns average alignment of neighboring fibers

    #arguments:
    #- props_pruned: array with calculated features of labeled fibers
    #- cell_fibs_neigh_index: array with indexes of neighboring fibers

    #function returns:
    #cell_fibs_neigh_alignment
    ###


    
    alignment_res = []

    for i in np.arange(len(cell_fibs_neigh_index)):

        fibs_cos = []

        alignment_res.append(fibs_cos)

        for j in cell_fibs_neigh_index[i]:
            
            fibs_alignment_mean = []

            fibs_cos.append(fibs_alignment_mean)
            
            for k in cell_fibs_neigh_index[i]:
                
                if j != k:

                    alignment_list = np.cos(abs(np.rad2deg(np.arctan2(
                    executed_fibs['props_pruned'][j].coords[np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[-1]][0] -
                    executed_fibs['props_pruned'][j].coords[np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[0]][0],
                    executed_fibs['props_pruned'][j].coords[np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[-1]][1] -
                    executed_fibs['props_pruned'][j].coords[np.lexsort(np.array(executed_fibs['props_pruned'][j].coords).T[:])[0]][1]))
                    - np.rad2deg(np.arctan2(
                    executed_fibs['props_pruned'][k].coords[np.lexsort(np.array(executed_fibs['props_pruned'][k].coords).T[:])[-1]][0] -
                    executed_fibs['props_pruned'][k].coords[np.lexsort(np.array(executed_fibs['props_pruned'][k].coords).T[:])[0]][0],
                    executed_fibs['props_pruned'][k].coords[np.lexsort(np.array(executed_fibs['props_pruned'][k].coords).T[:])[-1]][1] -
                    executed_fibs['props_pruned'][k].coords[np.lexsort(np.array(executed_fibs['props_pruned'][k].coords).T[:])[0]][1]))))

                    fibs_alignment_mean.append(abs(alignment_list))
                    
                else:
                    
                    pass
                
    
    mean_fibs_alignment = []
    
    for i in alignment_res:
        
        mean_lin = np.mean(i) 
        
        mean_fibs_alignment.append(mean_lin)
        
        

    return (mean_fibs_alignment)


def cell_fibs_spatial (executed_fibs, cells_coords_list, radius):
    
    
    ###
    #Function returns dataframe which contains information about neighboring fibers and their features

    #arguments:
    #- executed_fibs: dictionary executed with fiber_executer function that contains information about segmented fibers
    #- cells_coords_list: list with cells coordinates (y,x)
    #- radius: radius of neighborhood outline    

    #function returns:
    # - pd.DataFrame with calculatedd spatial information of neighboring fibers
    ###
    
    
    cell_fibs = cell_fibs_neigh (executed_fibs, cells_coords_list, radius)
    
    
    cell_length = cell_fibs_neigh_length (executed_fibs, cell_fibs['fibs_neigh_index'])
    
    
    cell_fibs_angle = cell_fibs_neigh_angle (executed_fibs, cell_fibs['fibs_neigh_index'])
    
    
    cell_fibs_strightness = cell_fibs_neigh_strightness (executed_fibs, cell_fibs['fibs_neigh_index'])
    
    
    cell_fibs_thikness = cell_fibs_neigh_thikness (executed_fibs, cell_fibs['fibs_neigh_index'])
    
    
    cell_fibs_alignment = cell_fibs_neigh_alignment (executed_fibs, cell_fibs['fibs_neigh_index'])
    
    
    return (pd.DataFrame({'fibs_neigh_num' : cell_fibs['fibs_neigh_num'],
                          'length' : cell_length,
                          'angle' : cell_fibs_angle,
                          'strightness' : cell_fibs_strightness,
                          'thickness' : cell_fibs_thikness,
                          'alignment' : cell_fibs_alignment}))
