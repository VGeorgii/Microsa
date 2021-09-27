def endPoints (skel):


    ##
    #Function returns boolen array that contains True values in the coordinates that fit as the endpoint of the fibers. 

    #arguments:
    #- ‘skel’: array with skeletonized fibers

    #function returns:
    #- 'endp': boolen array that contains True values in the coordinates that fit as the endpoint of the fibers
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


    endp =  endp1 + \
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



def endpoints_storage (skel_labels, nlabels):


    ###
    #Function returns array that include cooridnates of fiber endpoints. 

    #arguments:
    #- ‘skel_labels’: array that include labels of skeletonized fibers
    #- ‘nlabels’: total number of labels

    #function returns:
    #- 'endp': array that contains coordinates of each endpoint of the fiber
    ###


    endpoints_storage = np.where(endPoints(np.array(skel_labels)) == True)


    endpoints_list = []

    for i in np.arange(nlabels + 1):

        endpoints_list.append([])


    for i in np.arange(len(endpoints_storage[0])):

        endpoints_list[np.array(skel_labels[endpoints_storage[0][i]][endpoints_storage[1][i]])].append([endpoints_storage[0][i], endpoints_storage[1][i]])

    
    return (endpoints_list)



def array_transform (array):


    ###
    #Function transformed array where values that are equal to 0 transformed to big numbers to performe graph algorithm. 

    #arguments:
    #- ‘array’: array that include labels of skeletonized fibers

    #function returns:
    #- 'endp': transformed array where values that are equal to 0 transformed to big numbers
    ###

    
    array[array == 0] = 2*400000


    return (array)




def path_list (array, endpoints_list, nlabels):


    ###
    #Function returns array that include labels of pruned (without additional branches) fibers

    #arguments:
    #- ‘array’: array that include labels of skeletonized fibers
    #- ‘endpoints_list’: array that contains coordinates of each endpoint of the fiber
    #- ‘nlabels’: total number of labels

    #function returns:
    #- 'endp': transformed array where values that are equal to 0 transformed to big numbers
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
    
        pruned_array[path_list[i][:,0],path_list[i][:,1]] = i+1


    return pruned_array



def fibs_filter (segmented_fibs):

    filtered_fibs = []

    for i in np.arange(len(segmented_fibs)):

        if segmented_fibs[i].area > 3:

            filtered_fibs.append(segmented_fibs[i])
            
            
    return filtered_fibs



def fibers_executor (gray):

    ##
    #Function mainly performs segmentation of fibers and extract parameters that are required for other functions like number of labels, distance between centroids etc.

    #arguments:
    #image – gray scale image to process

    #function returns:
    #- 'skeleton': image with skeletonized objects, 
    #- 'distance': distance from skeleton to the edge of original object (radius), 
    #- 'skel_labels_pruned': pruned and labeled skeletonized objects, 
    #- 'props_pruned': properties of labeled objects, 
    #- 'nlabels_pruned': number of labeled objects
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
    
    props_pruned = fibs_filter(props_pruned)
    

    return {'skeleton': skeleton, 
            'distance': distance, 
            'skel_labels_pruned': skel_labels_pruned, 
            'props_pruned': props_pruned}
