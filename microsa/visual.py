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

    list_of_coords_types = np.concatenate((list_of_coords, [[i] for i in indexis]), axis=1)

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


def pict(array, size, colorlist=['black', 'red', 'green', 'yellow', 'orange', 'purple']):
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