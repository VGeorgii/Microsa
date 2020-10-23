# Mesa

Project was developed by G.Vasiukov and S.Novitskiy

**Mesa** (**M**icro**e**nvironment **s**patial **a**nalysis) is a package of useful functions for analysis of fibrous components of tissue microenvironment and combine that analysis in spatial dependent manner. The purpose of that library is to provide simplified tool  that gives an opportunity to combine the analysis of cellular and non-cellular components of tissue microenvironment.

**Utilization:**
-	Segmentation of fibril-like objects (tensor method) in 2D;
-	Calculation of fibers’ centroid, length, width, angle and linearity in 2D;
-	Estimation of fibers' features and providing a detailed report in pandas.DataFrame format;
-	Implementation of spatial analysis between cellular and non-cellular objects.

**Installation:**
- Using pip:<br />
    $ pip install mesa<br />
- Using Git repository:<br />
    $ git clone https://github.com/VGeorgii/Mesa.git<br />
    $ cd Mesa<br />
    $ python setup.py install<br />
    
**Requirements:**
-	Python 3.7;
-	Numpy >= 1.13;
-	SciPy >= 0.19;
-	scikit-image >= 0.12;
-	scikit-learn >= 0.18;
-	Pandas >= 0.19;
-	Matplotlib >= 2.0.

**Usage:**

Open a grayscale image, perform filtering, binarization and segmentation of fibers, calculate geometrical and spatial features of fibers, perform spatial analysis for other objects (cells).

**Fibers:**
  
def fibers_executor (binary_image):<br />
    return {'skeleton': skeleton, 'distance': distance, 'skel_labels_pruned': skel_labels_pruned, 'props_pruned': props_pruned, 'nlabels_pruned': nlabels_pruned}<br />
Function mainly performs segmentation of fibers and extract parameters that are required for other functions like number of labels, distance between centroids etc.

def fibs_cell_neigh (cells_coords_list, props_pruned, radius):<br />
    return ({'fibs_neigh_index': fibs_neigh_index, 'fibs_cell_neigh': fibs_cell_neigh})<br />
Function returns nested list that include fiber and indexes of cells’ (or other objects) neighbors. The distance is calculated between cells centroid and closest pixel of fibers skeleton.

def fibs_cell_type_neighbor(fibs_neigh_index, dataframe, cell_type):<br />
    return (cell_type_neigh)<br />
Function use the variable that was returned by fibs_cell_neigh function to classify neighbor objects (for example different type of cells) 

def fibs_total_number(nlabels_pruned):<br />
    return (len(nlabels_pruned))<br />
Function returns total number of fibers that were segmented by fibers_executor 

def fibs_fibs_neigh_num(nested_closest_neigbor_list):<br />
    return (number_of_fibs_neighbors)<br />
Function calculates number of other fibers neighbor in the given radius. The distance is calculated between fibers centroids.

def fibs_length (props_pruned):<br />
    return (fibs_length)<br />
Function returns list of fibers length  

def fibs_angle (props_pruned):<br />
    return (fibs_angle)<br />
Function returns list of fibers angle to against X axis  

def fibs_strightness(props_pruned):<br />
    return (fibs_strightness)<br />
Functions performs calculation of straightness of fibers. Str = line between ends of fibers/length of fibers 

def fibs_thikness(props_pruned, skeleton, dist):<br />
    return (labels_thikness)<br />
Function returns list of fibers thickness  

def fibs_linearity(fibers_coords, nested_closest_neigbor_list):<br />
    return (linearity)<br />
Function returns linearity between fiber of interest and other closest fiber neighbors. 

**Cells:**

def centroid_list_transform(dataframe, id, column_x, column_y):<br />
    return (coords_list)<br />
Function helps to put XY coordinates of objects as list of tuples.

def cells_cells_neigh(cells_coords_list, radius):<br />
    return ({'cell_neigh_index': cell_neigh_index, 'cell_cell_neigh': cell_cell_neigh})<br />
Function returns list of cell neighbors indexes

def cell_fibs_total_num(cell_cell_neigh):<br />
    return ([len(cell_cell_neigh[i]) for i in len(np.arange(cell_cell_neigh))])<br />
Function returns total number of fibs neighbors

def cell_cell_type_neigbor (dataframe, id, cells_neigh_index, cell_type_column, type):<br />
    return (type_neigbor)<br />
Function returns number of neighbors of particular type

def cell_fibs_neigh (props_pruned, cells_centroid):<br />
    return (fibs_neigh_index)<br />
Function returns number of fibs neighbors

def cell_fibs_neigh_num (fibs_neigh_index):<br />
    return (number_of_fibs_neighbors)<br />
Function returns number of fibs neighbors

def cell_fibs_neigh_length (props_pruned, fibs_neigh_index):<br />
    return (length)<br />
Function calculates average length of neighboring fibers

def cell_fibs_neigh_angle (props_pruned, fibs_neigh_index):<br />
    return (angle)<br />
Function returns average angle of neighboring fibers

def cell_fibs_neigh_strightness (props_pruned, fibs_neigh_index):<br />
    return (straightness)<br />
Function returns average straightness of neighboring fibers

def cell_fibs_neigh_strightness(skeleton, dist, props_pruned, fibs_neigh_index):<br />
    return (thikness)<br />
Function returns average thickness of neighboring fibers

def cell_fibs_neigh_strightness(fibers_coords, nested_closest_neigbor_list, fibs_neigh_index):<br />
    return (linearity)<br />
Function returns average linearity of neighboring fibers
