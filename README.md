# Microsa

Project was developed by G.Vasiukov and S.Novitskiy

**Microsa** (**Micro**environment **s**patial **a**nalysis) is a package of useful functions for analysis of fibrous components of tissue microenvironment and combine that analysis in spatial dependent manner. The purpose of that library is to provide simplified tool  that gives an opportunity to combine the analysis of cellular and non-cellular components of tissue microenvironment.

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
-	Numpy >= 1.14;
-	SciPy >= 1.5.3;
-	scikit-image >= 0.14;
-	Pandas >= 1.1.3;
-	Matplotlib >= 2.0.

**Usage:**

Open a grayscale image, perform filtering, binarization and segmentation of fibers, calculate geometrical and spatial features of fibers, perform spatial analysis for other objects (cells).

**Fibers:**
  
- fibers_executor:<br />
  Function mainly performs segmentation of fibers and extract parameters that are required for other functions like number of       labels, distance between centroids etc.

- fibs_cell_neigh:<br />
  Function returns nested list that include fiber and indexes of cells’ (or other objects) neighbors. The distance is calculated     between cells centroid and closest pixel of fibers skeleton.

- fibs_cell_type_neighbor:<br />
  Function use the variable that was returned by fibs_cell_neigh function to classify neighbor objects (for example different type   of cells) 

- fibs_total_number:<br />
  Function returns total number of fibers that were segmented by fibers_executor 

- fibs_fibs_neigh_num:<br />
  Function calculates number of other fibers neighbor in the given radius. The distance is calculated between fibers centroids.

- fibs_length:<br />
  Function returns list of fibers length  

- fibs_angle:<br />
  Function returns list of fibers angle to against X axis  

- fibs_strightness:<br />
  Functions performs calculation of straightness of fibers. 

- fibs_thikness:<br />
  Function returns list of fibers thickness  

- fibs_linearity:<br />
  Function returns linearity between fiber of interest and other closest fiber neighbors. 

**Cells:**

- centroid_list_transform:<br />
  Function helps to put XY coordinates of objects as list of tuples.

- cells_cells_neigh:<br />
  Function returns list of cell neighbors indexes

- cell_fibs_total_num:<br />
  Function returns total number of fibs neighbors

- cell_cell_type_neigbor:<br />
  Function returns number of neighbors of particular type

- cell_fibs_neigh:<br />
  Function returns number of fibs neighbors

- cell_fibs_neigh_num:<br />
  Function returns number of fibs neighbors

- cell_fibs_neigh_length:<br />
  Function calculates average length of neighboring fibers

- cell_fibs_neigh_angle:<br />
  Function returns average angle of neighboring fibers

- cell_fibs_neigh_strightness:<br />
  Function returns average straightness of neighboring fibers

- cell_fibs_neigh_strightness:<br />
  Function returns average thickness of neighboring fibers

- cell_fibs_neigh_strightness:<br />
  Function returns average linearity of neighboring fibers
