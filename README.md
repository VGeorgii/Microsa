# Microsa


Project was developed by G.Vasiukov and S.Novitskiy



**Microsa** (**Micro**environment **s**patial **a**nalysis) is a package of useful functions for analysis of fibrous components of tissue microenvironment and combine that analysis in spatial dependent manner. The purpose of presented package is to provide simplified tool which gives an opportunity to combine the analysis of cellular and non-cellular components of tissue microenvironment.



**Utilization:**
-	Segmentation of fibril-like objects (tensor method) in 2D;
-	Calculation of fibers’ centroid, length, width, angle and linearity in 2D;
-	Estimation of fibers' features and providing a detailed report in pandas.DataFrame format;
-	Implementation of spatial analysis between cellular and non-cellular objects.



**Installation:**
- Using pip:<br />
    $ pip install Microsa<br />
- Using Git repository:<br />
    $ git clone https://github.com/VGeorgii/Microsa.git<br />
    $ cd Microsa<br />
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



**Fiber segmentation module:**
  
- fibers_executor (image):<br />

  Function mainly performs segmentation of fibers and extract parameters that are required for other functions like number of labels, distance between centroids etc.

  arguments:
  image – gray scale image to process

  function returns:
  'skeleton': image with skeletonized objects, 
  'distance': distance from skeleton to the edge of original object (radius), 
  'skel_labels_pruned': pruned and labeled skeletonized objects, 
  'props_pruned': properties of labeled objects, 
  'nlabels_pruned': number of labeled objects



**Fiber geometrical feature calculation:**

- fibs_geom (executed_fibs, radius):<br />

  Function returns dataframe which contains information about fibers

  arguments:
  - executed_fibs: array with calculated features of labeled fibers
  - radius: radius of neighborhood outline    

  function returns:
  - pd.DataFrame with calculated fibers features (number, length, angle, strightness, thickness, linearity)



**Spatial:**

- fibs_spatial (cells_coords_list, executed_fibs, radius, cell_type = 'None', cell_type_list = 'None'):<br />

  Function returns dataframe which contains information about neighboring cells, their type, and features

  arguments:
  - cells_coords_list: list of cells coords
  - executed_fibs: array with calculated features of labeled fibers
  - radius: radius of neighborhood outline
  - cell_type: default 'None', list of cell types for spatial analysis. 'All' make function perform calculation for all types of     cell
  - cell_type_list: default 'None', list(column) with cell types

  function returns:
  - pd.DataFrame with calculatedd spatial information of neighboring cells
 

- cell_cell_spatial (cells_coords_list, radius, cell_type = 'None', cell_type_list = 'None', cell_feature_list = 'None'):<br />

  Function returns dataframe which contains information about neighboring cells, their type, and features

  arguments:
  - cells_coords_list: list of cells coords
  - radius: radius of neighborhood outline
  - cell_type: default 'None', list of cell types for spatial analysis. 'All' make function perform calculation for all types of     cell
  - cell_type_list: list(column) with cell types
  - cell_feature_list: list(colum) of feature of cells    

  function returns:
  - pd.DataFrame with calculatedd spatial information of neighboring cells
 

- cell_fibs_spatial (executed_fibs, cells_coords_list, radius):<br />
  
  Function returns dataframe which contains information about neighboring fibers and their features

  arguments:
  - executed_fibs: dictionary executed with fiber_executer function that contains information about segmented fibers
  - cells_coords_list: list with cells coordinates (y,x)
  - radius: radius of neighborhood outline    

  function returns:
  - pd.DataFrame with calculatedd spatial information of neighboring fibers



**Example of usage:**

import numpy as np<br />
from skimage import io<br />
from matplotlib import colors<br />
from skimage.morphology import medial_axis<br />
from skimage.measure import regionprops<br />
from skimage.measure import label<br />
from scipy import ndimage<br />
from skimage.graph import route_through_array<br />
from scipy.ndimage import binary_closing, binary_hit_or_miss<br />
from skimage.filters import frangi<br />
from scipy.spatial import distance<br />
import pandas as pd<br />
import matplotlib<br />
from matplotlib import pyplot as plt<br />



For example, we have dataframe with information about cells (localization (cells_coords), type of cell (cell_type), and features of cells (cell_feature_1, cell_feature_2))

![Table_1](https://user-images.githubusercontent.com/65576385/121410782-91fbb880-c928-11eb-97c7-ddf9229ab30d.PNG)



**Generating lists form dataframe columns**

cells_coords = list(dataframe.loc[:, 'cells_coords'])<br />
cell_type = list(dataframe.loc[:, 'cell_type'])<br />
cell_feature_1 = list(dataframe.loc[:, 'cell_feature_1'])<br />
cell_feature_2 = list(dataframe.loc[:, 'cell_feature_2'])<br />



**Importing image sample**

img = io.imread('C:/Users/test_image.tif')<br />
fig, ax = plt.subplots(figsize = (100,100))<br />
ax = plt.imshow(img)<br />

![set1_FIB_HEL](https://user-images.githubusercontent.com/65576385/121410856-9e801100-c928-11eb-8450-7878828fb97d.png)



**Image filtering (Frangi filter implementation)**

gray_frangi = np.array(frangi(img, sigmas=range(4, 6, 10), gamma = 25, black_ridges = False))<br />
gray_frangi_bit = (gray_frangi/255) > 0.000000001)<br />
fig, ax = plt.subplots(figsize = (50,50))<br />
ax = plt.imshow(gray_frangi_bit, cmap = 'gray')<br />

![frangi_fltr](https://user-images.githubusercontent.com/65576385/121409479-30871a00-c927-11eb-840a-d6dc910a96bc.png)



**Fibers skeletonization and labeling**

fibere_exe = fibers_executor (gray_frangi_bit)<br />
fig, ax = plt.subplots(figsize = (100,100))<br />
ax = plt.imshow(np.multiply(fibere_exe['skel_labels_pruned'] > 0, 1), cmap = 'gray')<br />
plt.axis('off')<br />
plt.show()<br />


![pruned_skeleton](https://user-images.githubusercontent.com/65576385/121409573-498fcb00-c927-11eb-897e-8a78424296e4.png)



**Fiber geometry**

fbs_mrph = fibs_geom (fibere_exe, 25)<br />

![Table_3](https://user-images.githubusercontent.com/65576385/121409713-688e5d00-c927-11eb-8513-0d08c429fade.PNG)



**Fiber spatial analysis**

fib_spat = fibs_spatial (cells_coords, fibere_exe, 25, cell_type = 'All', cell_type_list = cell_type)<br />

![Table_2](https://user-images.githubusercontent.com/65576385/121409949-ac816200-c927-11eb-8252-6be29b3e4a81.PNG)



**Cell spatial analysis**

cell_cell_spt = cell_cell_spatial (cells_coords, 25, cell_type = ['green', 'red'], cell_type_list = cell_type, cell_feature_list = [cell_feature_1, cell_feature_2])<br />

![Table_4](https://user-images.githubusercontent.com/65576385/121410054-ca4ec700-c927-11eb-972b-f341d89d1c2d.PNG)

cll_fbs_sptl = cell_fibs_spatial (fibere_exe, cells_coords, 25)<br />

![Table_5](https://user-images.githubusercontent.com/65576385/121410231-fe29ec80-c927-11eb-9a69-30a39c06b5d0.PNG)



**Visualization**

vis = visualization (gray_frangi_bit, cells_coords, cell_type, fibere_exe, 5)<br />

![overlaid_map](https://user-images.githubusercontent.com/65576385/121410307-113cbc80-c928-11eb-8a0b-51e2702ae168.png)







