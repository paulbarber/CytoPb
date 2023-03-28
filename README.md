# CytoPb: Cell Probabilities
For multiplexed imaging (IF, IMC, ...), this code uses marker probability densities to generate cell type maps and interaction signatures in a pixel-based, segmentation-free and threshold-free (gating-free) method.

The usual method for cell type determination in multiplexed images is to segment cells using a nuclear and/or membrane marker, and maybe use a random forest or deep learning classifier. Then to use channel gating (thresholding) to determine if the segmented cells are positive or negative for particular markers (like the cell differentiation (CD) markers). Making hard (irreversible) decisions about which pixels constitute particular cells and which levels of marker are positive is sometimes difficult and subjective. Especially if your marker staining did not turn out well.

We would rather use a probabilistic approach, maybe it could be called "fuzzy logic", to not make these decisions or to delay any hard decisions until the last moment in the analysis. Here, that approach is used to produce the usual cell abundance bar graphs that agree well with the usual method, and cell type maps that can be used for neighbourhood analysis etc.

This method copes with all the complex marker patterns usually used in immunology, including enforcing that a cell type should be negative for certain markers, in an easy to define matrix.

## Quick start instructions

Clone or download the repository files. Set your R working directory to the source code location. Copy "CytoPb_RUN ALL.R" to your data folder and modify it as required for the scripts to run and required options.

Check the folder of "example_input_files". panel.csv may be automatically created by one of the import scripts (Script 0). "cell_type_matrix.csv" defines the cell type of interest by positive (1) and negative (-1) marker channels. "cell_type_colours.txt" is an optional file to define the colours to use for the cell types.

The scripts work mainly with stacked tiff files which will and should be placed in an "img" folder. Use each script in turn. Read the comments at the top of each script about what the inputs and outputs are.

Script 0 is data import. It will expect a "raw" folder of raw data and will make tiff files in the img folder and a panel.csv file. 

Scripts 1, 2 and 3 run through producing the cell type probability maps based on a set of tiff stacks in the img folder and a panel.csv file. After script 2 you can check the channel png folder to see if reasonable positive marker areas have been identified. If not, the pos_value_table.csv file can be adjusted manually (decrease the value to assign more positive signal in this channel).

Script 2a is unfinished, but may allow for exploratory analysis by clustering pixels into cell type clusters.

Script 4 will measure colocalisation between pairs of cell types as defined within the script using a scale space approach.

## Outputs

* celltype_maps
Images of the most likely cell type at each pixel. Cell type colours will correspond to those in the Cell Total Plots pdf file.

* Cell Total Plots.pdf and csv
Satcked bar graphs of the total probability, probability density (normalised to image area) and area occupied in the cell type map (most likely cell type per pixel). Raw data behind the Cell Total Plots is in the csv file.

* celltype_png
Images of the cell type probability for every cell type and every image.

* channel_png
Images of the channel probability of positivity for every marker (channel) of every image.

* objects folder
A folder of R objects (cytomapper image lists) that contain the raw data behind the celltype_png and channel_png images. For follow on processing.

* Marker per CellType pdfs and csv
The raw marker expressions behind the pixels of the celltype_maps. Demonstrates which markers tend to be higher for each cell type compared to the other defined cell types (note, not compared to all pixels. There may be other cell types, not defined in the cell type matrix that have higher expressions).

* Positive Value Plot.pdf
A qualitative heatmap of the channel intensities defined to be marker positive (as defined in pos_value_plot.csv). Here is a visual check of which markers have high levels and those which may have failed or have no content.

## How to cite

Please cite the repository URL (https://github.com/paulbarber/CytoPb). There are no publications associated with the code yet. 


## Licence

CytoPb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

CytoPb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with CytoPb; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
