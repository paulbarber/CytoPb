# CytoPb: Cell Probabilities
For multiplexed imaging (IF, IMC, ...), this code uses marker probability densities to generate cell type maps and interaction signatures in a pixel-based, segmentation-free and threshold-free (gating-free) method.

The usual method for cell type determination in multiplexed images is to segment cells using a nuclear and/or membrane marker, and maybe use a random forest or deep learning classifier. Then to use channel gating (thresholding) to determine if the segmented cells are positive or negative for particular markers (like the cell differentiation (CD) markers). Making hard (irreversible) decisions about which pixels constitute particular cells and which levels of marker are positive is sometimes difficult and subjective. Especially if your marker staining did not turn out well.

We would rather use a probabilistic approach, maybe it could be called "fuzzy logic", to not make these decisions or to delay any hard decisions until the last moment in the analysis. Here, that approach is used to produce the usual cell abundance bar graphs that agree well with the usual method, and cell type maps that can be used for neighbourhood analysis etc.

This method copes with all the complex marker patterns usually used in immunology, including enforcing that a cell type should be negative for certain markers, in an easy to define matrix.

## Quick start instructions

Clone or download the repository files. Set your R working directory to the source code location. Copy "CytoPb_RUN ALL.R" to your data folder and modify it as required for the scripts to run and required options.

Check the folder of "example_input_files". panel.csv may be automatically created by one of the import scripts (Script 0). "cell_type_matrix.csv" defines the cell type of interest by positive (1) and negative (-1) marker channels. "cell_type_colours.txt" is an optional file to define the colours to use for the cell types.

Script 0 will expect a "raw" folder of raw data within that directory. The scripts work mainly with stacked tiff files which will and should be placed in an "img" folder. Use each script in turn. Read the comments at the top of each script about what the inputs and outputs are.

## How to cite

Please cite the repository URL (https://github.com/paulbarber/CytoPb). There are no publications associated with the code yet. 


## Licence

CytoPb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

CytoPb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with CytoPb; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
