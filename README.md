# CytoPb: Cell Probabilities
For multiplexed imaging (IF, IMC, ...), use of marker probability densities to generate cell type maps and interaction signatures in a pixel-based, segmentation-free and threshold-free (gating-free) method.

The usual method for cell type determination in multiplexed images is to segment cells using a nuclear and/or membrane marker, and maybe use a random forest or deep learning classifier. Then to use channel gating (thresholding) to determine if the segmented cells are positive or negative for particular markers (like the cell differentiation (CD) markers). Making hard (irreversible) decisions about which pixels constitute particular cells and which levels of marker are positive is sometimes difficult and subjective. Especially if your marker staining did not turn out well.

We would rather use a probabilistic approach, maybe it could be called "fuzzy logic", to not make these decisions or to delay any hard decisions until the last moment in the analysis. Here, that approach is used to produce the usual cell abundance bar graphs that agree well with the usual method, and cell type maps that can be used for neighbourhood analysis etc.

This method copes with all the complex marker patterns usually used in immunology, including enforcing that a cell type should be negative for certain markers, in an easy to define matrix.

## Quick start instructions

Clone or download the repository files. Set your R working directory to the location of your data. Script 0 will expect a "raw" folder of raw data within that directory. The scripts work mainly with stacked tiff files which will and should be placed in an "img" folder. Use each script in turn. read the comments at the top of each script about what the inputs and outputs are.

