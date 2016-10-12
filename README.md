# Plotting multiple principal components (PC) from classical multidimensional scaling (cmd) output and mark outlier clusters based on PC1 vs. PC2

**The function is mainly designed as integration to the [GenABEL workflow](https://zenodo.org/record/19738#.V_3eyvmLSHs)**

**More specifically, it points to the detection of genetic outliers using identity by state (IBS) procedure (see quality control, section 5.2 of the tutorial).**

## Description
Although the function is embedded in the GenABEL workflow, it can be used for every cmd- or principal component matrix.
It is important to have the sample names or ID names within the rownames of your matrix.
The function plots different mds plots to verify genetic outliers obtained in a PC1 vs. PC2 plot.
Further it returns a vector of IDs which represents the main sample-cluster (optional), the samples which should be kept in the analysis.

The input-arguments of the function are described within the code.
Since this function is not part of a official package one should read the comments written in the first lines.

## Installing
Just download the function itself which is stored in the folder `outlier_function`

Once stored on your hard drive use:
```
source("path/to/where/the/function/is/stored")
```
The function should be loaded to the global environment of your current R-session.

## Code comments
Within the function I extensively describe what the specific code blocks do.
