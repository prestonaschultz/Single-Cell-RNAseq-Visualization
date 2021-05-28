# Single-Cell-RNAseq-Visualization
This is a visualization tool for single cell RNA sequencing data. For directions on how to use and information on obtaining datasets, please run the two lines of code below.

# Installing Shiny package
If this is your first time using this site or are unsure if you have the Shiny package downloaded, copy and paste this and run it in your R console:

```
install.packages("shiny", repos="http://cran.us.r-project.org")
library(shiny)
```

# Running the single cell tool:
Run this code in your R console. The single cell tool should pop up once everything is downloaded. You can see any changes occuring/errors in your R console window.

```
library(shiny)
runGitHub("Single-Cell-RNAseq-Visualization", "prestonaschultz", launch.browser = TRUE, ref="main")
```

This should launch the site as a web browser. If you look in your R console, you can see all of the package installations, changes, and working parts of the site. 

# Wait times for Heatmaps and Data tables
Please give time when using the "Heatmaps" or "Results/Data Table" features. When you use these two features, you will see the clusters downloading in your R console. Please wait for these to download before clicking on anything else.

