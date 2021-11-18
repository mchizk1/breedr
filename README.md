# breedr

An R package for analyzing and visualizing pedigrees

### Package Contents:

* `nameR()` - A helper function that attempts to standardize names with pesky inconsistencies
* `plotigree()` - Quick assembly of pedigrees plots
* `breedr_f()` - Calculates Wright's F statistic from pedigree data
* `breedr_COP()` - Computes a Coefficient of Parentage matrix from pedigree data
* `replace_name()` - A function that replaces names using a key.  Helpful when going from selection numbers to cultivar names (similar to Excel's vlookup function).
* `habsburg` - A small set of pedigree data from an infamously inbred royal family...

### Installation Instructions

Simply run the following two lines of code in an R session

    install.packages("devtools")
    devtools::install_github("mchizk1/breedr")
    
### breedr input format - the Habsburg Dynasty example

Once upon a time, in a far away land, there lived a king who was very inbred. Run the following to see an example of the simple 3-column format required for breedr datasets.
Column 1 - individual ID
Column 2 - female parent ID
Column 3 - male parent ID

    library(breedr)
    head(habsburg)
       
### Plotigree() Demonstration

The `plotigree()` function uses the DiagrammeR package to make dynamic flowchart 
representations of family tree plots.

The "FULL" (default) method displays the entire tree with color coding for parental sex:

    plotigree(habsburg, "charles the bewitched", orientation = "LR")
    
![FULL](https://github.com/mchizk1/breedr/blob/main/FULL_method.png)
    
The "CA" method displays a similar tree, but is consolidated by common ancestry.
Common ancestors are highlighted in yellow:
    
    plotigree(habsburg, "charles the bewitched", method = "CA")
    
![CA](https://github.com/mchizk1/breedr/blob/main/CA_method.png)
       
### Inbreeding and Kinship Statistics

To calculate the inbreeding coefficients (F) of all unique individuals in a dataset, run this:

    breedr_F(habsburg)

And for a kinship matrix (coefficients of parentage), run this (this one can take a while): 

    breedr_COP(habsburg)
