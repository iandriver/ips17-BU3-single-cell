# Single Cell data for "Prospective isolation and single cell profiling of NKX2-1+/CD47+ human lung progenitors derived from pluripotent stem cells"
------
The data and associated scripts that were used to create the single cell figures
in the paper "Prospective isolation and single cell profiling of NKX2-1+/CD47+
human lung progenitors derived from pluripotent stem cells"

Interactive version of the monocle data can be found here: https://iandriver.shinyapps.io/ips17-bu3_single_cell_data/

Instructions
-----
Clone or download this repository and move working directory there:

```bash
    git clone https://github.com/iandriver/ips17-BU3-single-cell.git
    cd ips17-BU3-single-cell
```
or if working in RStudio:
```R
    setwd('path_to/ips17-BU3-single-cell')
```
Then run the script "monocle2_ips17_BU3.R".

```bash
    Rscript monocle2_ips17_BU3.R
```
or run in in RStudio.

This will rerun the analysis, but the pseudotime ordering is non-dimensional, so the resulting monocle object may be reversed. Rerun these two commands to flip the trajectory until it is correct:

```R
    ips17_BU3_data <- reduceDimension(ips17_BU3_data, max_components=3)
    ips17_BU3_data <- orderCells(ips17_BU3_data, reverse=TRUE)
```
Or simply load the existing monocle CellDataSet:

```R
    ips17_BU3_data <- readRDS('kmeans2_new_allgene200.RDS')
```
This can then be used for any plotting or analysis.
Simple interactive plotting of this data can also be found here: https://iandriver.shinyapps.io/ips17-bu3_single_cell_data/
