# Introduction to the rbsrp R Package

## Overview
R package "rbsrp" provides a procedure to select response variables for ultrahigh dimensional multi-outcome data, where both the dimensions of response and predictor outcomes are substantially greater than the sample size. It also provides the multi-test procedure including Benjamini-Hochberg and Bonferroni correction.

## Installation
You can install the rbsrp package in R using the following steps:
    install.packages("devtools")
    library(devtools)
    install_github("xliusufe/rbsrp")

## Usage
The detail of the usage of the package can be found in https://github.com/xliusufe/rbsrp/inst/rbsrp-manual.pdf.


# Replicating Table 1 using table1.r

## Overview
The `table1.r` script is used to replicate Table 1 in the main paper.

## Usage
1. Ensure that you have the necessary R packages installed, such as 'rbsrp', 'stats', 'MASS' and 'Matrix'.

2. Navigate to the directory where the `table1.r` script is located using the `setwd()` function if necessary.

3. Open your R terminal and run the script using the following commands:

    ```
    Rscript table1.r
    ```
    Or you can select the entire script and run it in the R environment.

4. After running the script, the output will be generated, typically in the form of a `.csv` file or printed directly to the terminal. Due to the number of seeds, the results may be different from those in Table 1.


# Replicating Table 2 using table2.r

## Overview
The `table2.r` script is used to replicate Table 2 in the main paper. 

## Usage
1. Ensure that you have the necessary R packages installed, such as 'rbsrp', 'stats', 'MASS' and 'Matrix'.

2. Navigate to the directory where the `table2.r` script is located using the `setwd()` function if necessary.

3. Open your R terminal and run the script using the following commands:

    ```
    Rscript table2.r
    ```
    Or you can select the entire script and run it in the R environment.

4. After running the script, the output will be generated, typically in the form of a `.csv` file or printed directly to the terminal. Due to the number of seeds, the results may be different from those in Table 2.


# Replicating Table 4 using table4.r

## Overview
The `table4.r` script is used to replicate Table 4 in the main paper.

## Usage
1. Ensure that you have the necessary R packages installed, such as 'rbsrp', 'stats', 'MASS' and 'Matrix'.

2. Navigate to the directory where the `table4.r` script is located using the `setwd()` function if necessary.

3. Open your R terminal and run the script using the following commands:

    ```
    Rscript table4.r
    ```
    Or you can select the entire script and run it in the R environment.

4. After running the script, the output will be generated, typically in the form of a `.csv` file. 

# Replicating all Tables about the real data analysis

## Overview
The `realdata.r` script is used to replicate Tables 2-23 in the Supplementary Material.

## Usage
1. Ensure that you have the necessary R packages installed, such as 'rbsrp', 'stats', 'MASS' and 'Matrix'.

2. Navigate to the directory where the `realdata.r` script is located using the `setwd()` function if necessary.

3. Open your R terminal and run the script using the following commands:

    ```
    Rscript realdata.r
    ```
    Or you can select the entire script and run it in the R environment.

4. After running the script, the output will be generated, typically in the form of `.csv` files. 
