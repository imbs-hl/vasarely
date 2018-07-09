# Vasarely Chart

Alexandra Banach


## Introduction
Vasarely is a implementation to produce a so called vasarely chart. It was developed by Manaster et al. in 1999 to check Hardy-Weinberg-Law for genetic data and to discover null alleles. Thus, for some genetic input data expected and real probabilty of allel frequencies can be calculated and plotted. For all possible combinations of allels the expected probabilty is described by a geom_raster plot. Each cell of the plot has a color which depends on the value of the expected probability. The real probability is described by a dot in the corresponding cell which also has a color that depends on the calculated probability. The chart and its colors for the two probabilities helps to check how expected and real probability differ from each other. The more the color differ from each other the more the both probabilities differ. In general the function can be used to represent expected and real probability for categorial variables but not only genetic data.


## Installation
To install the vasarely R package from GitHub using devtools, just run



	devtools::install_github("imbs-hl/vasarely")

## Usage


An example for the usage of the vasarely package can be found by running
	help(vasarely)
or 
?vasarely
. It is also explained which format for the input data is required and which optional parameters can be chosen for the function. If you have any questions or if you find any bugs, please contact us. 


## References


* Manaster,C.J., Nanthakumar, E., & Morin, P.A.(1999). Detecting null alleles with vasarely charts. In: IEEE Visualization, D.S. Ebert, M.Gross, &B. Hamann, ed., pp. 463-466. IEEE Computer Society: San Francisco
