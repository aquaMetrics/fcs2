# fcs2 - Fisheries classification scheme 2

This package provides a *limited* update to the original fcs2 packaged to allow it to run fish classification on R >= 3.4. This package was developed for use in Scotland, Northern Ireland and the Republic of Ireland as part of SNIFFER project WFD68c: Science Work. It carries out SNIFFER's implementation of the Environment Agency's Fisheries Classification Scheme 2 (FCS2). 
  
The original package provided by HR Wallingford ran on R <= 2.1.

## Future Developments

Currently, this package only runs the classification for Scotland. The original package also provides functions to rebuild the underlying classificaiton model based on new fish survey data. These functions no longer work due to changes to the dependency (INLA and openBUGS) in the intervening years. See github issues for outstanding tasks.

## Quick Start

In R: Download and install [INLA](http://www.r-inla.org/download) package.

Install the development fcs2 package from GitHub:

```
install.packages("devtools")
library(devtools)
install_github("aquaMetrics/fcs2")

```
Run the demo dataset using the `calcClassScot` function to get classification output based on predictive models for Scotland:

```
results <- calcClassScot(data = fcs2::demo_data)
```
Or to run your own data (save data as .csv file to working directory):

```
# read your survey results from working directory
data <- read.csv(file = "YOUR_FILE.csv")

# run classification using fish prediction model for Scotland
results <- calcClassScot(data)

# write the classification result to working directory
write.csv(results, file = "YOUR_RESULTS.csv")
```

## Optional Dependencies

If creating new models (rather than running the classification with existing models) install [openBUGS](http://www.openbugs.net/w/Downloads) on your computer. NOTE - using openBUGS with this package has not been tested -  see issues.

See `fcs2FitModel` function for more details.

## Contributing 

Please read the [Contributing guidelines](CONTRIBUTING.md) file for more details 

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
