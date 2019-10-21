# fcs2 - Fisheries classification scheme 2

Provides functions that carry out SNIFFER's implementation of the Environment Agency's Fisheries Classification Scheme 2 (FCS2). This package was developed for use in Scotland, Northern Ireland and the Republic of Ireland as part of SNIFFER project WFD68c: Science Work.

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

If creating new models (rather than running the classification with existing models) install [openBUGS](http://www.openbugs.net/w/Downloads) on your computer.

See `fcs2FitModel` function for more details.

## Contributing 

Please read the [Contributing guidelines](CONTRIBUTING.md) file for more details 

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
