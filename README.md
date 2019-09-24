# Work in progress 

This R package is in development. Please don't use in production.

<<<<<<< HEAD
# fcs2 - Fisheries classification scheme 2
=======
# fcs2 :fish: 
>>>>>>> af042f3d47ff323520e3e363b20abfc9738885f3

Provides functions that carry out SNIFFER's implementation of the Environment Agency's Fisheries Classification Scheme 2 (FCS2). This package was developed for use in Scotland, Northern Ireland and the Republic of Ireland as part of SNIFFER project WFD68c: Science Work.

## Quick Start

In R: Install the development version from GitHub:

```
install.packages("devtools")
library(devtools)
install_github("aquaMetrics/fcs2")
```

Run the demo dataset through the `calcClassScot` function to get classification output for Scotland:

```
results <- calcClassScot(data = fcs2::demo_data)
```
Or to run your own data (save as .csv to working directory):

```
data <- read.csv(file = "YOUR_FILE.csv")

results <- calcClassScot(data)

write.csv(results,file = "YOUR_REUSULTS.csv")
```

## Contributing 

Please read the [Contributing guidelines](CONTRIBUTING.md) file for more details 


## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
