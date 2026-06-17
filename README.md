CSRcalculator

CSRcalculator provides three functions for calculating plant CSR (Competitor–Stress-tolerator–Ruderal) strategies using three published models: Hodgson, MorphoPhys, and StrateFy.
Each function takes a species × traits data frame and returns C%, S%, R% values and a CSR strategy class appended to the original data frame.

Installation (local ZIP)

Unzip the folder CSRcalculator
in R:
install.packages("path/to/CSRcalculator", repos = NULL, type = "source")
Then load:
library(CSRcalculator)

Functions

hodgson(data, calcSLA = TRUE, calcLDMC = TRUE)
morphoPhys(data, calcLDMC = TRUE)
strateFy(data, calcSLA = TRUE, calcLDMC = TRUE)
e.g. to use strateFy for a data frame containing species, LFW, LDW and LA: 
# Run the StrateFy CSR model
csrData <- strateFy(myData)
# View results
head(csrData)
