## Example question:
## How does temperature (DP4.00001.001) predict
## Culex tarsalis occurrence (DP1.10043.001)?

###### STEP 1. Download Relevant NEON data
## Set up script
rm(list=ls())
options(stringsAsFactors = F)

## Install packages (if necessary)
# install.packages("devtools")
# install.packages("gtools")
# install.packages("neonUtilities")
# install_github(repo = 'NEONScience/NEON-geolocation/geoNEON')
# library(geoNEON)

## Load Libraries
library(devtools)
library(neonUtilities)

# Set working directory
setwd('~/GitHub/winter-journal-club-2018/data') # Path to repo

## Get the Data from NEON
## Mosquito Data = "DP1.10043.001"; Temperature Summary Data = "DP4.00001.001"

## MOSQUITOES --------
## Gives good summary information about what data are available for the product
dp = neonUtilities::getProductInfo(dpID = "DP1.10043.001") # what is available?

# Download all mosquito data as a zip file; one per month per site
for (site in dp$siteCodes$siteCode){
  # For each site with data, download the expanded package
  neonUtilities::zipsByProduct(dpID = "DP1.10043.001", # mosquito CO2 trapping
                               site = site, 
                               package = 'expanded', 
                               check.size = F)
}
rm(site) # cleanup

# Squash all the files together so that there is one file per table
neonUtilities::stackByTable(filepath = list.files(full.names = T, 
                                                  pattern = '10043'), folder = T)

# Read in each consolidated file into the session
# and name it according to the table
files = list.files("./filesToStack10043/stackedFiles/") # location of stacked files
for (x in files){
  # makes the name of each file
  nm = gsub('\\.csv', '_pub', x) # replace the '.csv' file ending with '_pub'
  # read in the csv file and assign the name accordingly
  assign(x = nm, 
         value = read.csv(file = list.files("./filesToStack10043/stackedFiles/", 
                                            full.names = T)[which(files%in%x)]))
}

## TEMPERATURE -------
## Gives good summary information about what data are available for the product
dp = neonUtilities::getProductInfo(dpID = "DP4.00001.001") # what is available?

# Download all temperature data as a zip file; one per month per site
for (site in dp$siteCodes$siteCode){
  # For each site with data, download the expanded package
  neonUtilities::zipsByProduct(dpID = "DP4.00001.001", # abiotic summary data
                               site = site, 
                               package = 'expanded', 
                               check.size = F)
}
rm(site) # cleanup

# Squash all the files together so that there is one file per table
neonUtilities::stackByTable(filepath = list.files(full.names = T, 
                                                  pattern = '00001'), 
                            folder = T)

# Read in each consolidated file into the session
# and name it according to the table
files = list.files("./filesToStack00001/stackedFiles") # location of stacked files
for (x in files){
  # makes the name of each file
  nm = gsub('\\.csv', '_pub', x) # replace the '.csv' file ending with '_pub'
  # read in the csv file and assign the name accordingly
  assign(x = nm, 
         value = read.csv(file = list.files("./filesToStack00001/stackedFiles/", 
                                            full.names = T)[which(files%in%x)]))
}

###### STEP 2. Make any necessary conversions 


###### STEP 3. Conduct the analysis
## Question: What is the relationship between 
## Culex tarsalis emergence and temperature?

