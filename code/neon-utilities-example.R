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
library(dplyr)

# Set working directory
setwd('~/GitHub/winter-journal-club-2018/data') # Path to repo

#NCDC temperature TMIN
temp = read.csv('temperature_NCDC.csv')


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
# Very little available for this summary product

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
## Gives good summary information about what data are available for the product
dp = neonUtilities::getProductInfo(dpID = "DP1.00003.001") # what is available?
# Very little available for this product

# Download all temperature data as a zip file; one per month per site
for (site in dp$siteCodes$siteCode){
  # For each site with data, download the expanded package
  neonUtilities::zipsByProduct(dpID = "DP1.00003.001", # abiotic summary data
                               site = site,
                               package = 'expanded',
                               check.size = F)
}
rm(site) # cleanup

# Squash all the files together so that there is one file per table
neonUtilities::stackByTable(filepath = list.files(full.names = T,
                                                  pattern = '00003'),
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

## SOIL MOISTURE -------
## Gives good summary information about what data are available for the product
dp = neonUtilities::getProductInfo(dpID = "DP1.00094.001") # what is available?

# Download all temperature data as a zip file; one per month per site
for (site in dp$siteCodes$siteCode){
  # For each site with data, download the expanded package
  neonUtilities::zipsByProduct(dpID = "DP1.00094.001", # abiotic summary data
                               site = site, 
                               package = 'expanded', 
                               check.size = F)
}
rm(site) # cleanup

# Squash all the files together so that there is one file per table
neonUtilities::stackByTable(filepath = list.files(full.names = T, 
                                                  pattern = '00094'), 
                            folder = T)

# Read in each consolidated file into the session
# and name it according to the table
files = list.files("./filesToStack00094/stackedFiles") # location of stacked files
for (x in files){
  # makes the name of each file
  nm = gsub('\\.csv', '_pub', x) # replace the '.csv' file ending with '_pub'
  # read in the csv file and assign the name accordingly
  assign(x = nm, 
         value = read.csv(file = list.files("./filesToStack00094/stackedFiles/", 
                                            full.names = T)[which(files%in%x)]))
}

## Gives good summary information about what data are available for the product
dp = neonUtilities::getProductInfo(dpID = "DP1.10086.001") # what is available?

# Download all temperature data as a zip file; one per month per site
for (site in dp$siteCodes$siteCode){
  # For each site with data, download the expanded package
  neonUtilities::zipsByProduct(dpID = "DP1.10086.001", # abiotic summary data
                               site = site, 
                               package = 'expanded', 
                               check.size = F)
}
rm(site) # cleanup

# Squash all the files together so that there is one file per table
neonUtilities::stackByTable(filepath = list.files(full.names = T, 
                                                  pattern = '10086'), 
                            folder = T)

# Read in each consolidated file into the session
# and name it according to the table
files = list.files("./filesToStack10086/stackedFiles") # location of stacked files
for (x in files){
  # makes the name of each file
  nm = gsub('\\.csv', '_pub', x) # replace the '.csv' file ending with '_pub'
  # read in the csv file and assign the name accordingly
  assign(x = nm, 
         value = read.csv(file = list.files("./filesToStack10086/stackedFiles/", 
                                            full.names = T)[which(files%in%x)]))
}

###### STEP 2. Make any necessary conversions 
# Join the sorting & trapping tables
mos = dplyr::left_join(x = mos_trapping_pub, 
                       y = mos_sorting_pub[,c("sampleID", "subsampleID", 
                                              "totalWeight", "subsampleWeight")], 
                       by = "sampleID")
# Add the individualCount of Culex tarsalis mosquitoes from the the ID table
for (sample in unique(mos_expertTaxonomistIDProcessed_pub$subsampleID
                      [mos_expertTaxonomistIDProcessed_pub$scientificName%in%"Culex tarsalis"])){
  mos$individualCount[mos$subsampleID%in%sample] = sum(mos_expertTaxonomistIDProcessed_pub$individualCount
                                                    [mos_expertTaxonomistIDProcessed_pub$subsampleID%in%sample], 
                                                    na.rm = T)
}
# Get an estimated count of mosquitoes (i.e. account for subsamples < totalWeight)
mos$estimatedAbundance[!is.na(mos$individualCount)] = ((mos$totalWeight[!is.na(mos$individualCount)]*
                                                         mos$individualCount[!is.na(mos$individualCount)])/
                                                         mos$subsampleWeight[!is.na(mos$individualCount)])
mos = mos[!is.na(mos$individualCount), ]
# Add in weather data
mos$date = substr(mos$collectDate,1,10)
mos = dplyr::left_join(mos,wss_daily_temp_pub)

# Take just the first date
mos$yearID = paste0(mos$siteID, substr(mos$collectDate,1,4))
# Sort by the date
mos = mos[order(mos$date), ]
# pick the first known occurrence each year
mos = mos[!duplicated(mos$yearID), ]
# Fix the date (change type from character to date)
for(d in unique(substr(mos$collectDate,1,10))){
  mos$juliandate[mos$date%in%d] = julian(as.Date(d), 
                                         origin = as.Date(paste0(format(as.Date(d), "%Y"),
                                                                 '-01-01')))
}

for(i in 1:nrow(mos)){
  mos$tempMin[i] = mean(temp[temp$siteID%in%mos$siteID[i]&
                               temp$date<=mos$date[i]&
                               (as.Date(mos$date[i])-14)<=temp$date,
                             'tempC'])
}

###### STEP 3. Conduct the analysis
## Question: What is the relationship between 
## Culex tarsalis emergence and temperature?

plot(mos$decimalLatitude, mos$juliandate)
summary(lm(mos$juliandate~mos$decimalLatitude))

