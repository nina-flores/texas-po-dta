#This script rund through an example of one of the autologisitc models

#set working directory
#read in data:
setwd("C:/Users/ninaf/Downloads/Power Outages Texas/R code/new-data-runs/data")
require(tidyverse)
require(dplyr)
require(tidycensus)
require(sf)
require(dplyr)
require(lubridate)
require(spdep)
require(spind)
require(gtools)
require(surveillance)
require(ngspatial)

dta<- read.csv("./po-data-texas-feb-2021.csv")

dta<-dta %>% mutate(RecordDateTime2 = ymd_hms(RecordDateTime))
#####
#this first step here generates a variable that identifies whether a county
#meets the 10,000+ threshold for 24+ hours:
###

###############################################################################
###############################################################################
# 10,000 out threshold
###############################################################################
CustomersOut_threshold = 10000 #setting threshold at 10,000 - for relative regressions, set this to .2, and use customers out rate instead of customers out

data_thres_24 <- dta %>%  arrange(CountyName, RecordDateTime2) %>% #arranging by county and time - this order allows us to find consecutive hours by county
  group_by(CountyName, grp = rleid(CustomersOut >= CustomersOut_threshold)) %>%   # grabs all counties meeting this threshold
  filter(n() >= 24, all(CustomersOut >= CustomersOut_threshold)) %>%              # keeps counties meeting threshold for 24 hours
  slice(1:24) %>% ungroup()

#Get unique counties that reach this:
data_thresh_counties_24 <- data_thres_24 %>% dplyr::select(CountyName) %>% unique()
counties_24 <- data_thresh_counties_24$CountyName                                   #creates vector of all counties meeting thresholds

#make the variable
texas_dta_thresh <- dta %>% mutate(dep = ifelse(dta$CountyName %in% counties_24, 1, 0)) #generating new variable with 1 indicating a county meeting the threshold-- this will be the depdendent variable

#only need dep and indep variables as unique observations
texas_dta_thresh <- texas_dta_thresh %>% dplyr::select(c("CountyName",
                                                         "Customers_Served", 
                                                         "overall_SVI", 
                                                         "DME_users_percent", 
                                                         "pct_black",
                                                         "urban_rural",
                                                         "POP_DENSITY",
                                                         "pct_hispanic",
                                                         "nursing_facilities_per_total_pop",
                                                         "hospitals_per_total_pop",
                                                         "dep")) %>% unique()


#####
#need to pull in a map for spatial information and join

###
setwd("C:/Users/ninaf/Downloads/Power Outages Texas")

texas <- read_sf(dsn = "./cb_2018_us_county_500k(2)/cb_2018_us_county_500k.shp") %>%# plugging a US census map and filtering to texas
  filter(STATEFP == "48")%>% arrange(NAME) 

#steps to ensure capitalization is correct so names match up
texas <- texas %>% dplyr::select(GEOID, NAME, geometry)
texas <- texas %>% mutate(CountyName = NAME)
#fixing mismatched capitalization of Mcculloch between df
texas$NAME[texas$NAME == "McCulloch"] <- "Mcculloch"
texas$CountyName[texas$CountyName == "McCulloch"] <- "Mcculloch"


###
#join spatial info to data
###

texas_dta <- full_join(texas, texas_dta_thresh)

#since some counties will not have 10,000 customers need to
#filter just to those with customers above 10,000
texas_dta<- texas_dta %>% filter(Customers_Served >= 10000) #remove this line for relative models

#################################################
###
##generate spatial weights matrix that will be used in the model:
tx_mat <- poly2adjmat(texas_dta)  # generates adjacency matrix from polygon data
###



##################################################
###Set up model
##################################################
pct_hisp_mod <- autologistic(dep ~pct_hispanic+ POP_DENSITY + urban_rural,   # formula
                          data = texas_dta,                                  #data
                          A = tx_mat,                                        #specifying adjacency matrix
                          control = list(confint = "bootstrap", bootit = 500, nodes = 6)) # specifies bootstrapped CI with 500 samples
summary(pct_hisp_mod) 

###link to rdrr.io: https://rdrr.io/cran/ngspatial/man/autologistic.html 
###link to vingette: https://cran.r-project.org/web/packages/ngspatial/ngspatial.pdf 


