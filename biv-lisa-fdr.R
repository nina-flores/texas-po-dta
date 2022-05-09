

#read in the data:

#-------------------------------------------------------------------------------
setwd("C:/Users/ninaf/Downloads/Power Outages Texas/R code/new-data-runs/data") # setting working directory to where my data is
require(tidyverse)
require(dplyr)
require(tidycensus)
require(ggplot2)
require(ggpubr)
require(sf)
require(tmap)
require(dplyr)
require(data.table)
require(spdep)
require(stringr)
require(spdep)
require(rgdal)
require(magrittr) # make sure to have this, or the matrix function wont work
require(gridExtra)
require(grid)
require(spatialreg)

#install.packages("rgeoda")
require(rgeoda)
require(stats)

#-------------------------------------------------------------------------------

###############################################################################
#Data Prep
#-------------------------------------------------------------------------------
#read in the SVI_SES data

dta <- read.csv("./LISA_dta-updated.csv") # this is a csv that includes all my independent variables of interest, as well
# as the average rate out per county, that I am going to standardize as a percentage
# and use as the dependent variable
dta$NAME <- dta$CountyName

dta <- dta %>% mutate(hosp_per_1000 = (count_hosp/tot_pop)*1000,
                      nf_per_1000 = (count_nf/tot_pop)*1000,
                      avg_CO_percentage = avg_CO_rate *100,
                      EDB_per_pop = (EDB / Beneficiaries)*100)

setwd("C:/Users/ninaf/Downloads/Power Outages Texas")
#join to spatial data:
texas_shapefile <- 
  read_sf(dsn = "./cb_2018_us_county_500k(2)/cb_2018_us_county_500k.shp") %>% # reading in US county shapefile from the census
  filter(STATEFP == "48") # filtering to just texas
texas_shapefile$NAME[texas_shapefile$NAME== "McCulloch"] <-"Mcculloch" # the csv and shapefile capitialize mcculloch differently,
# here I am adjusting for that

#join csv file to shapefile
texas_shapefile$GEOID <- as.integer(texas_shapefile$GEOID)
texas_dta <- full_join(texas_shapefile, dta)  #make sure when using full_join to list the shapefile/sf first, or spatial 
#remove sabine and sn aug

texas_dta <- texas_dta %>% filter(CountyName != "Sabine")
texas_dta <- texas_dta %>% filter(CountyName != "San Augustine")




#add ercot info:

non_ercot <- c("El Paso", 
               "Hudspeth", 
               "Gaines", 
               "Terry",
               "Yoakum",
               "Cochran", 
               "Hockley", 
               "Lamb", 
               "Bailey", 
               "Dallam", 
               "Sherman", 
               "Hansford", 
               "Ochiltree", 
               "Lipscomb", 
               "Hartley", 
               "Moore", 
               "Hutchinson", 
               "Hemphill", 
               "Bowie", 
               "Cass", 
               "Morris", 
               "Camp", 
               "Upshur", 
               "Gregg", 
               "Marion", 
               "Harrison", 
               "Shelby", 
               "San Augustine",
               "Panola",
               "Sabine", 
               "Trinity", 
               "Polk", 
               "Tyler", 
               "Jasper", 
               "Newton", 
               "Hardin", 
               "San Jacinto", 
               "Liberty", 
               "Jefferson", 
               "Orange")


#create dataset with only ercot data:
texas_dta <- texas_dta %>% mutate(ercot = ifelse(CountyName %in% non_ercot,"no", "yes"))
#-------------------------------------------------------------------------------

#create the queens matrix:
queen_w <- queen_weights(texas_dta)
summary(queen_w)

################################################################################
#apply this to biv lisa

#DME 

DME <- local_bimoran(queen_w, texas_dta[c('EDB_per_pop', 'avg_CO_percentage')], permutations  = 99999) #conducts the biv lisa
texas_dta$unadjusted_p.DME <- DME$p_vals #pulling out the determined p values

#FDR adjustment
texas_dta$adjusted_p.DME <- p.adjust(texas_dta$unadjusted_p.DME, method = "fdr") #adjusting pvales with fdr

#need the labels
texas_dta$labels.DME <- DME$c_vals # caution - I think the documentation of this code has HL and LH values flipped - I double checked this with past code and using GeoDa and I'm pretty sure 3 is LH and 4 is HL

#need to update the labels, so that those that are no longer significant with adjustment reflect that
texas_dta$labels.DME[texas_dta$adjusted_p.DME > .05] <- 0 #recoding the labels for anything now insignificant as not significant

################################################################################
#map it using the labels as our fill

DME_map<- ggplot(texas_dta) +
  geom_sf(data=texas_dta, aes(fill=as.character(labels.DME)), color="NA") +
  scale_fill_manual(values = c("grey80","red","dark blue",  "light blue", "pink", "grey50" ),
                    labels = c("Not significant", "High vulnerability - High outage %", "Low vulnerability - Low outage %", "Low vulnerability - High outage %", "High vulnerability - Low outage %","No data" ))+
  geom_sf(aes(fill = "grey0"), 
          data = texas_shapefile %>%filter(NAME == "Sabine" | NAME == "San Augustine"))+
  geom_sf(fill = "transparent", color = "black", size = 1,
          data = texas_dta %>% group_by(CountyName) %>% summarize())+
  geom_sf(fill = "transparent", aes(color = "#ffbf1f"), size = 1,
          data = texas_dta%>%dplyr::filter(ercot == "yes")%>% group_by(ercot) %>% dplyr::summarize())+  
  guides(fill = guide_legend(title="Bivariate LISA clusters"), color = guide_legend(title=""))+
  scale_color_manual(values = c( "#ffbf1f"), labels = c("ERCOT boundary"))+
  theme_minimal(base_size = 32) +
  
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank()) +
  
  theme(# Hide panel borders and remove grid lines
  panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) 

DME_map

#repeat with other variables of interest