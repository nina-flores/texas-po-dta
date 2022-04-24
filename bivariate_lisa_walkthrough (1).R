#Quick run through setting up bivariate lisa with texas SVI and POUS data


#read in the data:

setwd("C:/Users/ninaf/Downloads/Power Outages Texas") # setting working directory to where my data is
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



###############################################################################
#Data Prep
###############################################################################
#read in the SVI_SES data

dta <- read.csv("./LISA_dta.csv") # this is a csv that includes all my independent variables of interest, as well
# as the average rate of customers out per county, that I am going to standardize as a percentage
# and use as the dependent variable

dta$NAME <- dta$CountyName # for the future join, we need a common variable to merge the csv and the shapefile on. In the shapefile,
                           # all counties are under the column "NAME" so I am creating a "NAME" column with the county info
                           # in the csv data as well

dta <- dta %>% mutate(hosp_per_1000 = (count_hosp/tot_pop)*1000, #standardizing independent variables
                      nf_per_1000 = (count_nf/tot_pop)*1000,
                      avg_CO_percentage = avg_CO_rate *100) #standardizing dependent variables


#join to spatial data:
texas_shapefile <- 
  read_sf(dsn = "./cb_2018_us_county_500k(2)/cb_2018_us_county_500k.shp") %>% # reading in US county shapefile from the census as a shapefile
  filter(STATEFP == "48") # filtering to just texas
texas_shapefile$NAME[texas_shapefile$NAME== "McCulloch"] <-  "Mcculloch" # the csv and shapefile capitialize mcculloch differently,
# here I am adjusting for that

#join csv file to shapefile

texas_dta <- full_join(texas_shapefile, dta)  #make sure when using full_join to list the shapefile/sf first, or spatial 
                                              #information will not be retained



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


texas_dta <- texas_dta %>% mutate(ercot = ifelse(CountyName %in% non_ercot,"no", "yes")) # adding a column with a binary yes/no for ercot
                                                                                         # classification-- we add a layer for this later



##############################################################################
#This is what will be consistently used, so I am moving to the beginning:

#This is programming code from rafapereirabr at https://gist.github.com/rafapereirabr/5348193abf779625f5e8c5090776a228
#======================================================
# Adjacency Matrix (Queen)
#We will use his for each different comparison

nb <- poly2nb(texas_dta)
lw <- nb2listw(nb, style = "B", zero.policy = T)
W  <- as(lw, "symmetricMatrix")
W  <- as.matrix(W/ Matrix::rowSums(W))
W[which(is.na(W))] <- 0




#======================================================
# Programming some functions

# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  
  xp <- scale(x)[, 1]
  yp <- scale(y)[, 1]
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  
  list(global = global, local  = as.numeric(local))
}


# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 2000){
  
  if(is.null(y)) y = x
  
  n   = nrow(W)
  IDs = 1:n
  
  xp <- scale(x)[, 1]
  W[which(is.na(W))] <- 0
  
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  
  list(global_sims = global_sims,
       local_sims  = local_sims)
}

#======================================================




##############################################################################
#Now getting into creating our own lisa map:
#For this walk through I am using SVI:


x<- texas_dta$RPL_THEMES #independent variable
y<- texas_dta$avg_CO_percentage #dependent variable

# Calculating the index and its simulated distribution
# for global and local values

m <- moran_I(x, y, W)

# Global Moral
global_moran <- m[[1]][1]
#> 0.2218409

# Local values
m_i <- m[[2]] 

# local simulations
local_sims <- simula_moran(x, y, W)$local_sims


# global pseudo p-value  
# get all simulated global moran
global_sims <- simula_moran(x, y, W)$global_sims

# Proportion of simulated global values that are higher (in absolute terms) than the actual index 
moran_pvalue <- sum(abs(global_sims) > abs( global_moran )) / length(global_sims)
#> 0


# Identifying the significant values 
alpha <- .05  # for a 95% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig       <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )

#======================================================
# Preparing for plotting


# Convert shape file into sf object
map_sf     <- st_as_sf(texas_dta) # we already read our data in as a sf, so this line might be useless
map_sf$sig <- sig


# Identifying the LISA clusters
xp <- scale(x)[,1]
yp <- scale(y)[,1]


patterns <- as.character( interaction(xp > 0, W%*%yp > 0) )
patterns <- patterns %>% 
  str_replace_all("TRUE","High") %>% 
  str_replace_all("FALSE","Low")

patterns[map_sf$sig==0] <- "Not significant"
map_sf$patterns <- patterns


# Rename LISA clusters
map_sf$patterns2 <- factor(map_sf$patterns, levels=c("High.High", "High.Low", "Low.High", "Low.Low", "Not significant"),
                           labels=c("High Vulnerability (CDC SVI) - High Outage Percentage ", "High Vulnerability (CDC SVI) - Low Outage Percentage", "Low Vulnerability (CDC SVI) - High Outage Percentage","Low Vulnerability (CDC SVI) - Low Outage Percentage", "Not significant"))
                           



### PLOT

SVI_out_rate<- ggplot(map_sf) +
  geom_sf(data=map_sf, aes(fill=patterns2), color="NA") +
  scale_fill_manual(values = c("red", "pink", "light blue", "dark blue", "grey80")) + 
  guides(fill = guide_legend(title="LISA Clusters"))+
  labs(title = "SVI and Average Percent Out",
       subtitle = "During Texas Power Outages, Feb 10-24, 2021")+
  
  geom_sf(fill = "transparent", color = "black", size = 1,
          data = map_sf %>% group_by(CountyName) %>% summarize())+  #This is a line I added to get the county boundaries to show up
  
  geom_sf(fill = "transparent", color = "#ffbf1f", size = 1.5,
          data = map_sf%>%dplyr::filter(ercot == "yes")%>% group_by(ercot) %>% dplyr::summarize())+  #this is a line I added to add the ercot boundaries-- I kept it in in case you want to have thicker lines for state boundaries-- this is one way to add that-- ex: group_by(state_name) %>% summarize()
  guides(fill = guide_legend(title="LISA Clusters"))+
  theme_minimal(base_size = 14) +
  
  theme(axis.text.x = element_blank(),
        
        axis.text.y = element_blank(),
        
        axis.ticks = element_blank(),
        
        rect = element_blank()) +
  
  theme(
    
    # Hide panel borders and remove grid lines
    
    panel.grid.major = element_blank(),
    
    panel.grid.minor = element_blank(),
    
  ) 




