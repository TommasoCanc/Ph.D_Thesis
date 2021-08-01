# R scripts are stored in https://github.com/TommasoCanc/Ph.D_Thesis.git

# Load libraries
library(rgdal) # <- GIS
library(raster) # <- GIS
library(dplyr) # <- Data managing
library(tidyr) # <- Data managing
library(data.table) # <- Data managing
library(vegan) # <- Ecology pack
library(ggplot2) # <- Plot
library(ggridges) # <- Plot
library(osfr) # <- Connect to OSF


#####################
# Connection to OSF #
#####################

main.path <- "/Users/tommasocancellario/Desktop/OSF_prove"
# Create folder where you can download the files. 
# For more information: ?dir.create
dir.create(main.path)
setwd(main.path)

# Connect to Open Science Framework
cr_project <- osf_retrieve_node("4p63t")
osf_files <- osf_ls_files(cr_project, type = "folder")


##################
# MZNA MACRO-CHE #
##################
#--------#
# Plot 1 #
#--------#

# Download file Gbif_MZNA_MACRO-CHE.txt ----------------------------------------
# path: Path pointing to a local directory where the downloaded 
# files will be saved. Default is to use the current working directory.
file.number <- which(osf_ls_files(osf_files[1, ])[ ,1] == "Gbif_MZNA_MACRO-CHE.txt")
osf_download(osf_ls_files(osf_files[1, ])[file.number, ], path = main.path); rm(file.number)

# Load the occurrence file. 
taxa <- fread("./Gbif_MZNA_MACRO-CHE.txt")

# Filter taxa object
taxa.1 <- taxa[ ,c("kingdom", "phylum", 'class', 'order', 'family', 'genus', 'organismQuantity')]
head(taxa.1); rm(taxa)

# Convert taxa.1 to dataframe and select only the record with abundance > 0
taxa.1 <- as.data.frame(taxa.1)
taxa.1 <- taxa.1[taxa.1$organismQuantity > 0, ]

# Create an object with the taxonomy
taxonomy <- taxa.1[ ,1:6]
taxonomy <- taxonomy[!duplicated(taxonomy), ]

# Calculate the number of record for each order
taxa.1.df <- as.data.frame(table(taxa.1$order))
colnames(taxa.1.df)[1] <- "order"
taxa.1.df$order <- as.character(taxa.1.df$order)

# Merge the dataframe taxa.1.df with the taxonomy
taxa.1.merge <- merge(taxa.1.df, taxonomy[,c(2,4)], by="order", all = F)
taxa.1.merge <- taxa.1.merge[!duplicated(taxa.1.merge), ]
taxa.1.merge <- taxa.1.merge[taxa.1.merge$order != "", ]

# Plot
ggplot(taxa.1.merge, aes(x= reorder(order, Freq), y= Freq, fill = phylum)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(breaks = seq(0, 9000, 1000)) +
  coord_flip()

# Delete unnecessary objects
rm(taxa.1.df, taxa.1.merge, taxonomy)

#-----------------#
# Basic Statistic #
#-----------------#

# How many phylum?
unique(taxa.1$phylum[which(taxa.1$phylum != "")])

# How many classes?
unique(taxa.1$class[which(taxa.1$class != "")])

# How many orders?
unique(taxa.1$order[which(taxa.1$order != "")])

# How many families?
unique(taxa.1$family[which(taxa.1$family != "")])

# How many genus?
unique(taxa.1$genus[which(taxa.1$genus != "")])

##################################
# Odonata dataset (MZNA-Odonata) #
##################################

#--------#
# Plot 2 #
#--------#

# Download files ---------------------------------------------------------------
file_name <- c("Gbif_Spain.txt", "Gbif_Italy.txt", "Gbif_Germany.txt",
               "Gbif_France.txt", "Gbif_Netherlands.txt", "Gbif_UK.txt")
for(i in file_name){
file.number <- which(osf_ls_files(osf_files[1, ])[ ,1] == i)
osf_download(osf_ls_files(osf_files[1, ])[file.number, ], path = main.path); rm(file.number)
}; rm(i, file_name)

# Load files
spain <- fread("./Gbif_Spain.txt")
italy <- fread("./Gbif_Italy.txt")
germany <- fread("./Gbif_Germany.txt")
france <- fread("./Gbif_France.txt")
nether <- fread("./Gbif_Netherlands.txt")
uk <- fread("/Gbif_UK.txt")

# Filter record after 2000 and remove unused columns
spain.2000 <- spain[spain$year >= 2000, ]
spain.2000 <- spain.2000[ ,c("level0Gid", "level0Name","year")]
italy.2000 <- italy[italy$year >= 2000, ]
italy.2000 <- italy.2000[ ,c("level0Gid", "level0Name","year")]
germany.2000 <- germany[germany$year >= 2000, ]
germany.2000 <- germany.2000[ ,c("level0Gid", "level0Name","year")]
france.2000 <- france[france$year >= 2000, ]
france.2000 <- france.2000[ ,c("level0Gid", "level0Name","year")]
nether.2000 <- nether[nether$year >= 2000, ]
nether.2000 <- nether.2000[ ,c("level0Gid", "level0Name","year")]
uk.2000 <- uk[uk$year >= 2000, ]
uk.2000 <- uk.2000[ ,c("level0Gid", "level0Name","year")]

# Convert to dataframe
spain.2000 <- as.data.frame(spain.2000)
italy.2000 <- as.data.frame(italy.2000)
germany.2000 <- as.data.frame(germany.2000)
france.2000 <- as.data.frame(france.2000)
nether.2000 <- as.data.frame(nether.2000)
uk.2000 <- as.data.frame(uk.2000)

# Delete unnecessary objects
rm(spain, italy, germany, france, nether, uk)

# Create a total dataframe and calculate the number of records per year
total <- rbind(spain.2000, italy.2000, germany.2000, france.2000, nether.2000, uk.2000)
total$count <- 1
total_ends <- total %>% 
  group_by(level0Gid, year) %>% 
  summarise(Frequency = sum(count))

# Plot
ggplot(total_ends, aes(year, Frequency, color = level0Gid)) +
  geom_line() +
  scale_colour_viridis_d()

# Delete unnecessary objects
rm(spain.2000, italy.2000, germany.2000, france.2000, nether.2000, uk.2000, total)

##########################
# MZNA OdoLarvae Dataset #
##########################

# Download files ---------------------------------------------------------------
file_name <- c("Elevation_30_sec.tif", "Ebro_Basin.zip")
for(i in file_name){
  file.number <- which(osf_ls_files(osf_files[1, ])[ ,1] == i)
  osf_download(osf_ls_files(osf_files[1, ])[file.number, ], path = main.path); rm(file.number)
}; rm(i, file_name)

# Load raster and shape files
elev <- raster("./Elevation_30_sec.tif") # <- Elevation
unzip("./Ebro_Basin.zip", exdir=main.path)
basin <- shapefile("./cuenca/Limite_Cuenca_Ebro.shp") # <- Ebro Basin

# Assign the same projection to raster and shapefile
basin <- spTransform(basin, crs(elev))

# Crop elevation raster with EU extension
elev.sub <- crop(elev, extent(basin))
# Mask
elev.sub <- mask(elev.sub, basin); rm(elev, basin)

# Load table containing the Odonata point occurrences
file.number <- which(osf_ls_files(osf_files[1, ])[ ,1] == "odonata_occurrence_sample_postgis.csv")
osf_download(osf_ls_files(osf_files[1, ])[file.number, ], path = main.path); rm(file.number)
odo <- read.csv("./odonata_occurrence_sample_postgis.csv")

# Plot and quickly check
plot(elev.sub)
points(odo$decimallongitude, odo$decimallatitude, pch = 20, cex = .5)

# Retrieve coordinate from Odonata occurrences sample points 
odo.coord <- odo[,c("locality","decimallongitude", "decimallatitude")]
odo.coord <- odo.coord[!duplicated(odo.coord), ]

# Calculate site abundance
odo.abb.tot <-odo.abb %>%
  group_by(locality) %>%
  summarise(abundance = sum(organismquantity))

# Merge site abundance with geographic coordinate
odo.merge <- merge(x = odo.abb.tot, y =  odo.coord)

# Convert raster to dataframe to use in ggplot
elev.sub.df <- raster::as.data.frame(elev.sub, xy=TRUE)
colnames(elev.sub.df)[3] <- "elevation"
elev.sub.df <- elev.sub.df[!is.na(elev.sub.df$elevation),]

# Plot
ggplot() +
  geom_raster(data=elev.sub.df, aes(x = x, y = y), fill = "lightgray") +
  geom_point(data=odo.merge, aes(x=decimallongitude, y=decimallatitude, size=abundance, color=abundance)) +
  theme_void() 

#----------#
# Richness #
#----------#

# Abundance site : genus
odo.abb <- odo.abb %>%
  group_by(locality, genericname) %>%
  summarise(abundance = sum(organismquantity))

# Transformation from long format to wide format
odo.abb.wide <- spread(odo.abb, genericname, abundance)
odo.abb.wide[is.na(odo.abb.wide)] <- 0
odo.abb.wide <- as.data.frame(odo.abb.wide)
rownames(odo.abb.wide) <- odo.abb.wide[ ,1]
odo.abb.wide <- odo.abb.wide[ ,-1]

# Observed family richness
richness <- as.data.frame(specnumber(odo.abb.wide)) # <- Calcualte family richness with vegan pack 
richness$locality <- rownames(richness)
rownames(richness) <- 1:nrow(richness)

# Add coordinates to richness object
odo.coord <- merge(x = richness, y = odo.coord, by = "locality")
colnames(odo.coord)[2] <- "richness"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Does a Correlation Exist between family richness and altitude? #
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# Extract odonata sample point altitude
odo.coord.spatial <- odo.coord
coordinates(odo.coord.spatial) = ~ decimallongitude + decimallatitude
odo.alt <- raster::extract(elev.sub, odo.coord.spatial)

odo.coord$altitude <- odo.alt; rm(odo.coord.spatial, odo.alt)

# To check the possible correlation we add the table with every sample points.
file.number <- which(osf_ls_files(osf_files[1, ])[ ,1] == "odonata_occurrence_sample_postgis.csv")
osf_download(osf_ls_files(osf_files[1, ])[file.number, ], path = main.path); rm(file.number)
all.sample <- read.csv("./sample_points.csv")

# We filter only the sample points without Odonata
'%ni%' <- Negate("%in%") # <- Create th efunction %ni%
no.odo <- all.sample[all.sample$locality %ni% odo.coord$locality, ]; rm(all.sample)

# Extract altitude for sample points without Odonata
no.odo.spatial <- no.odo
coordinates(no.odo.spatial) = ~ decimallongitude + decimallatitude
no.odo.alt <- raster::extract(elev.sub, no.odo.spatial)

# Create a dataframe
no.odo$altitude <- no.odo.alt; rm(no.odo.spatial, no.odo.alt)
no.odo <- data.frame(locality = no.odo$locality,
                     richness = 0,
                     decimallongitude = no.odo$decimallongitude,
                     decimallatitude = no.odo$decimallatitude,
                     altitude = no.odo$altitude)

# Merge dataset
odo.coord <- odo.coord[ ,c("locality", "richness", "decimallongitude",
                           "decimallatitude", "altitude")]
odo.coord <- rbind(odo.coord, no.odo)

# Linear model
summary(lm(richness ~ altitude, odo.coord)) # <- 0.0424 * 

# Density plot family richness : altitude
ggplot(odo.coord, aes(x = altitude, y = as.factor(richness))) +
  geom_density_ridges(aes(fill = as.factor(richness)), alpha = 0.8) +
  scale_fill_viridis_d(option = "viridis", begin = 0.2, end = 0.8) +
  theme_bw()

