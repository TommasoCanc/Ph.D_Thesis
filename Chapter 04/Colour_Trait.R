##################
# Abdomen Colour #
##################

# Load library
library(countcolors) # <- Manage color

# Main path 
main.path <- "YOUR_PATH"

# Connect to Open Science Framework
cr_project <- osf_retrieve_node("4p63t")
osf_files <- osf_ls_files(cr_project, type = "folder")

# Download abdomen images
file.number <- which(osf_ls_files(
  osf_files[2, ], pattern = "Abdomen.zip")[ ,1] == "Abdomen.zip")
osf_download(osf_ls_files(osf_files[2, ], pattern = "Abdomen.zip")
             [file.number, ], path = main.path)
rm(file.number)

# Unzip file
unzip("./Abdomen.zip", exdir = main.path); file.remove("./Abdomen.zip")

# Create folder to save results
dir.create("./Color")
dir.create("./Color/Plot3d")
dir.create("./Color/Histogram")
dir.create("./Color/Pixel_coverage")

# Load abdomen images
img.name.1 <- "anax_ephippiger_01_abdomen.jpg"
img.name.2 <- "anax_ephippiger_02_abdomen.jpg"
img.name.3 <- "anax_ephippiger_03_abdomen.jpg"

# Remove withe color
lower <- c(0.9, 0.9, 0.9)
upper <- c(1, 1, 1)

# Load images
img.1 <- colordistance::loadImage(paste0("./Abdomen/", img.name.1), lower=lower, upper=upper)
img.2 <- colordistance::loadImage(paste0("./Abdomen/", img.name.2), lower=lower, upper=upper)
img.3 <- colordistance::loadImage(paste0("./Abdomen/", img.name.3), lower=lower, upper=upper)

jpeg(paste0("./Color/Plot3d/", img.name.1), 
     width = 3000, height = 1600, quality = 100)
par(mfrow = c(1,3))
colordistance::plotPixels(img.1)
colordistance::plotPixels(img.2)
colordistance::plotPixels(img.3)
dev.off()

jpeg(paste0("./Color/Histogram/", img.name.1), 
     width = 800, height = 800, quality = 100)
par(mfrow = c(1,3))
color.1 <- colordistance::getImageHist(paste0("./Abdomen/", img.name.1), 
                                 bins=c(2, 2, 2), lower=lower, upper=upper)
color.2 <- colordistance::getImageHist(paste0("./Abdomen/", img.name.2), 
                                 bins=c(2, 2, 2), lower=lower, upper=upper)
color.3 <- colordistance::getImageHist(paste0("./Abdomen/", img.name.3), 
                                 bins=c(2, 2, 2), lower=lower, upper=upper)
dev.off()


color.1 <- round(color.1, digits = 3)
color.2 <- round(color.2, digits = 3)
color.3 <- round(color.3, digits = 3)

color.sort.1 <- sort(color.1$Pct, decreasing = T)[1:2]
color.sort.2 <- sort(color.2$Pct, decreasing = T)[1:2]
color.sort.3 <- sort(color.3$Pct, decreasing = T)[1:2]

col.1.1 <- as.numeric(color.1[which(color.sort.1[1] == color.1$Pct), ][,1:3])
col.1.2 <- as.numeric(color.1[which(color.sort.1[2] == color.1$Pct), ][,1:3])

col.2.1 <- as.numeric(color.2[which(color.sort.2[1] == color.2$Pct), ][,1:3])
col.2.2 <- as.numeric(color.2[which(color.sort.2[2] == color.2$Pct), ][,1:3])

col.3.1 <- as.numeric(color.3[which(color.sort.3[1] == color.3$Pct), ][,1:3])
col.3.2 <- as.numeric(color.3[which(color.sort.3[2] == color.3$Pct), ][,1:3])

# Color mean
mean.col.1 <- colMeans(rbind(col.1.1, col.2.1, col.3.1))*255
mean.col.2 <- colMeans(rbind(col.1.2, col.2.2, col.3.2))*255
per.col.1 <- colMeans(rbind(color.sort.1, color.sort.2, color.sort.3))
mean.col.1;mean.col.2;per.col.1

jpeg(paste0("./Color/Pixel_coverage/", img.name.1), 
     width = 1500, height = 800, quality = 100)
par(mfrow = c(1,3))
countcolors::countColors(paste0("./Abdomen/", img.name.1), color.range="spherical", 
                         center = c(col.1.1, col.1.2), radius = c(0.1, 0.1),
                         bg.lower=NULL, bg.upper=NULL, plotting = TRUE,
                         target.color=c("magenta", "cyan"))
countcolors::countColors(paste0("./Abdomen/", img.name.2), color.range="spherical", 
                         center = c(col.2.1, col.2.2), radius = c(0.1, 0.1),
                         bg.lower=NULL, bg.upper=NULL, plotting = TRUE,
                         target.color=c("magenta", "cyan"))
countcolors::countColors(paste0("./Abdomen/", img.name.3), color.range="spherical", 
                         center = c(col.3.1, col.3.2), radius = c(0.1, 0.1),
                         bg.lower=NULL, bg.upper=NULL, plotting = TRUE,
                         target.color=c("magenta", "cyan"))
dev.off()

# Clean Workspace
rm(list = ls())
