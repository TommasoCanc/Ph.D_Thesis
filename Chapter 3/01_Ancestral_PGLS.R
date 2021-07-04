###########################
# Traits and phylogenesis #
###########################

# Load libraries
library(ape) # <- Genetic
library(splits) # <- Genetic
library(caper) # <- Genetic
library(phytools) # <- Genetic
library(MCMCglmm) # <- Multi-response Generalized Linear Mixed Models
library(xlsx) # <- Open .xlsx files
library(viridis) # <- Manage color
library(osfr) # <- Connect to OSF

#####################
# Connection to OSF #
#####################

main.path <- "/Users/tommasocancellario/Desktop/OSF"
# Create folder where you can download the files. 
# For more information: ?dir.create
setwd(main.path)

# Connect to Open Science Framework
cr_project <- osf_retrieve_node("4p63t")
osf_files <- osf_ls_files(cr_project, type = "folder")

# Download main variables
# Download phylogenetic tree 
file.number <- which(osf_ls_files(osf_files[3, ])[ ,1] == "tree_treeannotator.nex")
osf_download(osf_ls_files(osf_files[3, ])[file.number, ], path = main.path)
rm(file.number)
# Download .xlsx with traits
file.number <- which(osf_ls_files(osf_files[3, ])[ ,1] == "Odonata_traits_SDM.xlsx")
osf_download(osf_ls_files(osf_files[3, ])[file.number, ], path = main.path)
rm(file.number)
# Download Bcc_2050 and Bcc_2070 
file.number <- which(osf_ls_files(osf_files[3, ], pattern = "Bcc_2050")[ ,1] == "Bcc_2050.csv")
osf_download(osf_ls_files(osf_files[3, ], pattern = "Bcc_2050")[file.number, ], path = main.path)
rm(file.number)
file.number <- which(osf_ls_files(osf_files[3, ], pattern = "Bcc_2070")[ ,1] == "Bcc_2070.csv")
osf_download(osf_ls_files(osf_files[3, ], pattern = "Bcc_2070")[file.number, ], path = main.path)
rm(file.number)
# Download Miroc_2050 and Miroc_2070 
file.number <- which(osf_ls_files(osf_files[3, ], pattern = "Miroc_2050")[ ,1] == "Miroc_2050.csv")
osf_download(osf_ls_files(osf_files[3, ], pattern = "Miroc_2050")[file.number, ], path = main.path)
rm(file.number)
file.number <- which(osf_ls_files(osf_files[3, ], pattern = "Miroc_2070")[ ,1] == "Miroc_2070.csv")
osf_download(osf_ls_files(osf_files[3, ], pattern = "Miroc_2070")[file.number, ], path = main.path)
rm(file.number)
# Download Nor_2050 and Nor_2070 
file.number <- which(osf_ls_files(osf_files[3, ], pattern = "Nor_2050")[ ,1] == "Nor_2050.csv")
osf_download(osf_ls_files(osf_files[3, ], pattern = "Nor_2050")[file.number, ], path = main.path)
rm(file.number)
file.number <- which(osf_ls_files(osf_files[3, ], pattern = "Nor_2070")[ ,1] == "Nor_2070.csv")
osf_download(osf_ls_files(osf_files[3, ], pattern = "Nor_2070")[file.number, ], path = main.path)
rm(file.number)

# Functions
'%ni%' <- Negate('%in%')

# Create folder to save results
dir.create("./Ancestor_recostriction")


# Load phylogenetic tree
fict.phylo <- read.nexus("./tree_treeannotator.nex")
#inspect the resulting object:
str(fict.phylo)
fict.phylo$tip.label

# Load the traits table
traits <- read.xlsx(file="./Odonata_traits_SDM.xlsx", sheetIndex = 1)
head(traits)

# Load model tables
miroc2050 <- read.csv("./Miroc_2050.csv")
colnames(miroc2050)[10] <- "Species"
miroc2070 <- read.csv("./Miroc_2070.csv")
colnames(miroc2070)[10] <- "Species"
nor2050 <- read.csv("./Nor_2050.csv")
colnames(nor2050)[10] <- "Species"
nor2070 <- read.csv("./Miroc_2070.csv")
colnames(nor2070)[10] <- "Species"
bcc2050 <- read.csv("./Bcc_2050.csv")
colnames(bcc2050)[10] <- "Species"
bcc2070 <- read.csv("./Nor_2070.csv")
colnames(bcc2070)[10] <- "Species"

# Merge traits and variables
miroc2050.mg <- merge(traits, miroc2050, by = "Species", all.x = TRUE )
miroc2070.mg <- merge(traits, miroc2070, by = "Species", all.x = TRUE )
nor2050.mg <- merge(traits, nor2050, by = "Species", all.x = TRUE )
nor2070.mg <- merge(traits, nor2070, by = "Species", all.x = TRUE )
bcc2050.mg <- merge(traits, bcc2050, by = "Species", all.x = TRUE )
bcc2070.mg <- merge(traits, bcc2070, by = "Species", all.x = TRUE )


# Remove outliers
outliers <- c("Trithemis arteriosa", "Sympecma paedisca", 
              "Coenagrion lunulatum", "Calopteryx haemorrhoidalis")
miroc2050.noOut <- miroc2050.mg[miroc2050.mg$Species %ni% outliers, ]
miroc2070.noOut <- miroc2070.mg[miroc2070.mg$Species %ni% outliers, ]
nor2050.noOut <- nor2050.mg[nor2050.mg$Species %ni% outliers, ]
nor2070.noOut <- nor2070.mg[nor2070.mg$Species %ni% outliers, ]
bcc2050.noOut <- bcc2050.mg[bcc2050.mg$Species %ni% outliers, ]
bcc2070.noOut <- bcc2070.mg[bcc2070.mg$Species %ni% outliers, ]

rm(miroc2050.mg, miroc2070.mg, nor2050.mg, nor2070.mg, bcc2050.mg, bcc2070.mg,
   miroc2050, miroc2070, nor2050, nor2070, bcc2050, bcc2070, traits, outliers)

# Drop tip
tree <- drop.tip(fict.phylo, "Trithemis_arteriosa")
tree <- drop.tip(tree, "Sympecma_paedisca")
tree <- drop.tip(tree, "Coenagrion_lunulatum")
tree <- drop.tip(tree, "Calopteryx_haemorrhoidalis")
tree$tip.label # 101 labels

plot(tree, cex=.3)

#-----------------------------------#
# Ancestor character reconstruction #
#-----------------------------------#

# Latitude difference MIROC ----------------------------------------------------
lat_miroc2050 <- miroc2050.noOut$LatDif
sp <- gsub(" ", "_", miroc2050.noOut$Species)
names(lat_miroc2050) <- sp
lat_miroc2050_1 <- lat_miroc2050[tree$tip.label]

tree_miroc_2050 <- contMap(tree, lat_miroc2050_1, fsize = 1, scale = 1, lwd = 6, 
                           tip.labels = TRUE, outline = TRUE, plot = FALSE)

lat_miroc2070 <- miroc2070.noOut$LatDif
names(lat_miroc2070) <- sp
lat_miroc2070_1 <- lat_miroc2070[tree$tip.label]

tree_miroc_2070 <- contMap(tree, lat_miroc2070_1, fsize = 1, scale = 1, 
                           lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

# Save Plot
pdf("./Ancestor_recostriction/Latitude_miroc.pdf")
plot(tree_miroc_2050)
title("Latitude Miroc_2050")
plot(tree_miroc_2070)
title("Latitude Miroc_2070")
dev.off()

rm(tree_miroc_2050, tree_miroc_2070, lat_miroc2050, lat_miroc2050_1,
   lat_miroc2070, lat_miroc2070_1)

# Altitude difference MIROC
alt_miroc2050 <- miroc2050.noOut$DifavgAltitude
names(alt_miroc2050) <- sp
alt_miroc2050_1 <- alt_miroc2050[tree$tip.label]

tree_miroc_2050 <-contMap(tree, alt_miroc2050_1, fsize = 1, scale = 1, 
                          lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

alt_miroc2070 <- miroc2070.noOut$DifavgAltitude
names(alt_miroc2070) <- sp
alt_miroc2070_1 <- alt_miroc2070[tree$tip.label]

tree_miroc_2070 <-contMap(tree, alt_miroc2070_1, fsize = 1, scale = 1, 
                          lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

# Save Plot
pdf("./Ancestor_recostriction/Altitude_miroc.pdf")
plot(tree_miroc_2050)
title("Altitude Miroc_2050")
plot(tree_miroc_2070)
title("Altitude Miroc_2070")
dev.off()

rm(tree_miroc_2050, tree_miroc_2070, alt_miroc2050, alt_miroc2050_1,
   alt_miroc2070, alt_miroc2070_1)

# Area difference MIROC
area_miroc2050 <- miroc2050.noOut$cellDifArea
names(area_miroc2050) <- sp
area_miroc2050_1 <- area_miroc2050[tree$tip.label]

tree_miroc_2050 <-contMap(tree, area_miroc2050_1, fsize = 1, scale = 1, 
                          lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

area_miroc2070 <- miroc2070.noOut$cellDifArea
names(area_miroc2070) <- sp
area_miroc2070_1 <- area_miroc2070[tree$tip.label]

tree_miroc_2070 <-contMap(tree, area_miroc2070_1, fsize = 1, scale = 1, 
                          lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

# Save Plot
pdf("./Ancestor_recostriction/Area_miroc.pdf")
plot(tree_miroc_2050)
title("Area Miroc_2050")
plot(tree_miroc_2070)
title("Area Miroc_2070")
dev.off()

rm(tree_miroc_2050, tree_miroc_2070, area_miroc2050, area_miroc2050_1,
   area_miroc2070, area_miroc2070_1)

# Latitude difference NOR ------------------------------------------------------
lat_nor2050 <- nor2050.noOut$LatDif
sp <- gsub(" ", "_", nor2050.noOut$Species)
names(lat_nor2050) <- sp
lat_nor2050_1 <- lat_nor2050[tree$tip.label]

tree_nor_2050 <-contMap(tree, lat_nor2050_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

lat_nor2070 <- nor2070.noOut$LatDif
names(lat_nor2070) <- sp
lat_nor2070_1 <- lat_nor2070[tree$tip.label]

tree_nor_2070 <-contMap(tree, lat_nor2070_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

# Save Plot
pdf("./Ancestor_recostriction/Latitude_nor.pdf")
plot(tree_nor_2050)
title("Latitude Nor_2050")
plot(tree_nor_2070)
title("Latitude Nor_2070")
dev.off()

rm(tree_nor_2050, tree_nor_2070, lat_nor2050, lat_nor2050_1,
   lat_nor2070, lat_nor2070_1)

# Altitude difference NOR
alt_nor2050 <- nor2050.noOut$DifavgAltitude
names(alt_nor2050) <- sp
alt_nor2050_1 <- alt_nor2050[tree$tip.label]

tree_nor_2050 <-contMap(tree, alt_nor2050_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

alt_nor2070 <- nor2070.noOut$DifavgAltitude
names(alt_nor2070) <- sp
alt_nor2070_1 <- alt_nor2070[tree$tip.label]

tree_nor_2070 <-contMap(tree, alt_nor2070_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

# Save Plot
pdf("./Ancestor_recostriction/Altitude_nor.pdf")
plot(tree_nor_2050)
title("Altitude Nor_2050")
plot(tree_nor_2070)
title("Altitude Nor_2070")
dev.off()

rm(tree_nor_2050, tree_nor_2070, alt_nor2050, alt_nor2050_1,
   alt_nor2070, alt_nor2070_1)

# Area difference NOR
area_nor2050 <- nor2050.noOut$cellDifArea
names(area_nor2050) <- sp
area_nor2050_1 <- area_nor2050[tree$tip.label]

tree_nor_2050 <-contMap(tree, area_nor2050_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

area_nor2070 <- nor2070.noOut$cellDifArea
names(area_nor2070) <- sp
area_nor2070_1 <- area_nor2070[tree$tip.label]

tree_nor_2070 <-contMap(tree, area_nor2070_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

# Save Plot
pdf("./Ancestor_recostriction/Area_nor.pdf")
plot(tree_nor_2050)
title("Area Nor_2050")
plot(tree_nor_2070)
title("Area Nor_2070")
dev.off()

rm(tree_nor_2050, tree_nor_2070, area_nor2050, area_nor2050_1,
   area_nor2070, area_nor2070_1)

# Latitude difference BCC ------------------------------------------------------
lat_bcc2050 <- bcc2050.noOut$LatDif
sp <- gsub(" ", "_", bcc2050.noOut$Species)
names(lat_bcc2050) <- sp
lat_bcc2050_1 <- lat_bcc2050[tree$tip.label]

tree_bcc_2050 <-contMap(tree, lat_bcc2050_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

lat_bcc2070 <- bcc2070.noOut$LatDif
names(lat_bcc2070) <- sp
lat_bcc2070_1 <- lat_bcc2070[tree$tip.label]

tree_bcc_2070 <-contMap(tree, lat_bcc2070_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

# Save Plot
pdf("./Ancestor_recostriction/Latitude_bcc.pdf")
plot(tree_bcc_2050)
title("Latitude Bcc_2050")
plot(tree_bcc_2070)
title("Latitude Bcc_2070")
dev.off()

rm(tree_bcc_2050, tree_bcc_2070, lat_bcc2050, lat_bcc2050_1,
   lat_bcc2070, lat_bcc2070_1)

# Altitude difference BCC
alt_bcc2050 <- bcc2050.noOut$DifavgAltitude
names(alt_bcc2050) <- sp
alt_bcc2050_1 <- alt_bcc2050[tree$tip.label]

tree_bcc_2050 <-contMap(tree, alt_bcc2050_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

alt_bcc2070 <- bcc2070.noOut$DifavgAltitude
names(alt_bcc2070) <- sp
alt_bcc2070_1 <- alt_bcc2070[tree$tip.label]

tree_bcc_2070 <-contMap(tree, alt_bcc2070_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

# Save Plot
pdf("./Ancestor_recostriction/Altitude_bcc.pdf")
plot(tree_bcc_2050)
title("Altitude Bcc_2050")
plot(tree_bcc_2070)
title("Altitude Bcc_2070")
dev.off()

rm(tree_bcc_2050, tree_bcc_2070, 
   alt_bcc2050, alt_bcc2050_1,
   alt_bcc2070, alt_bcc2070_1)

# Area difference NOR
area_bcc2050 <- bcc2050.noOut$cellDifArea
names(area_bcc2050) <- sp
area_bcc2050_1 <- area_bcc2050[tree$tip.label]

tree_bcc_2050 <-contMap(tree, area_bcc2050_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

area_bcc2070 <- bcc2070.noOut$cellDifArea
names(area_bcc2070) <- sp
area_bcc2070_1 <- area_bcc2070[tree$tip.label]

tree_bcc_2070 <-contMap(tree, area_bcc2070_1, fsize = 1, scale = 1, 
                        lwd = 6, tip.labels = TRUE, outline = TRUE, plot = FALSE)

# Save Plot
pdf("./Ancestor_recostriction/Area_bcc.pdf")
plot(tree_bcc_2050)
title("Area Bcc_2050")
plot(tree_bcc_2070)
title("Area Bcc_2070")
dev.off()

rm(tree_bcc_2050, tree_bcc_2070, area_bcc2050, area_bcc2050_1,
   area_bcc2070, area_bcc2070_1)

# Habitat_2
habitat2 <- miroc2050.noOut$Habitat_2
names(habitat2)<- sp
habitat2.1 <- habitat2[tree$tip.label]

est.habitat <- ace(habitat2.1, tree, type = "discrete", 
                   method = "ML", CI = T, marginal = TRUE)

pdf("./Ancestor_recostriction/Habitat.pdf")
plot.phylo(tree, type="phylogram", use.edge.length = TRUE, 
           label.offset = 0.1, cex = 0.5)
co <- c("yellow", "blue")
tiplabels(pch = 22, bg = co[as.factor(habitat2.1)], cex = 0.9, adj = 0.55)
nodelabels(pie = est.habitat$lik.anc, piecol = co, cex = 0.2)
title("Habitats")
dev.off()

rm(co, habitat2, habitat2.1)


#---------------------#
# Phylogenetic signal #
#---------------------#

# MIROC
miroc2050_k_lat <-  phylosig(tree, miroc2050.noOut$LatDif, "K", test = TRUE)
miroc2050_lam_lat <- phylosig(tree, miroc2050.noOut$LatDif, "lambda", test = TRUE)
miroc2070_k_lat <- phylosig(tree, miroc2070.noOut$LatDif, "K", test = TRUE)
miroc2070_lam_lat <- phylosig(tree, miroc2070.noOut$LatDif, "lambda", test = TRUE)

miroc2050_k_alt <- phylosig(tree, miroc2050.noOut$DifavgAltitude, "K", test = TRUE)
miroc2050_lam_alt <- phylosig(tree, miroc2050.noOut$DifavgAltitude, "lambda", test = TRUE)
miroc2070_k_alt <- phylosig(tree, miroc2070.noOut$DifavgAltitude, "K", test = TRUE)
miroc2070_lam_alt <- phylosig(tree, miroc2070.noOut$DifavgAltitude, "lambda", test = TRUE)

miroc2050_k_area <- phylosig(tree, miroc2050.noOut$cellDifArea, "K", test = TRUE)
miroc2050_lam_area <- phylosig(tree, miroc2050.noOut$cellDifArea, "lambda", test = TRUE)
miroc2070_k_area <- phylosig(tree, miroc2070.noOut$cellDifArea, "K", test = TRUE)
miroc2070_lam_area <- phylosig(tree, miroc2070.noOut$cellDifArea, "lambda", test = TRUE)

phylosig.k <- c(miroc2050_k_lat$K, miroc2070_k_lat$K, 
                miroc2050_k_alt$K, miroc2070_k_alt$K, 
                miroc2050_k_area$K,  miroc2070_k_area$K)

phylosig.lam <- c(miroc2050_lam_lat$lambda, miroc2070_lam_lat$lambda,
                  miroc2050_lam_alt$lambda, miroc2070_lam_alt$lambda,
                  miroc2050_lam_area$lambda, miroc2070_lam_area$lambda)

phylo_miroc <- data.frame(phylosig = c(phylosig.k, phylosig.lam),
                          test = rep(c("k", "lambda"),  each =6),
                          variable = rep(c("lat", "alt", "area"), 
                                         each = 2, times = 2),
                          scenario = rep(c("miroc_2050", "miroc_2070"), 6))

rm(miroc2050_k_lat, miroc2070_k_lat, miroc2050_k_alt, miroc2070_k_alt, 
   miroc2050_k_area,  miroc2070_k_area, miroc2050_lam_lat, miroc2070_lam_lat,
   miroc2050_lam_alt, miroc2070_lam_alt, miroc2050_lam_area, miroc2070_lam_area,
   phylosig.k, phylosig.lam)

# NOR
nor2050_k_lat <- phylosig(tree, nor2050.noOut$LatDif, "K", test = TRUE)
nor2050_lam_lat <- phylosig(tree, nor2050.noOut$LatDif, "lambda", test = TRUE)
nor2070_k_lat <- phylosig(tree, nor2070.noOut$LatDif, "K", test = TRUE)
nor2070_lam_lat <- phylosig(tree, nor2070.noOut$LatDif, "lambda", test = TRUE)

nor2050_k_alt <- phylosig(tree, nor2050.noOut$DifavgAltitude, "K", test = TRUE)
nor2050_lam_alt <- phylosig(tree, nor2050.noOut$DifavgAltitude, "lambda", test = TRUE)
nor2070_k_alt <- phylosig(tree, nor2070.noOut$DifavgAltitude, "K", test = TRUE)
nor2070_lam_alt <- phylosig(tree, nor2070.noOut$DifavgAltitude, "lambda", test = TRUE)

nor2050_k_area <- phylosig(tree, nor2050.noOut$cellDifArea, "K", test = TRUE)
nor2050_lam_area <- phylosig(tree, nor2050.noOut$cellDifArea, "lambda", test = TRUE)
nor2070_k_area <- phylosig(tree, nor2070.noOut$cellDifArea, "K", test = TRUE)
nor2070_lam_area <- phylosig(tree, nor2070.noOut$cellDifArea, "lambda", test = TRUE)

phylosig.k <- c(nor2050_k_lat$K, nor2070_k_lat$K, 
                nor2050_k_alt$K, nor2070_k_alt$K, 
                nor2050_k_area$K,  nor2070_k_area$K)

phylosig.lam <- c(nor2050_lam_lat$lambda, nor2070_lam_lat$lambda,
                  nor2050_lam_alt$lambda, nor2070_lam_alt$lambda,
                  nor2050_lam_area$lambda, nor2070_lam_area$lambda)

phylo_nor <- data.frame(phylosig = c(phylosig.k, phylosig.lam),
                          test = rep(c("k", "lambda"),  each =6),
                          variable = rep(c("lat", "alt", "area"), 
                                         each = 2, times = 2),
                          scenario = rep(c("nor_2050", "nor_2070"), 6))

rm(nor2050_k_lat, nor2070_k_lat, nor2050_k_alt, nor2070_k_alt, 
   nor2050_k_area,  nor2070_k_area, nor2050_lam_lat, nor2070_lam_lat,
   nor2050_lam_alt, nor2070_lam_alt, nor2050_lam_area, nor2070_lam_area,
   phylosig.k, phylosig.lam)

# BCC
bcc2050_k_lat <- phylosig(tree, bcc2050.noOut$LatDif, "K", test = TRUE)
bcc2050_lam_lat <- phylosig(tree, bcc2050.noOut$LatDif, "lambda", test = TRUE)
bcc2070_k_lat <- phylosig(tree, bcc2070.noOut$LatDif, "K", test = TRUE)
bcc2070_lam_lat <- phylosig(tree, bcc2070.noOut$LatDif, "lambda", test = TRUE)

bcc2050_k_alt <- phylosig(tree, bcc2050.noOut$DifavgAltitude, "K", test = TRUE)
bcc2050_lam_alt <- phylosig(tree, bcc2050.noOut$DifavgAltitude, "lambda", test = TRUE)
bcc2070_k_alt <- phylosig(tree, bcc2070.noOut$DifavgAltitude, "K", test = TRUE)
bcc2070_lam_alt <- phylosig(tree, bcc2070.noOut$DifavgAltitude, "lambda", test = TRUE)

bcc2050_k_area <- phylosig(tree, bcc2050.noOut$cellDifArea, "K", test = TRUE)
bcc2050_lam_area <- phylosig(tree, bcc2050.noOut$cellDifArea, "lambda", test = TRUE)
bcc2070_k_area <- phylosig(tree, bcc2070.noOut$cellDifArea, "K", test = TRUE)
bcc2070_lam_area <- phylosig(tree, bcc2070.noOut$cellDifArea, "lambda", test = TRUE)

phylosig.k <- c(bcc2050_k_lat$K, bcc2070_k_lat$K, 
                bcc2050_k_alt$K, bcc2070_k_alt$K, 
                bcc2050_k_area$K,  bcc2070_k_area$K)

phylosig.lam <- c(bcc2050_lam_lat$lambda, bcc2070_lam_lat$lambda,
                  bcc2050_lam_alt$lambda, bcc2070_lam_alt$lambda,
                  bcc2050_lam_area$lambda, bcc2070_lam_area$lambda)

phylo_bcc <- data.frame(phylosig = c(phylosig.k, phylosig.lam),
                        test = rep(c("k", "lambda"),  each =6),
                        variable = rep(c("lat", "alt", "area"), 
                                       each = 2, times = 2),
                        scenario = rep(c("bcc_2050", "bcc_2070"), 6))

rm(bcc2050_k_lat, bcc2070_k_lat, bcc2050_k_alt, bcc2070_k_alt, 
   bcc2050_k_area,  bcc2070_k_area, bcc2050_lam_lat, bcc2070_lam_lat,
   bcc2050_lam_alt, bcc2070_lam_alt, bcc2050_lam_area, bcc2070_lam_area,
   phylosig.k, phylosig.lam)

# Join phylo signal
phylo_signal <- rbind(phylo_miroc, phylo_nor, phylo_bcc)

write.csv(phylo_signal, "./Phylogenetic_signal.csv", row.names = FALSE)

rm(phylo_signal, phylo_miroc, phylo_nor, phylo_bcc)

#------#
# PGLS #
#------#

# Check variables correlation
str(miroc2050.noOut)
psych::pairs.panels(miroc2050.noOut[,c("TotMax..mm.", "AbMax..mm.", "HwMax..mm.", 
                                       "Habitat_2","FlightSeason..month.")])

# PGLS MIROC 2050
str(miroc2050.noOut)
miroc2050.noOut$Habitat_2 <- as.factor(miroc2050.noOut$Habitat_2)
miroc2050.noOut$flight_mode <- as.factor(miroc2050.noOut$flight_mode)

miroc2050.noOut$Species <- gsub(" ", "_", miroc2050.noOut$Species)
test.data  <-  comparative.data(phy = tree, data = miroc2050.noOut, 
                                names.col = Species, vcv = TRUE, vcv.dim = 3)
# PGLS model can be fitted:
odo_lat <- pgls(LatDif ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_lat$phyres)
qqnorm(odo_lat$phyres)
qqline(odo_lat$phyres)
plot(x = fitted(odo_lat), y = odo_lat$phyres, pch = 19)
miroc2050.latitude <- summary(odo_lat)
miroc2050.latitude <- as.data.frame(miroc2050.latitude$coefficients)
miroc2050.latitude$y <- c("Latitude_MIROC_2050")

odo_alt <- pgls(DifavgAltitude ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
               data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_alt$phyres)
qqnorm(odo_alt$phyres)
qqline(odo_alt$phyres)
plot(x = fitted(odo_alt), y = odo_alt$phyres, pch = 19)
miroc2050.altitude <- summary(odo_alt)
miroc2050.altitude <- as.data.frame(miroc2050.altitude$coefficients)
miroc2050.altitude$y <- c("Altitude_MIROC_2050")

odo_area <- pgls(cellDifArea ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                 data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_area$phyres)
qqnorm(odo_area$phyres)
qqline(odo_area$phyres)
plot(x = fitted(odo_area), y = odo_area$phyres, pch = 19)
miroc2050.area <- summary(odo_area)
miroc2050.area <- as.data.frame(miroc2050.area$coefficients)
miroc2050.area$y <- c("Area_MIROC_2050")

# PGLS table MIROC 2050
pgls_miroc_2050 <- rbind(miroc2050.latitude, miroc2050.altitude, miroc2050.area)
rm(miroc2050.latitude, miroc2050.altitude, miroc2050.area,
   odo_lat, odo_alt, odo_area)


# PGLS MIROC 2070
str(miroc2070.noOut)
miroc2070.noOut$Habitat_2 <- as.factor(miroc2070.noOut$Habitat_2)
miroc2070.noOut$flight_mode <- as.factor(miroc2070.noOut$flight_mode)

miroc2070.noOut$Species <- gsub(" ","_",miroc2070.noOut$Species)
test.data= comparative.data(phy = tree, data = miroc2070.noOut, 
                            names.col = Species, vcv = TRUE, vcv.dim = 3)

# PGLS model can be fitted:
odo_lat <- pgls(LatDif ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_lat$phyres)
qqnorm(odo_lat$phyres)
qqline(odo_lat$phyres)
plot(x = fitted(odo_lat), y = odo_lat$phyres, pch = 19)
miroc2070.latitude <- summary(odo_lat)
miroc2070.latitude <- as.data.frame(miroc2070.latitude$coefficients)
miroc2070.latitude$y <- c("Latitude_MIROC_2070")

odo_alt <-  pgls(DifavgAltitude ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                 data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_alt$phyres)
qqnorm(odo_alt$phyres)
qqline(odo_alt$phyres)
plot(x = fitted(odo_alt), y = odo_alt$phyres, pch = 19)
miroc2070.altitude <- summary(odo_alt)
miroc2070.altitude <- as.data.frame(miroc2070.altitude$coefficients)
miroc2070.altitude$y <- c("Altitude_MIROC_2070")

odo_area <-  pgls(cellDifArea ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                  data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_area$phyres)
qqnorm(odo_area$phyres)
qqline(odo_area$phyres)
plot(x = fitted(odo_area), y = odo_area$phyres, pch = 19)
miroc2070.area <- summary(odo_area)
miroc2070.area <- as.data.frame(miroc2070.area$coefficients)
miroc2070.area$y <- c("Area_MIROC_2070")

# PGLS table MIROC 2070
pgls_miroc_2070 <- rbind(miroc2070.latitude, miroc2070.altitude, miroc2070.area)
rm(miroc2070.latitude, miroc2070.altitude, miroc2070.area,
   odo_lat, odo_alt, odo_area)

# PGLS NOR 2050
str(nor2050.noOut)
nor2050.noOut$Habitat_2 <- as.factor(nor2050.noOut$Habitat_2)
nor2050.noOut$flight_mode <- as.factor(nor2050.noOut$flight_mode)

nor2050.noOut$Species <- gsub(" ","_",nor2050.noOut$Species)
test.data <- comparative.data(phy = tree, data = nor2050.noOut, 
                              names.col = Species, vcv = TRUE, vcv.dim = 3)
# PGLS model can be fitted:
odo_lat <-  pgls(LatDif ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                 data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_lat$phyres)
qqnorm(odo_lat$phyres)
qqline(odo_lat$phyres)
plot(x = fitted(odo_lat), y = odo_lat$phyres, pch = 19)
nor2050.latitude <- summary(odo_lat)
nor2050.latitude <- as.data.frame(nor2050.latitude$coefficients)
nor2050.latitude$y <- c("Latitude_NOR_2050")

odo_alt <-  pgls(DifavgAltitude ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                 data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_alt$phyres)
qqnorm(odo_alt$phyres)
qqline(odo_alt$phyres)
plot(x = fitted(odo_alt), y = odo_alt$phyres, pch = 19)
nor2050.altitude <- summary(odo_alt)
nor2050.altitude <- as.data.frame(nor2050.altitude$coefficients)
nor2050.altitude$y <- c("Altitude_NOR_2050")

odo_area <- pgls(cellDifArea ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                 data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_area$phyres)
qqnorm(odo_area$phyres)
qqline(odo_area$phyres)
plot(x = fitted(odo_area), y = odo_area$phyres, pch = 19)
nor2050.area <- summary(odo_area)
nor2050.area <- as.data.frame(nor2050.area$coefficients)
nor2050.area$y <- c("Area_NOR_2050")

# PGLS table NOR 2050
pgls_nor_2050 <- rbind(nor2050.latitude, nor2050.altitude, nor2050.area)
rm(nor2050.latitude, nor2050.altitude, nor2050.area,
   odo_lat, odo_alt, odo_area)

# PGLS NOR 2070
str(nor2070.noOut)
nor2070.noOut$Habitat_2 <- as.factor(nor2070.noOut$Habitat_2)
nor2070.noOut$flight_mode <- as.factor(nor2070.noOut$flight_mode)

nor2070.noOut$Species <- gsub(" ","_",nor2070.noOut$Species)
test.data <-  comparative.data(phy = tree, data = nor2070.noOut, 
                               names.col = Species, vcv = TRUE, vcv.dim = 3)
# PGLS model can be fitted:
odo_lat <- pgls(LatDif ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_lat$phyres)
qqnorm(odo_lat$phyres)
qqline(odo_lat$phyres)
plot(x = fitted(odo_lat), y = odo_lat$phyres, pch = 19)
nor2070.latitude <- summary(odo_lat)
nor2070.latitude <- as.data.frame(nor2070.latitude$coefficients)
nor2070.latitude$y <- c("Latitude_NOR_2070")

odo_alt <- pgls(DifavgAltitude ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_alt$phyres)
qqnorm(odo_alt$phyres)
qqline(odo_alt$phyres)
plot(x = fitted(odo_alt), y = odo_alt$phyres, pch = 19)
nor2070.altitude <- summary(odo_alt)
nor2070.altitude <- as.data.frame(nor2070.altitude$coefficients)
nor2070.altitude$y <- c("Altitude_NOR_2070")

odo_area <- pgls(cellDifArea ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_area$phyres)
qqnorm(odo_area$phyres)
qqline(odo_area$phyres)
plot(x = fitted(odo_area), y = odo_area$phyres, pch = 19)
nor2070.area <- summary(odo_area)
nor2070.area <- as.data.frame(nor2070.area$coefficients)
nor2070.area$y <- c("Area_NOR_2070")

# PGLS table NOR 2070
pgls_nor_2070 <- rbind(nor2070.latitude, nor2070.altitude, nor2070.area)
rm(nor2070.latitude, nor2070.altitude, nor2070.area,
   odo_lat, odo_alt, odo_area)

# PGLS BCC 2050
str(bcc2050.noOut)
bcc2050.noOut$Habitat_2 <- as.factor(bcc2050.noOut$Habitat_2)
bcc2050.noOut$flight_mode <- as.factor(bcc2050.noOut$flight_mode)

bcc2050.noOut$Species <- gsub(" ","_",bcc2050.noOut$Species)
test.data <- comparative.data(phy = tree, data = bcc2050.noOut, 
                              names.col = Species, vcv = TRUE, vcv.dim = 3)

#now the PGLS model can be fitted:
odo_lat <-  pgls(LatDif ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                 data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_lat$phyres)
qqnorm(odo_lat$phyres)
qqline(odo_lat$phyres)
plot(x = fitted(odo_lat), y = odo_lat$phyres, pch = 19)
bcc2050.latitude <- summary(odo_lat)
bcc2050.latitude <- as.data.frame(bcc2050.latitude$coefficients)
bcc2050.latitude$y <- c("Latitude_BCC_2050")

odo_alt <-  pgls(DifavgAltitude ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                 data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_alt$phyres)
qqnorm(odo_alt$phyres)
qqline(odo_alt$phyres)
plot(x = fitted(odo_alt), y = odo_alt$phyres, pch = 19)
bcc2050.altitude <- summary(odo_alt)
bcc2050.altitude <- as.data.frame(bcc2050.altitude$coefficients)
bcc2050.altitude$y <- c("Altitude_BCC_2050")

odo_area <- pgls(cellDifArea ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                 data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_area$phyres)
qqnorm(odo_area$phyres)
qqline(odo_area$phyres)
plot(x = fitted(odo_area), y = odo_area$phyres, pch = 19)
bcc2050.area <- summary(odo_area)
bcc2050.area <- as.data.frame(bcc2050.area$coefficients)
bcc2050.area$y <- c("Area_BCC_2050")

# PGLS table BCC 2050
pgls_bcc_2050 <- rbind(bcc2050.latitude, bcc2050.altitude, bcc2050.area)
rm(bcc2050.latitude, bcc2050.altitude, bcc2050.area,
   odo_lat, odo_alt, odo_area)

# PGLS BCC 2070
str(bcc2070.noOut)
bcc2070.noOut$Habitat_2 <- as.factor(bcc2070.noOut$Habitat_2)
bcc2070.noOut$flight_mode <- as.factor(bcc2070.noOut$flight_mode)

bcc2070.noOut$Species <- gsub(" ","_",bcc2070.noOut$Species)
test.data <- comparative.data(phy = tree, data = bcc2070.noOut, 
                              names.col = Species, vcv = TRUE, vcv.dim = 3)

#now the PGLS model can be fitted:
odo_lat <- pgls(LatDif ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_lat$phyres)
qqnorm(odo_lat$phyres)
qqline(odo_lat$phyres)
plot(x = fitted(odo_lat), y = odo_lat$phyres, pch = 19)
bcc2070.latitude <- summary(odo_lat)
bcc2070.latitude <- as.data.frame(bcc2070.latitude$coefficients)
bcc2070.latitude$y <- c("Latitude_BCC_2070")

odo_alt <- pgls(DifavgAltitude ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_alt$phyres)
qqnorm(odo_alt$phyres)
qqline(odo_alt$phyres)
plot(x = fitted(odo_alt), y = odo_alt$phyres, pch = 19)
bcc2070.altitude <- summary(odo_alt)
bcc2070.altitude <- as.data.frame(bcc2070.altitude$coefficients)
bcc2070.altitude$y <- c("Altitude_BCC_2070")

odo_area <- pgls(cellDifArea ~ TotMax..mm. + FlightSeason..month. + Habitat_2, 
                 data = test.data, kappa = "ML", lambda = "ML", delta = "ML")
hist(odo_area$phyres)
qqnorm(odo_area$phyres)
qqline(odo_area$phyres)
plot(x = fitted(odo_area), y = odo_area$phyres, pch = 19)
bcc2070.area <- summary(odo_area)
bcc2070.area <- as.data.frame(bcc2070.area$coefficients)
bcc2070.area$y <- c("Area_BCC_2070")

# PGLS table BCC 2070
pgls_bcc_2070 <- rbind(bcc2070.latitude, bcc2070.altitude, bcc2070.area)
rm(bcc2070.latitude, bcc2070.altitude, bcc2070.area,
   odo_lat, odo_alt, odo_area)


# Coefficient table total
odo.table.sm <- rbind(pgls_miroc_2050, pgls_miroc_2070,
                      pgls_nor_2050, pgls_nor_2070,
                      pgls_bcc_2050, pgls_bcc_2070)
write.csv(odo.table.sm, "./Odonata_pgls_summary.csv")

rm(pgls_miroc_2050, pgls_miroc_2070,
   pgls_nor_2050, pgls_nor_2070,
   pgls_bcc_2050, pgls_bcc_2070)

# Clean workspace
rm(list = ls())