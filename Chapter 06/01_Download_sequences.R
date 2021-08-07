######################
# Download Sequences #
######################

# Load libraries
library(bold) # <- Bold connection
library(osfr) # <- Connect to OSF
library(xlsx) # <- Read xlsx

#####################
# Connection to OSF #
#####################

main.path <- "YOUR_PATH"
# Create folder where you can download the files. 
# For more information: ?dir.create
dir.create(main.path)
setwd(main.path)

# Connect to Open Science Framework
cr_project <- osf_retrieve_node("4p63t")
osf_files <- osf_ls_files(cr_project, type = "folder")

file.number <- which(osf_ls_files(osf_files[4, ])[ ,1] == "Complessive_distibution.xlsx")
osf_download(osf_ls_files(osf_files[4, ])[file.number, ], path = main.path)
rm(file.number)

file.number <- which(osf_ls_files(osf_files[4, ])[ ,1] == "NCBI_download.R")
osf_download(osf_ls_files(osf_files[4, ])[file.number, ], path = main.path)
rm(file.number)

# Load function
source("./NCBI_download.R")


# LOAD DATA ####
taxa <- read.xlsx("./Complessive_distibution.xlsx", sheetIndex = 1, header = TRUE)

unique(taxa$Class)

ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Ephemeroptera"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Odonata"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Plecoptera"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Bivalvia"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Branchiobdellida"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Bryozoa"]) 
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Diptera"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Crustacea"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Polychaeta"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Trichoptera"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Megaloptera"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Coelenterata"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Heteroptera"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Hydrachnidia"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Coleoptera"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Hymenoptera"])
# #ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Platyhelminthes"])
# #ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Planipennia"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Porifera"])
# #ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Araneae"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Lepidoptera"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Oligochaeta"])
# ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Nematoda"])
# #ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Nemertea"])
# #ephe <- as.character(taxa$species_Tax_Review[taxa$Class == "Namatomorpha"])

# Download record from NCBI ----------------------------------------------------
ephe.ncbi <- seqNcbi(ephe[1:length(ephe)], tr = "[Organism] AND COI[Gene]", 
                     coord = T, dataBase = "nucleotide")
# ephe.ncbi <- seqNcbi(ephe[1:length(ephe)], 
# tr = "[Organism] AND (CO1 [GENE] OR COI [GENE] OR COX1 [GENE] OR COXI [GENE])", 
# coord = T, dataBase = "nucleotide")

write.csv(ephe.ncbi, "...", row.names=F)

# Download record from BOLD ----------------------------------------------------
bold <- data.frame()

for(i in 1:length(ephe)){ 
  a <- bold_seqspec(taxon = ephe[i], sepfasta = TRUE)
  
  if(!is.na(a)){
    b <- a$data
    c <- as.data.frame(unlist(a$fasta))
    colnames(c) <- "nucleotide"
    
    bold.1 <- cbind(b, c)
    bold <- rbind(bold, bold.1)
    
    #print(i/length(ephe))
  }
print(i)
  }

write.csv(bold, "...", row.names=F)

rm(list = ls())