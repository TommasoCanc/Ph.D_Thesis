# Function to download record from NCBI
seqNcbi <- function(x, tr = "[Organism]", coord = F, dataBase = "nucleotide"){
  
  # Packages
  if (!require("rentrez")) install.packages("rentrez")
  require(rentrez)
  if (!require("readtext")) install.packages("readtext")
  require(readtext)
  if (!require("stringr")) install.packages("stringr")
  require(stringr)
  if (!require("progress")) install.packages("progress")
  require(progress)
  
  if(coord == F){
    ncbi <- data.frame()
    for (i in 1:length(x)) {
      taxon <- x[i]
      print(taxon)
      a <- entrez_search(db=dataBase, term=paste0(taxon,tr),use_history = T)
      a <- entrez_search(db=dataBase, term=paste0(taxon,tr),retmax = a$count, use_history = T)
      
      ncbi.1 <- data.frame()
      
      pb <- progress_bar$new(total = length(a$ids))
      if (length(a$ids) != 0) {
        for (j in 1:length(a$ids)) {
          ncbi.2 <- as.data.frame(matrix(NA, ncol=12))
          colnames(ncbi.2) <- c("sampleid", "species_name", "original_name","markercode", 
                                "nucleotides", "nucleotides_bp", "definition", "voucher", "pubmed",
                                "country", "collection_date", "INV")
          
          gbank <- entrez_fetch(db="nucleotide", id=a$ids[j], rettype="gbwithparts", retmode = "text")
          
          ncbi.2$sampleid[1] <- as.character(gsub("\"","",gsub("^.*ACCESSION\\s*|\\s*\n.*$", "", gbank)))
          taxonomy <- gsub("\"", "", gsub("^.*ORGANISM\\s*|\\s*.\nREFERENCE.*$", "", gbank))
          taxonomy <- unlist(strsplit(taxonomy, "\n"))
          ncbi.2$species_name[1] <- taxonomy[1]
          
          ncbi.2$original_name[1] <- taxon
          
          ncbi.2$markercode <- ifelse(grepl("product=", gbank) == T,gsub("\"","",gsub("^.*product=\\s*|\\s*\n.*$", "", gbank)), NA)
          
          # VOUCHER
          ncbi.2$voucher[1] <- ifelse(grepl("specimen_voucher=", gbank) == T, as.character(gsub("\"","",gsub("^.*specimen_voucher=\\s*|\\s*\n.*$", "", gbank))), NA)
          
          # DEFINITION
          ncbi.2$definition[1] <- gsub("\"", "", gsub("^.*DEFINITION\\s*|\\s*.\n.*$", "", gbank))
          
          # PUBMED
          ncbi.2$pubmed[1] <- ifelse(grepl("PUBMED", gbank) == T, as.character(gsub("\"", "", gsub("^.*PUBMED\\s*|\\s*\n.*$", "", gbank))), NA)
          
          # COUNTRY
          ncbi.2$country[1] <- ifelse(grepl("country", gbank) == T, as.character(gsub("\"", "", gsub("^.*country=\\s*|\\s*\n.*$", "", gbank))), NA)
      
          # Collection date
          ncbi.2$collection_date[1] <- ifelse(grepl("collection_date", gbank) == T, as.character(gsub("\"", "", gsub("^.*collection_date=\\s*|\\s*\n.*$", "", gbank))), NA)
          
          # INV
          ncbi.2$INV[1] <- ifelse(grepl("INV", gbank) == T, as.character(gsub("\"", "", gsub("^.*INV\\s*|\\s*\n.*$", "", gbank))), NA)
          
          
          seq <- gsub("\n", "", gsub("^.*ORIGIN\\s*|\\//.*$", "", gbank))
          ncbi.2$nucleotides[1] <- gsub(" ", "", gsub("[[:digit:]]+", "", seq))
          ncbi.2$nucleotides_bp[1] <- nchar(as.character(ncbi.2$nucleotides[1]))
          
          ncbi.1 <- rbind(ncbi.1, ncbi.2)
          
          pb$tick()
          Sys.sleep(i / length(a$ids))
        }
      }
      
      ncbi <- rbind(ncbi, ncbi.1)
      rm(ncbi.1,ncbi.2)
      print(paste(i, "of", length(x)))
      
      
    }
  }
  
  if(coord == T){
    ncbi <- data.frame()
    for (i in 1:length(x)) {
      taxon <- x[i]
      print(taxon)
      a <- entrez_search(db=dataBase, term=paste0(taxon,tr),use_history = T)
      a <- entrez_search(db=dataBase, term=paste0(taxon,tr),retmax = a$count, use_history = T)
      
      ncbi.1 <- data.frame()
      
      pb <- progress_bar$new(total = length(a$ids))
      if (length(a$ids) != 0) {
        for (j in 1:length(a$ids)) {
          
          ncbi.2 <- as.data.frame(matrix(NA, ncol=15))
          colnames(ncbi.2) <- c("sampleid", "species_name", "original_name","country",
                                "lat", "lon", "markercode", "nucleotides", "nucleotides_bp",
                                "definition", "voucher", "pubmed", "country", "collection_date",
                                "INV")
          
          gbank <- entrez_fetch(db="nucleotide", id=a$ids[j], rettype="gbwithparts", retmode = "text")
          
          ncbi.2$sampleid[1] <- as.character(gsub("\"","",gsub("^.*ACCESSION\\s*|\\s*\n.*$", "", gbank)))
          taxonomy <- gsub("\"", "", gsub("^.*ORGANISM\\s*|\\s*.\nREFERENCE.*$", "", gbank))
          taxonomy <- unlist(strsplit(taxonomy, "\n"))
          ncbi.2$species_name[1] <- taxonomy[1]
          
          ncbi.2$original_name[1] <- taxon
          
          ncbi.2$country <- as.character(ifelse(grepl("country", gbank) == T,gsub("\"","",gsub("^.*country=\\s*|\\s*\n.*$", "", gbank)), NA))
          
          lat_lon <- ifelse(grepl("lat_lon", gbank) == T,gsub("\"","",gsub("^.*lat_lon=\\s*|\\s*\n.*$", "", gbank)), NA)
          lon <- as.numeric(word(lat_lon,-2))
          ncbi.2$lon <- ifelse(word(lat_lon,-1) == "W", -abs(lon), lon)
          lat <- as.numeric(gsub("\"", "", as.character(word(lat_lon,1))))
          ncbi.2$lat <- ifelse(word(lat_lon,2) == "S", -abs(lat), lat)
          
          ncbi.2$markercode <- ifelse(grepl("product=", gbank) == T,gsub("\"","",gsub("^.*product=\\s*|\\s*\n.*$", "", gbank)), NA)
          
          seq <- gsub("\n", "", gsub("^.*ORIGIN\\s*|\\//.*$", "", gbank))
          ncbi.2$nucleotides[1] <- gsub(" ", "", gsub("[[:digit:]]+", "", seq))
          ncbi.2$nucleotides_bp[1] <- nchar(as.character(ncbi.2$nucleotides[1]))
          
          # VOUCHER
          ncbi.2$voucher[1] <- ifelse(grepl("specimen_voucher=", gbank) == T, as.character(gsub("\"","",gsub("^.*specimen_voucher=\\s*|\\s*\n.*$", "", gbank))), NA)
          
          # DEFINITION
          ncbi.2$definition[1] <- gsub("\"", "", gsub("^.*DEFINITION\\s*|\\s*.\n.*$", "", gbank))
          
          # PUBMED
          ncbi.2$pubmed[1] <- ifelse(grepl("PUBMED", gbank) == T, as.character(gsub("\"", "", gsub("^.*PUBMED\\s*|\\s*\n.*$", "", gbank))), NA)
          
          # COUNTRY
          ncbi.2$country[1] <- ifelse(grepl("country", gbank) == T, as.character(gsub("\"", "", gsub("^.*country=\\s*|\\s*\n.*$", "", gbank))), NA)
          
          # Collection date
          ncbi.2$collection_date[1] <- ifelse(grepl("collection_date", gbank) == T, as.character(gsub("\"", "", gsub("^.*collection_date=\\s*|\\s*\n.*$", "", gbank))), NA)
          
          # INV
          ncbi.2$INV[1] <- ifelse(grepl("INV", gbank) == T, as.character(gsub("\"", "", gsub("^.*INV\\s*|\\s*\n.*$", "", gbank))), NA)
          
          ncbi.1 <- rbind(ncbi.1, ncbi.2)
          
          pb$tick()
          Sys.sleep(i / length(a$ids))
        }
      }
      
      ncbi <- rbind(ncbi, ncbi.1)
      rm(ncbi.1,ncbi.2)
      print(paste(i, "of", length(x)))
      
      
    }
  }
  
  ncbi
}