# x: DNAStringSet object
# Detect how many cores have your pc detectCores()
# reads <- readDNAStringSet(path, format = "fasta")

haplo <- function(x, clus = F, cpu = 1, fasta = F){
  # Packages
  if (!require("QSutils")) install.packages("QSutils")
  require(QSutils)
  if (!require("doParallel")) install.packages("doParallel")
  require(doParallel)
  
  ht <- Collapse(x)
  
  if(clus == F){
    
    ht.seq <- data.frame()
    ht.report <- data.frame()
    
    for(i in 1:length(ht$nr)){
      haplo.seq <- data.frame(haplotype =  ht$hseqs[i],
                           haplotype_lab = i,
                           haplotype_num = ht$nr[i],
                           label = names(x[x %in% ht$hseqs[[i]]])[1])
      
      haplo.report <- data.frame(haplotype_number = rep(i, ht$nr[i]),
                                 label = names(x[x %in% ht$hseqs[[i]]]))
      
      ht.seq <- rbind(ht.seq, haplo.seq) 
      ht.report <- rbind(ht.report, haplo.report)
    }
  }
  
  if(clus == T){
    
    ht.seq <- data.frame()
    ht.report <- data.frame()
    
    cl <- makeCluster(cpu) # Numer of CPU
    registerDoParallel(cl)
    
    ht.seq <- foreach(i=1:length(ht$nr), .combine=rbind, .packages = "QSutils") %dopar% {
      haplo.seq = data.frame(haplotype =  ht$hseqs[i],
                             haplotype_lab = i,
                             haplotype_num = ht$nr[i],
                             label = names(x[x %in% ht$hseqs[[i]]])[1])
      haplo.seq
    }
    
    
    ht.report <- foreach(i=1:length(ht$nr), .combine=rbind, .packages = "QSutils") %dopar% {
      haplo.report = data.frame(haplotype_number = rep(i, ht$nr[i]),
                                label = names(x[x %in% ht$hseqs[[i]]]))
      haplo.report
    }
    #stop cluster
    stopCluster(cl)
    
  }
  
  if(fasta == F) {
  
    list(haplotype_seq = ht.seq, 
       haplotype_report = ht.report)
  }
  
  if(fasta == T) {
    y <- t(sapply(strsplit(ht.seq$haplotype,""), tolower))
    rownames(y) <- ht.seq$label
    y.fasta <- as.DNAbin(y)
    
    list(haplotype_seq = ht.seq, 
         haplotype_report = ht.report,
         fasta = y.fasta)
  }
}
