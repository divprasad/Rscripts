#install.packages("RCurl")

library(RCurl)


url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
#url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"

filenames <- getURL(url, ftp.use.epsv = FALSE,dirlistonly = TRUE) 

#for (i in )

  #{
  
   
    if i % 50 ==0
    
    {
      # Sys.sleep(1)
    }
    
    

 # }




