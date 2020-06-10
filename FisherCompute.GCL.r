FisherCompute.GCL <- function(freq, loci, prec = 4){
	
  if(!require("pacman")) install.packages("pacman"); library(pacman); pacman::p_load(tidyverse) #Install packages, if not in library and then load them.

  #Getting pvalues for each locus.
  pval <- sapply(loci, function(locus){
    
   my.freq <- freq %>% 
      filter(locus == !!locus) %>% 
     select(silly, allele_no, freq) 
   
   freq.mat <- matrix(my.freq$freq, ncol = max(my.freq$allele_no), dimnames = list(unique(my.freq$silly), unique(my.freq$allele_no)), byrow = TRUE)
   
   fisher.test(freq.mat, workspace = 8000000, hybrid = TRUE)$p.value
    
  })
  
  pval[pval==0] = min(pval[pval!=0])#I think this is something Jim added so you can calculate overall p-value. I'm not sure why he decided to use the minimum pvalue to replace zeros.
  
  #Using Fisher's method to get overall pvalue
  overall <- pchisq(q =- 2*sum(log(pval)), df = 2*nloci, lower.tail = FALSE) %>% 
    round(prec)
  
  bylocus <- tibble::tibble(locus = names(pval), pval = pval %>% round(prec))
  
  test <- paste0(unique(freq$silly), collapse = ".")
  
  output <- tibble::tibble(silly = test,  overall = !!overall, bylocus = !!bylocus %>% tidyr::nest({{test}} := dplyr::everything()))

  return(output)
  
}
 


