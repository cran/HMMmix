.First.lib <- function(lib, pkg){
    library.dynam("HMMmix", pkg, lib)
    cat("HMMmix Loaded \n")
   
    
}
