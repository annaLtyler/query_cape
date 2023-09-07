#' Generates a null distribution for the pairscan
#' 
#' This function generates a null distribution
#' for the pairscan. For each permutation,
#' it runs a single scan and selects markers
#' in the same manner as for the true test.
#'
#' @param data_obj a \code{\link{Cape}} object
#' @param geno_obj a genotype object
#' @param kin_obj a kinship object
#' @param scan_what A character string uniquely identifying whether eigentraits
#'   or raw traits should be scanned. Options are "eigentraits", "raw_traits"
#' @param pairscan_null_size The total size of the null distribution.
#' This is DIFFERENT than the number of permutations to run. Each permutation
#' generates n choose 2 elements for the pairscan. So for example, a permutation
#' that tests 100 pairs of markers will generate a null distribution of size 4950.
#' This process is repeated until the total null size is reached. If the null size
#' is set to 5000, two permutations of 100 markers would be done to get to a null
#' distribution size of 5000.
#' @param max_pair_cor A numeric value between 0 and 1 indicating the maximum
#'   Pearson correlation that two markers are allowed. If the correlation
#'   between a pair of markers exceeds this threshold, the pair is not tested.
#'   If this value is set to NULL, min_per_genotype must have a numeric value.
#' @param min_per_geno The minimum number of individuals allowable per
#'   genotype. If for a given marker pair, one of the genotypes is
#'   underrepresented, the marker pair is not tested. If this value is NULL,
#'   max_pair_cor must have a numeric value.
#' @param model_family Indicates the model family of the phenotypes. This can be 
#'   either "gaussian" or "binomial".
#' @param marker_selection_method options are "top_effects", "uniform", "effects_dist", "by_gene"
#' @param run_parallel Whether to run the analysis on multiple CPUs.
#' @param n_cores The number of CPUs to use if run_parallel is TRUE.
#' @param verbose Whether to print progress to the screen. Defaults to FALSE.
#' 
#' @return This function returns a list with two elements, one containing
#' the results of the permutations, and the other containing the markers
#' that were tested in the individual permutations.
#' @keywords internal
#' 
pairscan_null_kin_query <- function(data_obj, geno_obj = NULL, kin_obj = NULL, 
  query_genotype, scan_what = c("eigentraits", "raw_traits"), pairscan_null_size = NULL, 
  max_pair_cor = NULL, min_per_geno = NULL, model_family = "gaussian", 
  marker_selection_method = c("top_effects", "uniform", "effects_dist", "by_gene"), 
  run_parallel = FALSE, n_cores = 4, verbose = FALSE){

  marker_selection_method = "from_list"
  
  ref_allele <- data_obj$ref_allele
    
  if(is.null(pairscan_null_size)){
    stop("The total number of permutations must be specified.")
  }  
  
  #If the user does not specify a scan_what, 
  #default to eigentraits, basically, if eigen,
  #et, or ET are anywhere in the string, use the
  #eigentraits, otherwise, use raw phenotypes
  type_choice <- c(grep("eig", scan_what, ignore.case = TRUE), grep("ET", scan_what, ignore.case = TRUE)) #look for any version of eigen or eigentrait, the user might use.
  if(length(type_choice) > 0){ #if we find any, use the eigentrait matrix
    pheno <- data_obj$ET
  }else{
    pheno <- data_obj$pheno #otherwise, use the raw phenotype matrix
  }
  

  num_pheno <- dim(pheno)[2]
    
  #make a list to hold the results. 
  #Tables from each of the phenotypes will be
  #put into this list
  results_perm_list <- vector(mode = "list", length = num_pheno)
  
  #generate the null distribution for the pairscan 
  #do the pairscan on these markers for all traits
  
  if(verbose){
    cat("\nGenerating null distribution...\n")
  }
  
  all_pairs_tested <- NULL
  
  final_perm <- 1
  while(final_perm < pairscan_null_size){
    
    # TODO what is happening here? should this be a deep copy?
    perm_data_obj <- data_obj
    perm_data_obj$pheno <- perm_data_obj$pheno[sample(nrow(perm_data_obj$pheno)),]
    
    #we don't need to run the singlescans in query cape because we test
    #pairs exhaustively, so there is no marker selection step  

    #run the pairscan for each permuted phenotype and the pairs we just found
    if(verbose){cat("Performing marker pair scans of permuted traits with kinship correction...\n")}
    
    #run a pairscan on these markers and each permuted phenotype
    pairscan_results <- pairscan_query(data_obj = perm_data_obj, geno_obj = geno_obj, 
        query_genotype = query_genotype, 
        kin_obj = kin_obj, scan_what = scan_what, run_parallel = run_parallel, 
        n_cores = n_cores, verbose = verbose, overwrite_alert = FALSE, 
        pairscan_null_size = 0, max_pair_cor = data_obj$max_pair_cor, 
        min_per_genotype = data_obj$min_per_genotype)
    
    #integrate the results into the permutation object
    one_perm <- pairscan_results$pairscan_results
    #because there will be different numbers of markers each time, just take 
    #the marker names, the intercept, and the effects for marker1 marker2 and 
    #their interaction
    # last.col = dim(one_perm[[1]])[2]
    # take.col <- c(1:3, (last.col-2):last.col)
    if(final_perm == 1){ #if this is the first time through, 
      #just copy the results into the results_perm_list
      for(p in 1:num_pheno){
        results_perm_list[[p]] <- one_perm[[p]]
      }
        all_pairs_tested <- one_perm[[1]][[1]][,1:2]
    }else{
      if(!is.null(one_perm)){
        for(p in 1:num_pheno){
            for(i in 1:length(one_perm[[p]]))
                results_perm_list[[p]][[i]] <- rbind(results_perm_list[[p]][[i]], one_perm[[p]][[i]])
          }
        all_pairs_tested <- rbind(all_pairs_tested, one_perm[[1]][[1]][,1:2])
        }
      }
    
    final_perm <- nrow(results_perm_list[[1]][[1]]) #end looping through phenotypes
    if(verbose){cat("\t", final_perm, " null tests: ", round((final_perm/pairscan_null_size)*100), "%...\n", sep = "")} 
    }#end when we have enough permutations
    
  
  names(results_perm_list) <- colnames(pheno)
  results_list <- list("pairscan_perm" = results_perm_list, "pairs_tested_perm" = all_pairs_tested)
  return(results_list)

} 
  
  

