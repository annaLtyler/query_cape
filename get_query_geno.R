get_query_geno <- function(data_obj, geno_obj, query_genotype){

    covar_info <- get_covar(data_obj)
    covar_names <- covar_info$covar_names
    covar_table <- covar_info$covar_table
    colnames(covar_table) <- covar_names

    #create a large 2D matrix with all alleles
    geno <- get_geno(data_obj, geno_obj)

    gene <- matrix(NA, nrow = nrow(geno), ncol = dim(geno_obj)[2]*(dim(geno)[3]-1)) #take off the query allele, which is the last marker
    start_idx <- 1
    n_alleles <- dim(geno)[2]
    n_markers <- dim(geno)[3] - 1
    for(i in 1:n_markers){
        gene[,start_idx:(start_idx+n_alleles-1)] <- geno[,,i]
        start_idx = start_idx + n_alleles
    }
    fill_markers <- rep(dimnames(geno_obj)[[3]][-c(dim(geno_obj)[3])], each = n_alleles)
    fill_alleles <- rep(dimnames(geno_obj)[[2]], n_markers)
    colnames(gene) <- paste(fill_markers, fill_alleles, sep = "_")

    #add covariates and query genotype. Put the query genotype
    #in the last position
    #make sure the query genotype has the same individuals as geno
    common.ind <- intersect(rownames(geno), rownames(query_genotype))
    query_geno_idx <- match(common.ind, rownames(query_genotype))
    geno_idx <- match(common.ind, rownames(geno))


    gene <- cbind(gene[geno_idx,], covar_table, query_genotype[query_geno_idx,])
    colnames(gene)[ncol(gene)] <- "query"

    #if we are only looking at two alleles, A = 1-B, so we don't
    #need to test both. If this is the case, remove all alleles
    #that match the reference allele
    if(ncol(geno_obj) == 2){
        ref.allele <- data_obj$ref_allele
        keep.idx <- which(fill_alleles != ref.allele)
        gene <- gene[,keep.idx]
    }

    return(gene)
}