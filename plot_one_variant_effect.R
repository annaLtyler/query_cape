plot_one_variant_effect <- function(data_obj, geno_obj, marker_name,
    allele_name = NULL, pheno_type = c("pheno", "norm_pheno", "ET"), p_or_q = 0.05, 
    verbose = FALSE){

    query_genotype <- data_obj$query_genotype
    trait_cols <- categorical_pal(8)
    pheno_type <- pheno_type[1]

    if(pheno_type == "pheno"){
        pheno <- data_obj$pheno
    }
    if(pheno_type == "norm_pheno"){
        pheno <- apply(data_obj$pheno, 2, rankZ)
    }
    if(pheno_type == "ET"){
        pheno <- data_obj$ET
    }

    common.ind <- Reduce("intersect", list(rownames(pheno), rownames(geno_obj), rownames(query_genotype)))
    common.pheno.locale <- match(common.ind, rownames(pheno))
    common.geno.locale <- match(common.ind, rownames(geno_obj))
    common.query.locale <- match(common.ind, rownames(query_genotype))
    matched.pheno <- pheno[common.pheno.locale,]
    matched.geno <- geno_obj[common.geno.locale,,]
    matched.query <- query_genotype[common.query.locale,]

    var_int <- write_variant_influences(data_obj, p_or_q = p_or_q, 
        include_main_effects = FALSE, mark_covar = FALSE, 
        write_file = FALSE)

    test.source.idx <- grep(marker_name, var_int[,"Source"]) 
    test.target.idx <- grep(marker_name, var_int[,"Target"])

    if(length(test.source.idx) == 0 && length(test.target.idx) == 0){
        return("No interactions found for this marker")
    }

    if(!is.null(allele_name)){
        source.allele.idx <- sapply(allele_name, 
            function(x) grep(paste0("_", x), var_int[test.source.idx,"Source"]))
        test.source.idx <- test.source.idx[unlist(source.allele.idx)]
        target.allele.idx <- sapply(allele_name, 
            function(x) grep(paste0("_", x), var_int[test.target.idx,"Target"]))
        test.target.idx <- test.target.idx[unlist(target.allele.idx)]
    }
        
    get_int <- function(phenoV, geno1, geno2){
        geno_pairs <- pair.matrix(c(0,1), self.pairs = TRUE, ordered = TRUE)
        pheno.groups <- apply(geno_pairs, 1, function(x) phenoV[intersect(which(geno1 == x[1]), which(geno2 == x[2]))])
        #boxplot(pheno.groups, names = apply(geno_pairs, 1, function(x) paste(x, collapse = "_")))
        group.means <- sapply(pheno.groups, function(x) mean(x, na.rm = TRUE))
        #centered.groups <- lapply(pheno.groups, function(x) x-group.means[1])
        #boxplot(centered.groups, names = apply(geno_pairs, 1, function(x) paste(x, collapse = "_")))
        group.se <- sapply(pheno.groups, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
        centered.means <- group.means - group.means[1]
        add.pred <- centered.means[2] + centered.means[4]

        result <- rbind(c(centered.means[4], centered.means[2], add.pred, centered.means[3]),
              c(group.se[4], group.se[2], group.se[2] + group.se[4], group.se[3]))  
        rownames(result) <- c("effect", "se")
        colnames(result) <- c("main1", "main2", "add", "int")
        return(result)
    }

    plot.effects <- function(marker.idx, source.target = c("Source", "Target")){
        marker.label <- var_int[marker.idx,source.target]
        split.marker <- strsplit(marker.label, "_")
        marker.names <- sapply(split.marker, function(x) x[1])
        allele.names <- sapply(split.marker, function(x) x[2])
        marker.geno <- sapply(1:length(marker.names), function(x) matched.geno[,allele.names[x], marker.names[x]])
        colnames(marker.geno) <- marker.names
        
        marker.int <- lapply(1:ncol(matched.pheno), function(y) lapply(1:ncol(marker.geno), function(x) get_int(matched.pheno[,y], marker.geno[,x], matched.query)))
        names(marker.int) <- colnames(matched.pheno)

        for(p in 1:length(marker.int)){
            mean.vals <- marker.int[[p]][[1]][1,]
            se.vals <- marker.int[[p]][[1]][2,]
            max.val <- max(mean.vals+se.vals)
            min.val <- min(mean.vals-se.vals)
            if(max.val < 0){max.val <- 0}
            if(min.val > 0){min.val <- 0}

            a <- barplot(mean.vals, 
            main = paste(colnames(matched.pheno)[p], marker.label, sep = "\n"), 
            names = c(marker.label, "query", "Additive", "Actual"),
            col = trait_cols[as.numeric(allele.names)], ylim = c(min.val, max.val))
            segments(a[,1], mean.vals-se.vals, a[,1], mean.vals+se.vals)
            segments(a[,1]-0.2, mean.vals-se.vals, a[,1]+0.2, mean.vals-se.vals)
            segments(a[,1]-0.2, mean.vals+se.vals, a[,1]+0.2, mean.vals+se.vals)
            abline(h = 0)
        }

        add.dev <- sapply(marker.int, function(x) x[[1]][1,"int"] - x[[1]][1,"add"])
        if(!is.null(nrow(add.dev))){
            colnames(add.dev) <- colnames(pheno)
        }else{
            names(add.dev) <- colnames(pheno)
        }
        return(add.dev)
    }

    source.deviation <- vector(mode = "list", length = length(test.source.idx))
    names(source.deviation) <- var_int[test.source.idx,"Source"]
    
    if(length(test.source.idx) > 0){
        for(s in 1:length(test.source.idx)){
            #par(mfrow = c(1,2))
            source.deviation[[s]] <- plot.effects(test.source.idx[s], "Source")
        }
    }

    target.deviation <- vector(mode = "list", length = length(test.target.idx))
    names(target.deviation) <- var_int[test.target.idx,"Target"]

    if(length(test.target.idx) > 0){
        for(s in 1:length(test.target.idx)){
            target.deviation[[s]] <- plot.effects(test.target.idx[s], "Target")
        }
    }

    results <- list("deviation_with_query_as_source" = target.deviation,
        "deviation_with_query_as_target" = source.deviation)
    invisible(results)

}
