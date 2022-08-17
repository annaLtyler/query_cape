plot_variant_effects_query <- function(data_obj, geno_obj, 
    pheno_type = c("pheno", "norm_pheno", "ET"), 
    p_or_q = 0.05, scale_coord = 1, verbose = FALSE){

    query_genotype <- data_obj$query_genotype
    trait_cols <- categorical_pal(ncol(geno_obj))
    pheno_type <- pheno_type[1]

    coord_label = paste0("Position (bp/", scale_coord, ")")
    if(scale_coord == 1){
        coord_label = "Position (bp)"
    }
    if(scale_coord == 1000){
        coord_label = "Position (Kb)"
    }
    if(scale_coord == 1e6){
        coord_label = "Position (Mb)"
    }

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

    #each interaction in query_cape has the query genotype as
    #the source or the target.
    #plot these types of interactions separately.
    source.chr <- as.numeric(var_int[,2])
    source.pos <- as.numeric(var_int[,3])/scale_coord
    target.chr <- as.numeric(var_int[,5])
    target.pos <- as.numeric(var_int[,6])/scale_coord

    u_chr <- sort(unique(c(source.chr, target.chr)))
    u_chr <- u_chr[-1] #take off the 0 chromosome. this is the query marker
        
    get_int <- function(phenoV, geno1, geno2){
        geno_pairs <- pair.matrix(c(0,1), self.pairs = TRUE, ordered = TRUE)
        pheno.groups <- apply(geno_pairs, 1, function(x) phenoV[intersect(which(geno1 == x[1]), which(geno2 == x[2]))])
        #boxplot(pheno.groups);abline(h = 0)
        group.means <- sapply(pheno.groups, function(x) mean(x, na.rm = TRUE))
        #group.se <- sapply(pheno.groups, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
        ref.effect <- group.means[1]
        centered.means <- group.means - ref.effect
        source.effect <- centered.means[4]
        target.effect <- centered.means[2]
        
        add.pred <- source.effect + target.effect

        result <- c("source.effect" = source.effect, "target.effect" = target.effect, 
            "additive" = add.pred, "actual" = centered.means[3])
        #barplot(result)

    
        return(result)
    }


    source.deviation <- target.deviation <- vector(mode = "list", length = length(u_chr))
    names(source.deviation) <- names(target.deviation) <- u_chr

    source.ch.pos <- target.ch.pos <- vector(mode = "list", length = length(u_chr))
    names(source.ch.pos) <- names(target.ch.pos) <- u_chr
    
    source.allele <- target.allele <- vector(mode = "list", length = length(u_chr))
    names(source.allele) <- names(target.allele) <- u_chr

    for(ch in u_chr){
        if(verbose){report.progress(ch, length(u_chr))}
        target.chr.locale <- which(target.chr == ch)
        source.chr.locale <- which(source.chr == ch)

        if(length(target.chr.locale) > 0){
            target.markers <- var_int[target.chr.locale,"Target",drop=FALSE]
            split.target <- strsplit(target.markers, "_")
            target.marker.names <- sapply(split.target, function(x) x[1])
            target.allele.names <- sapply(split.target, function(x) x[2])
            target.allele[[ch]] <- target.allele.names
            target.geno <- sapply(1:length(target.marker.names), function(x) matched.geno[,target.allele.names[x], target.marker.names[x],drop=FALSE])
            colnames(target.geno) <- target.markers
            target.int <- lapply(1:ncol(matched.pheno), function(y) t(apply(target.geno, 2, function(x) get_int(matched.pheno[,y], matched.query, x))))
            target.ch.pos[[ch]] <- target.pos[target.chr.locale]
            target.effects <- sapply(target.int, function(x) x[,"actual"] - x[,"additive"])
            if(is.null(dim(target.effects))){
                target.effects <- matrix(target.effects, nrow = 1)
            }
            colnames(target.effects) <- colnames(pheno)
            target.deviation[[ch]] <- target.effects

        }

        if(length(source.chr.locale) > 0){
            source.markers <- var_int[source.chr.locale,"Source",drop=FALSE]        
            split.source <- strsplit(source.markers, "_")
            source.marker.names <- sapply(split.source, function(x) x[1])
            source.allele.names <- sapply(split.source, function(x) x[2])

            source.allele[[ch]] <- source.allele.names

            source.geno <- sapply(1:length(source.marker.names), function(x) matched.geno[,source.allele.names[x], source.marker.names[x],drop=FALSE])
            colnames(source.geno) <- source.markers    
            source.int <- lapply(1:ncol(matched.pheno), function(y) t(apply(source.geno, 2, function(x) get_int(matched.pheno[,y], x, matched.query))))

            source.ch.pos[[ch]] <- source.pos[source.chr.locale]

            source.effects <- sapply(source.int, function(x) x[,"actual"] - x[,"additive"])
            if(is.null(dim(source.effects))){
                source.effects <- matrix(source.effects, nrow = 1)
            }
            colnames(source.effects) <- colnames(pheno)
            source.deviation[[ch]] <- source.effects

        }
    }

    #find ylim across all chromosomes
    ylim = c(min(c(min(unlist(target.deviation), na.rm = TRUE),min(unlist(source.deviation), na.rm = TRUE))),
           max(c(max(unlist(target.deviation), na.rm = TRUE),max(unlist(source.deviation), na.rm = TRUE))))


    for(ch in 1:length(source.deviation)){
        #quartz(width = 10, height = 6)
        xlim <- c(1, max(c(source.ch.pos[[ch]], target.ch.pos[[ch]])))

        layout.mat <- matrix(c(1,3,2,4), ncol = 2, byrow = FALSE)
        layout(layout.mat, widths = c(1,0.3))
        par(mar = c(2,4,2,2))

        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
        for(p in 1:ncol(pheno)){
            #plot effect with query as target
            points(source.ch.pos[[ch]], source.deviation[[ch]][,p], type = "h",
                col = trait_cols[p])
            #add allele colors
            points(source.ch.pos[[ch]], source.deviation[[ch]][,p], type = "p",
                col = trait_cols[as.numeric(source.allele[[ch]])], pch = 16,
                cex = 0.5)

        } 
        axis(1); axis(2)
        abline(h = 0)
        mtext("Deviation of Trait from Additive", side = 2, line = 2.5)
        mtext("Query as Target", side = 4)
        
        plot.new()
        plot.window(xlim = c(0,1), ylim = c(0,1))        
        par(mar = c(0,0,0,0))
        legend(0.5,0.5, col = trait_cols[1:ncol(pheno)], lty = 1, lwd = 2,
            legend = colnames(pheno))

        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
        par(mar = c(4,4,0,2))
        for(p in 1:ncol(pheno)){
            points(target.ch.pos[[ch]], target.deviation[[ch]][,p], type = "h",
                col = trait_cols[p])
            points(target.ch.pos[[ch]], target.deviation[[ch]][,p], type = "p",
                col = trait_cols[as.numeric(target.allele[[ch]])], pch = 16,
                cex = 0.5)

        } 
        axis(1); axis(2)
        abline(h = 0)
        mtext("Deviation of Trait from Additive", side = 2, line = 2.5)
        mtext(coord_label, , side = 1, line = 2.5)
        mtext("Query as Source", side = 4)

        plot.new()
        plot.window(xlim = c(0,1), ylim = c(0,1))        
        par(mar = c(0,0,0,0))
        legend(x = 0.5, y = 1, col = trait_cols, pch = 16,
            legend = 1:ncol(geno_obj))

        mtext(paste("Chr", ch), side = 3, outer = TRUE, line = -1)
    }

    results <- list("deviation_with_query_as_source" = target.deviation,
        "deviation_with_query_as_target" = source.deviation)
    invisible(results)

}
