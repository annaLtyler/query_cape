plot_clinical_effects_query <- function(data_obj, geno_obj, 
    pheno_type = c("pheno", "norm_pheno", "ET"), p_or_q = 0.05, 
    path = ".", verbose = FALSE){

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

    #each interaction in query_cape has the query genotype as
    #the source or the target.
    #plot these types of interactions separately.
    source.chr <- as.numeric(var_int[,2])
    source.pos <- as.numeric(var_int[,3])
    target.chr <- as.numeric(var_int[,5])
    target.pos <- as.numeric(var_int[,6])

    u_chr <- sort(unique(c(source.chr, target.chr)))
    u_chr <- u_chr[-1] #take off the 0 chromosome. this is the query marker
        
    get_effect <- function(phenoV, geno1, geno2){
        geno_pairs <- pair.matrix(c(0,1), self.pairs = TRUE, ordered = TRUE)
        pheno.groups <- apply(geno_pairs, 1, function(x) phenoV[intersect(which(geno1 == x[1]), which(geno2 == x[2]))])
        #boxplot(pheno.groups);abline(h = 0)
        group.means <- sapply(pheno.groups, function(x) mean(x, na.rm = TRUE))
        #group.se <- sapply(pheno.groups, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
        ref.effect <- group.means[1]
        centered.means <- group.means - ref.effect
        source.effect <- centered.means[4]
        target.effect <- centered.means[2]
        actual.effect <- centered.means[3]

        add.pred <- source.effect + target.effect

        if(is.finite(add.pred)){
            #are the main effects pushing in the same direction 
            #or opposite directions
            if(sign(source.effect) == sign(target.effect)){
                main.effect <- "coherent"
            }else{
                main.effect <- "incoherent"
            }

            sign.flipped = FALSE
            #if(main.effect == "coherent"){
            #    if(sign(actual.effect) != sign(source.effect)){
            #        sign.flipped = TRUE
            #    }
            #}
            if(sign(add.pred) != sign(actual.effect)){
                sign.flipped = TRUE
            }

            #is the actual trait closer to the reference
            #than predicted
            if(abs(actual.effect) < abs(add.pred)){
                int.effect <- "alleviating"
                }else{
                int.effect <- "aggravating"
                }
            }else{
                int.effect <- "none"
                main.effect <- "none"
            }
		
		pheno.result <- c("main1" = source.effect, "main2" = target.effect, 
        "additive" = add.pred, "actual" = actual.effect, "main.effects" = main.effect, 
        "effect" = int.effect, "sign.flipped" = sign.flipped)
		return(pheno.result)
    }

    source.effects <- target.effects <- vector(mode = "list", length = length(u_chr))
    names(source.effects) <- names(target.effects) <- u_chr

    for(ch in u_chr){
        if(verbose){report.progress(ch, length(u_chr))}
        target.chr.locale <- which(target.chr == ch)
        target.markers <- var_int[target.chr.locale,"Target"]

        source.chr.locale <- which(source.chr == ch)
        source.markers <- var_int[source.chr.locale,"Source"]

        if(length(target.markers)){        
            split.target <- strsplit(target.markers, "_")
            target.marker.names <- sapply(split.target, function(x) x[1])
            target.allele.names <- sapply(split.target, function(x) x[2])

            target.geno <- sapply(1:length(target.marker.names), function(x) matched.geno[,target.allele.names[x], target.marker.names[x]])
            colnames(target.geno) <- target.markers
            target.int <- lapply(1:ncol(matched.pheno), function(y) t(apply(target.geno, 2, function(x) get_effect(matched.pheno[,y], matched.query, x))))
            names(target.int) <- colnames(matched.pheno)
            for(ph in 1:length(target.int)){
                target.int[[ph]] <- cbind(target.int[[ph]], target.allele.names)
            }
            target.effects[[ch]] <- target.int
        }
        
        if(length(source.markers) > 0){
            split.source <- strsplit(source.markers, "_")
            source.marker.names <- sapply(split.source, function(x) x[1])
            source.allele.names <- sapply(split.source, function(x) x[2])

            source.geno <- sapply(1:length(source.marker.names), function(x) matched.geno[,source.allele.names[x], source.marker.names[x]])
            colnames(source.geno) <- source.markers
            source.int <- lapply(1:ncol(matched.pheno), function(y) t(apply(source.geno, 2, function(x) get_effect(matched.pheno[,y], x, matched.query))))
            names(source.int) <- colnames(matched.pheno)
            for(ph in 1:length(source.int)){
                source.int[[ph]] <- cbind(source.int[[ph]], source.allele.names)
            }

            source.effects[[ch]] <- source.int
        }
        
    }

    plot_motif_effects <- function(motif.table, query_position = c("Source", "Target"), 
        trait_name = ""){

        all.num.effects <- apply(motif.table[,c("main1", "main2", "additive", "actual")], 
            2, as.numeric)
        int.types <- motif.table[,c("main.effects", "effect", "sign.flipped")]
        allele.col <- grep("allele.names", colnames(motif.table))
        alleles <- as.numeric(motif.table[,allele.col])

        #get all combinations of coherent/incoherent and aggravating/alleviating
        #do not include those with flipped signs
        type.idx <- apply(type.sets, 1, 
            function(x) Reduce("intersect", list(which(int.types[,1] == x[1]), 
            which(int.types[,2] == x[2]), which(int.types[,3] == x[3]))))

        ylim <- c(min(all.num.effects), max(all.num.effects))
        plot.height <- ylim[2]-ylim[1]
        
        #dim(all.num.effects)
        pdf(file.path(path, paste0("motif_effects_", query_position, "_", trait_name, ".pdf")), 
            width = 12, height = 6)
        #layout.mat <- matrix(c(3,1,4,2,6,5), byrow = FALSE, nrow = 2)
        #layout(layout.mat)
        #par(mfrow = c(2,4))
        layout.mat <- matrix(c(1:8, 0,9,9,0), nrow = 3, byrow = TRUE)
        layout(layout.mat, heights = c(1,1,0.1))
        par(mar = c(2,4,4,2))
        for(it in 1:nrow(type.sets)){
            if(length(type.idx[[it]]) > 0){
                plot.new()
                plot.window(xlim = c(1,4), ylim = ylim)
                null <- lapply(1:length(type.idx[[it]]), function(x) points(x = 1:4, 
                    y = all.num.effects[type.idx[[it]][x],], type = "b", 
                    col = trait_cols[alleles[type.idx[[it]]][x]]))
                axis(2)
                text(x = 1:4, y = (ylim[1]-(plot.height*0)), 
                labels = c("Source", "Target", "Additive", "Actual"))
                abline(h = 0)
            }else{
                plot.text("No motifs in this category")
            }
            if(type.sets[it,3] == "FALSE"){
                mtext(paste(type.sets[it,1:2], collapse = " "), side = 3)
            }else{
                mtext(paste(c(type.sets[it,1:2], "\nsign flipped"), collapse = " "), side = 3)
            }
            
        }
        mtext(trait_name, side = 3, line = -1.5, outer = TRUE)
        par(mar = c(0,0,0,0))
        plot.new()
        plot.window(xlim = c(0,1), ylim = c(0,1))
        legend(0.1,1, legend = colnames(geno_obj), horiz = TRUE,
            col = trait_cols[1:ncol(geno_obj)], lty = 1, lwd = 3)
        dev.off()
    }

    #do interactions have similar or different effects across multiple conditions?
    all.target.effects <- lapply(1:ncol(pheno), function(y) Reduce("rbind", sapply(target.effects, function(x) x[[y]])))
    all.source.effects <- lapply(1:ncol(pheno), function(y) Reduce("rbind", sapply(source.effects, function(x) x[[y]])))
    names(all.target.effects) <- names(all.source.effects) <- colnames(matched.pheno)

    #target.main <- lapply(all.target.effects, function(x) as.factor(x[,"main.effects"]))
    #names(target.main) <- colnames(matched.pheno)
    #target.main.counts <- table(target.main)
    #source.main <- lapply(all.source.effects, function(x) as.factor(x[,"main.effects"]))
    #names(source.main) <- colnames(matched.pheno)
    #source.main.counts <- table(source.main)

    main.types <- c("coherent", "incoherent")
    int.effect <- c("aggravating", "alleviating")
    int.main.sets <- rbind(cbind(rep(main.types, length(int.effect)), 
        rep(int.effect, each = length(main.types))))
    type.sets <- rbind(cbind(int.main.sets, rep("FALSE", nrow(int.main.sets))),
                       cbind(int.main.sets, rep("TRUE", nrow(int.main.sets))))
    
    for(p in 1:length(all.source.effects)){
        plot_motif_effects(motif.table = all.target.effects[[p]], 
            query_position = "Source", trait_name = colnames(matched.pheno)[p])
        plot_motif_effects(motif.table = all.source.effects[[p]], 
            query_position = "Target", trait_name = colnames(matched.pheno)[p])
    }

    result <- list("Query_as_Source" = all.target.effects, 
        "Query_as_Target" = all.source.effects)
    invisible(result)
}
