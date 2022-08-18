plot_variant_influences_query <- function(data_obj, geno_obj, p_or_q = 0.05){

    query_genotype <- data_obj$query_genotype
    trait_cols <- categorical_pal(ncol(data_obj$pheno))
    allele_cols <- categorical_pal(ncol(geno_obj))

    var_inf <- write_variant_influences(data_obj, p_or_q = p_or_q, 
        include_main_effects = TRUE, mark_covar = FALSE, 
        write_file = FALSE)
    main_locale <- which(var_inf[,"Target"] %in% colnames(data_obj$pheno))
    just_main <- var_inf[main_locale,]
    
    main_chr <- as.numeric(just_main[,2])
    main_pos <- as.numeric(just_main[,3])
    main_trait <- just_main[,4]
    u_trait <- unique(main_trait)
    main_col <- rep(NA, nrow(just_main))
    for(p in 1:length(u_trait)){
        main_col[which(main_trait == u_trait[p])] <- trait_cols[p]
    }
    main_effect <- as.numeric(just_main[,"Effect"])
    main.lim <- c(min(main_effect), max(main_effect))
    
    var_int <- write_variant_influences(data_obj, p_or_q = p_or_q, 
        include_main_effects = FALSE, mark_covar = FALSE, 
        write_file = FALSE)

    main_allele <- sapply(strsplit(just_main[,"Source"], "_"), function(x) x[2])
    
    #each interaction in query_cape has the query genotype as
    #the source or the target.
    #plot these types of interactions separately.
    source.chr <- as.numeric(var_int[,2])
    source.pos <- as.numeric(var_int[,3])

    target.chr <- as.numeric(var_int[,5])
    target.pos <- as.numeric(var_int[,6])

    all.effect <- as.numeric(var_int[,"Effect"])
    int.lim <- c(min(all.effect), max(all.effect))

    u_chr <- sort(unique(c(source.chr, target.chr)))
    u_chr <- u_chr[-1] #take off the 0 chromosome. this is the query marker
    
    #quartz(width = 12, height = 6)
    
    #query as source
    for(ch in u_chr){
        #quartz(width = 8, height = 6)
        layout.mat <- matrix(c(1,1,1,3,3,3,2,2,4,4,5,5), ncol = 2, byrow = FALSE)
        layout(layout.mat, widths = c(1,0.3))
        par(mar = c(2,4,2,2))

        target.chr.locale <- which(target.chr == ch)
        source.chr.locale <- which(source.chr == ch)

        target.chr.pos <- target.pos[target.chr.locale]
        source.chr.pos <- source.pos[source.chr.locale]
        target.chr.effect <- all.effect[target.chr.locale]
        source.chr.effect <- all.effect[source.chr.locale]

        #plot main effects
        main.chr.locale <- which(main_chr == ch)
        main.chr.pos <- main_pos[main.chr.locale]
        main.chr.effect <- main_effect[main.chr.locale]
        main.chr.col <- main_col[main.chr.locale]
        main.chr.alleles <- main_allele[main.chr.locale]

        chr.boundaries <- c(0, max(c(target.chr.pos, source.chr.pos, main.chr.pos)))
        
        plot.new()
        plot.window(xlim = chr.boundaries, ylim = main.lim)
        points(main.chr.pos, main.chr.effect, col = main.chr.col, type = "h") 
        pt.cols <- allele_cols[match(main.chr.alleles, colnames(geno_obj))]
        points(main.chr.pos, main.chr.effect, col = pt.cols, type = "p",
            pch = 16, cex = 0.5)
        abline(h = 0)       
        axis(2)
        mtext(side = 2, "Test Marker Main Effect", line = 2.5)
        mtext(side = 3, paste("Chr", u_chr[ch]))        

        par(mar = c(0,0,0,0))
        plot.new()
        plot.window(xlim = c(0,1), ylim = c(0,1))
        legend(0,0.2, col = trait_cols[1:length(u_trait)], lty = 1, lwd = 3,
            legend = u_trait, cex = 0.7)
        
        par(mar = c(4,4,0,4))
        plot.new()
        plot.window(xlim = chr.boundaries, ylim = int.lim)
        target.allele <- sapply(strsplit(var_int[target.chr.locale,"Target"], "_"), function(x) x[2])
        target.col <- allele_cols[match(target.allele, colnames(geno_obj))]
        points(target.chr.pos, target.chr.effect, type = "h", col = "#a6611a")
        points(target.chr.pos, target.chr.effect, type = "p", col = target.col, cex = 0.5, pch = 16)

        source.allele <- sapply(strsplit(var_int[source.chr.locale,"Source"], "_"), function(x) x[2])
        source.col <- allele_cols[match(source.allele, colnames(geno_obj))]
        points(source.chr.pos, source.chr.effect, type = "h", col = "#018571")
        points(source.chr.pos, source.chr.effect, type = "p", col = source.col, cex = 0.5, pch = 16)
        abline(h = 0)
        axis(2);axis(1)
        mtext(side = 2, "Interaction Effect Size", line = 2.5)
        mtext(side = 1, "Position", line = 2.5)
        
        par(mar = c(0,0,0,0))
        plot.new()
        plot.window(xlim = c(0,1), ylim = c(0,1))
        legend(x = 0, y = 0.8, legend = colnames(geno_obj), col = trait_cols, pch = 16)

        plot.new()
        plot.window(xlim = c(0,1), ylim = c(0,1))
        legend(0, 0.8, legend = c("Query as Source", "Query as Target"), 
            col = c("#a6611a", "#018571"), lty = 1, lwd = 3, cex = 0.7)
    }

}