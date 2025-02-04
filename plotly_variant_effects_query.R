plotly_variant_effects_query <- function(data_obj, geno_obj, 
    pheno_type = c("pheno", "norm_pheno", "ET"), 
    p_or_q = 0.05, scale_coord = 1, gene.table = NULL, 
    gene.bp.window = 5000, verbose = FALSE){

    query_genotype <- data_obj$query_genotype
    allele_cols <- categorical_pal(ncol(geno_obj))
    names(allele_cols) <- 1:ncol(geno_obj)
    pheno_cols <- categorical_pal(ncol(data_obj$pheno))
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

    num.pheno <- ncol(pheno)
    pheno.names <- colnames(pheno)
    names(pheno_cols) <- pheno.names
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
        geno1 <- round(geno1, 2)
        geno2 <- round(geno2, 2)
        geno_pairs <- pair.matrix(c(0,1), self.pairs = TRUE, ordered = TRUE)
        pheno.groups <- apply(geno_pairs, 1, function(x) phenoV[intersect(which(geno1 == x[1]), which(geno2 == x[2]))])
        #boxplot(pheno.groups);abline(h = 0)
        group.means <- sapply(pheno.groups, function(x) mean(x, na.rm = TRUE))
        group.se <- sapply(pheno.groups, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
        ref.effect <- group.means[1]
        centered.means <- group.means - ref.effect
        source.effect <- centered.means[4] #difference between 0/0 and 1/0
        target.effect <- centered.means[2] #difference between 0/0 and 0/1
        
        add.pred <- source.effect + target.effect
        add.error <- sum(group.se[2:3])
        add.range <- c((add.pred - add.error), (add.pred + add.error))

        result <- c("source.effect" = source.effect, "target.effect" = target.effect, 
            "additive" = add.pred, "additive_range" = add.range, "actual" = centered.means[3])
        #barplot(result)

    
        return(result)
    }


    get_nearby_genes <- function(chr, pos, bp.window = 1000){
        chr.locale <- which(as.numeric(gene.table[,1]) == chr)
        chr.table <- gene.table[chr.locale,]
        above.min <- which(as.numeric(chr.table[,"start"]) >= pos-bp.window)
        below.max <- which(as.numeric(chr.table[,"stop"]) <= pos+bp.window)
        gene.idx <- intersect(above.min, below.max)
        window.genes <- chr.table[gene.idx,,drop=FALSE]
        gene.names <- window.genes[,"feature"]
        name.locale <- which(window.genes[,"gene"] != "")
        if(length(name.locale) > 0){
            gene.names[name.locale] <- window.genes[name.locale,"gene"]
        }
        return(gene.names)
    }

    source.deviation <- target.deviation <- vector(mode = "list", length = length(u_chr))
    names(source.deviation) <- names(target.deviation) <- u_chr

    source.ch.pos <- target.ch.pos <- vector(mode = "list", length = length(u_chr))
    names(source.ch.pos) <- names(target.ch.pos) <- u_chr
    
    source.allele <- target.allele <- vector(mode = "list", length = length(u_chr))
    names(source.allele) <- names(target.allele) <- u_chr

    source.gene.names <- target.gene.names <- vector(mode = "list", length = length(u_chr))
    names(source.gene.names) <- names(target.gene.names) <- u_chr

    fig.list <- vector(mode = "list", length = length(u_chr))

    for(ch in u_chr){
        if(verbose){report.progress(ch, length(u_chr))}
        target.chr.locale <- which(target.chr == ch)
        source.chr.locale <- which(source.chr == ch)

        if(length(target.chr.locale) > 0){
            target.idx <- which(u_chr == ch)
            target.markers <- var_int[target.chr.locale,"Target",drop=FALSE]
            split.target <- strsplit(target.markers, "_")
            target.marker.names <- sapply(split.target, function(x) x[1])
            target.allele.names <- sapply(split.target, function(x) x[2])
            target.allele[[target.idx]] <- target.allele.names
            target.geno <- sapply(1:length(target.marker.names), function(x) matched.geno[,target.allele.names[x], target.marker.names[x],drop=FALSE])
            colnames(target.geno) <- target.markers
            target.int <- lapply(1:ncol(matched.pheno), 
                function(y) t(apply(target.geno, 2, 
                function(x) get_int(phenoV = matched.pheno[,y], geno1 = matched.query, geno2 = x))))
            target.ch.pos[[target.idx]] <- target.pos[target.chr.locale]
            target.add.min <- sapply(target.int, function(x) x[,"additive_range1"])
            target.add.max <- sapply(target.int, function(x) x[,"additive_range2"])
            target.actual <- sapply(target.int, function(x) x[,"actual"])
            
            #idx = 1
            #test.mat <- rbind(target.add.min[,idx], target.actual[,idx], target.add.max[,idx])
            #barplot(test.mat, beside = TRUE, col = c("lightblue", "orange", "lightgreen"), main = colnames(matched.pheno)[idx])

            target.in.range <- intersect(which(target.actual > target.add.min), which(target.actual < target.add.max))

            target.effects <- sapply(target.int, function(x) x[,"actual"] - x[,"additive"])
            if(is.null(dim(target.effects))){
                target.effects <- matrix(target.effects, nrow = 1)
            }
                
            target.effects[target.in.range] <- 0 #zero out effects that are within the additive error
            
            colnames(target.effects) <- colnames(pheno)
            target.deviation[[target.idx]] <- target.effects

            if(!is.null(gene.table)){
            target.nearest.gene <- sapply(1:length(target.chr.locale), 
            function(x) get_nearby_genes(ch, (target.ch.pos[[target.idx]][x]*scale_coord), 
            gene.bp.window))
            target.gene.names[[target.idx]] <- target.nearest.gene
        }
        }

        if(length(source.chr.locale) > 0){
            source.idx <- which(u_chr == ch)
            source.markers <- var_int[source.chr.locale,"Source",drop=FALSE]        
            split.source <- strsplit(source.markers, "_")
            source.marker.names <- sapply(split.source, function(x) x[1])
            source.allele.names <- sapply(split.source, function(x) x[2])

            source.allele[[source.idx]] <- source.allele.names

            source.geno <- sapply(1:length(source.marker.names), function(x) matched.geno[,source.allele.names[x], source.marker.names[x],drop=FALSE])
            colnames(source.geno) <- source.markers    
            source.int <- lapply(1:ncol(matched.pheno), 
                function(y) t(apply(source.geno, 2, 
                function(x) get_int(matched.pheno[,y], x, matched.query))))

            source.ch.pos[[source.idx]] <- source.pos[source.chr.locale]
            source.add.min <- sapply(source.int, function(x) x[,"additive_range1"])
            source.add.max <- sapply(source.int, function(x) x[,"additive_range2"])
            source.actual <- sapply(source.int, function(x) x[,"actual"])
            
            #idx = 1
            #test.mat <- rbind(source.add.min[,idx], source.actual[,idx], source.add.max[,idx])
            #barplot(test.mat, beside = TRUE, col = c("lightblue", "orange", "lightgreen"), main = colnames(matched.pheno)[idx])

            source.in.range <- intersect(which(source.actual > source.add.min), which(source.actual < source.add.max))

            source.effects <- sapply(source.int, function(x) x[,"actual"] - x[,"additive"])
            if(is.null(dim(source.effects))){
                source.effects <- matrix(source.effects, nrow = 1)
            }

            source.effects[source.in.range] <- 0 #zero out effects that are within the additive error

            colnames(source.effects) <- colnames(pheno)
            source.deviation[[source.idx]] <- source.effects

            if(!is.null(gene.table)){
                source.nearest.gene <- sapply(1:length(source.chr.locale), 
                    function(x) get_nearby_genes(ch, (source.ch.pos[[source.idx]][x]*scale_coord), 
                    gene.bp.window))
                source.gene.names[[source.idx]] <- source.nearest.gene
            }

        }
    }

    generate_fig <- function(plot.pos, plot.dev, plot.allele, plot.gene.names, plot.label,
        xlim, ylim){
        
        if(length(plot.pos) == 0){
            fig <- plot_ly(x = segment_region(xlim[1], xlim[2], 100), 
            y = segment_region(ylim[1], ylim[2], 100), colors = allele_cols, 
            mode = "markers", type = "scatter", visible = FALSE)

        }else{
            fig <- plot_ly() %>%
            add_segments(x = plot.pos, y = plot.dev[,pheno.names[1]], 
            xend = plot.pos, yend = 0, colors = allele_cols, 
            name = pheno.names[1], color = I(pheno_cols[1])) %>%
            add_markers(x = plot.pos, y = plot.dev[,pheno.names[1]], 
                type = "scatter", color = ~plot.allele, mode = "markers",
                text = plot.gene.names, size = 1.5)

            for(ph in 2:num.pheno){
                fig <- add_segments(fig, x = plot.pos, y = plot.dev[,pheno.names[ph]], 
                    xend = plot.pos, yend = 0, name = pheno.names[ph], color = I(pheno_cols[ph]))
                fig <- add_trace(fig, type = "scatter", mode = "markers", 
                    x = plot.pos, y = plot.dev[,pheno.names[ph]], size = 1.5,
                    color = plot.allele, text = plot.gene.names, showlegend = FALSE)
                }
            }
        
        #add labels
        fig <- layout(fig, title = plot.label, 
            plot_bgcolor = "white", 
            xaxis = list(title = coord_label), 
            yaxis = list(title = 'Deviation of Trait from Additive'))
        
        return(fig)

    }

    #find ylim across all chromosomes
    ylim = c(min(c(min(unlist(target.deviation), na.rm = TRUE),min(unlist(source.deviation), na.rm = TRUE))),
           max(c(max(unlist(target.deviation), na.rm = TRUE),max(unlist(source.deviation), na.rm = TRUE))))


    for(ch in 1:length(source.deviation)){
        #quartz(width = 10, height = 6)
        if(is.null(source.ch.pos[[ch]])){next()}
        xlim <- c(1, max(c(source.ch.pos[[ch]], target.ch.pos[[ch]])))

        #create the scatter plot object with gene names as the hover text
        source.x <- source.ch.pos[[ch]]
        source.y <- source.deviation[[ch]]
        source.allele.labels <- as.factor(source.allele[[ch]])
        source.gene.labels <- sapply(source.gene.names[[ch]], function(x) paste(unlist(x), collapse = ", "))

        source.fig <- generate_fig(plot.pos = source.x, plot.dev = source.y, 
            plot.allele = source.allele.labels, 
            plot.gene.names = source.gene.labels, plot.label = paste0("Chr", ch),
            xlim, ylim)        

        #create the scatter plot object for the query as target with gene names as the hover text
        target.x <- target.ch.pos[[ch]]
        target.y <- target.deviation[[ch]]
        target.allele.labels <- target.allele[[ch]]
        target.gene.labels <- sapply(target.gene.names[[ch]], function(x) paste(unlist(x), collapse = ", "))

        target.fig <- generate_fig(target.x, target.y, target.allele.labels, 
            target.gene.labels, paste0("Chr", ch),
            xlim, ylim)
        

        annotations <- list(list( 
            x = 0.5,  
            y = 1.0,  
            text = "Query as Source",  
            xref = "paper",  
            yref = "paper",  
            xanchor = "center",  
            yanchor = "top",  
            showarrow = FALSE 
        ),  
        list( 
            x = 0.5,  
            y = 0.4,  
            text = "Query as Target",  
            xref = "paper",  
            yref = "paper",  
            xanchor = "center",  
            yanchor = "top",  
            showarrow = FALSE 
        ))

        full.fig <- subplot(source.fig, style(target.fig, showlegend = FALSE), 
        nrows = 2, shareX = TRUE, shareY = TRUE, titleX = TRUE, titleY = TRUE) %>%
        layout(annotations = annotations)

        fig.list[[ch]] <- full.fig
        

    }
    
    invisible(fig.list)

}
