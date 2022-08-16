plot_pairscan_query <- function(data_obj, pairscan_obj){

    pairscan_results <- pairscan_obj$pairscan_results
    pairscan_perm <- pairscan_obj$pairscan_perm

    all_chr <- data_obj$chromosome[which(data_obj$chromosome != 0)]
    chr <- unique(all_chr)
    all_markers <- data_obj$geno_names[[3]]

    trait_markers <- sapply(strsplit(pairscan_results[[1]][[1]][,2], "_"), function(y) y[1])
    all_query_effects <- lapply(pairscan_results, function(x) as.numeric(x[[1]][,3]))
    all_test_effects <- lapply(pairscan_results, function(x) as.numeric(x[[1]][,4]))
    all_int_effects <- lapply(pairscan_results, function(x) as.numeric(x[[1]][,5]))
    
    all_perm_query_effects <- lapply(pairscan_perm, function(x) as.numeric(x[[1]][,3]))
    all_perm_test_effects <- lapply(pairscan_perm, function(x) as.numeric(x[[1]][,4]))
    all_perm_int_effects <- lapply(pairscan_perm, function(x) as.numeric(x[[1]][,5]))

    query_conf <- lapply(all_perm_query_effects, 
        function(x) c(get.percentile(x, 1), get.percentile(x, 99)))

    test_conf <- lapply(all_perm_test_effects, 
        function(x) c(get.percentile(x, 1), get.percentile(x, 99)))

    int_conf <- lapply(all_perm_int_effects, 
        function(x) c(get.percentile(x, 1), get.percentile(x, 99)))

    #boxplot(all_test_effects);abline(h = c(test_conf[[1]]))
    #boxplot(all_query_effects);abline(h = c(query_conf[[2]]))

    effect_lim <- c(min(c(unlist(all_test_effects)), unlist(all_query_effects)), 
        max(c(unlist(all_test_effects), unlist(all_query_effects))))
    int_effect_lim <- c(min(unlist(all_int_effects)), max(unlist(all_int_effects)))
    for(ch in chr){
        #quartz(width = 8, height = 6)
        chr_locale <- which(all_chr == ch)
        chr_markers <- all_markers[chr_locale]
        marker_effect_locale <- match(chr_markers, trait_markers)
        marker_pos_locale <- match(chr_markers, all_markers)
        marker_pos <- data_obj$marker_location[marker_pos_locale]

        for(p in 1:length(pairscan_results)){
            layout(matrix(c(1,1,1,2,2,2,3,4,5), nrow = 3, byrow = TRUE))
            plot(marker_pos, all_query_effects[[p]][marker_effect_locale], type = "h", 
                ylim = effect_lim, ylab = "Effect Size", xlab = "Genomic Position",
                main = paste("Main Effects Chr", ch, names(pairscan_results)[p]))
            points(marker_pos, all_test_effects[[p]][marker_effect_locale], 
                col = "red", type = "h")
            abline(h = c(0, test_conf[[p]]))
            abline(h = query_conf[[p]], col = "red")

            plot(marker_pos, all_int_effects[[p]][marker_effect_locale], 
                type = "h", ylab = "Interaction Effect Size", 
                xlab = "Genomic Position", ylim = int_effect_lim,
                main = paste("Interaction Effects\nChr", ch, names(pairscan_results)[p]))
            abline(h = 0)
            abline(h = int_conf[[p]])
            
            plot.hexbin.as.plot(all_query_effects[[p]], all_int_effects[[p]], 
                xlab = "Query Main Effects", ylab = "Interaction Effects", 
                main = "Query vs. Interaction", min.cex = 1, max.cex = 3,
                legend.pos = "topright", round.legend = 500)

            plot.hexbin.as.plot(all_test_effects[[p]], all_int_effects[[p]], 
                xlab = "Test Main Effects", ylab = "Interaction Effects", 
                main = "Test vs. Interaction", min.cex = 1, max.cex = 3,
                legend.pos = "topright", round.legend = 500)

            plot.hexbin.as.plot(all_query_effects[[p]], all_test_effects[[p]], 
                xlab = "Query Main Effects", ylab = "Test Main Effects", 
                main = "Query vs. Test", min.cex = 1, max.cex = 3,
                legend.pos = "topright", round.legend = 500)

        #    boxplot(list(all_query_effects[[p]][marker_effect_locale], 
        #        all_test_effects[[p]][marker_effect_locale], 
        #        all_int_effects[[p]][marker_effect_locale]), 
        #        names = c("Query", "Test", "Interaction"))
        #    abline(h = 0)
        }
    }

}