#This function looks at whether interactions are aggravating or
#alleviating. Then based on whether high/low values of each
#phenotype is good/bad, the function decides which interactions are 
#beneficial for a patient, and which are harmful
#clinical.benefit is a vector that says whether each phenotype is
#beneficial if it is high (h), low (l), or doesn't matter (*)

clinical.benefits.query <- function(data.obj, clinical.benefit, p.or.q = 0.05, 
	covar = NULL, scan.what = c("normalized.traits", "raw.traits"), 
	geno.coding = c("Additive", "Dominant", "Recessive")){
	
    query_genotype <- data.obj$query_genotype
	geno.coding = geno.coding[1]

	  if(pheno_type == "pheno"){
        pheno <- data_obj$pheno
    }
    if(pheno_type == "norm_pheno"){
        pheno <- apply(data_obj$pheno, 2, rankZ)
    }
    if(pheno_type == "ET"){
        pheno <- data_obj$ET
    }

	covar.info <- get_covar(data.obj)
	
	if(geno.coding == "Dominant"){
		geno.obj[which(geno.obj >= 0.5)] <- 1
		}
	if(geno.coding == "Recessive"){
		geno.obj[which(geno.obj <= 0.5)] <- 0
		}

	var.inf <- write_variant_influences(data.obj, p.or.q, include_main_effects = FALSE, 
		write_file = FALSE)	
	
	if(nrow(var.inf) == 0){
		stop("There are no interactions.")
		}
		
	
	#=========================================================
	# internal functions
	#=========================================================
		
	get.genotype <- function(marker.name){
		marker.locale <- which(colnames(geno) == marker.name)
		if(length(marker.locale) == 1){
			return(geno[,marker.locale])
			}else{
			marker.locale <- which(data.obj$p.covar == marker.name)	
			if(length(marker.locale) == 1){
				return(data.obj$p.covar.table[,marker.locale])
				}
			}
		}	
		
		
	get.pheno.vals <- function(marker1, marker2, phenotype){		
		marker1.geno <- get.genotype(marker1)
		marker2.geno <- get.genotype(marker2)
	
		max1 <- which(marker1.geno == max(marker1.geno))
		max2 <- which(marker2.geno == max(marker2.geno))
		min1 <- which(marker1.geno == min(marker1.geno))
		min2 <- which(marker2.geno == min(marker2.geno))

		baseline <- mean(phenotype[intersect(min1, min2)])
		just1.pheno <- mean(phenotype[intersect(max1, min2)]) - baseline
		just2.pheno <- mean(phenotype[intersect(min1, max2)]) - baseline
		both.pheno <- mean(phenotype[intersect(max1, max2)]) - baseline
		add.exp <- just1.pheno + just2.pheno
		
		if(is.finite(add.exp)){
			if(add.exp < 0){
				int.expect <- "expect.negative"
				}else{
				int.expect <- "expect.positive"
				}
			
			if(abs(both.pheno) < abs(add.exp)){
				int.effect <- "alleviating"
				}else{
				int.effect <- "aggravating"
				}
			}else{
				int.expect <- "none"
				int.effect <- "none"
				}
		
		pheno.result <- c(just1.pheno, just2.pheno, add.exp, both.pheno, int.expect, int.effect)
		return(pheno.result)
		}


	#find which effects on phenotypes are beneficial
	#and which are harmful
	evaluate.effect <- function(result.mat){
		clin.effects <- rep(NA, nrow(result.mat))
		
		high.locale <- which(clinical.benefit == "h")
		low.locale <- which(clinical.benefit == "l")

		diff.expect <- as.numeric(result.mat[,"actual"]) - as.numeric(result.mat[,"expected.additive"])
		higher.than.expected <- which(diff.expect > 0)
		lower.than.expected <- which(diff.expect < 0)		
				
		#==============
		#beneficial
		#==============
		#for traits that we want high, interactions are
		#beneficial if the actual value is higher than 
		#the expected value
		ben.high <- intersect(high.locale, higher.than.expected)
		if(length(ben.high) > 0){
			clin.effects[ben.high] <- "beneficial"
			}
	
		#for traits that we want low, interactions
		#are beneficial if the actual value is lower
		#than the expected value
		ben.low <- intersect(low.locale, lower.than.expected)
		if(length(ben.low) > 0){
			clin.effects[ben.low] <- "beneficial"
			}
		
		#==============
		#harmful
		#==============
		#for traits that we want high, interactions are
		#harmful if the actual value is lower than 
		#the expected value
		harm.high <- intersect(high.locale, lower.than.expected)
		if(length(harm.high) > 0){
			clin.effects[harm.high] <- "harmful"
			}
	
		#for traits that we want low, interactions
		#are harmful if the actual value is higher than
		#expected
		harm.low <- intersect(low.locale, higher.than.expected)
		if(length(harm.low) > 0){
			clin.effects[harm.low] <- "harmful"
			}
		
		return(clin.effects)
		
		}
	
	scan.effects <- function(full.effects.mat){
		all.effects <- full.effects.mat[,"clinical"]
		all.effects <- all.effects[which(!is.na(all.effects))]
		u_effects <- unique(all.effects)
		if(length(u_effects) == 1){
			return(u_effects)
			}else{
			return("mixed")
			}
		}
	#=========================================================

	all.int.effects <- vector(mode = "list", length = nrow(var.inf))
	int.names <- apply(var.inf, 1, function(x) paste(x[1], x[4], sep = "_"))
	names(all.int.effects) <- int.names
	for(i in 1:nrow(var.inf)){
		marker1 <- var.inf[i,1]
		marker2 <- var.inf[i,4]		

		# plot.effects(data.obj, marker1, marker2, plot.type = "b", error.bars = "se", covar = covar)

		#find phenotype values for individuals with different genotype combinations
		#and 
		all.pheno.effects <- matrix(NA, nrow = ncol(pheno), ncol = 6)
		for(ph in 1:ncol(pheno)){
			all.pheno.effects[ph,] <- get.pheno.vals(marker1, marker2, phenotype = pheno[,ph])	
			}
		rownames(all.pheno.effects) <- colnames(pheno)
		colnames(all.pheno.effects) <- c("marker1", "marker2", "expected.additive", "actual", "expected.effect", "interaction.effect")
		clin <- evaluate.effect(all.pheno.effects)

		all.pheno.effects <- cbind(all.pheno.effects, clin)
		colnames(all.pheno.effects)[7] <- "clinical"
		
		all.int.effects[[i]] <- all.pheno.effects
	
		}	

	clin.int <- unlist(lapply(all.int.effects, scan.effects))
	final.clin <- cbind(var.inf[,1], var.inf[,4], clin.int)
	colnames(final.clin) <- c("Source", "Target", "Clinical_Effect")
	rownames(final.clin) <- NULL

	final.results <- list("clinical_summary" = final.clin, "details" = all.int.effects)
	return(final.results)
	}