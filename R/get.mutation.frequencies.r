get.mutation.frequencies <- function(xmut.ids, tcga.cancer.type=NULL, tcga.cancer.type.2=NULL,reference.data=NULL,reference.genes=NULL,combine.with.TCGA=FALSE ) {
  # library(data.table)
  # library(variantprobs)
  # library(dplyr)
  # library(dtplyr)
  # library(tibble)
  
  # print("Running test")
  library(grDevices)
  library(graphics)
  library(stats)
  library(utils)
  library(data.table)
  library(variantprobs)
  library(dplyr)
  library(dtplyr)
  library(tibble)
  
  if(!is.null(reference.data))
    reference.data <- as_tibble(reference.data)
  if(!is.null(reference.genes))
    reference.genes <- as_tibble(reference.genes)
  
  xmut.ids <- gsub("X",23,xmut.ids)
  plus1<-TRUE
  # if (! (all(substr(xmut.ids,1,1) %in% c(c(1:9),"X","Y") | substr(xmut.ids,1,2) %in% c(c(10:22)) )))
  #   stop("xmut.ids should be of the following format: {Chromosome Location RefAllele AltAllele}, each entry separated by space, \n where chromosome is a number 1-22 or X or Y; \n location is genomic location in  GRCh37 build; \n RefAllele is a reference allele and AltAllele is Alternative allele/Tumor_Seq_Allele2. \n For example '10 100003849 G A', is the mutation at chromosome 10, genomic location 100003849, where reference allele G is substituted with A, or \n '10 100011448 - CCGCTGCAAT' is the insertion of 'CCGCTGCAAT' at chromosome 10, location 100011448. \n The ref and alt alleles follow standard TCGA maf file notations.")
  
  if (is.null(tcga.cancer.type) & is.null(reference.data)) stop("You either need to specify 'tcga.cancer.type' to use TCGA based frequnecies or specify 'reference.data' to get reference dataset for frequency calculation")
  if(!(tcga.cancer.type %in% names(GT_probs)) && is.null(reference.data))
    stop(paste0("The cancer type you selected is not among the ones currently supported, please see ",names(GT_probs)))
  if (combine.with.TCGA & (is.null(tcga.cancer.type) | is.null(reference.data))) stop("If you choose combine.with.TCGA=TRUE, i.e. to combine reference dataset with TCGA data, then you  need to  specify both'tcga.cancer.type'  and 'reference.data' ")
  # if (!combine.with.TCGA & (!is.null(tcga.cancer.type) & !is.null(reference.data))) warning("You choose combine.with.TCGA=FALSE, but specified both TCGA and reference data - only TCGA data will be used unless you choose combine.with.TCGA=TRUE")
  
  # match variants to their genes #
  if( is.null(reference.data) || (!is.null(reference.data) && is.null(reference.genes)) || (!is.null(reference.data) && !is.null(reference.genes) && ncol(reference.genes) == 2)){
    matched_genes <- as.data.table(do.call('rbind',lapply(xmut.ids, function(x){
      # print(x)
      temp <- strsplit(x, " ")[[1]]
      chr <- as.numeric(temp[1])
      pos <- as.numeric(temp[2])
      # print("test")
      gene <- d.gene$Hugo_Symbol[d.gene$Chromosome == chr & pos > d.gene$Start & pos < d.gene$End]
      if(length(gene)==0)
        gene <- "MissingGeneName"
      return(c(x,gene[1],NA,NA))
    })))
    colnames(matched_genes) <- c("variant","gene","freq","aggregated_estimate")
  }
  
  
  if(is.null(reference.data)){
    # get estimates for GT at that sites #
    site_spe_gt <- as.data.table(GT_probs[[tcga.cancer.type]])
    site_spe_gt <- site_spe_gt[site_spe_gt$warning ==0,]
    site_spe_agg_gt <- as.data.table(GT_agg_probs[[tcga.cancer.type]])
    # get estimates per genes found in the matching #
    genes <- unique(matched_genes$gene)
    
    for(i in 1:length(genes)){
      gene_temp <- genes[i]
      variants <- matched_genes$variant[matched_genes$gene == gene_temp]
      site_spe_gt_temp <- site_spe_gt[site_spe_gt$Gene == gene_temp,]
      if(nrow(site_spe_gt_temp) > 0){
        variants_temp <- ifelse(variants %in% site_spe_gt_temp$Variant, variants, "each_unseen")
        matched_genes[match(variants,matched_genes$variant) , 3] <- as.numeric(site_spe_gt_temp$GT[match(variants_temp,site_spe_gt_temp$Variant)])
        matched_genes[match(variants,matched_genes$variant) , 4] <- FALSE
      }
      else{
        variants_temp <- ifelse(variants %in% site_spe_agg_gt$Variant, variants, "each_unseen")
        matched_genes[match(variants,matched_genes$variant) , 3] <- as.numeric(site_spe_agg_gt$GT_agg[match(variants_temp,site_spe_agg_gt$Variant)])
        matched_genes[match(variants,matched_genes$variant) , 4] <- TRUE
      }
    }
    
    if(!is.null(tcga.cancer.type.2)){
      # get new estimates of tcga #
      site_spe_gt.2 <- as.data.table(GT_probs[[tcga.cancer.type.2]])
      site_spe_gt.2 <- site_spe_gt.2[site_spe_gt.2$warning ==0,]
      site_spe_agg_gt.2 <- as.data.table(GT_agg_probs[[tcga.cancer.type.2]])
      
      matched_genes.2 <- data.frame(matrix(nrow = nrow(matched_genes), ncol = 2))
      for(i in 1:length(genes)){
        gene_temp <- genes[i]
        variants <- matched_genes$variant[matched_genes$gene == gene_temp]
        site_spe_gt_temp.2 <- site_spe_gt.2[site_spe_gt.2$Gene == gene_temp,]
        if(nrow(site_spe_gt_temp.2) > 0){
          variants_temp <- ifelse(variants %in% site_spe_gt_temp.2$Variant, variants, "each_unseen")
          matched_genes.2[match(variants,matched_genes$variant) , 1] <- as.numeric(site_spe_gt_temp.2$GT[match(variants_temp,site_spe_gt_temp.2$Variant)])
          matched_genes.2[match(variants,matched_genes$variant) , 2] <- FALSE
        }
        else{
          variants_temp <- ifelse(variants %in% site_spe_agg_gt.2$Variant, variants, "each_unseen")
          matched_genes.2[match(variants,matched_genes$variant) , 1] <- as.numeric(site_spe_agg_gt.2$GT_agg[match(variants_temp,site_spe_agg_gt.2$Variant)])
          matched_genes.2[match(variants,matched_genes$variant) , 2] <- TRUE
        }
      }
      colnames(matched_genes.2) <- c("freq.2","aggregated_estimate.2")
      out <- cbind(matched_genes,matched_genes.2)
      
      ### fix aggregation mismatch ###
      for(i in 1:nrow(out)){
        if(sum(as.logical(out$aggregated_estimate[i]), out$aggregated_estimate.2[i]) == 1){
          # if first site has true --> check how many variants #
          if(!as.logical(out$aggregated_estimate[i])){
            temp <- site_spe_gt[site_spe_gt$Gene == out$gene[i],]
            temp <- temp[-which(temp$Variant %in% c("atleast_1new","each_unseen")),]
            
            agg_vals_temp <- sort(unique(site_spe_agg_gt$GT_agg[-which(site_spe_agg_gt$Variant %in% c("atleast_1new"))]))
            vals_temp <- sort(unique(site_spe_gt$GT[site_spe_gt$Variant != "atleast_1new" & site_spe_gt$Gene == out$gene[i]]))
            
            if(is.na(sort(unique(temp$GT))[3] < as.numeric(out$freq[i]))){
              rank_temp <- which.min(abs(vals_temp - as.numeric(out$freq[i])))
              out$freq[i] <- agg_vals_temp[rank_temp]
              out$aggregated_estimate[i] <- TRUE
            }
            
            if(!is.na(sort(unique(temp$GT))[3] <= as.numeric(out$freq[i]))){
              if(sort(unique(temp$GT))[3] < as.numeric(out$freq[i])){
                rank_temp <- which.min(abs(vals_temp - as.numeric(out$freq[i])))
                out$freq[i] <- agg_vals_temp[rank_temp]
                out$aggregated_estimate[i] <- TRUE
              }
            }
          }
          
          # if second site has true #
          if(!as.logical(out$aggregated_estimate.2[i])){
            temp <- site_spe_gt.2[site_spe_gt.2$Gene == out$gene[i],]
            temp <- temp[-which(temp$Variant %in% c("atleast_1new","each_unseen")),]
            
            agg_vals_temp <- sort(unique(site_spe_agg_gt.2$GT_agg[-which(site_spe_agg_gt.2$Variant %in% c("atleast_1new"))]))
            vals_temp <- sort(unique(site_spe_gt.2$GT[site_spe_gt.2$Variant != "atleast_1new" & site_spe_gt.2$Gene == out$gene[i]]))
            
            if(is.na(sort(unique(temp$GT))[3] < as.numeric(out$freq.2[i]))){
              rank_temp <- which.min(abs(vals_temp - as.numeric(out$freq.2[i])))
              out$freq.2[i] <- agg_vals_temp[rank_temp]
              out$aggregated_estimate.2[i] <- TRUE
            }
            
            if(!is.na(sort(unique(temp$GT))[3] <= as.numeric(out$freq.2[i]))){
              if(sort(unique(temp$GT))[3] > as.numeric(out$freq.2[i])){
                rank_temp <- which.min(abs(vals_temp - as.numeric(out$freq.2[i])))
                out$freq.2[i] <- agg_vals_temp[rank_temp]
                out$aggregated_estimate.2[i] <- TRUE
              }
            }
          }
        }
      }
    }
  }
  
  
  
  ##############################
  ##### REFERENCE MAF FILE #####
  ##############################
  
  if(!is.null(reference.data)){
    
    
    # if (all(colnames(reference.data)!="Chromosome")) stop("reference.data should include column 'Chromosome' that will be used to create mutation IDs")
    # if (all(colnames(reference.data)!="PatientID")) stop("reference.data should include column 'PatientID' - could be the same as sample IDs Tumor_Sample_Barcode")
    # if (all(colnames(reference.data)!="Start_Position")) stop("reference.data should include column 'Start_Position' that will be used to create mutation IDs")
    # if (all(colnames(reference.data)!="Reference_Allele")) stop("reference.data should include column 'Reference_Allele' that will be used to create mutation IDs")
    # if (all(colnames(reference.data)!="Tumor_Seq_Allele2")) stop("reference.data should include column 'Tumor_Seq_Allele2' that will be used to create mutation IDs")
    # if (all(colnames(reference.data)!="Tumor_Sample_Barcode")) stop("reference.data should include column 'Tumor_Sample_Barcode' that will be used in frequency calculation")
    # if (all(colnames(reference.data)!="Cancer_Code")) stop("reference.data should include column 'Cancer_Code' that will be used in frequency calculation")
    
    
    # no gene size file or baits #
    if(is.null(reference.genes)){
      reference.data.temp <- as_tibble(reference.data) %>%
        filter(Cancer_Code %in% c(tcga.cancer.type, tcga.cancer.type.2))
      
      sites <- unique(reference.data.temp$Cancer_Code)
      genes <- unique(reference.data.temp$Hugo_Symbol)
      GT_probs <- lapply(sites, function(y){
        
        # print(y)
        
        out <- lapply(genes, function(x){
          # TCGA estimates #
          # get frequencies for each gene #
          var_freq <- as.data.table(reference.data.temp %>% 
                                      filter(Hugo_Symbol == x, Cancer_Code == y) %>% 
                                      group_by(Variant) %>% 
                                      summarise(v_f = length(unique(patient_id))))
          
          v_f <- var_freq$v_f
          names(v_f) <-var_freq$Variant
          # sample size #
          m <- length(unique(reference.data.temp$patient_id[reference.data.temp$Cancer_Code == y]))
          N <- sum(3*d.gene$Exome_Size[d.gene$Hugo_Symbol == x])
          N0 = N -length(v_f)
          # run Good-Turing #
          GT <- goodturing_probs(counts = v_f,m=m,N=N, N0=N0)
          GT_test <- tryCatch(goodturing_probs(counts = v_f,m=m,N=N, N0=N0),
                              error=function(e) e, warning=function(w) w)
          # keep warned genes #
          warning_found_tcga = 0
          if(is(GT_test,"warning") )
            warning_found_tcga = 1
          
          out <- as.data.frame(GT) %>%
            rownames_to_column("Variant") %>%
            mutate(Cancer_Code = y,
                   warning = warning_found_tcga#,
                   #Gene = x
            )
          out$Gene <- x
          return(out)
        })
        
        estimates <- data.table::rbindlist(out, fill = TRUE)
        return(estimates)
      })
      names(GT_probs) <- sites
      # print(colnames(GT_probs[[1]]))
      
      ################################################
      # for each site combine genes that were warned #
      GT_agg_probs <- list()
      GT_agg_N0 <- list()
      for(i in 1:length(GT_probs)){
        # print(i)
        temp <- GT_probs[[i]]
        # find all genes that were warned #
        genes_agg <- as.character(unlist(unique(temp[temp$warning == 1, "Gene"])))
        
        # get frequencies for each gene #
        var_freq <- as.data.table(reference.data.temp %>% 
                                    filter(Hugo_Symbol %in% genes_agg, Cancer_Code == names(GT_probs)[i]) %>% 
                                    group_by(Variant) %>% 
                                    summarise(v_f = length(unique(patient_id))))
        
        v_f <- var_freq$v_f
        names(v_f) <-var_freq$Variant
        # sample size #
        m <- length(unique(reference.data.temp$patient_id[reference.data.temp$Cancer_Code == names(GT_probs)[i]]))
        N <- sum(3*d.gene$Exome_Size[d.gene$Hugo_Symbol %in% genes_agg])
        N0 = N -length(v_f)
        # run Good-Turing #
        GT_agg <- goodturing_probs(counts = v_f,m=m,N=N, N0=N0)
        GT_agg <- as.data.table(as.data.frame(GT_agg) %>% 
                                  rownames_to_column("Variant") %>% 
                                  mutate(Cancer_Code = names(GT_probs)[i]
                                  ))
        GT_agg_probs[[i]] <- GT_agg
        GT_probs[[i]] <- GT_probs[[i]][GT_probs[[i]]$warning == 0,]
        GT_agg_N0[[i]] <- c(N0,N)
      }
      names(GT_agg_N0) <- names(GT_probs)
      names(GT_agg_probs) <- names(GT_probs)
      # print(colnames(GT_agg_probs[[1]]))
    }
    # if file is provided #
    else{
      
      # gene sizes #
      if(ncol(reference.genes) == 2){
        
        
        
        reference.data.temp <- reference.data %>%
          filter(Cancer_Code %in% c(tcga.cancer.type, tcga.cancer.type.2))
        
        sites <- unique(reference.data.temp$Cancer_Code)
        genes <- unique(reference.data.temp$Hugo_Symbol)
        GT_probs <- lapply(sites, function(y){
          
          # print(y)
          
          out <- lapply(genes, function(x){
            # TCGA estimates #
            # get frequencies for each gene #
            var_freq <- as.data.table(reference.data.temp %>% 
                                        filter(Hugo_Symbol == x, Cancer_Code == y) %>% 
                                        group_by(Variant) %>% 
                                        summarise(v_f = length(unique(patient_id))))
            
            v_f <- var_freq$v_f
            names(v_f) <-var_freq$Variant
            # sample size #
            m <- length(unique(reference.data.temp$patient_id[reference.data.temp$Cancer_Code == y]))
            N <- 3*as.numeric(unlist(reference.genes[reference.genes$Hugo_Symbol == x, 2])) 
            N0 = N -length(v_f)
            # run Good-Turing #
            GT <- goodturing_probs(counts = v_f,m=m,N=N, N0=N0)
            GT_test <- tryCatch(goodturing_probs(counts = v_f,m=m,N=N, N0=N0),
                                error=function(e) e, warning=function(w) w)
            # keep warned genes #
            warning_found_tcga = 0
            if(is(GT_test,"warning") )
              warning_found_tcga = 1
            
            out <- as.data.frame(GT) %>%
              rownames_to_column("Variant") %>%
              mutate(Cancer_Code = y,
                     warning = warning_found_tcga
              )
            out$Gene <- x
            return(out)
          })
          
          estimates <- data.table::rbindlist(out, fill = TRUE)
          return(estimates)
        })
        names(GT_probs) <- sites
        # print(colnames(GT_probs[[1]]))
        
        ################################################
        # for each site combine genes that were warned #
        GT_agg_probs <- list()
        GT_agg_N0 <- list()
        for(i in 1:length(GT_probs)){
          # print(i)
          temp <- GT_probs[[i]]
          # find all genes that were warned #
          genes_agg <- as.character(unlist(unique(temp[temp$warning == 1, "Gene"])))
          
          # get frequencies for each gene #
          var_freq <- as.data.table(reference.data.temp %>% 
                                      filter(Hugo_Symbol %in% genes_agg, Cancer_Code == names(GT_probs)[i]) %>% 
                                      group_by(Variant) %>% 
                                      summarise(v_f = length(unique(patient_id))))
          
          v_f <- var_freq$v_f
          names(v_f) <-var_freq$Variant
          # sample size #
          m <- length(unique(reference.data.temp$patient_id[reference.data.temp$Cancer_Code == names(GT_probs)[i]]))
          N <- N <- 3*sum(as.numeric(unlist(reference.genes[reference.genes$Hugo_Symbol %in% genes_agg, 2])))
          N0 = N -length(v_f)
          # run Good-Turing #
          GT_agg <- goodturing_probs(counts = v_f,m=m,N=N, N0=N0)
          GT_agg <- as.data.table(as.data.frame(GT_agg) %>% 
                                    rownames_to_column("Variant") %>% 
                                    mutate(Cancer_Code = names(GT_probs)[i]
                                    ))
          GT_agg_probs[[i]] <- GT_agg
          GT_probs[[i]] <- GT_probs[[i]][GT_probs[[i]]$warning == 0,]
          GT_agg_N0[[i]] <- c(N0,N)
        }
        names(GT_agg_N0) <- names(GT_probs)
        names(GT_agg_probs) <- names(GT_probs)
        # print(colnames(GT_agg_probs[[1]]))
      }
      
      
      # baits #
      if(ncol(reference.genes) > 2){
        
        # add sizes to all genes #
        reference.genes <- as.data.table(reference.genes %>% 
                                           group_by(Hugo_Symbol) %>% 
                                           mutate(Exome_Size = sum(End - Start + 1)) %>% 
                                           ungroup())
        
        
        matched_genes <- as.data.table(do.call('rbind',lapply(xmut.ids, function(x){
          temp <- strsplit(x, " ")[[1]]
          chr <- as.numeric(temp[1])
          pos <- as.numeric(temp[2])
          gene <- reference.genes$Hugo_Symbol[reference.genes$Chromosome == chr & pos > reference.genes$Start & pos < reference.genes$End]
          if(length(gene)==0)
            gene <- "MissingGeneName"
          return(c(x,gene[1],NA,NA))
        })))
        colnames(matched_genes) <- c("variant","gene","freq","aggregated_estimate")
        
        ###########################
        
        reference.data.temp <- as_tibble(reference.data) %>%
          filter(Cancer_Code %in% c(tcga.cancer.type, tcga.cancer.type.2))
        
        sites <- unique(reference.data.temp$Cancer_Code)
        genes <- unique(reference.data.temp$Hugo_Symbol)
        GT_probs <- lapply(sites, function(y){
          
          # print(y)
          
          out <- lapply(genes, function(x){
            
            # TCGA estimates #
            # get frequencies for each gene #
            var_freq <- as.data.table(reference.data.temp %>% 
                                        filter(Hugo_Symbol == x, Cancer_Code == y) %>% 
                                        group_by(Variant) %>% 
                                        summarise(v_f = length(unique(patient_id))))
            
            v_f <- var_freq$v_f
            names(v_f) <-var_freq$Variant
            # sample size #
            m <- length(unique(reference.data.temp$patient_id[reference.data.temp$Cancer_Code == y]))
            N <- 3*as.numeric(unlist(reference.genes[reference.genes$Hugo_Symbol == x,"Exome_Size"])) 
            N0 = N -length(v_f)
            # run Good-Turing #
            GT <- goodturing_probs(counts = v_f,m=m,N=N, N0=N0)
            GT_test <- tryCatch(goodturing_probs(counts = v_f,m=m,N=N, N0=N0),
                                error=function(e) e, warning=function(w) w)
            # keep warned genes #
            warning_found_tcga = 0
            if(is(GT_test,"warning") )
              warning_found_tcga = 1
            
            out <- as.data.frame(GT) %>%
              rownames_to_column("Variant") %>%
              mutate(Cancer_Code = y,
                     warning = warning_found_tcga
              )
            out$Gene = x
            return(out)
          })
          
          estimates <- data.table::rbindlist(out, fill = TRUE)
          return(estimates)
        })
        names(GT_probs) <- sites
        # print(colnames(GT_probs[[1]]))
        
        
        ################################################
        # for each site combine genes that were warned #
        GT_agg_probs <- list()
        GT_agg_N0 <- list()
        for(i in 1:length(GT_probs)){
          # print(i)
          temp <- GT_probs[[i]]
          # find all genes that were warned #
          genes_agg <- as.character(unlist(unique(temp[temp$warning == 1, "Gene"])))
          # print(genes_agg)
          # get frequencies for each gene #
          var_freq <- as.data.table(reference.data.temp %>% 
                                      filter(Hugo_Symbol %in% genes_agg, Cancer_Code == sites[i]) %>% 
                                      group_by(Variant) %>% 
                                      summarise(v_f = length(unique(patient_id))))
          
          var_freq <- as.data.table(reference.data.temp %>% 
                                      filter(Hugo_Symbol %in% genes_agg, Cancer_Code == names(GT_probs)[i]) %>% 
                                      group_by(Variant) %>% 
                                      summarise(v_f = length(unique(patient_id))))
          
          # print(nrow(var_freq))
          v_f <- var_freq$v_f
          names(v_f) <-var_freq$Variant
          # sample size #
          m <- length(unique(reference.data.temp$patient_id[reference.data.temp$Cancer_Code == names(GT_probs)[i]]))
          N <- 3*sum(as.numeric(unlist(reference.genes[reference.genes$Hugo_Symbol %in% genes_agg, "Exome_Size"])))
          # print(N)
          N0 = N -length(v_f)
          # run Good-Turing #
          GT_agg <- goodturing_probs(counts = v_f,m=m,N=N, N0=N0)
          GT_agg <- as.data.table(as.data.frame(GT_agg) %>% 
                                    rownames_to_column("Variant") %>% 
                                    mutate(Cancer_Code = names(GT_probs)[i]
                                    ))
          GT_agg_probs[[i]] <- GT_agg
          GT_probs[[i]] <- GT_probs[[i]][GT_probs[[i]]$warning == 0,]
          GT_agg_N0[[i]] <- c(N0,N)
        }
        names(GT_agg_N0) <- names(GT_probs)
        names(GT_agg_probs) <- names(GT_probs)
        # print(colnames(GT_agg_probs[[1]]))
      }
      
    }
    
    
    ###################################################
    # get estimates for GT at that sites #
    site_spe_gt <- as.data.table(GT_probs[[tcga.cancer.type]])
    site_spe_gt <- site_spe_gt[site_spe_gt$warning ==0,]
    site_spe_agg_gt <- as.data.table(GT_agg_probs[[tcga.cancer.type]])
    # get estimates per genes found in the matching #
    genes <- unique(matched_genes$gene)
    # print(colnames(site_spe_gt))
    # print(colnames(site_spe_agg_gt))
    for(i in 1:length(genes)){
      gene_temp <- genes[i]
      variants <- matched_genes$variant[matched_genes$gene == gene_temp]
      site_spe_gt_temp <- site_spe_gt[site_spe_gt$Gene == gene_temp,]
      if(nrow(site_spe_gt_temp) > 0){
        variants_temp <- ifelse(variants %in% site_spe_gt_temp$Variant, variants, "each_unseen")
        matched_genes[match(variants,matched_genes$variant) , 3] <- as.numeric(site_spe_gt_temp$GT[match(variants_temp,site_spe_gt_temp$Variant)])
        matched_genes[match(variants,matched_genes$variant) , 4] <- FALSE
      }
      else{
        variants_temp <- ifelse(variants %in% site_spe_agg_gt$Variant, variants, "each_unseen")
        matched_genes[match(variants,matched_genes$variant) , 3] <- as.numeric(site_spe_agg_gt$GT_agg[match(variants_temp,site_spe_agg_gt$Variant)])
        matched_genes[match(variants,matched_genes$variant) , 4] <- TRUE
      }
    }
    
    if(!is.null(tcga.cancer.type.2)){
      # get new estimates of tcga #
      site_spe_gt.2 <- as.data.table(GT_probs[[tcga.cancer.type.2]])
      site_spe_gt.2 <- site_spe_gt.2[site_spe_gt.2$warning ==0,]
      site_spe_agg_gt.2 <- as.data.table(GT_agg_probs[[tcga.cancer.type.2]])
      
      matched_genes.2 <- data.frame(matrix(nrow = nrow(matched_genes), ncol = 2))
      for(i in 1:length(genes)){
        gene_temp <- genes[i]
        variants <- matched_genes$variant[matched_genes$gene == gene_temp]
        site_spe_gt_temp.2 <- site_spe_gt.2[site_spe_gt.2$Gene == gene_temp,]
        if(nrow(site_spe_gt_temp.2) > 0){
          variants_temp <- ifelse(variants %in% site_spe_gt_temp.2$Variant, variants, "each_unseen")
          matched_genes.2[match(variants,matched_genes$variant) , 1] <- as.numeric(site_spe_gt_temp.2$GT[match(variants_temp,site_spe_gt_temp.2$Variant)])
          matched_genes.2[match(variants,matched_genes$variant) , 2] <- FALSE
        }
        else{
          variants_temp <- ifelse(variants %in% site_spe_agg_gt.2$Variant, variants, "each_unseen")
          matched_genes.2[match(variants,matched_genes$variant) , 1] <- as.numeric(site_spe_agg_gt.2$GT_agg[match(variants_temp,site_spe_agg_gt.2$Variant)])
          matched_genes.2[match(variants,matched_genes$variant) , 2] <- TRUE
        }
      }
      colnames(matched_genes.2) <- c("freq.2","aggregated_estimate.2")
      out <- cbind(matched_genes,matched_genes.2)
      
      ### fix aggregation mismatch ###
      for(i in 1:nrow(out)){
        if(sum(as.logical(out$aggregated_estimate[i]), out$aggregated_estimate.2[i]) == 1){
          # if first site has true --> check how many variants #
          if(!as.logical(out$aggregated_estimate[i])){
            temp <- site_spe_gt[site_spe_gt$Gene == out$gene[i],]
            temp <- temp[-which(temp$Variant %in% c("atleast_1new","each_unseen")),]
            
            agg_vals_temp <- sort(unique(site_spe_agg_gt$GT_agg[-which(site_spe_agg_gt$Variant %in% c("atleast_1new"))]))
            vals_temp <- sort(unique(site_spe_gt$GT[site_spe_gt$Variant != "atleast_1new" & site_spe_gt$Gene == out$gene[i]]))
            
            if(is.na(sort(unique(temp$GT))[3] < as.numeric(out$freq[i]))){
              rank_temp <- which.min(abs(vals_temp - as.numeric(out$freq[i])))
              out$freq[i] <- agg_vals_temp[rank_temp]
              out$aggregated_estimate[i] <- TRUE
            }
            
            if(!is.na(sort(unique(temp$GT))[3] <= as.numeric(out$freq[i]))){
              if(sort(unique(temp$GT))[3] < as.numeric(out$freq[i])){
                rank_temp <- which.min(abs(vals_temp - as.numeric(out$freq[i])))
                out$freq[i] <- agg_vals_temp[rank_temp]
                out$aggregated_estimate[i] <- TRUE
              }
            }
          }
          
          # if second site has true #
          if(!as.logical(out$aggregated_estimate.2[i])){
            temp <- site_spe_gt.2[site_spe_gt.2$Gene == out$gene[i],]
            temp <- temp[-which(temp$Variant %in% c("atleast_1new","each_unseen")),]
            
            agg_vals_temp <- sort(unique(site_spe_agg_gt.2$GT_agg[-which(site_spe_agg_gt.2$Variant %in% c("atleast_1new"))]))
            vals_temp <- sort(unique(site_spe_gt.2$GT[site_spe_gt.2$Variant != "atleast_1new" & site_spe_gt.2$Gene == out$gene[i]]))
            
            if(is.na(sort(unique(temp$GT))[3] < as.numeric(out$freq.2[i]))){
              rank_temp <- which.min(abs(vals_temp - as.numeric(out$freq.2[i])))
              out$freq.2[i] <- agg_vals_temp[rank_temp]
              out$aggregated_estimate.2[i] <- TRUE
            }
            
            if(!is.na(sort(unique(temp$GT))[3] <= as.numeric(out$freq.2[i]))){
              if(sort(unique(temp$GT))[3] > as.numeric(out$freq.2[i])){
                rank_temp <- which.min(abs(vals_temp - as.numeric(out$freq.2[i])))
                out$freq.2[i] <- agg_vals_temp[rank_temp]
                out$aggregated_estimate.2[i] <- TRUE
              }
            }
          }
        }
      }
    }
    
  }
  
  
  if(!is.null(tcga.cancer.type.2))
    return(out)
  else
    return(matched_genes)
}



