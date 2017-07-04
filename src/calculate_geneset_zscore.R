
calculate_geneset_zscore <- function(cs,annotation_matrix_gs, ncores=8, bienrich_pval=0.05,shrinking_geneset=FALSE){
	# cs: candidate genes and loci;
	# annotation_matrix_gs: annotation matrix;
	# ncores: number of cores for parallel computing
	# bienrich_pval: p-value cutoff to perform greedy search for the potention multiple mixed gene sets.
	# shrinking_geneset: shrinking the candidate gene set by full overlap, the larger one will be remove to accelerate computing 
	
	f_zcs2gs_geneset_mc <- function(csv,mtxg,geneset,annoted,ncores=8,rm.annoted=T){
	  annotedgene <- as.numeric(apply(data.frame(mtxg[,annoted]),1,max))
	  if(rm.annoted){
		mtxg <- mtxg[,-match(annoted,colnames(mtxg))]
	  }
	  require(parallel)
	  z <- mclapply(1:ncol(mtxg),function(i){
		x<-mtxg[,i]
		y<-csv
		y[which((y+annotedgene-x)==2)]<-0
		x=factor(x,levels = c(0,1));
		y=factor(y,levels = c(0,1));
		f <- fisher.test(table(y,x))
		p <- f$p.value
		s <- ifelse(f$estimate>1,1,-1)
		s*qnorm(p/2,lower.tail = F)
	  },mc.cores=ncores)
	  z <- data.frame(z=unlist(z),anno=colnames(mtxg),stringsAsFactors = F)
	  z
	}

	f_zcs2gs <- function(cs,mtxg){
	  apply(mtxg,2,function(x){
	  x=factor(x,levels = c(0,1));
	  y=factor(cs,levels = c(0,1));
	  f <- fisher.test(table(cs,x))
	  p <- f$p.value
	  s <- ifelse(f$estimate>1,1,-1)
	  s*qnorm(p/2,lower.tail = F)
	  })
	}
	f_zcs2gs_mc <- function(cs,mtxg,ncores=ncores){
	  require(parallel)
	  z <- mclapply(1:ncol(mtxg),function(x){
		x<-mtxg[,x]
		x=factor(x,levels = c(0,1));
		y=factor(cs,levels = c(0,1));
		f <- fisher.test(table(cs,x))
		p <- f$p.value
		s <- ifelse(f$estimate>1,1,-1)
		s*qnorm(p/2,lower.tail = F)
	  },mc.cores=ncores)
	  unlist(z)
	}
	f_zcs2gs_geneset_mc <- function(cs,mtxg,geneset,ncores=8){
	  require(parallel)
	  z <- mclapply(1:ncol(mtxg),function(x){
		x<-mtxg[,x]
		x=factor(x,levels = c(0,1));
		y=factor(cs,levels = c(0,1));
		f <- fisher.test(table(cs,x))
		p <- f$p.value
		s <- ifelse(f$estimate>1,1,-1)
		s*qnorm(p/2,lower.tail = F)
	  },mc.cores=ncores)
	  unlist(z)
	}
	f_zg2cs <- function(ga,zcs2gs){
	  a <- cor.test(zcs2gs,ga) # logistic regression, anova, cor.tes, bayeslm, and t.test. anova/cor.test are identical and is the best.
	  p <- a$p.value
	  s <- ifelse(a$estimate>0,1,-1)
	  c(z=s*qnorm(p/2,lower.tail = F),p=p)
	}

	
	annotation_matrix <- annotation_matrix_gs
	anno <- colnames(annotation_matrix)
	genes <- rownames(annotation_matrix)

	table(cs[[2]])
	
	
	geneset <- cs[[1]]
	geneset <- intersect(geneset,genes)
	cs <- cs[!is.na(match(cs[[1]],geneset)),]

	csv <- match(genes,geneset)
	csv[!is.na(csv)]<- 1
	csv[is.na(csv)] <- 0
	ncores=12
	zcs2gs <- f_zcs2gs_mc(csv, annotation_matrix,ncores ) 
	anno[which(zcs2gs>qnorm(1e-4,lower.tail=F))]
	anno[order(zcs2gs,decreasing = T)][1:10]

	x <-  cs[[1]][cs[[2]]=="l11"]

	#annotation_matrix_weigth <- 1/apply(annotation_matrix,2,sum)^2


	# restrict each loci two genes

	maxgene_perloci <- 2

	lz <- list()
	csn <- names(sort(table(cs[[2]]),decreasing=T))
	csn <- names(which(sort(table(cs[[2]]),decreasing=T)>maxgene_perloci))
	geneset_tmp <- geneset
	genes_tmp <- genes
	removed_tmp <- character()


	cat('preprocessing step to balance gene density;\n')
	while(max(table(cs[[2]]))>maxgene_perloci){
	  i <- names(sort(table(cs[[2]]),decreasing=T)[1])
	  print (i)
	  x <- cs[[1]][cs[[2]]==i]
	  geneset_tmp <- setdiff(intersect(geneset_tmp,genes_tmp),c(x,removed_tmp))
	  csv <- match(genes_tmp,geneset_tmp)
	  csv[!is.na(csv)]<- 1
	  csv[is.na(csv)] <- 0
	  ncores=12
	  zcs2gs <- f_zcs2gs_mc(csv, annotation_matrix,ncores ) 
	  if(length(x)==1){
		  zg <- annotation_matrix[x,]#*annotation_matrix_weigth
	  }else{
		  zg <- t(data.frame(annotation_matrix[x,]))#*annotation_matrix_weigth
	  }
	  lz[[i]] <- cor(zg,zcs2gs)
	  if(length(lz[[i]])>maxgene_perloci){
		removed_tmp <- c(removed_tmp,rownames(lz[[i]])[order((lz[[i]]))][1])
	  }
	  print(lz[[i]])
	  print(removed_tmp)
	  cs1 <- setdiff(cs[[1]],removed_tmp)
	  cs <- cs[!is.na(match(cs[[1]],cs1)),]
	}

	geneset <- setdiff(geneset,removed_tmp)

	## maximum enrichment with allowing to remove enriched gene, cutoff

	geneset <- cs[[1]]
	geneset <- intersect(geneset,genes)
	geneset <- setdiff(geneset,removed_tmp)
	csv <- match(genes,geneset)
	csv[!is.na(csv)]<- 1
	csv[is.na(csv)] <- 0
	ncores=ncores
    
	##### shrinking geneset
	if(shrinking_geneset){
		geneset_sum <- apply(annotation_matrix,2,function(x){
		  csv+x==2
		})
		geneset_size <- apply(annotation_matrix,2,sum)
		geneset_keep <- mclapply(1:ncol(geneset_sum),function(i){
		  x <- geneset_sum[,i]
		  apply(geneset_sum,2,function(y){
			sum(x!=y)==0
		  })
		},mc.cores=ncores)
		geneset_keep_mx <- sapply(geneset_keep,function(x){x})

		geneset_keep_id <- apply(geneset_keep_mx,1,function(x){
		  size <- as.numeric(x)*geneset_size
		  rownames(geneset_keep_mx)[which(size==min(size[size>0]))]
		})
		annotation_matrix <- annotation_matrix[,unique(unlist(geneset_keep_id))]
	}
	


	anno <- colnames(annotation_matrix)
	zcs2gs <- f_zcs2gs_mc(csv, annotation_matrix,ncores ) 
	topanno <- order(zcs2gs,decreasing = T)[1]
	anno[topanno]
	mappedgene <- intersect(genes[which(apply(data.frame(annotation_matrix[,anno[topanno]]),1,sum)>0)],geneset)
	mappedgenecsv <- match(genes,mappedgene)
	mappedgenecsv[!is.na(mappedgenecsv)]<- 1
	mappedgenecsv[is.na(mappedgenecsv)] <- 0

	annoted <- character()
	annoted <- unique(c(annoted,anno[topanno]))

	ztop <- zcs2gs[topanno]

	while(ztop > qnorm(bienrich_pval,lower.tail=F) ){
	  z <- f_zcs2gs_geneset_mc(csv,annotation_matrix,geneset,annoted,ncores=ncores)
	  topanno <- order(z[,1],decreasing = T)[1]
	  print(z[topanno,])
	  annoted <- unique(c(annoted,z[topanno,2])) 
	  print(intersect(genes[which(annotation_matrix[,z[topanno,2]]>=1)],geneset))
	  ztop <- z[topanno,1]
	}

	annoted_res <- annoted

	zcs2gs_res <- f_zcs2gs_geneset_mc(csv,annotation_matrix,geneset,annoted_res,ncores=ncores,rm.annoted=F)

	
	if(shrinking_geneset){
		list(zcs2gs_res,annotation_matrix)
	}else{
		zcs2gs_res
	}
}