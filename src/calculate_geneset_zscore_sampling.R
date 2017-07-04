
calculate_geneset_zscore_sampling <- function(times=10000,cs,annotation_matrix, ncores=8, bienrich_pval=0.05,shrinking_geneset=FALSE){

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
	f_zg2cs <- function(ga,zcs2gs){
	  a <- cor.test(zcs2gs,ga) # logistic regression, anova, cor.tes, bayeslm, and t.test. anova/cor.test are identical and is the best.
	  p <- a$p.value
	  s <- ifelse(a$estimate>0,1,-1)
	  c(z=s*qnorm(p/2,lower.tail = F),p=p)
	}

	
	f_permute_gs <- function(permutn,geneset,genes,annotation_matrix,ncores=ncores,bienrich_pval=0.05){
	  print (permutn);
	  geneset <- sample(genes,length(geneset))
	  csv <- match(genes,geneset)
	  csv[!is.na(csv)]<- 1
	  csv[is.na(csv)] <- 0
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
	  ztop
	  
	  while(ztop > qnorm(bienrich_pval,lower.tail=F) ){
		z <- f_zcs2gs_geneset_mc(csv,annotation_matrix,geneset,annoted,ncores=ncores)
		topanno <- order(z[,1],decreasing = T)[1]
		print(z[topanno,])
		annoted <- unique(c(annoted,z[topanno,2])) 
		print(intersect(genes[which(annotation_matrix[,z[topanno,2]]>=1)],geneset))
		ztop <- z[topanno,1]
	  }
	  zcs2gsres <- f_zcs2gs_geneset_mc(csv,annotation_matrix,geneset,annoted,ncores=ncores,rm.annoted=F)
	  zcs2gsres
	}
	geneset = unique(cs[[1]])
	genes = rownames(annotation_matrix)
	permute_gs <- lapply(1:times,f_permute_gs,geneset,genes,annotation_matrix,ncores=ncores,bienrich_pval=bienrich_pval)
	permute_gs
}
