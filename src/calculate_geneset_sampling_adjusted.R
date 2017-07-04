
calculate_geneset_sampling_adjusted <- function(zcs2gs_res,permute_gs){
	permute_gs_mx <- sapply(permute_gs,function(x){x[,1]})
	rownames(zcs2gs_res) <- zcs2gs_res[[2]]
	pvalue <- sapply(1:length(zcs2gs_res[[2]]),function(x){
	  z <- zcs2gs_res[x,1]
	  if(z==0){1}else
		{
		d <- permute_gs_mx[x,]
		d <- d[d>0]
		if(length(d)<100){
		  d <- as.numeric(permute_gs_mx)
		  l <- length(d)
		  d <- d[d>0]
		  mean <- mean(d)
		  var <- var(d)
		  shape = mean^2/var   ## 0.056
		  scale = var/mean
		  pgamma(z,shape,scale = scale,lower.tail = F)*length(d)/length(permute_gs_mx)
		}else{
			  mean <- mean(d)
			  var <- var(d)
			  shape = mean^2/var   ## 0.056
			  scale = var/mean
			  pgamma(z,shape,scale = scale,lower.tail = F)*length(d)/1000
		}
	  }
	})

	names(pvalue) <- zcs2gs_res[[2]]
	
	pvalueall <- sapply(1:length(zcs2gs_res[[2]]),function(x){
		d <- permute_gs_mx[x,]
		res <- d
		res[res<=0] <- -1
		d <- d[d>0]
		if(length(d)<100){
		  d <- as.numeric(permute_gs_mx)
		  l <- length(d)
		  d <- d[d>0]
		  mean <- mean(d)
		  var <- var(d)
		  shape = mean^2/var   ## 0.056
		  scale = var/mean
		  d <- permute_gs_mx[x,]
		  res[res>0] <- pgamma(d[d>0],shape,scale = scale,lower.tail = F)*length(d)/length(permute_gs_mx)
		}else{
			  mean <- mean(d)
			  var <- var(d)
			  shape = mean^2/var   ## 0.056
			  scale = var/mean
			  d <- permute_gs_mx[x,]
			  res[res>0] <-  pgamma(d[d>0],shape,scale = scale,lower.tail = F)*length(d)/1000
		}
		res[res<=0] <- 1
		res
	})
	colnames(pvalueall) <- zcs2gs_res[[2]]
	fdrall <- apply(pvalueall,1,function(x){
	  sapply(pvalue,function(y){
		sum(y>x)
	  })
	})

	fdr <- apply(fdrall,1,sum)/1000
	fdr[fdr>1] <- 1
	# build table

	# pvalue fdr geneset raw z-score totalgene include-gene genelist


	cnt<-t(sapply(names(pvalue),function(x){
	  c(length(rownames(annotation_matrix)[which(annotation_matrix[,x]==1)]),
		length(intersect(cs[[1]],rownames(annotation_matrix)[which(annotation_matrix[,x]==1)])),
		paste(intersect(cs[[1]],rownames(annotation_matrix)[which(annotation_matrix[,x]==1)]),collapse = " "))
	}))

	fun_enrich <- data.frame(Geneset=names(pvalue),pvalue=pvalue,fdr=fdr,zraw=zcs2gs_res[,1],cnt)
	fun_enrich <- fun_enrich[cnt[,2]!="0",]
	
	fun_enrich <- fun_enrich[,c(1,2,3,4,7)]
	colnames(fun_enrich)[5]<-"genes"
	fun_enrich
}

