# This is a fork of Adatiss containing a fix of a arithmetic bug preventing the processing of GTEx data.
# For details and changes compared to the base version, see: https://github.com/mwgrassgreen/AdaTiSS/issues/2.

#' @title AdaTiSS
#' @description To obtain tissue specificity scores.
#' @author Meng Wang
#' \email{mengw1@stanford.edu}


#---------------------------------
# preprocessing: pre-filtering
# Input
# exp.mx: raw expression matrix with rows of genes and columns of samples
# dat.type: 'TPM or RPKM' - for RNA-seq data
#           'intensity' - for intensity data generated from mass spectrometry or microarray
# proc.zero: to adjust zero expression before taking the logarithm
#           'ceiled to 1' - max(TPM, 1) to take expression < 1 to be 1
#           'added 1' - (TPM + 1) to add 1 to the raw expression
#           'perturbed by a small value' - to add a random small amount to the low expression 
# filter.col.prp: to filter genes from non-expressed (NA) samples in a large proportion especially for intensity data
#                 default = 1 - to filter genes with NA proportion >= sample number * filter.col.prp
# exp.thres: the threshold for small expression especially for RNA-seq data
#            default = 1 - to filter genes with proportion of (TPM <= exp.thres) >= filter.col.prp
# Output
# exp.mx.f.log: filtered expression maxtrix in log scale

preproc.filter.fn = function (exp.mx, dat.type = "TPM or RPKM", proc.zero = 'ceiled to 1', filter.col.prp = 1, exp.thres=1) {
	if (dat.type == "intensity") {
		filter.row = rowSums(is.na(exp.mx)) >= (ncol(exp.mx)*filter.col.prp)
		exp.mx.f = exp.mx[!filter.row, ]
		exp.mx.f.log = log2(exp.mx.f)
	} else {
		if (dat.type == "TPM or RPKM") {
		   filter.row = rowSums(exp.mx <= exp.thres) >= (ncol(exp.mx)*filter.col.prp)
		   exp.mx.f = exp.mx[!filter.row, ]
		   if ( proc.zero == 'ceiled to 1') {
		   	   exp.mx.f[exp.mx.f < 1] = 1
		   }
		   if (proc.zero == 'added 1') {
		   	  exp.mx.f = exp.mx.f + 1
		   }
		   if (proc.zero == 'perturbed by a small value') {
		      exp.mx.f[exp.mx.f < 0.01] = runif(sum(exp.mx.f < 0.01), 0.001, 0.01)		   	
		   }
		   exp.mx.f.log = log2(exp.mx.f)		   
		} else {
			#exp.mx.f.log = NULL
			stop("to reset dat.type as 'intensity' or 'TPM or RPKM'")
		}
	}
	
	return(exp.mx.f.log)
}

#--------------------------------------------
# to summary sample expression in tissue level 
# by taking the median value of the sample expression from the same tissue
# Input
# X: preprocessed expression matrix in log scale
# p.dat: pheonotype info (eg. tissue type) for each sample
# Output
# tiss.abd: tissue level expression matrix with rows of genes and columns of tissue types
tiss.abd.fn = function(X, p.dat) {	
	tiss.nm.ls = sort(unique(p.dat[,2]))
	tiss.abd = matrix(NA, nrow(X), length(tiss.nm.ls))
	for (i in 1:length(tiss.nm.ls)) {
	  tiss.nm = tiss.nm.ls[i]
	  id.col = p.dat[p.dat[,2] == tiss.nm, 1]
	  X.sub = matrix(X[, id.col], nrow=nrow(X))
	  tiss.abd[,i] = apply(X.sub, 1, median, na.rm=TRUE)
	}
	rownames(tiss.abd) = rownames(X)
	colnames(tiss.abd) = tiss.nm.ls
	return(tiss.abd)
}


#----------------
# to obtain tissue specificity scores
# Input
# X: expression matrix in log scale
# tiss.abd: summarized tissue level expression (default: NULL)
#           if not providing tiss.ada (NULL), only output normalized expressio in sample level
# adjust: whether to adjust zero expression (default: FALSE)
#        to set TRUE if working on RNA-seq data with many zero expression
# adjust.opt: zero expression adjustment option (default: NULL)
#             adjust.opt = 0 - consider two cases to determine whether zeroes contribute to the population estimation
#             adjust.opt = other value - consider that all the zeroes contribute to the population estimation
# Output
# ada.s: score matrix in sample level
# ada.z: score matrix in tissue level
# pop.fit.mx: population fitting info including
#             n.observed - number of observed sample size
#             gam.sel - selected gamma parameter for fitting the population
#             mu0.hat - estimated population mean
#             sd0.hat - estimated population standard deviation
#             pi0.hat - estimated sample proportion in population
#             crt - criterion for selecting the gamma
#             pop.adjust - whether adjust population fitting at presence of zeroes
#             (note: to take another care on the genes with 'pi0.hat' <= 0.5)
AdaTiSS = function(X, tiss.abd=NULL, adjust=FALSE, adjust.opt=NULL) {

	pop.fit.mx = matrix(NA, nrow(X), 7)
	rownames(pop.fit.mx) = rownames(X)
	
	#------ to do: to remove 'crt' argument
	colnames(pop.fit.mx) = c("n.observed", "gam.sel", "mu0.hat", "sd0.hat" , "pi0.hat", "crt", "pop.adjust")
	id.ls.1 = rownames(X)[rowSums(!is.na(X)) >= 20]
	length(id.ls.1)
	for (i in 1:length(id.ls.1)) { 
	  if(i %% 500 == 0) print(i)
	  id = id.ls.1[i]
	  x.0 = X[id, ] 
	  
	   # to estimate population parameters
	  if (adjust == FALSE) {
	    gam.limit = ifelse(sum(!is.na(x.0)) <= 100, 1, 3)
	    result.x = adapt.gam.rob.fit.fn(x.0, gam.seq=seq(0,gam.limit,by=0.1), bin.num=round(length(x.0)/10))	  	
	  } else {
	  	# adjustment for zero expression especially for RNA-seq data
	  	result.x = adapt.gam.rob.fit.fn.adjust(x.0, adjust.opt = adjust.opt)
	  }
	  pop.fit.mx[id, names(result.x[["est.hat"]])] = result.x[["est.hat"]]
	}
	
	id.ls.2 = setdiff(rownames(X), id.ls.1)
	length(id.ls.2)
	pop.info.2 = apply(X[id.ls.2, ], 1, function(x) c(sum(!is.na(x)), median(x, na.rm=TRUE), mad(x, na.rm=TRUE), sum(abs(x-median(x, na.rm=TRUE)) <= 2*mad(x, na.rm=TRUE), na.rm=TRUE )/sum(!is.na(x)) ))
    pop.info.2 = t(pop.info.2)
    pop.fit.mx[id.ls.2, c("n.observed", "mu0.hat", "sd0.hat" , "pi0.hat")] = pop.info.2
    
    pop.fit.mx[, 'sd0.hat'] = pmax(pop.fit.mx[, 'sd0.hat'], 0.01)
    ada.s = (X - outer(pop.fit.mx[rownames(X), "mu0.hat"], rep(1, ncol(X))))/outer(pop.fit.mx[rownames(X), "sd0.hat"], rep(1, ncol(X)))
    if (!is.null(tiss.abd)) {
    	ada.z = (tiss.abd - outer(pop.fit.mx[rownames(tiss.abd), "mu0.hat"], rep(1, ncol(tiss.abd))))/outer(pop.fit.mx[rownames(tiss.abd), "sd0.hat"], rep(1, ncol(tiss.abd)))
    	ada.z[ada.z > 10 & !is.na(ada.z)] = 10
    	ada.z[ada.z < -10 & !is.na(ada.z)] = -10
    } else {
    	ada.z = NULL
    }


    return(list(ada.s = ada.s, ada.z=ada.z, pop.fit.mx=pop.fit.mx))
}



# ------------------------------------------------------------------------------
# data-adaptive selection procedure without population adjustment
adapt.gam.rob.fit.fn = function (x.00, gam.seq, step=50, mu.fix=NULL, var.fix=NULL, bin.num=NULL) {

	x.0 = x.00[!is.na(x.00)]	
	nm = c("mu0.hat", "sd0.hat", "pi0.hat", 'crt')
	par.hat = matrix(NA, length(gam.seq), length(nm))
	rownames(par.hat) = gam.seq
	colnames(par.hat) = nm
	
	for (i in 1:length(gam.seq)) {
		  	gam = gam.seq[i]
			x.mu = ifelse(is.null(mu.fix), mean(x.0), mu.fix)
		  	x.var = ifelse(is.null(var.fix), var(x.0), var.fix)
			result = est.fn(x.0, x.mu, x.var, gam, fix.mu=!is.null(mu.fix), fix.var=!is.null(var.fix), step=step)		
			mu.hat = result$mu.est
			var.hat = result$var.est
			if (!is.na(var.hat)) {
				est.result = efdr.0.fn(x.0, mu.hat, var.hat, gam, bin.num)
				par.hat[i, ] = est.result[colnames(par.hat)]
			}

	}
	crt.hat.0 = abs(pmin(par.hat[,'crt'], 10) - 1)
	ind.comp = !is.na(crt.hat.0)
	gam.comp = gam.seq[ind.comp]
	if (length(gam.comp)  == 0 ) {
		 est.hat=NA
		 x.n=NA
		 w=NA
		 gam.comp=NA
		 crt.hat.0=NA
		 para.hat.mx=NA
	} else {
		crt.hat.0 = crt.hat.0[ind.comp]
		par.hat = matrix(par.hat[ind.comp, ], ncol=length(nm))
		colnames(par.hat) = nm
		rownames(par.hat) = gam.comp
		
		crt.est.0 = pmin(par.hat[,'crt'], 10)
		gam.sel = gam.comp[which.min(crt.hat.0)]
		
		gam.sel.char = as.character(gam.sel)
		est.hat = c(length(x.0), gam.sel, par.hat[gam.sel.char, c("mu0.hat", "sd0.hat", "pi0.hat", "crt")])
		names(est.hat)[1:2] = c("n.observed", "gam.sel")
		
		w.nu = dnorm(x.00, est.hat["mu0.hat"], est.hat["sd0.hat"])^est.hat["gam.sel"]
		w.nu[is.na(x.00)] = NA
		w = w.nu/sum(w.nu, na.rm=TRUE)
	}
    est.hat['pi0.hat'] = min(1, est.hat['pi0.hat'])
    est.hat['crt'] = min(10, est.hat['crt'])                                                 
	return(list(est.hat=est.hat, x.w=w, gam.comp=gam.comp, crt.hat.0=crt.hat.0,  para.hat.mx=par.hat )) 
	
}


# ------------------------------------------------------------------------------
# data-adaptive selection procedure with population adjustment at presence of zero expression especially in RNA-seq data
adapt.gam.rob.fit.fn.adjust = function (x.00, gam.seq=NULL, adjust.opt = 0, step=50, mu.fix=NULL, var.fix=NULL, bin.num=NULL) {

	  y = x.00[!is.na(x.00)]
	  y.0 = y[y <= 0]
	  y.1 = y[y > 0]	
	  mu0 = 0
	  sd0 = ifelse(length(y.0) <= 1, 0, sd(y.0))
	  p0 = length(y.0)/length(y)
	  if (length(y.1) >= 20 && length(unique(y.1)) > 1) {
	  	gam.limit = ifelse(sum(!is.na(y.1)) <= 100, 1, 3)	  	
	  	result.x = adapt.gam.rob.fit.fn(y.1, gam.seq=seq(0,gam.limit,by=0.1), step=step, mu.fix=mu.fix, var.fix=var.fix, bin.num=round(length(y.1)/10))
	    x.est = result.x[["est.hat"]]
	    mu1 = x.est["mu0.hat"]
	    sd1 = x.est["sd0.hat"]
	    p1 = (1-p0)*min(1, x.est["pi0.hat"])
	    gam.sel = x.est['gam.sel']
	    ada.crt = x.est['crt']
	  } else {
	  	mu1 = median(y.1)
	  	sd1 = ifelse(length(y.1) == 1, 0, mad(y.1))
	  	p1 = 1-p0
	  	gam.sel = NA
	  	ada.crt = NA
	  }

      if (p0 > 0) {
      	pop.adjust = TRUE
		  if (adjust.opt == 0) {
		  	### consider two cases to determine whether zeroes contribute to the population estimation
			  if (p1 >= 0.7 & mu1 > 3*sd1) {
			  	result.x.adj = c(mu1, sd1, p1)
			  } else {
			  	thres = mu1 + 3*sd1
		        y.in = y[y < thres]
		        result.x.adj = c( mean(y.in), sd(y.in), length(y.in)/length(y))  
			  }	  	
		  } else {
		  	  ### or consider that all the zeroes contribute to the population estimation
			  pop.prp = p0 + (1-p0)*min(1, x.est["pi0.hat"])
			  pop.mean = p0*mu0 + (1-p0)*mu1
			  pop.sd = sqrt(p0*sd0^2 + (1-p0)*sd1^2 + (p0*mu0^2 + (1-p0)*mu1^2 - (p0*mu0 + (1-p0)*mu1)^2) ) 
		      result.x.adj = c(pop.mean, pop.sd, pop.prp)	  	
		  }      	
      } else {
      	pop.adjust = FALSE
      	result.x.adj = c(mu1, sd1, p1)
      }

      result = c(length(y), gam.sel, result.x.adj, ada.crt, pop.adjust)
      names(result) = c('n.observed', 'gam.sel', 'mu0.hat', 'sd0.hat', 'pi0.hat', 'crt', "pop.adjust")
      result['pi0.hat'] = min(1, result['pi0.hat'])  
      result['crt'] = min(10, result['crt'])                                
      return(list(est.hat=result))
  
}


# ------------------------------------------------------------------------------
# expected of fdr criterion
efdr.0.fn = function (x, mu.hat, var.hat, gam, bin.num=NULL) {
	    x = x[!is.na(x)]
	    den.fit = dnorm(x, mu.hat, sqrt(var.hat))
		frac.hat= mean(den.fit^gam)*sqrt(2*pi*var.hat)^gam * sqrt(1 + gam)
			
		my.hist = bk.cnt.fn(x, bin.num)
		bin.bk = my.hist[[1]]
		cnt.bk = my.hist[[2]]

        p0.hat.bin  = numeric(length(bin.bk)-1)
        p0.hat.bin[1] = pnorm(bin.bk[2], mu.hat, sqrt(var.hat))
        p0.hat.bin[length(bin.bk)-1] = 1 - pnorm(bin.bk[length(bin.bk)-1], mu.hat, sqrt(var.hat))
	    for ( j in 3:(length(bin.bk)-1)) {
	      	   p0.hat.bin[j-1] = pnorm(bin.bk[j], mu.hat, sqrt(var.hat)) - pnorm(bin.bk[j-1], mu.hat, sqrt(var.hat))
	     }
	    p.hat.bin =  cnt.bk/sum(cnt.bk)
	    null.efdr.hat = min(1,frac.hat) * sum(p0.hat.bin^2/p.hat.bin)	
	    est.sum = c(gam, mu.hat, sqrt(var.hat),  frac.hat, null.efdr.hat)
	    names(est.sum) = c("gamma", "mu0.hat", "sd0.hat", "pi0.hat",  'crt')
        return(est.sum)
}



# ------------------------------------------------------------------------------
# estimation under a fixed gamma
est.fn = function(x, mu.0, var.0, gam,  tol=10^(-4), step=step, fix.mu=FALSE, fix.var=FALSE) {
	          x = x[!is.na(x)]
	          n = length(x)
	   	      dum = dnorm(x, mu.0, sqrt(var.0))^gam
	          w.0 = dum/sum(dum)
	         
	       	  int = 1
	          flag = FALSE
	    	  diff.par.int = c()
	          mu.int = mu.0
	          var.int = var.0
	          while ( flag == FALSE) {
	          		     mu.1 = ifelse(fix.mu, mu.0, sum(w.0 * x))
	          	         var.1 =  ifelse(fix.var, var.0, (1+gam) * sum(w.0 * (x-mu.1)^2)) 
	          	         if (var.1 < 10^(-4)) {
							  mu.0 = NA
							  var.0 = NA
							  diff.par = NA
							  diff.par.int = c( diff.par.int, diff.par)
   							  break;
						 }
     	          	     	 
	          	         diff.par = abs(mu.1 - mu.0) + abs(sqrt(var.1) - sqrt(var.0))
	          	         diff.par.int = c( diff.par.int, diff.par)
	          	         if (diff.par < tol | int > step ){  
	          	         	flag = TRUE
	          	         	break;
	          	         } else {
	          	       	 	mu.0 = mu.1
	          	         	var.0 = var.1
	          	         	dum = dnorm(x, mu.0, sqrt(var.0))^gam
	                        w.0 = dum/sum(dum)
	                        mu.int = c(mu.int, mu.0)
	                        var.int = c(var.int, var.0)
	                        int = int + 1
	          	         }
	          }
	          return(list(mu.est=mu.0, var.est=var.0, w=w.0, diff.par.est=diff.par.int))
}


# ------------------------------------------------------------------------------
# to merge intervals s.t. each bin has positive number of data points
bk.cnt.fn = function (x, bin.num=NULL) {
	if (is.null(bin.num)) {
		if (length(x) > 1000) {
			bk.num = 20
		} 
		if (length(x) <= 1000 & length(x) > 500) {
			bk.num = 10
		} 
		if (length(x) <= 500){
			bk.num= 5
		}
	} else {
		bk.num = bin.num
	}
	h = hist(x, breaks=bk.num, plot=FALSE)
	# h$counts
	# h$breaks
	ind.zero = (1:length(h$counts))[h$counts == 0]
	if (length(ind.zero) != 0) {
	   bk.start = h$breaks[ind.zero]
		bk.end = h$breaks[ind.zero + 2]
		bk.update = c(bk.start[1])
		if (length(bk.start) > 1) {
				for (i in 1: (length(bk.start)-1) ) {
					  if (bk.end[i] < bk.start[i+1]) bk.update = c(bk.update, bk.end[i])					      
				}
	    	}
	   bk.update = c(bk.update,  bk.end[length(bk.end)] )
		bk.pts = sort(unique(c(h$breaks[-c(ind.zero, ind.zero+1)], bk.update)))
		cnt = unname(table(cut(x, breaks=bk.pts)))
	} else {
		bk.pts = h$breaks
		cnt = h$counts
	}
    return( list( bk.pts, cnt))
}

