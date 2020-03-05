args <- commandArgs(TRUE)
B <- as.integer(args[1L])
L <- as.integer(args[2L])
U <- as.integer(args[3L])
X <- as.character(args[4L])

N <- 168

preprocess <- function(x){
	if(length(x)>B+1L){
		a <- do.call(rbind,lapply(x,function(u)c(substr(u,1L,1L),substr(u,nchar(u),nchar(u)))))
		p <- grepl('[^:0-9:]',a[-length(x),2L])
		q <- grepl('[^:0-9:]',a[-1L,1L])
		r <- p&q
		r <- c(TRUE,!r)
		out <- tapply(x,cumsum(r),paste0,collapse='\n',simplify=FALSE)
		return(out)
	} else{
		return(x)
	}
}
process <- function(i){
	x <- readLines(paste0('path',i,'.csv')) # change path
	x <- preprocess(x)
	#lengths(strsplit(x,','))
	y <- gregexpr('\"',x)
	z <- gregexpr(',',x)
	f <- function(z,y){
		if(length(y)){
			a <- matrix(y,ncol=2,byrow=TRUE)
			b <- lapply(seq_len(nrow(a)),function(i) a[i,1L]<=z&z<=a[i,2L])
			out <- z[rowSums(do.call(cbind,b))==0L]
			attr(out,'match.length') <- rep(1,length(out))
			return(out)
		} else{
			return(z)
		}
	}
	w <- mapply(f,z,y,SIMPLIFY=FALSE)
	g <- function(i){
		a <- cbind(c(1,w[[i]]+1L),c(w[[i]]-1L,nchar(x[[i]])))
		b <- a[,2L]-a[,1L]+1L
		m <- a[,1L]
		attr(m,'match.length') <- b
		out <- unlist(regmatches(x[[i]],list(m)))
		if(length(out)>N){
			w <- rev(which(out==""))[1L]
			out <- out[-w]
		}
		return(out)
	}
	y <- lapply(seq_len(length(x)),g)
	z <- as.data.frame(do.call(rbind,y[-1L]))
	colnames(z) <- y[[1L]]
	return(z)
}
fn <- function(i){
	out <- process(i)
	cat('\r',i); flush.console()
	return(out)
}
a <- dir('path/sim_out_SGE') # change path
b <- L:U
tmp <- lapply(b,fn)
results <- do.call(rbind,tmp)
save(results,file=paste0('path',X,'.RData')) # change path
