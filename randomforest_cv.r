##########################################################################################
## the most of following codes were stolen shamefully from the ref's method.
## Reference: Feng Q, Liang S, Jia H, Stadlmayr A, Tang L, Lan Z, et al. 
##            Gut microbiome development along the colorectal adenoma-carcinoma sequence. 
##            Nat Commun. 2015;6:6528.
##########################################################################################

library(optparse)
option_list = list(
	make_option("--prof", type="character", metavar="character", default=NULL, help="profiling file, samples in column and features in row"),
	make_option("--outpre", type="character", metavar="character", default=NULL, help="output prefix. default: %default [required]"),
	make_option("--group", type="character", metavar="character", default=NULL, help="group file, 2 columns, table header must be: Sample<tab>Group"),
	make_option("--folder", type="integer", metavar="integer", default=10, help="cv folder. default: %default"),
	make_option("--trial", type="integer", metavar="integer", default=5, help="n trials. default: %default"),
	make_option("--step", type="double", metavar="double", default="0.9", help="the fraction of variables to remove at each step"),
	make_option("--width", type="double", metavar="double", default="8.0", help="plot width. default: %default"),
        make_option("--height", type="double", metavar="double", default="6.0", help="height width. default: %default")
)
opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="randomforest_cv.r"), print_help_and_exit=TRUE, positional_arguments=FALSE)
file1 = opt$prof
outp = opt$outpre
file2 = opt$group
folder = 10
trial = 5
step = 0.9
width = 8.0
height = 6.0

if (is.null(file1)) {cat("input prof is not specified\n"); quit(save="no");}
if (is.null(outp)) {cat("output prefix is not specified\n"); quit(save="no");}
if (is.null(file2)) {cat("group is not specified\n"); quit(save="no");}
if (!is.null(opt$folder)) {folder <- opt$folder}
if (!is.null(opt$trial)) {trial <- opt$trial}
if (!is.null(opt$step)) {step <- opt$step}
if (!is.null(opt$width)) {width <- opt$width}
if (!is.null(opt$height)) {height <- opt$height}


rfcv1 <- function (trainx, trainy, cv.fold = 5, scale = "log", step=0.5, mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE, ...){
	meancol <- nlevels(as.factor(trainy)) + 1
	classRF <- is.factor(trainy)
	n <- nrow(trainx)
	p <- ncol(trainx)
	if (scale == "log") {
		k <- floor(log(p, base = 1/step))
		n.var <- round(p * step^(0:(k - 1)))
		same <- diff(n.var) == 0
		if (any(same))
			n.var <- n.var[-which(same)]
		if (!1 %in% n.var)
			n.var <- c(n.var, 1)
	}
	else {
		n.var <- seq(from = p, to = 1, by = step)
	}
	k <- length(n.var)
	cv.pred <- vector(k, mode = "list")
	for (i in 1:k) cv.pred[[i]] <- trainy
	if (classRF) {
		f <- trainy
	}
	else {
		f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
	}
	nlvl <- table(f)
	idx <- numeric(n)
	for (i in 1:length(nlvl)) {
		idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, length = nlvl[i]))
	}
	res=list()
	for (i in 1:cv.fold) {
		all.rf <- randomForest(trainx[idx != i, , drop = FALSE], trainy[idx != i], 
					trainx[idx == i, , drop = FALSE], trainy[idx == i], 
					mtry = mtry(p), importance = TRUE, ...)
		cv.pred[[1]][idx == i] <- all.rf$test$predicted
		impvar <- (1:p)[order(all.rf$importance[,1], decreasing = TRUE)]
		res[[i]]=impvar
		for (j in 2:k) {
			imp.idx <- impvar[1:n.var[j]]
			sub.rf <- randomForest(trainx[idx != i, imp.idx, drop = FALSE], trainy[idx != i], 
						trainx[idx ==i, imp.idx, drop = FALSE], trainy[idx == i],
						mtry = mtry(n.var[j]), importance = recursive, ...)
			cv.pred[[j]][idx == i] <- sub.rf$test$predicted
			if (recursive) {
				impvar <- (1:length(imp.idx))[order(sub.rf$importance[,1], decreasing = TRUE)]
			}
			NULL
		}
		NULL
	}
	if (classRF) {
		error.cv <- sapply(cv.pred, function(x) mean(trainy != x))
	}
	else {
		error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
	}
	names(error.cv) <- names(cv.pred) <- n.var
	list(n.var = n.var, error.cv = error.cv, predicted = cv.pred,res=res)
}

library('randomForest')
library('plyr')
set.seed(999)

data <- read.table(file1)
data <- data[rowSums(data)!=0,]
colnames(data) <- gsub("^X", "", colnames(data))
colnames(data) <- gsub("\\.", "-", colnames(data))
colName <- colnames(data)
phe <- read.table(file2, header=T,sep="\t")
rownames(phe) <- phe[,1]
phe <- phe[colName,]
rowName <- rownames(phe)
gid <- intersect(colName ,rowName)
data <- data[,pmatch(gid,colName)]
phe <- phe[pmatch(gid,rowName),]
phe$Group=as.factor(as.character(phe$Group))
outcome=phe$Group
X <- as.data.frame(t(data))
nvar = ncol(X)
X$outcome = outcome


######5*10_crossvalidation####
result <- replicate(trial, rfcv1(X[,-ncol(X)], X$outcome, cv.fold=folder, step=step), simplify=FALSE)
error.cv <- sapply(result, "[[", "error.cv")
error.cv.cbm<-cbind(rowMeans(error.cv), error.cv)
write.table(error.cv.cbm,file=paste(outp,".cm_error.txt", sep=""), sep="\t",quote=F, col.names=F)

cutoff<-min (error.cv.cbm[,1])+sd(error.cv.cbm[,1])
error.cv.cut<-error.cv.cbm[error.cv.cbm[,1]<cutoff,]
cutn <- which.min(rownames(error.cv.cut))
cutv <- rownames(error.cv.cut)[cutn]

cat("cutoff", cutoff, "\n", sep=" ")
cat("cutn", cutn, "\n", sep=" ")
cat("cutv", cutv, "\n", sep=" ")

library('ggplot2')
library('reshape2')

cv_d <- cbind(result[[1]]$n.var, error.cv, rowMeans(error.cv))
cv_d <- as.data.frame(cv_d)
cv_meancol <- paste('V', ncol(cv_d), sep="")
plot_d <- melt(cv_d,id.var="V1")
colnames(plot_d) = c('Name','Group','Mean')
plot_s <- c()
plot_g <- c()
for(i in 1:nrow(plot_d)){
	if(plot_d$Group[i] == cv_meancol){
		plot_s[i] <- 1
		plot_g[i] <- c('Mean')
	}else{
		plot_s[i] <- 0.5
		plot_g[i] <- c('Ntrial')
	}
}
plot_d$size <- plot_s
plot_d$group <- plot_g
write.table(plot_d,file=paste(outp,".cv_error.matrix", sep=""), sep="\t",quote=F, row.names=F)


pdf(file=paste(outp,".cv.pdf",sep=""), width, height)
matplot(result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type="l",
        lwd=c(2, rep(1, ncol(error.cv))), col=1, lty=1, log="x",
        xlab="Number of variables", ylab="CV Error")
abline(v=cutv,col="pink",lwd=2)

dev.off()


#####pick N marker by corossvalidation#######
ntrial = trial * folder
k=1
b <- matrix(0,ncol=nvar,nrow=ntrial)
for(i in 1:trial){
	for(j in 1:folder){
		b[k,]<-result[[i]]$res[[j]]
		k=k+1
	}
}
Nmarker.list<-b[,1:cutv]
list<-c()
k=1
for(i in 1:cutv){
	for(j in 1:ntrial){
		list[k]<-Nmarker.list[j,i]
		k=k+1
	}
}
Nmarker.sort<-as.matrix(table(list))
Nmarker.sort<-Nmarker.sort[rev(order(Nmarker.sort[,1])),]
cutv <- as.numeric(cutv)
cat("cutv", cutv, "\n", sep=" ")
pick<- as.numeric(names(head(Nmarker.sort,cutv)))
tmp=X[,-ncol(X)]
Nmarker.pick<-colnames(tmp)[pick]

write.table(Nmarker.pick,file=paste(outp,".cross_validation_pick.txt", sep=""), sep="\t",quote=F, row.names=F, col.names=F)
