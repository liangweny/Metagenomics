##########################################################################################
## the most of following codes were stolen shamefully from the ref's method.
## Reference: Feng Q, Liang S, Jia H, Stadlmayr A, Tang L, Lan Z, et al.
##            Gut microbiome development along the colorectal adenoma-carcinoma sequence.
##            Nat Commun. 2015;6:6528.
##########################################################################################

library(optparse)
option_list = list(
	make_option("--prof1", type="character", metavar="character", default=NULL, help="profiling file for training, samples in column and features in row"),
	make_option("--prof2", type="character", metavar="character", default=NULL, help="profiling file for predict, samples in column and features in row"),
	make_option("--outpre", type="character", metavar="character", default=NULL, help="output prefix. default: %default [required]"),
	make_option("--list", type="character", metavar="character", default=NULL, help="1 column, no header. filter features in prof by the list"),
	make_option("--group1", type="character", metavar="character", default=NULL, help="group file, 2 columns, table header must be: Sample<tab>Group"),
	make_option("--group2", type="character", metavar="character", default=NULL, help="group file for predict, 2 columns, table header must be: Sample<tab>Group"),
	make_option("--name", type="character", metavar="character", default=NULL, help="group name for prediction. default: %default [required]"),
	make_option("--roc", type="logical", action="store_true", help="plot roc curve"),
	make_option("--width", type="double", metavar="double", default="8.0", help="plot width. default: %default"),
        make_option("--height", type="double", metavar="double", default="8.0", help="height width. default: %default")
)
opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="randomforest_cv.r"), print_help_and_exit=TRUE, positional_arguments=FALSE)
file1 = opt$prof1
outp = opt$outpre
file2 = opt$group1
file3 = opt$list
file4 = opt$prof2
file5 = opt$group2
name = opt$name
width = 8.0
height = 8.0

if (is.null(file1)) {cat("input prof is not specified\n"); quit(save="no");}
if (is.null(outp)) {cat("output prefix is not specified\n"); quit(save="no");}
if (is.null(file2)) {cat("group is not specified\n"); quit(save="no");}
if (is.null(name)) {cat("group name is not specified\n"); quit(save="no");}
if (!is.null(opt$width)) {width <- opt$width}
if (!is.null(opt$height)) {height <- opt$height}



library('randomForest')
library('plyr')
set.seed(999)

data <- read.table(file1)
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
if(!is.null(file3)){
	marker_f <- read.table(file3, header=F)
	marker <- marker_f$V1
	marker <- as.vector(marker)
	all_marker <- colnames(X)
	X <- X[,pmatch(marker,all_marker)]
}
X$outcome = outcome
name_index <- match(name, levels(as.factor(outcome)))
cat(name_index)

train1 <-data.frame(X)
train1.rf <- randomForest(outcome ~ ., data =train1, importance = TRUE)
if(!is.null(file4)){
	dpre <- read.table(file4)
	colnames(dpre) <- gsub("^X", "", colnames(dpre))
	colnames(dpre) <- gsub("\\.", "-", colnames(dpre))
	colName2 = colnames(dpre)
	if(!is.null(file5)){
		phe2 <- read.table(file5, header=T,sep="\t")
		rownames(phe2) <- phe2[,1]
		phe2 <- phe2[colName2,]
		rowName2 <- rownames(phe2)
		gid2 <- intersect(colName2 ,rowName2)
		phe2 <- phe2[pmatch(gid2,rowName2),]
		phe2$Group=as.factor(as.character(phe2$Group))
		dpre <- dpre[,pmatch(gid2,colnames(dpre))]
	}
	test <- as.data.frame(t(dpre))
	if(!is.null(file3)){
		index_marker <- pmatch(marker,colnames(test),nomatch = 0)
		na_marker <- NULL
		flag = 1
		for(i in 1:length(index_marker)){
			if(index_marker[i] == 0){
				na_marker[i] <- marker[i]
				flag = 0
			}
		}
		if(flag == 0){
			na_marker <- na_marker[!is.na(na_marker)]
			rep_t <- nrow(test)*length(na_marker)
			na_c <- rep(0, times=rep_t)
			na_d <- matrix(na_c, nrow=nrow(test), ncol=length(na_marker))
			na_d <- as.data.frame(na_d)
			colnames(na_d) <- na_marker
			rownames(na_d) <- rownames(test)
			test <- cbind(test, na_d)
		}
		test <- test[,pmatch(marker,colnames(test))] 
	}
        colnames(test) <-  gsub("-", ".", colnames(test))
	train1.pre <- predict(train1.rf,test,type="prob")
	p.train<-train1.pre[,name_index]
	write.table(p.train,file=paste(outp,".cross_validation.predict.in.test.txt",sep=""), sep="\t",quote=F,col.names=F)
}else{
	train1.pre <- predict(train1.rf,type="prob")
	p.train<-train1.pre[,name_index]
	write.table(p.train,file=paste(outp,".cross_validation.predict.in.train.txt",sep=""), sep="\t",quote=F,col.names=F)
}


if(!is.null(opt$roc)){
	library(pROC)

	roc_flag = 0
	if(!is.null(file4) && !is.null(file5)){
		if(nlevels(phe2$Group) > 1){
		pdf(file=paste(outp,".ROC.in.test.pdf",sep=""), width, height)
		roc1 <- roc(phe2$Group,p.train,plot=T,percent=T,ci=T)
		roc_flag = 1
		}
	}else{	
		pdf(file=paste(outp,".ROC.in.train.pdf",sep=""), width, height)
		roc1 <- roc(outcome,p.train,plot=T,percent=T,ci=T)
		roc_flag = 1
	}
	if(roc_flag == 1){
		sens.ci <- ci.se(roc1, specificities=seq(0, 100, 5))
		write.table(sens.ci,file=paste(outp, ".sens.ci.txt",sep=""), sep="\t",quote=F,col.names=F)
		roc_df <- cbind(Specificity=roc1$specificities, Sensitivity=roc1$sensitivities)
		write.table(roc_df,file=paste(outp, ".ROC.txt",sep=""), sep="\t",quote=F,col.names=T, row.names=F)

		plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
		plot(sens.ci, type="bars")
		plot(roc1,col=2,add=T)
		legend("bottomright",c(paste("AUC=",round(roc1$ci[2],2),"%"),
		paste("95% CI:",round(roc1$ci[1],2),"%-",round(roc1$ci[3],2),"%")))

		dev.off()
	}
}
