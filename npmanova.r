########################################################################################
## Most of codes were originally created by qinyouwen
## Reference: Anderson, M.J. 2001. A new method for non-parametric multivariate analysis
##            of variance. _Austral Ecology_, *26*: 32-46. 
########################################################################################

library(optparse)
option_list = list(
	make_option("--prof", type="character", metavar="character", default=NULL, help="profiling file, samples in column and features in row"),
	make_option("--phe", type="character", metavar="character", default=NULL, help="phenotype file, samples in row and features in column"),
	make_option("--dist", type="character", metavar="character", default="bray", help="distance method used in 'vegdist'. default: %default"),
	make_option("--perm", type="integer", metavar="integer", default=10000, help="number of permutations. default: %default"),
	make_option("--padjust", type="character", metavar="character", default="BH", help="padjust methon. default: %default"),
        make_option("--outfile", type="character", metavar="character", default=NULL, help="output file. default: %default [required]")
)
opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="npmanova.r"), print_help_and_exit=TRUE, positional_arguments=FALSE)

file1 = opt$prof
outfile  = opt$outfile
file2 = opt$phe
perm_times  = 1000
dist_method  = "bray"
adjust_method = "BH"

if (is.null(file1)) {cat("input prof is not specified\n"); quit(save="no");}
if (is.null(outfile)) {cat("output prefix is not specified\n"); quit(save="no");}
if (is.null(file2)) {cat("phenotype file is not specified\n"); quit(save="no");}
if (!is.null(opt$perm)) {perm_times <- opt$perm}
if (!is.null(opt$dist)) {dist_method <- opt$dist}
if (!is.null(opt$padjust)) {adjust_method <- opt$padjust}

dat <- read.table(file1, header = T, sep = "\t", row.names = 1)
dat <- t(dat)
rownames(dat) <- gsub("X", "", rownames(dat))
phe <- read.table(file2, header = T, sep = "\t", row.names = 1)
phe <- as.data.frame(phe)
label <- colnames(phe)
dat <- dat[as.vector(rownames(phe)), ]

res_colnames <- c("phenotype", "Df", "SumsOfSqs", "MeanSqs", "F.Model", "R2", "Pvalue", "Padjust")
write.table(t(res_colnames), outfile, quote = F, sep = "\t", append = F, col.names = F, row.names = F)
outmatrix <- matrix(0, nrow = ncol(phe), ncol = length(res_colnames))

library("vegan")
dist_m <- vegdist(dat, method = dist_method) ## calculating distance
dist_m <- as.matrix(dist_m)

for(i in 1:ncol(phe)) {
	set.seed(0)
	subset = which(!is.na(phe[, i]))
	dist_tmp <- dist_m[subset, subset]
	dist_tmp <- as.dist(dist_tmp)
	if(length(levels(as.factor(phe[subset, i]) == 1))) { next }
	if(length(levels(as.factor(phe[subset, i]) <= 5))) { phe[, i] = as.factor(phe[, i]) }
	cat(i, "\t", label[i], "\tphenotype\n")

	run <- paste("phe_tmp <- data.frame(", label[i], " = phe[subset, i])")
	eval(parse(text = run))

	run <- paste("res <- adonis(dist_tmp~", label[i], ", data = phe_tmp, permutations = ", perm_times, ")", sep = "")
	eval(parse(text = run))
	tab <- res$aov.tab
	out <- t(c(label[i], tab$Df[1], tab$SumsOfSqs[1], tab$MeanSqs[1], tab$F.Model[1], tab$R2[1], tab$"Pr(>F)"[1]))
	outmatrix[i, 1:length(out)] <- out
}
outmatrix[, ncol(outmatrix)] <- p.adjust(outmatrix[, ncol(outmatrix)-1], method = adjust_method)
write.table(outmatrix, outfile, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
