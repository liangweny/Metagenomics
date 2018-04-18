##########################################################################################
## Reference: Feng Q, Liang S, Jia H, Stadlmayr A, Tang L, Lan Z, et al. Gut microbiome
##            development along the colorectal adenoma-carcinoma sequence. Nat Commun.
##            2015;6:6528
##########################################################################################

library(optparse)
option_list = list(
        make_option("--infile", type="character", metavar="character", default=NULL, help="input file, 2 columns, table header must be: KO\\tPadjust"),
        make_option("--outpre", type="character", metavar="character", default=NULL, help="output prefix. default: %default [required]"),
        make_option("--feature", type="character", metavar="character", default=NULL, help="kegg map feature file, 2 columns, table header must be: Path\\tKO"),
        make_option("--seed", type="integer", metavar="number", default=9999, help="random seed number")
)

opt = parse_args(OptionParser(usage="%prog[option]", option_list=option_list, add_help_option=TRUE, prog="KEGGreporter.r"), print_help_and_exit=TRUE, positional_arguments=FALSE)

file = opt$infile
file_f = opt$feature
outp = opt$outpre
seed = 9999
monteCarloTimes = 100000

if (is.null(file)) {cat("input file is not specified\n"); quit(save="no");}
if (is.null(outp)) {cat("output prefix is not specified\n"); quit(save="no");}
if (is.null(file_f)) {cat("feature is not specified\n"); quit(save="no");}
if (!is.null(opt$seed)) {seed <- opt$seed}


library('hash')
library('plyr')
set.seed(seed)


dat <- read.table(file, header=T)
rownames(dat) <- dat[,1]
rowName = rownames(dat)

dat0 <- dat
dat0$Zscore = qnorm( 1 - dat0$Padjust/2, lower.tail=TRUE)
numberOfRow = nrow(dat0)
dathash = hash()
for(i in 1:numberOfRow){
	dathash[dat0$KO[i]] = dat0$Zscore[i]
}

feature <- read.table(file_f, header=T)

##################################################################
## the following codes were originally created by songli
####################################################$$$$##########

featurePool = ddply(feature, .(Path), nrow)
names(featurePool) = c("Feature", "Count")
numberOfFeatureInPool = nrow( featurePool )

# sample K( number of genes in this pathway ) genes from whole gene set, using those Zscore to calculate Z_feature.
MCmethodToEstimateMeanAndSD <- function( geneZscoreInWholeSet, numberOfGenesInPathway, permutationTime ){
	monteCarloZscoreSet = c()
	for ( iter in 1:permutationTime ) {
		randomSubsetZscoreFromWholeSet = sample( geneZscoreInWholeSet, numberOfGenesInPathway, replace = FALSE )
		monteCarloZscoreSet[ iter ] = mean( randomSubsetZscoreFromWholeSet )
	}
	monteCarloMean = mean( monteCarloZscoreSet )
	monteCarloSD = sd( monteCarloZscoreSet )
	rm( monteCarloZscoreSet )
	return( c(monteCarloMean, monteCarloSD) )
}

# circle all features to find out significantly feature
outputDataframe = data.frame()

for ( featureID in 1:numberOfFeatureInPool ){
	featureName = featurePool$Feature[ featureID ]
	NumberOfGenesInFeature = featurePool$Count[ featureID ]
	cat( featureName, NumberOfGenesInFeature, "\n", sep = "\t" )
	genesListForCurrentPathway = feature$KO[ which(feature$Path == featureName) ]
	observedZfeatureForPathway = mean( values( dathash[ genesListForCurrentPathway ] ) )
	monteCarloEstimate = MCmethodToEstimateMeanAndSD( dat0$Zscore, NumberOfGenesInFeature, monteCarloTimes )
	estimateMean = monteCarloEstimate[ 1 ]
	estimateSD = monteCarloEstimate[ 2 ]
	adjustZfeature = (observedZfeatureForPathway - estimateMean) / estimateSD
	tmpDataframe = data.frame( featureName, NumberOfGenesInFeature, observedZfeatureForPathway, adjustZfeature, estimateMean, estimateSD)
        outputDataframe = rbind( outputDataframe, tmpDataframe )
}

names( outputDataframe ) = c( "Feature", "#KOs", "RawZfeature", "AdjustZfeature", "Mean", "SD")
write.table( file = paste(outp,".reporterZscore.txt",sep=""), outputDataframe, sep="\t", quote=F, row.names=F )
