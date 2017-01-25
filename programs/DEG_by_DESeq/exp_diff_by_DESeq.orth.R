#!/bin/env Rscript

# This program is the same as exp_diff_by_DESeq.R, except that it
# compares gene expression among species.

##########################################################
################### Functions ############################
# now let us define some functions
 # read the htseq_count account
read_htseq_result<-function(f) # give the filename
{
	if(length(f) < 1)
	{
		stop("Nothing received in read_htseq_result")
	}
	if(! file.exists(f))
	{
		cat(f, " does not exist")
		return(NULL)
	}
	cnt<-read.delim(pipe(paste("grep -v '^_'",f,sep=' ')))
	return(cnt[c("gene_id","total")])
}

run_DESeq<-function(mat,sampleCondition=colnames(mat),cmpOrder=unique(sampleCondition),alpha=0.01,method="per-condition",sharingMode="maximum",...)
{
	if(isTRUE(as.logical(get_param('quiet'))))
	{
		suppressMessages(require(DESeq))
	}else
	{
		require(DESeq)
	}
	colnames(mat)<-NULL
	cds<-newCountDataSet(mat, sampleCondition)
	normStr<-'DESeq'
	if(!is.null(geneSetForNorm)) # use specified set to normalize counts
	{
		cat("[DESeq] Estimating size factors using autosomal genes only\n");
		autoMat<-mat[rownames(mat) %in% geneSetForNorm,]
		stopifnot(nrow(autoMat) > 100)
		autoCds<-newCountDataSet(autoMat,sampleCondition)
		autoCds<-estimateSizeFactors(autoCds)
		sizeFactors(cds)<-sizeFactors(autoCds)
		normStr<-paste(normStr,"autosomal genes only",sep=" + ")
	}else
	{
		cds<-estimateSizeFactors(cds)
	}
	# we need report the correlation between replicates. Some are so
	# bad for a reliable test
	report_cor(cds,sampleCondition,normMethod=normStr)
	cds<-estimateDispersions(cds, method=method, sharingMode=sharingMode,...)
	res<-nbinomTest(cds, cmpOrder[1],cmpOrder[2]) # ref level first
	res$change[with(res, baseMean > 0)]<-'no'
	res$change[with(res, baseMean > 0 & padj < alpha & log2FoldChange > 0)]<-'up'
	res$change[with(res, baseMean > 0 & padj < alpha & log2FoldChange < 0)]<-'down'
	return(res)
}

run_DESeq2<-function(mat,sampleCondition=colnames(mat),
	cmpOrder=unique(sampleCondition),alpha=0.01,testType='Wald',
	fitType="parametric",reducedModel=NULL,...)
{
	# require's 'quietly' option not working
	if(isTRUE(as.logical(get_param('quiet'))))
	{
		suppressMessages(require(DESeq2))
	}else
	{
		require(DESeq2)
	}
	dds<-DESeqDataSetFromMatrix(countData=mat, colData =
	data.frame(condition=sampleCondition), design = ~condition);
	# DESeq is a wrapper of 3 functions for size factors, dispersions,
	# and test
	#if(testType == "LRT")
	#{
	#	dds<-DESeq(dds, test=testType,fitType=fitType,reduced=reducedModel)
	#}else
	#{
	#	dds<-DESeq(dds, test=testType,fitType=fitType)
	#}
	
	# to limit normalization using a subset of genes, I will run each
	# step of DESeq separately
	normStr<-"DESeq2"
	if(!is.null(geneSetForNorm)) # use specified set to normalize counts
	{
		cat("[DESeq2] Estimating size factors using autosomal genes only\n");
		dds<-estimateSizeFactors(dds, controlGenes = rownames(dds)
		%in% geneSetForNorm);
		normStr<-paste(normStr,"autosomal genes only",sep=" + ")
	}else
	{
		dds<-estimateSizeFactors(dds)
	}
	report_cor(dds,sampleCondition,normMethod=normStr)
	dds<-estimateDispersions(dds, fitType=fitType)
	if(testType == "LRT")
	{
		dds<-nbinomLRT(dds, reduced=reducedModel, betaPrior=F)
	}else
	{
		dds<-nbinomWaldTest(dds, betaPrior=T)
	}

	tmp.contrast<-c("condition",cmpOrder[2],cmpOrder[1]) # ref come later
	res<-results(dds,contrast=tmp.contrast,alpha=alpha,...)
	res<-as.data.frame(res)
	res<-cbind(id=rownames(res),res); # consistent with DESeq1 result
	res$change[with(res, baseMean > 0)]<-'no'
	res$change[with(res, baseMean > 0 & padj < alpha & log2FoldChange > 0)]<-'up'
	res$change[with(res, baseMean > 0 & padj < alpha & log2FoldChange < 0)]<-'down'
	# add the baseMeanA and baseMeanB
	baseMeans<-sapply(cmpOrder, function(x)
	rowMeans(counts(dds,norm=T)[,colData(dds)$condition == x]))
	colnames(baseMeans)<-paste('baseMean',colnames(baseMeans),sep="")
	res<-cbind(res[1:2],baseMeans,res[-(1:2)])
	return(res)
}

# compare the results from two calculations
compare_DESeqs<-function(deseq,deseq2)
{
	stopifnot(is.data.frame(deseq) & is.data.frame(deseq2))
	cmp<-merge(deseq[c('id','change')],deseq2[c('id','change')],by=1,all=T)
	names(cmp)[-1]<-c("DESeq","DESeq2")
	return(with(cmp, addmargins(table(DESeq, DESeq2))))
}

report_cor<-function(cds,cond,meth='pearson',normMethod="DESeq2 standard")
{
	threshold<-0.95
	cnt<-counts(cds,norm=T)
	corr<-cor(cnt,method=meth)
	rownames(corr)<-cond
	colnames(corr)<-cond
	for(i in seq_len(nrow(corr)-1))
	{
		for(j in seq(i+1,nrow(corr)))
		{
			if(corr[j,i] < threshold)
			{
				corr[i,j]<--10
			}else
			{
				corr[i,j]<-10
			}
		}
	}
	cat("## Correlation among sample read counts by",normMethod,"\n",
	file=sumConn)
	capture.output(corr,file=sumConn) # write to the summary file
}

read_orthologs<-function(f,...)
{
	orthogroups<-read.delim(f,...)
	colnames(orthogroups)<-tolower(colnames(orthogroups))
	colnames(orthogroups)[1]<-'orthogroup'
	return(orthogroups)
}

read_gene_loc<-function(inFile) # only need gene location file
{
	if(!file.exists(inFile)) {stop("Can not find file ",inFile)}
	tmp1<-read.delim(inFile)[1:2]
	names(tmp1)<-c("gene_id","chr")
	tmp1$xlinked[tmp1$chr %in% c("X","XH","XHET","XL","XR")]<-1
	tmp1$xlinked[is.na(tmp1$xlinked)]<-0
	tmp1$xlinked[tmp1$chr %in% c("U","UH")]<-NA # put unknown chromosomes as NA
	return(tmp1)
}

# get only autosomal genes, in ortholog-mode, the genes being
# autosomal in all informative species
get_normalize_geneset<-function()
{
	geneLocRows<-grep('^gene-loc-file',names(env.config),value=T);
	stopifnot(length(geneLocRows) > 0);
	geneSetForNorm<-NULL
	if(orthoMode)
	{
		orthoLoc<-orthologs;
		sps<-c()
		for(row in geneLocRows)
		{
			m<-regexec("gene-loc-file:(\\w+)$",row);
			sp<-regmatches(row,m)[[1]][2];
			sps<-c(sps,sp)
			geneLoc<-read_gene_loc(get_param(row))
			colnames(geneLoc)<-c(sp,paste(c('chr','xlinked'),sp,sep='.'))
			orthoLoc<-merge(orthoLoc,geneLoc,by=sp,all.x=T);
		}
		# find the orthogroups being autosomal in all species
		selectedCols<-paste('xlinked',sps,sep='.')
		goodRows<-apply(orthoLoc[selectedCols],1, function(x)
		{ auto<-all(x==0); return(ifelse(is.na(auto),F,auto))} );
		geneSetForNorm<-orthoLoc[goodRows,'orthogroup']
	}else
	{
		if(length(geneLocRows) > 1)
		{
			stop("Only one gene-loc-file is expected in non-ortholog",
			" mode\n");
		}
		geneLoc<-read_gene_loc(get_param(geneLocRows))
		geneSetForNorm<-geneLoc[geneLoc$xlinked == 0,1]
	}
	return(geneSetForNorm);
}

# read gene lengths from files specified in env.config
read_gene_length<-function()
{
	geneLenRows<-grep('^gene-len-file',names(env.config),value=T);
	stopifnot(length(geneLenRows) > 0);
	result<-NULL;
	if(orthoMode)
	{
		result<-list();
		for(row in geneLenRows)
		{
			m<-regexec("gene-len-file:(\\w+)$",row);
			sp<-regmatches(row,m)[[1]][2];
			geneLen<-read.delim(get_param(row))[1:2]
			colnames(geneLen)<-c(sp,paste('len',sp,sep='.'))
			result[[sp]]<-geneLen # store result for each sp
		}
	}else
	{
		if(length(geneLenRows) > 1)
		{
			stop("Only one gene-len-file is expected in non-ortholog",
			" mode\n");
		}
		result<-read.delim(get_param(geneLenRows))[1:2]
	}
	return(result);
}

# for genes whose length shorter than lenThres, their read counts will
# be set to NA
filter_genes_by_length<-function(d,lenThres)
{
	stopifnot(is.numeric(lenThres));
	geneLen<-read_gene_length();
	if(orthoMode) # set NA for each species
	{
		for(sp in names(geneLen)) # geneLen is a list
		{
			# lenCol<-paste('len',sp,sep=".")
			spGeneLen<-geneLen[[sp]];
			# here genes without length information are kept
			badRows<-d[[sp]] %in% spGeneLen[spGeneLen[,2]<lenThres,1]
			pat<-paste(':',sp,'\\.\\d+$',sep='')
			selectedCols<-grep(pat,colnames(d),value=F)
			d[badRows,selectedCols]<-NA;
		}
	}else
	{
		# here geneLen is a data.frame
		badRows<-d[[1]] %in% geneLen[geneLen[,2] < lenThres,1];
		d[badRows,-1]<-NA;
	}
	return(d);
}

# read into configuration
read_config<-function(f)
{
	d<-read.delim(f,head=T,stringsAsFactors=F,row.names=1,comment.char='#') # format: parameter	value
	rownames(d)<-tolower(rownames(d)) # always using lower case
	# use an environment to store the configuration, which has no name
	# partial matching issue
	newEnv<-new.env(parent=emptyenv())
	tmp<-sapply(rownames(d),function(x) assign(x, d[x,1], envir=newEnv))
	return(newEnv)
}

# get configuration for a parameter
get_param<-function(param)
{
	if(!exists(param, envir=env.config, inherits=F)) return(NA);
	return(get(param,envir=env.config, inherits=F))
}

mywrite.table<-function(...)
{
	write.table(..., sep = "\t", row = F, quote = F)
}


################### End of Functions ########################
#############################################################

# read into arguments
args<-commandArgs(TRUE)

if(length(args) < 3)
{
	cat("
	Usage: exp_diff_by_DESeq.R <htseq-file-info> <cmp-samples> <config-file>

	<htseq-file-info> contains the file and sample names for read
	counts. It should also has a third column for species name when
	the comparison is between species. one sample per row.

	<cmp-samples> lists the comparisons of sample pairs. In each row,
	reference sample first, followed by alternative sample

	<config-file> contains other information such as ortholog-file,
	parameter settings for DESeq2, etc.
	")
	stop();
#	q("no", status = 1, runLast = FALSE’)
}

cat("Step 1: processing arguments\n")
# Get input files of htseq read counts, the sample name must match
# those in sample comparison file, it also has one column for species
# name
fileInfo<-read.delim(args[1],head=T,stringsAsFactors=F,comment.char='#') # format:
# filepath	sample_name	(species)

# the file containing sample pairs to compare
sampleCmp<-read.delim(args[2],head=T,stringsAsFactors=F,comment.char='#') # format: ref_group	target_group
# read into config information, including
# 1. gene location files
# 2. gene length files
# 3. orthologous relationship file
# 4. outfile basename
env.config<-read_config(args[3])

# save.image("tmp.RData");
# stop();
## Global variables #######
spSep<-':'; # the seprator between tissue and species names in sample label
fdr<-ifelse(is.na(get_param('fdr-threshold')), 0.01,
as.numeric(get_param('fdr-threshold')))
exe<-'exp_diff_by_DESeq.orth.R';
## End global variables ####
orthologs<-NA;
orthoMode<-FALSE;
if(is.na(get_param('ortholog-file')))
{
	cat("[INFO] No 'ortholog-file' found in config-file\n");
	cat("[INFO] so all comparisons are within species\n");
}else
{
	orthologs<-read_orthologs(get_param('ortholog-file'));
	orthoMode<-T;
}

# now let us read into the counts and store them separately for each
# species

cat("Step 2: Read read counts\n");
htseqCounts<-data.frame()
if(orthoMode) # if orthologs, start with it
{
	htseqCounts<-orthologs;
}

for(i in seq_len(nrow(fileInfo)))
{
	htseqFile<-fileInfo[i,1]
	sampleName<-fileInfo[i,2]
	counts<-read_htseq_result(htseqFile)
	stopifnot(nrow(counts) > 0)
	# modify the matrix's colname
	if(orthoMode) # with othologs
	{
		if(is.null(fileInfo[i,3]))
		{
			stop("[ERROR] A third column in ", args[1],
			" is expected for species name\n")
		}
		
		# add species information
		sampleName<-paste(sampleName,fileInfo[i,3],sep=spSep)
		mergeCol<-fileInfo[i,3]
		all.y<-F
	}else
	{
		mergeCol<-1
		all.y<-T
	}
	sampleName<-paste(sampleName,i,sep='.') # add number suffix
	colnames(counts)[2]<-sampleName

	if(nrow(htseqCounts) > 0)
	{
		htseqCounts<-merge(htseqCounts,counts,sort=F,by.x=mergeCol,by.y=1,
		all.x=T,all.y=all.y)
	}else # this occurs only for no-ortholog case
	{
		htseqCounts<-counts
	}
}

# remove shorter genes if requested
if(!is.na(get_param('len-threshold')))
{
	htseqCounts<-filter_genes_by_length(htseqCounts,as.numeric(get_param('len-threshold')))
}

# convert the dataframe into a matrix
if(orthoMode)
{
	numInfoCol<-ncol(orthologs)
	dataRowNames<-htseqCounts[['orthogroup']]
}else
{
	numInfoCol<-1
	dataRowNames<-htseqCounts[[1]]
}
tmp1<-as.matrix(htseqCounts[-seq_len(numInfoCol)])
rownames(tmp1)<-dataRowNames
# use sample conditions such as tissue types as matrix columns
condition<-sub('\\.\\d+$','',colnames(tmp1))
colnames(tmp1)<-condition
tmp1->htseqCounts

cat("Step 3: Preprocess data\n");
# normalize read counts (i.e., estimate factor sizes) using autosomal
# genes only, as X-linked genes often evolve faster and differ among
# tissues (e.g., in testes and ovary)
geneSetForNorm<-NULL;
if(isTRUE(as.logical(get_param('normalize-by-autosome'))))
{
	geneSetForNorm<-get_normalize_geneset();
	cat("Using",length(geneSetForNorm),"autosomal genes/orthogroups",
	"for normalizing counts\n");
}

# now compare two types of samples each time
cat("Step 4: Find DEGs with DESeq/DESeq2\n");
outbase<-get_param('outbase') # get the basename of the outputs
fitType<-ifelse(is.na(get_param('fittype')),"parametric",get_param("fittype"))
summaryFile<-paste(outbase,"summary.tsv",sep='.')
 # open connection to summary file firstly
sumConn<-file(summaryFile,'w')
cat("# Produced by ", exe, args, file=sumConn,sep=' ');

for(i in seq_len(nrow(sampleCmp)))
{
	groups<-unlist(sampleCmp[i,])
	grpName1<-gsub(',','_',groups[1])
	grpName2<-gsub(',','_',groups[2])
	outFile<-paste("DEG",grpName1,'vs',grpName2,sep="_")
	outFile<-gsub(spSep,'.',outFile)
	# the program allows one group may contain multiple types of
	# samples, separated by ','
	group1<-unlist(strsplit(groups[1],','))
	group2<-unlist(strsplit(groups[2],','))
	cat('Analysing',grpName1,'and',grpName2,"\n")
	# submat<-htseqCounts[,colnames(htseqCounts) %in% c(group1,group2)]
	submat<-htseqCounts[,colnames(htseqCounts) %in% c(group1,group2)]
	# also modify the colnames into standard one
	colnames(submat)[colnames(submat) %in% group1] <- 'A' # reference group
	colnames(submat)[colnames(submat) %in% group2] <- 'B' # target group
	# remove orthogroups with NA values
	submat<-submat[apply(submat,1,function(x) !any(is.na(x))),]

	# one can also remove genes with too few reads if requested

	# now run DEG analysis using DESeq or DESeq2
	if(isTRUE(as.logical(get_param('deseq'))))
	{
		cat(">>>\n# ",i,':',grpName1,'(A)', 'versus',grpName2,"(B)",
		"by DESeq\n\n", file=sumConn,sep=' ');
		res.deseq<-run_DESeq(submat,colnames(submat),c('A','B'),
		fitType=fitType, alpha=fdr);
		# write out the result
		out=paste(outFile,'DESeq','tsv',sep=".");
		cat("Writing out",out,"\n",sep=" ");
		mywrite.table(res.deseq,out)
		cat("<<<End\n",file=sumConn,sep='')
	}
	if(isTRUE(as.logical(get_param('deseq2'))))
	{
		cat(">>>\n# ",i,':',grpName1,'(A)', 'versus',grpName2,
		"(B) by DESeq2\n\n", file=sumConn,sep=' ');
		res.deseq2<-run_DESeq2(submat,colnames(submat),c('A','B'),
		fitType=fitType, alpha=fdr,
		altHypothesis="greaterAbs", independentFiltering = TRUE);
		out=paste(outFile,'DESeq2','tsv',sep=".");
		cat("Writing out",out,"\n",sep=" ");
		mywrite.table(res.deseq2,out)
		cat("<<<End\n",file=sumConn,sep='')
	}
	if(isTRUE(as.logical(get_param('cmp_deseqs'))))
	{
		res.cmp<-compare_DESeqs(res.deseq,res.deseq2);
		#out=paste(outFile,'DESeq_1_vs_2','tsv',sep=".");
		#mywrite.table(res.cmp,out)
		cat(">>>\n* ",i,':', "Comparison between DESeq and DESeq2\n",
		file=sumConn,sep=' ');
		capture.output(res.cmp,file=sumConn);
	}
}

close(sumConn)

cat("Whole work is done\n");

