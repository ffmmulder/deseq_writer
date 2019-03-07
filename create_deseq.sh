#input should be config file containing
#PROJECT_ID
#OUTPUT	(optional, working directory if not specified)
#REFERENCE	 (optional, human if not specified)
#REF_CONDITION (optional, first entry if not specified)
#REF_GROUP (optional, first entry if not specified)
#ANNOTATION_SRC (optional, BIOMART if not specified)
#check if settings config file was passed as (only parameter)

debug=1

if [ $# -ne 1 ]; then
    echo "Invalid parameter, usage: create_deseq <path/to/settings.cfg>"
    exit 1
else
	cfg=`realpath ${1}`
	#assume location of cfg file as working directory
	cd `dirname $cfg`

	debugLine $(pwd)
	ANNOTATION_SRC="BIOMART"
	OUTPUT_DIR=$(pwd)
	#load settings
	source $cfg	
fi

function debugLine {
	if [ $debug -eq 1 ] ; then
		echo $1
	fi
}

#set necessary files
META_DATA=${PROJECT_ID}_MetaData.txt
COUNTS=${PROJECT_ID}_readCounts_raw.txt

#debug code
debugLine ${META_DATA}
debugLine ${COUNTS}
debugLine $OUTPUT_DIR

#set annotation db to use, default to human if none specified
#for bioconductor db
function getBioconductorAnnotationDB {
	ref=$1
	#default to human if none specified, make uppercase for easy comparison
	if [ "$ref" == "" ] ; then
		echo org.Hs.eg.db
	else
		ref=$(echo ${ref} | awk '{print toupper($0)}')
	fi

	if [ "$ref" == "HUMAN" ] ; then
		echo org.Hs.eg.db
	elif [ "$ref" == "RAT" ] ; then
		echo org.Rn.eg.db
	elif [ "$ref" == "PIG" ] ; then
		echo org.Ss.eg.db
	elif [ "$ref" == "CHICKEN" ] ; then
		echo org.Gg.eg.db
	elif [ "$ref" == "ZEBRAFISH" ] ; then
		echo org.Dr.eg.db
	elif [ "$ref" == "DOG" ] || [ "$ref" == "CANINE" ]; then
		echo org.Cf.eg.db
	elif [ "$ref" == "FRUITFLY" ] || [ "$ref" == "FLY" ]; then
		echo org.Dm.eg.db
	elif [ "$ref" == "MOUSE" ] ; then
		echo org.Mm.eg.db
	fi	
}

#for biomart db
function getBiomartAnnotationDB {
	ref=$1
	#default to human if none specified, make uppercase for easy comparison
	if [ "$ref" == "" ] ; then
		echo hsapiens_gene_ensembl
	else
		ref=$(echo ${ref} | awk '{print toupper($0)}')
	fi
	
	if [ "$ref" == "HUMAN" ] ; then
		echo hsapiens_gene_ensembl
	elif [ "$ref" == "RAT" ] ; then
		echo rnorvegicus_gene_ensembl
	elif [ "$ref" == "PIG" ] ; then
		echo sscrofa_gene_ensembl
	elif [ "$ref" == "CHICKEN" ] ; then
		echo ggallus_gene_ensembl
	elif [ "$ref" == "ZEBRAFISH" ] ; then
		echo drerio_gene_ensembl
	elif [ "$ref" == "DOG" ] || [ "$ref" == "CANINE" ]; then
		echo cfamiliaris_gene_ensembl
	elif [ "$ref" == "FRUITFLY" ] || [ "$ref" == "FLY" ]; then
		echo dmelanogaster_gene_ensembl
	elif [ "$ref" == "MOUSE" ] ; then
		echo mmusculus_gene_ensembl
	fi	
}

#for biomart db
function getBiomartAnnotationSymbol {
	ref=$1
	#default to human if none specified, make uppercase for easy comparison
	if [ "$ref" == "" ] ; then
		echo hgnc_symbol
	else
		ref=$(echo ${ref} | awk '{print toupper($0)}')
	fi
	
	if [ "$ref" == "HUMAN" ] ; then
		echo hgnc_symbol
	elif [ "$ref" == "RAT" ] ; then
		echo rgd_symbol #mgi_symbol
	elif [ "$ref" == "PIG" ] ; then
		echo hgnc_symbol
	elif [ "$ref" == "CHICKEN" ] ; then
		echo hgnc_symbol
	elif [ "$ref" == "ZEBRAFISH" ] ; then
		echo zfin_id_symbol #hgnc_symbol
	elif [ "$ref" == "DOG" ] || [ "$ref" == "CANINE" ]; then
		echo hgnc_symbol
	elif [ "$ref" == "FRUITFLY" ] || [ "$ref" == "FLY" ]; then
		echo external_gene_name #?NO symbol version found?!
	elif [ "$ref" == "MOUSE" ] ; then
		echo mgi_symbol #hgnc_symbol
	fi	
}

#Generate comparison code to use in the R file
function constructComparisonList {
	#comparisons <- list(c("condition","P53KO","WT"),c("condition","P53KO.OncoRNF43","WT"), c("condition","P53KO.OncoRNF43","P53KO"), c("group","2","1"))

	#First for the conditions
	result="comparisons <- list("
	
	conditions=$1
	size=$(echo ${#conditions[@]})
	#echo $size
	for ((i=0; i<size; i++)); do
		#echo "Cond ${size}: ${conditions[${i}]}"
		for((j=i+1; j<=size; j++)); do			
			if [ "${conditions[${i}]}" != "" ] && [ "${conditions[${j}]}" != "" ]; then
				result="${result} c(\"condition\",\"${conditions[${i}]}\",\"${conditions[${j}]}\"),"
			fi
		done
		#echo "${conditions[${i}]}"
		
	done
	
	#Next for the groups
	groups=$2
	size=$(echo ${#groups[@]})
	#echo $size
	for ((i=0; i<size; i++)); do
		#echo "Grp ${size}: ${groups[${i}]}"
		for((j=i+1; j<=size; j++)); do
			if [ "${groups[${i}]}" != "" ] && [ "${groups[${j}]}" != "" ]; then
				result="${result} c(\"group\",\"${groups[${i}]}\",\"${groups[${j}]}\"),"
			fi
		done	
		#echo "${groups[${i}]}"		
	done	
	#remove last ,
	result=${result%?}
	result="${result} )"
	#result="${result} )"
	echo $result;	
}

#if Reference condition not set, use first value from metadata file?
if [ "$REF_CONDITION" == "" ]; then 
	REF_CONDITION=$(cut -f3 ${META_DATA} | head -2 | tail -1)
fi

#if Reference group not set, use first value from metadata file?
if [ "$REF_GROUP" == "" ]; then
	REF_GROUP=$(cut -f4 ${META_DATA} | head -2 | tail -1)	
fi

debugLine "ref_cond: ${REF_CONDITION}"
debugLine "ref_group: ${REF_GROUP}"

#create array of conditions and groups
conditions_1=($(tail -n+2 ${META_DATA} | cut -f3 | sort | uniq))
groups_1=($(tail -n+2 ${META_DATA} | cut -f4 | sort | uniq))

size=$(echo ${#conditions_1[@]})
for ((i=0; i<=size; i++)); do
	echo "c_1 ${groups_1[${i}]}"
	if [ "${conditions_1[${i}]}" != "" ]; then
		conditions+=(${conditions_1[${i}]})
		echo "c ${groups_1[${i}]}"
	fi
done

size=$(echo ${#groups_1[@]})
for ((i=0; i<=size; i++)); do
	echo "g_1 ${groups_1[${i}]}"
	if [ "${groups_1[${i}]}" != "" ]; then
		groups+=(${groups_1[${i}]})
		echo "g ${groups_1[${i}]}"
	fi
done
debugLine "conditions_1: ${conditions_1[*]}"
debugLine "groups_1: ${groups_1[*]}"


debugLine "conditions: ${conditions[*]}"
debugLine "groups: ${groups[*]}"

#generate comparisons code-line
comparisons=$(constructComparisonList $conditions $groups)
debugLine "comparison_line: $comparisons"

#create R file
touch RNASeqAnalysis_create.R

#echo $annotation_db

#First echo first part of R script to output file
#in case stuff needs to be installed
#if (!requireNamespace("BiocManager"))
#	install.packages("BiocManager")
#BiocManager::install()

#install.packages("DESeq2")
#install.packages("gplots")
#install.packages("heatmap3")
#BiocManager::install("DESeq2", version = "3.8")
echo 'library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("heatmap3")
library("BiocParallel")
library("biomaRt")
register(MulticoreParam(2))

# TODO CHANGE THIS (project)
project <- "'${PROJECT_ID}'"
setwd("'${OUTPUT_DIR}'")

#added reference
reference <- "'${REFERENCE}'"

countfile <- paste0(project,"_readCounts_raw.txt")
rawcounts <- data.frame(read.table(countfile, sep='\''\t'\'', header=T, row.names=1))
md <- data.frame(read.table(paste0(project,"_MetaData.txt"), sep='\''\t'\'', header=T, row.names=1, comment.char=""))
' > RNASeqAnalysis_create.R

if [ "$GENES_FILE" != "" ]; then
	echo 'genesofintrest <- "'${GENES_FILE}'"
	' >> RNASeqAnalysis_create.R
fi

debugLine "Source: ${ANNOTATION_SRC}"
#depending on annotation source used, echo next part
if [ "${ANNOTATION_SRC}" == "BIOCONDUCTOR" ]; then
	debugLine "bioconductor"
	annotation_db=$(getBioconductorAnnotationDB $REFERENCE)

	debugLine "bioconductor: '${annotation_db}''"

	echo '# --------------------------------------------------------------------------------
# ################################################################################
# --------------------------------------------------------------------------------
library("'${annotation_db}'")
genesymbols <- select('${annotation_db}', rownames(rawcounts), "SYMBOL", "ENSEMBL")
symbols <- genesymbols[!duplicated(genesymbols$ENSEMBL),]
' >> RNASeqAnalysis_create.R
elif [ "${ANNOTATION_SRC}" == "BIOMART" ]; then
	debugLine "biomart"
	dataset=$(getBiomartAnnotationDB ${REFERENCE})
	ens_symbol=$(getBiomartAnnotationSymbol ${REFERENCE})
	debugLine "biomart: ${dataset}"
	echo '# --------------------------------------------------------------------------------
# ################################################################################
# --------------------------------------------------------------------------------
library(biomaRt)
ensembl.symbol <- "'${ens_symbol}'"
ensembl <- useMart("ensembl", dataset="'${dataset}'")
annot<-getBM(c("ensembl_gene_id", ensembl.symbol), filters = "ensembl_gene_id", values = rownames(rawcounts),mart=ensembl)
symbols <- annot[!duplicated(annot$ensembl_gene_id),]
' >> RNASeqAnalysis_create.R
fi

debugLine "Past annotation part"
#and continue with the rest
echo '
# --------------------------------------------------------------------------------
# ################################################################################
# --------------------------------------------------------------------------------
# TODO CHANGE THIS (Group and ref label)
md$condition <- relevel(as.factor(md$Condition), ref="'${REF_CONDITION}'")
md$group <- relevel(as.factor(md$Group), ref="'${REF_GROUP}'")

dds <- DESeqDataSetFromMatrix(countData=rawcounts[,rownames(md)], colData=md, design= ~ group + condition)
dds <- DESeq(dds, parallel=F)
idx <- rowSums( counts(dds, normalized=T) >= 5 ) >= nrow(md)/2

dds <- DESeq(dds[idx,], parallel=F)
symbols <- symbols[idx,]
resultsNames(dds)

group.colors <- brewer.pal(length(unique(md$condition)), "Dark2")
names(group.colors) <- as.character(unique(md$condition))
# --------------------------------------------------------------------------------
# ################################################################################
# --------------------------------------------------------------------------------
# QA
rld <- rlog(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

#HEATMAPS
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=T)), decreasing=T)[1:1000]
hmcol <- colorRampPalette(c("darkgoldenrod1","white", "purple"))(n = 100)

# TODO CHANGE THIS (colors)
ColSideColors <- apply(md[,2:3], 2, as.character)
ColSideColors[,1] <- group.colors[ColSideColors[,1]]
ColSideColors[,2] <- group.colors[ColSideColors[,2]]
rownames(ColSideColors) <- rownames(md)


#SAMPLE DISTANCES
library("RColorBrewer")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rownames(sampleDistMatrix), sep=" | ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#PCA
pcadata <- plotPCA(rld, intgroup=colnames(md), returnData=T)

#QC
pdf(file=paste("QCplots_",project,".pdf",sep=""), width=12, height=10, pointsize=12)
  # PCA
  print(ggplot(pcadata, aes(PC1, PC2)) + geom_point(aes(color=condition, shape=Group), size=3, alpha=.8) + geom_text(aes(label=SampleName), vjust=2, size=2) + scale_color_manual(values=group.colors))
  # Sample distances
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
  # Heatmap of top expressed genes
  heatmap3(assay(vsd)[select,rownames(md)], col=hmcol, showRowDendro=F, scale="row", balanceColor=T, cexCol=1, cexRow=.01, ColSideColors=ColSideColors, labCol=paste(md$condition, rownames(md), sep=" | "), margins=c(10,3) )
  legend("topright", fill=group.colors, legend=names(group.colors))
dev.off()
' >> RNASeqAnalysis_create.R

if [ "$GENES_FILE" != "" ]; then
echo '
# --------------------------------------------------------------------------------
# ################################################################################
# --------------------------------------------------------------------------------
#
# A PRIORI GENES OF INTREST
gn <- data.frame(read.table(genesofintrest, sep='\''\t'\'', header=T))
pdf(file=paste("Selected_",project,".pdf",sep=""), width=10, height=10, pointsize=10)
for (i in c(1:nrow(gn))) {
  dat <- counts(dds,normalized=T)[match(gn$ENSEMBL[i], rownames(dds)),]
  toplot <- data.frame(count=as.numeric(dat))
  rownames(toplot) <- names(dat)
  total <- merge(toplot, md, by="row.names")
  colnames(total)[1] <- "Sample"
  genename <- gn$SYMBOL[i]
  if (is.na(genename) | genename=="NA") {
    genename <- gn$ENSEMBL[i]
  }
  plotty <- ggplot(total, aes(condition, count)) + geom_boxplot( aes(fill=condition), colour='\''black'\'', alpha=.8)  +  scale_fill_manual(values=group.colors)
  plotty <- plotty + geom_point(aes(fill=condition), shape=21, position ="dodge") 
  plotty <- plotty + ggtitle(genename) + geom_text(aes(label=SampleName), vjust=2, size=2) + expand_limits(y=0)
  print(plotty)
}
dev.off()
' >> RNASeqAnalysis_create.R
fi

echo '# --------------------------------------------------------------------------------
# ################################################################################
# --------------------------------------------------------------------------------

df <- data.frame(row.names=rownames(dds), MeanExpression=rowMeans(counts(dds,normalized=T)), SYMBOL=symbols$'${ens_symbol}')
# TODO CHANGE THIS (comparisons)
'${comparisons}'

pdf(file=paste("DE_",project,".pdf",sep=""), width=10, height=10, pointsize=10)
for (test in comparisons) {
  res <- results(dds, contrast=c(test[1], test[2], test[3]))
  label <- paste0(test[2]," vs. ",test[3]," in ", test[1])
  plotMA(res, main=label)

  colnames(res) <- paste(label, colnames(res), sep="_")
  df <- cbind(df, res[,c(2,5,6)])
}
dev.off()

write.csv(df, file=paste0("DEseq_",project,".csv" ))


# --------------------------------------------------------------------------------
# ################################################################################
# --------------------------------------------------------------------------------
# POSTERIOR GENES OF INTREST" 
df <- read.csv(file=paste0("DEseq_",project,".csv" ), header = T)
lab.row <- symbols[select,ensembl.symbol]
' >> RNASeqAnalysis_create.R

#loop through conditions and groups, generating code to generate and print corresponding plots
size=$(echo ${#conditions[@]})

for ((i=0; i<=size; i++)); do
	for((j=i+1; j<=size; j++)); do
		if [ "${conditions[${i}]}" != "" ] && [ "${conditions[${j}]}" != "" ]; then
			echo 'comparison <- paste0(project,"_'${conditions[${i}]}'_vs_'${conditions[${j}]}'_padj0-01")
# TODO CHANGE THIS (df$)
select <- which(df$'${conditions[${i}]}'.vs..'${conditions[${j}]}'.in.condition_padj <= 0.01 & abs(df$'${conditions[${i}]}'.vs..'${conditions[${j}]}'.in.condition_log2FoldChange)>=2)

write.csv( rownames(rld)[select], file=paste0(comparison,".csv"))
pdf(file=paste0(comparison,".pdf",sep=""), width=12, height=14, pointsize=12)
#heatmap3(assay(rld)[select,rownames(md)], col=hmcol, showRowDendro=F, scale="row", balanceColor=T, cexCol=1, cexRow=.4, ColSideColors=ColSideColors, labCol=paste(md$condition, rownames(md), sep=" | "), labRow=symbols[select,"SYMBOL"], main="'${conditions[${i}]}' vs. '${condition[${j}]}'", margins=c(10,7))
heatmap3(assay(rld)[select,rownames(md)], col=hmcol, showRowDendro=F, scale="row", balanceColor=T, cexCol=1, cexRow=.4, ColSideColors=ColSideColors, labCol=paste(md$condition, rownames(md), sep=" | "), labRow=lab.row, main=comparison, margins=c(10,7))
legend("topright", fill=group.colors, legend=names(group.colors))
dev.off()
' >>  RNASeqAnalysis_create.R
		fi
	done
done

size=$(echo ${#groups[@]})
for ((i=0; i<=size; i++)); do
	for((j=i+1; j<=size; j++)); do
		if [ "${groups[${i}]}" != "" ] && [ "${groups[${j}]}" != "" ]; then
			echo 'comparison <- paste0(project,"_'${groups[${i}]}'_vs_'${groups[${j}]}'_padj0-01")
# TODO CHANGE THIS (df$)
select <- which(df$'${groups[${i}]}'.vs..'${groups[${j}]}'.in.group_padj <= 0.01 & abs(df$'${groups[${i}]}'.vs..'${groups[${j}]}'.in.group_log2FoldChange)>=2)

write.csv( rownames(rld)[select], file=paste0(comparison,".csv"))
pdf(file=paste0(comparison,".pdf",sep=""), width=12, height=14, pointsize=12)
#heatmap3(assay(rld)[select,rownames(md)], col=hmcol, showRowDendro=F, scale="row", balanceColor=T, cexCol=1, cexRow=.4, ColSideColors=ColSideColors, labCol=paste(md$group, rownames(md), sep=" | "), labRow=symbols[select,"SYMBOL"], main="Group '${groups[${i}]}' vs. '${groups[${j}]}'", margins=c(10,7))
heatmap3(assay(rld)[select,rownames(md)], col=hmcol, showRowDendro=F, scale="row", balanceColor=T, cexCol=1, cexRow=.4, ColSideColors=ColSideColors, labCol=paste(md$condition, rownames(md), sep=" | "), labRow=lab.row, main=comparison, margins=c(10,7))
legend("topright", fill=group.colors, legend=names(group.colors))
dev.off()
' >>  RNASeqAnalysis_create.R
		fi
	done	
done

echo "
# --------------------------------------------------------------------------------
# ################################################################################
# --------------------------------------------------------------------------------
" >> RNASeqAnalysis_create.R
