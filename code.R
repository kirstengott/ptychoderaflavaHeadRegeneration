## ----setup, include = FALSE----------------------------------------------

library(knitr)
opts_chunk$set(message=FALSE, results = 'hide', fig.keep="none", warning=FALSE, error=FALSE, fig.width = 15, fig.height = 15, fig.path="./analysis_figures/")


## ----libraries-----------------------------------------------------------
library(tidyverse)
library(pheatmap)
library(GO.db)
library(RColorBrewer)
library(kiRsten)

options(dplyr.width = Inf)



## ----functions-----------------------------------------------------------

simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep = "", collapse = " ")
}

get.anc.terms <- function(terms, ontol){
    ## Return the Ancestor terms for a list of terms in a particular ontology
    ## terms: a vector of GO terms specific to "ontol"
    ## ontol: one of either "CC", "BP", "MF", for which "terms" came from.
    require(GO.db)
    require(parallel)
    require(dplyr)
    ## uppercase the ontology
    ontol <- toupper(ontol)
    ## Convert the object to a list
    ont.name <- paste0("GO", ontol, "ANCESTOR")
    ont      <- base::get(ont.name)
    anc      <- as.list(ont)
    ## Remove GO IDs that do not have any ancestor
    anc       <- anc[!is.na(anc)]
    anc.terms <- anc[terms]
    ## find the GO term description for each parent term
    terms <- lapply(anc.terms, function(x){
        y  <- Term(GOTERM[x])
        dx <- data.frame(parent = names(y), parent_term = y)
        invisible(dx)
    })        
    ## attach the child term to the data frame of parent terms
    terms1 <- lapply(names(terms), function(x){
        df            <- terms[[x]]
        df$child      <- x
        df$child_term <- Term(GOTERM[x])
        invisible(df)
    })    
    ## bind all of the terms together
    terms.df <- dplyr::bind_rows(terms1)
    invisible(terms.df)
}

fix.ontol <- function(x){
    if (x =="MF"){
        new.x <- "Molecular Function"}
    if (x =="CC"){
        new.x <- "Cellular Component"}
    if (x =="BP"){
        new.x <- "Biological Process"}
    invisible(new.x)
}


calc.height <- function(go.df){
    if (nrow(go.df) < 5){
        height <- 170}
    else{
        height <- (nrow(go.df)/6)*100}
    invisible(height)
}


## ----code----------------------------------------------------------------
## The goal of this section is to find which differentially expressed genes are annotated to various different processes in development.



## First read in the enriched GO terms for upregulated and downregulated genes
up   <- read.table("allGO_upreg_GOterms.txt", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
down <- read.table("allGO_downreg_GOterms.txt", header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)


## select the unique GO ids from the enrichment analysis
up.terms.BP   <- filter(up, Ontol == "BP")$go_id
down.terms.BP <- filter(down, Ontol == "BP")$go_id


## map the ancestral terms
up.anc.BP   <- get.anc.terms(up.terms.BP, "BP")
down.anc.BP <- get.anc.terms(down.terms.BP, "BP")



## list the descriptors we are interested in looking at
tags <- c("signaling", 
          "proliferation", 
          "differentiation", 
          "dedifferentiation", 
          "organogenesis", 
          "morphogenesis", 
          "regeneration", 
          "ectoderm", 
          "mesoderm", 
          "endoderm", 
          "division")

## identify the GO term mappings to the 'tag' descriptors we would like to look at
terms <- lapply(tags, function(x){
    BP.parents <- unique(c(up.anc.BP$parent_term, down.anc.BP$parent_term))
    unique(grep(x, BP.parents, value=TRUE))
})

names(terms) <- tags


## pull out the descriptions that match our tags
desc <- c("cell-cell signaling",
          "BMP signaling pathway",
          "Notch signaling pathway",
          "Wnt signaling pathway",
          "fibroblast growth factor receptor signaling pathway",
          "cell proliferation",
          "stem cell proliferation",
          "cell differentiation",
          "endodermal cell differentiation",
          "mesodermal cell differentiation",
          "epidermal cell differentiation",
          "neuron differentiation", 
          "organ morphogenesis")



## munge the biological process dfs to get unique terms

bp <- rbind(up.anc.BP, down.anc.BP)

bp.sub <- filter(bp, parent_term %in% desc) %>% 
  distinct(parent) %>%
  .$parent



bp.sub.df           <- bp[which(bp$parent %in% bp.sub), ]
colnames(bp.sub.df) <- c("parent", "parent_term", "go_id", "Term")


bp.sub.up   <- bp.sub.df %>% filter(go_id %in% up$go_id)
bp.sub.down <- bp.sub.df %>% filter(go_id %in% down$go_id)

## join in the GO term enrichment
up.sub   <- left_join(up, bp.sub.up)
down.sub <- left_join(down, bp.sub.down)


## count number of genes in each parent
up.sub.1 <- up.sub %>% 
  dplyr::select(parent, parent_term, go_id, Term, Genes) %>% 
  arrange(parent) %>% 
  group_by(parent) %>% 
  mutate(parentNumGenes = length(unique(c(unlist(strsplit(Genes, split=", ")))))) %>%
  ungroup() %>%
  dplyr::select(parent, parent_term, parentNumGenes) %>%
  distinct()


## format the data frame for the picture
up.sub.1[is.na(up.sub.1)] <- "Other"
up.sub.1$parent_term      <- sapply(as.character(up.sub.1$parent_term), simpleCap)

num.genes        <- up.sub.1$parentNumGenes
names(num.genes) <- up.sub.1$parent_term

## relevel the factor for the picture
up.sub.1         <- within(up.sub.1, 
                           parent_term <- factor(parent_term, levels=names(sort(num.genes, decreasing=FALSE))))
png(filename="analysis_figures/UpDEAncestors.png", width=1000, height=1000)
ggplot(up.sub.1, aes(y=parentNumGenes, x=parent_term)) +
  geom_bar(stat="identity") +
  coord_flip() +
  theme_bw() +
  ylab("Number of Genes") + xlab("Gene Ontology Ancestral Term") +
  ggtitle("Ancestral GO terms of Up Regulated Genes") +
  theme(axis.title.y = element_text(size = 25, vjust = 1),
        axis.title.x = element_text(size = 25, vjust = 0.01),
        plot.title = element_text(size = 30),
        text=element_text(size = 15))
dev.off()

## Now do the same for the down regulated genes

down.sub.1 <- down.sub %>%
  dplyr::select(parent, parent_term, go_id, Term, Genes) %>%
  arrange(parent) %>%
  group_by(parent) %>%
  mutate(parentNumGenes = length(unique(c(unlist(strsplit(Genes, split=", ")))))) %>%
  ungroup() %>%
  dplyr::select(parent, parent_term, parentNumGenes) %>%
  distinct()


down.sub.1[is.na(down.sub.1)] <- "Other"
down.sub.1$parent_term        <- sapply(as.character(down.sub.1$parent_term), simpleCap)

num.genes        <- down.sub.1$parentNumGenes
names(num.genes) <- down.sub.1$parent_term

down.sub.1               <- within(down.sub.1, parent_term <- 
                                     factor(parent_term, levels=names(sort(num.genes, decreasing=FALSE))))

png(filename="analysis_figures/DownDEAncestors.png", width=1000, height=1000)
ggplot(down.sub.1, aes(y=parentNumGenes, x=parent_term)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() + 
  ylab("Number of Genes") + 
  xlab("Gene Ontology Ancestral Term") + 
  ggtitle("Ancestral GO terms of Down Regulated Genes") + 
  theme(axis.title.y = element_text(size = 25, vjust = 1),
        axis.title.x = element_text(size = 25, vjust = 0.01), 
        plot.title = element_text(size = 30), 
        text=element_text(size = 15))
dev.off()

## Now to the same except combining the enriched terms for the up and downregulated genes

up.sub.small   <- up.sub %>% dplyr::select(go_id, Term, Genes, parent, parent_term)
down.sub.small <- down.sub %>% dplyr::select(go_id, Term, Genes, parent, parent_term)

all.sub <- bind_rows(up.sub.small, down.sub.small)

all.sub.1 <- all.sub %>% 
  dplyr::select(parent, parent_term, go_id, Term, Genes) %>%
  arrange(parent) %>%
  group_by(parent) %>%
  mutate(parentNumGenes = length(unique(c(unlist(strsplit(Genes, split=", ")))))) %>%
  ungroup() %>%
  dplyr::select(parent, parent_term, parentNumGenes) %>%
  distinct()


sum.tot        <- sum(all.sub.1$parentNumGenes)
all.sub.1$perc <- (all.sub.1$parentNumGenes/sum.tot)*100
all.sub.1      <- na.omit(all.sub.1)
 
all.sub.1$parent_term       <- sapply(as.character(all.sub.1$parent_term), simpleCap)
all.sub.1[4, "parent_term"] <- "FGFR Signaling Pathway"

new.levels <- c("cell-cell signaling",
                "BMP signaling pathway",
                "Notch signaling pathway",
                "Wnt signaling pathway",
                "FGFR signaling pathway",
                "cell proliferation",
                "stem cell proliferation",
                "organ morphogenesis",
                "cell differentiation",
                "neuron differentiation",
                "endodermal cell differentiation",
                "mesodermal cell differentiation",
                "epidermal cell differentiation")

new.levels <- sapply(new.levels, simpleCap)
all.sub.1  <- within(all.sub.1, parent_term <- factor(parent_term, levels=new.levels))
                                  

num.genes        <- all.sub.1$parentNumGenes
names(num.genes) <- all.sub.1$parent_term


png(filename="analysis_figures/AllDEAncestors.png", width=1000, height=1000)
ggplot(all.sub.1, aes(y=perc, x=parent_term)) +
  geom_bar(stat="identity", fill="purple") +
  coord_flip() + 
  theme_bw() + 
  ylab("Percentage of Differentially Expressed Genes") + 
  xlab("Gene Ontology Term") + 
  ggtitle("Select GO terms of Differentially Expressed Genes") +
  theme(axis.title.y = element_text(size = 25, vjust = 1),
        axis.title.x = element_text(size = 25, vjust = 0.01), 
        plot.title = element_text(size = 30), 
        text=element_text(size = 20))
dev.off()


## Now we will make a more generalized graph for GO slim enrichment

all.up        <- read.table("goSlim_upreg_GOterms.txt", sep="\t", stringsAsFactors=FALSE, quote="",header=TRUE)
all.up$Parent <- sapply(all.up$Parent, function(x){is.na(x) <- "Offspring"})



all.down        <- read.table("goSlim_downreg_GOterms.txt", sep="\t", stringsAsFactors=FALSE, quote="",header=TRUE)
all.down$Parent <- sapply(all.down$Parent, function(x){is.na(x) <- "Offspring"})



## seperate out the time points to make the main plots


## 96 hour df up
go.df                     <- dplyr::filter(all.up, Num_Genes > 0, Parent == "Offspring", X0h.96h_BH.correction <= 0.01)
go.df.final               <- go.df[order(go.df$Num_Genes), ]
num.genes                 <- go.df.final$Num_Genes
names(num.genes)          <- go.df.final$Term
go.df.final               <- within(go.df.final, Term <- factor(Term, levels=names(sort(num.genes, decreasing=FALSE))))
go.df.final$BH.correction <- sapply(go.df.final$X0h.96h_BH.correction, function(x){paste0("   p \u2264 ", format(signif(x, digits = 1), scientific=TRUE))}) 
go.df.final$Ontology      <- sapply(go.df.final$Ontol, FUN=fix.ontol)
go.df.final$Time          <- "96h"
go.df.final$set           <- "Up"
go.df.final.96.up         <- go.df.final[, c("go_id", "Term", "Ontol", "Num_Genes", "Ontology", "BH.correction", "Time", "set")]


## 72 hour df up
go.df                     <- dplyr::filter(all.up, Num_Genes > 0, Parent == "Offspring", X0h.72h_BH.correction <= 0.01)
go.df.final               <- go.df[order(go.df$Num_Genes), ]
num.genes                 <- go.df.final$Num_Genes
names(num.genes)          <- go.df.final$Term
go.df.final               <- within(go.df.final, Term <- factor(Term, levels=names(sort(num.genes, decreasing=FALSE))))
go.df.final$BH.correction <- sapply(go.df.final$X0h.72h_BH.correction, function(x){paste0("   p \u2264 ", format(signif(x, digits = 1), scientific=TRUE))}) 
go.df.final$Ontology      <- sapply(go.df.final$Ontol, FUN=fix.ontol)
go.df.final$Time          <- "72h"
go.df.final$set           <- "Up"
go.df.final.72.up         <- go.df.final[, c("go_id", "Term", "Ontol", "Num_Genes", "Ontology", "BH.correction", "Time", "set")]




## 96 hour df down
go.df                     <- dplyr::filter(all.down, Num_Genes > 0, Parent == "Offspring", X0h.96h_BH.correction <= 0.01)
go.df.final               <- go.df[order(go.df$Num_Genes), ]
num.genes                 <- go.df.final$Num_Genes
names(num.genes)          <- go.df.final$Term
go.df.final               <- within(go.df.final, Term <- factor(Term, levels=names(sort(num.genes, decreasing=FALSE))))
go.df.final$BH.correction <- sapply(go.df.final$X0h.96h_BH.correction, function(x){paste0("   p \u2264 ", format(signif(x, digits = 1), scientific=TRUE))}) 
go.df.final$Ontology      <- sapply(go.df.final$Ontol, FUN=fix.ontol)
go.df.final$Time          <- "96h"
go.df.final$set           <- "Down"
go.df.final.96.down       <- go.df.final[, c("go_id", "Term", "Ontol", "Num_Genes", "Ontology", "BH.correction", "Time", "set")]


## 72 hour df down


go.df                     <- dplyr::filter(all.down, Num_Genes > 0, Parent == "Offspring", X0h.72h_BH.correction <= 0.01)
go.df.final               <- go.df[order(go.df$Num_Genes), ]
num.genes                 <- go.df.final$Num_Genes
names(num.genes)          <- go.df.final$Term
go.df.final               <- within(go.df.final, Term <- factor(Term, levels=names(sort(num.genes, decreasing=FALSE))))
go.df.final$BH.correction <- sapply(go.df.final$X0h.72h_BH.correction, function(x){paste0("   p \u2264 ", format(signif(x, digits = 1), scientific=TRUE))}) 
go.df.final$Ontology      <- sapply(go.df.final$Ontol, FUN=fix.ontol)
go.df.final$Time          <- "72h"
go.df.final$set           <- "Down"
go.df.final.72.down       <- go.df.final[, c("go_id", "Term", "Ontol", "Num_Genes", "Ontology", "BH.correction", "Time", "set")]





go.df      <- rbind(go.df.final.72.up, go.df.final.96.up, go.df.final.72.down, go.df.final.96.down)
go.df$set  <- factor(go.df$set, levels=c("Up", "Down"))
go.df$Term <- sapply(as.character(go.df$Term), simpleCap)
y.max      <- ceiling(max(go.df$Num_Genes)) + (ceiling(max(go.df$Num_Genes))/3)   



png(filename="analysis_figures/GOslimGraph.png", width=1000, height=1000)
ggplot(go.df, aes(y=Num_Genes, x=Term, fill=Ontology, label = BH.correction)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  theme_bw() +
  geom_text(hjust=0, size=4.3) + 
  ylab("Number of Genes") + xlab("Gene Ontology Term") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,y.max)) +
  scale_fill_manual(values=c("magenta","black", "aquamarine3")) + 
  facet_grid(set ~ Time) + 
  theme(axis.title.y = element_text(size = 25, vjust = 1.75), 
        axis.title.x = element_text(size = 25, vjust = 0.01), 
        legend.title=element_text(size=15), 
        text=element_text(size = 20),
        legend.justification=c(1,0),
        legend.position=c(1,0),
        strip.text.x = element_text(face="bold", size=15), 
        strip.text.y = element_text(face="bold", size=15))
dev.off()




## Now we will make a heatmap for all of the genes annotated from transcription factor.org

sig <- read.table("./significant.genes.txt", header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
tf  <- scan("/home/ejr/projects/Pflava/5_transcription_factors/tfs.out", what="character") ## read in the tf annotation
tf  <- sig[which(sig$Genes %in% tf), ]


## get the log2foldchanges of the significant transcription factors
tf.lfc           <- tf[ ,grep("log2FoldChange", colnames(tf))]
init             <- tf$"sp_desc"
colnames(tf.lfc) <- unlist(lapply(colnames(tf.lfc), function(x){
  strsplit(x, split="log2FoldChange.")[[1]][2]}))



## set the rownames of the lfc data frame.
rownames(tf.lfc) <- init

tf.lfc <- data.matrix(tf.lfc)

hr <- hclust(as.dist(1-cor(t(tf.lfc), method ="pearson")), method = "ward") # Clusters rows by Pearson correlation as distance method.
## cut the dendrogram into ten groups and assign each group a color
mycl          <- cutree(hr, k = 6)
clusters      <- c("Cluster A", "Cluster B", "Cluster C", "Cluster D", "Cluster E", "Cluster F")
mycols        <- clusters[as.vector(mycl)]
names(mycols) <- names(mycl)



annotation           <- data.frame(Cluster=mycols)
rownames(annotation) <- rownames(tf.lfc)
annotation$Cluster   <- as.factor(annotation$Cluster)
ann.cols             <- list(Cluster = c("Cluster A" = "green", "Cluster B" = "red", "Cluster C" = "blue", "Cluster D" = "pink", "Cluster E" = "black", "Cluster F" = "yellow"))

## map the colors used in the heatmap the the range of value in the data.frame
col.breaks    <- seq(-3, 3, length.out=101) ## needs to be one element longer than color vector
prettycolors  <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

tf.lfc.t <- make_pheatmap_df(tf.lfc)

                                 

pheatmap(tf.lfc.t, color = prettycolors, breaks = col.breaks, 
         clustering_distance_rows = as.dist(1-cor(t(tf.lfc), method = "pearson")), 
         clustering_method = "ward", cluster_cols=FALSE, 
         main = paste0("Clustering of  ", nrow(tf.lfc)," Genes"), 
         cellwidth = 15, cellheight=10, border_color=NA, fontsize=8, 
         filename = 'analysis_figures/TF_heatmap.png')



## Now we will make heatmaps for up and down regulated genes, as well as heatmaps for each individual cluster.


## define how many clusters to use
num.up.clusters   <- 4
num.down.clusters <- 4


## read in the analysis data
load("full_analysis.RData")

if('dx' %in% ls()){
    contrasts <- dx$contrasts
}

## read in the table of significant genes
sig.table <- read.table("significant.genes.txt", quote="", sep="\t", header=TRUE, stringsAsFactors=FALSE)
sp        <- filter(sig.table, sp_id != ".")
orf       <- filter(sig.table, length > 1)


## find the subset of genes that have swissprot hits and ORFS
sp.orf           <- filter(orf, Genes %in% sp$Genes)
rownames(sp.orf) <- sp.orf$Genes


## create a dataframe to reference later for swissprot hits
namer           <- data.frame(cbind(sp.orf$Genes, sp.orf$sp_desc))
colnames(namer) <- c("Genes", "sp")


## Create log fold change tables for heatmap creation for all genes within the pvalue cutoff specified in the config.file.
lfc.table  <- sp.orf[ ,grep("log2FoldChange", colnames(sp.orf), value=TRUE)]
pval.table <- sp.orf[ ,grep("padj", colnames(sp.orf), value=TRUE)]



## find the significant genes
genes.below.cutoff <- na.omit(rownames(pval.table[rowSums(pval.table < pval) > 0,,drop=FALSE]))



########## MAKE HEATMAP OF ALL UPREGULATED GENES AND CLUSTER THEM ##########

## first find the subset of genes that are upregulated to a fold change > 1 in at least one time point.
genes.up.cutoff    <- na.omit(rownames(lfc.table[rowSums(lfc.table > 1) > 0,,drop=FALSE]))
up.sub             <- genes.up.cutoff[which(genes.up.cutoff %in% genes.below.cutoff)]


p                  <- lfc.table[up.sub,,drop=FALSE]
colnames(p)        <- contrasts
rownames(p)        <- make.unique(as.character(namer[which(namer$Genes %in% rownames(p)), "sp"]))
p                  <- data.matrix(p)


## Clusters rows by Pearson correlation as distance method.
hr <- hclust(as.dist(1-cor(t(p), method = "pearson")), method = "complete") 

## cut the dendrogram into ten groups and assign each group a color
mycl            <- cutree(hr, k = num.up.clusters)
cluster.letters <- LETTERS[seq( from = 1, to = num.up.clusters )]
clusters        <- paste0("Cluster ", cluster.letters)
mycols          <- clusters[as.vector(mycl)]
names(mycols)   <- names(mycl)


## write out a table of the gene-cluster information
up.table        <- sig.table[which(sig.table$Genes %in% up.sub), c("Genes", "sp_desc")]
up.df           <- data.frame(up.table, cluster=mycols)
colnames(up.df) <- c("gene","sp_desc", "cluster")
write.table(up.df, file="upCluster_Genes.tsv", row.names=FALSE, quote=FALSE, sep="\t")


annotation           <- data.frame(Cluster=mycols)
rownames(annotation) <- rownames(p)
annotation$Cluster   <- as.factor(annotation$Cluster)
up.colors            <- rainbow(num.up.clusters)
names(up.colors)     <- clusters
ann.cols             <- list(Cluster = up.colors)


## map the colors used in the heatmap the the range of value in the data.frame
col.breaks    <- seq(-3, 3, length.out=101) ## needs to be one element longer than color vector
prettycolors  <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)


## change the values in the data.frame so they fit on the color break scale.
p.t <- make_pheatmap_df(p)


if (ncol(p) <= 7){
    cell.width <- 18
} else {
    cell.width <- 15
}

pheatmap(p.t, color = prettycolors, breaks = col.breaks, 
         clustering_distance_rows = as.dist(1-cor(t(p), method = "pearson")), 
         clustering_method = "complete", cluster_cols=FALSE, 
         main = paste0("Clustering of ", nrow(p)," Genes"), 
         cellwidth = cell.width, cellheight=10, border_color=NA, 
         annotation_row=annotation, annotation_colors = ann.cols, 
         fontsize=8, filename="analysis_figures/upregulated_heatmap.png")


########## CREATE HEATMAPS FOR UPREGULATED SUB CLUSTERS ##########


## pull out heatmap information for the different clusters
clusters <- unique(mycols)

cluster.dfs <- lapply(clusters, function(x){
    sub     <- names(mycols)[mycols == x]
    dx      <- data.matrix(p[which(rownames(p)%in% sub), ])
})
names(cluster.dfs) <- clusters

cluster.dfs.t        <- lapply(cluster.dfs, FUN=make_pheatmap_df)
names(cluster.dfs.t) <- clusters


## make heatmaps for the different clusters
lapply(clusters, function(x){
    orig.df   <- cluster.dfs[[x]]
    df.edit   <- cluster.dfs.t[[x]]
    if(nrow(orig.df) < 10){
        new.height = 5
    } else {
        new.height = NA}
    file.name <- gsub(" ", "", x, fixed = TRUE)
    pheatmap(df.edit, color = prettycolors, breaks = col.breaks, 
             clustering_distance_rows = as.dist(1-cor(t(orig.df), method = "pearson")), 
             clustering_method = "complete", cluster_cols=FALSE, 
             main = paste0(x, ": ", nrow(df.edit)," Genes"), 
             cellwidth = cell.width,  cellheight=10, border_color=NA, 
             fontsize=8, heigth = new.height, 
             filename = paste0("analysis_figures/up_", file.name, "_heatmap.png"))
})


########## MAKE HEATMAP OF ALL DOWNREGULATED GENES AND CLUSTER THEM ##########


genes.down.cutoff    <- na.omit(rownames(lfc.table[rowSums(lfc.table < -1) > 0,,drop=FALSE]))
down.sub             <- genes.down.cutoff[which(genes.down.cutoff %in% genes.below.cutoff)]


p                  <- lfc.table[down.sub,,drop=FALSE]
colnames(p)        <- contrasts
rownames(p)        <- make.unique(as.character(namer[which(namer$Genes %in% rownames(p)), "sp"]))
p                  <- data.matrix(p)


## Clusters rows by Pearson correlation as distance method.
hr <- hclust(as.dist(1-cor(t(p), method = "pearson")), method = "complete") 

## cut the dendrogram into ten groups and assign each group a color
mycl            <- cutree(hr, k = num.down.clusters)
cluster.letters <- LETTERS[seq( from = 1, to = num.down.clusters )]
clusters        <- paste0("Cluster ", cluster.letters)
mycols          <- clusters[as.vector(mycl)]
names(mycols)   <- names(mycl)



## write out a table of the gene cluster information
down.table        <- sig.table[which(sig.table$Genes %in% down.sub), c("Genes", "sp_desc")]
down.df           <- data.frame(down.table, cluster=mycols)
colnames(down.df) <- c("gene","sp_desc", "cluster")
write.table(down.df, file="downCluster_Genes.tsv", row.names=FALSE, quote=FALSE, sep="\t")


annotation           <- data.frame(Cluster=mycols)
rownames(annotation) <- rownames(p)
annotation$Cluster   <- as.factor(annotation$Cluster)
down.colors          <- rainbow(num.down.clusters)
names(down.colors)   <- clusters
ann.cols             <- list(Cluster = down.colors)



p.t <- make_pheatmap_df(p)

pheatmap(p.t, color = prettycolors, breaks = col.breaks, 
         clustering_distance_rows = as.dist(1-cor(t(p), method = "pearson")), 
         clustering_method = "complete", cluster_cols=FALSE, 
         main = paste0("Clustering of ", nrow(p)," Genes"),
         cellwidth = cell.width, cellheight= 10, border_color=NA,
         annotation_row=annotation, annotation_colors = ann.cols, 
         fontsize=8, filename="analysis_figures/downregulated_heatmap.png")


########## CREATE HEATMAPS FOR DOWNREGULATED SUB CLUSTERS ##########


## pull out heatmap information for the different clusters
clusters <- unique(mycols)

cluster.dfs <- lapply(clusters, function(x){
    sub     <- names(mycols)[mycols == x]
    dx      <- data.matrix(p[which(rownames(p)%in% sub), ])
})
names(cluster.dfs) <- clusters

cluster.dfs.t        <- lapply(cluster.dfs, FUN=make_pheatmap_df)
names(cluster.dfs.t) <- clusters


## make heatmaps for the different clusters
lapply(clusters, function(x){
    orig.df   <- cluster.dfs[[x]]
    df.edit   <- cluster.dfs.t[[x]]
    if(nrow(orig.df) < 10){
        new.height = 5
    } else {
        new.height = NA}
    file.name <- gsub(" ", "", x, fixed = TRUE)
    pheatmap(df.edit, color = prettycolors, breaks = col.breaks, 
             clustering_distance_rows = as.dist(1-cor(t(orig.df), method = "pearson")), 
             clustering_method = "complete", cluster_cols=FALSE, 
             main = paste0(x, ": ", nrow(df.edit)," Genes"), 
             cellwidth = cell.width,  cellheight=10, border_color=NA,
             fontsize=8, height = new.height, filename = paste0("analysis_figures/down_", file.name, "_heatmap.png"))
})



########## MAKE TREND LINE PLOTS OF CLUSTERS ##########

sub.sig           <- sig.table[,grep("log2FoldChange", colnames(sig.table), value=TRUE)]
colnames(sub.sig) <- lapply(colnames(sub.sig), function(x){strsplit(x, split="log2FoldChange.")[[1]][2]})
sub.sig$gene      <- sig.table$Gene


up   <- up.df[,c("gene", "cluster")]
down <- down.df[,c("gene", "cluster")]

up.df        <- inner_join(up, sub.sig, by="gene")
up.df$gene   <- NULL
down.df      <- inner_join(down, sub.sig, by="gene")
down.df$gene <- NULL

up.df.t   <- up.df %>% gather(Timepoint, Log2FoldChange, -cluster)
down.df.t <- down.df %>% gather(Timepoint, Log2FoldChange, -cluster)


up.means   <- arrange(up.df.t, cluster) %>% 
  group_by(cluster, Timepoint) %>% 
  summarise(aveFoldChange = mean(Log2FoldChange), error = (1.9645*se(Log2FoldChange)))  

down.means <- arrange(down.df.t, cluster) %>% 
  group_by(cluster, Timepoint) %>% 
  summarise(aveFoldChange = mean(Log2FoldChange), error = (1.9645*se(Log2FoldChange)))  

## In general, I recommend 95% CI interval for the error, which shows the range of your estimated mean.
### 95% CI is  mean +/- 1.9645*standard error.
up.means$Cluster   <- paste0("Up ", up.means$cluster)
down.means$Cluster <- paste0("Down ", down.means$cluster)


limits <- aes(ymax = aveFoldChange + error, ymin = aveFoldChange - error)


## grab the colors of the clusters from before
names(up.colors)   <- NULL
names(down.colors) <- NULL


png(filename = "analysis_figures/upCluster_line.png", width = 1000, height = 1000)
ggplot(up.means, aes(x = Timepoint, y = aveFoldChange, group = Cluster, colour = Cluster)) + 
  geom_line(size = 2.5) + 
  theme_bw() +
  ggtitle("Upregulated Gene Cluster Expression Patterns") + 
  ylab("Mean(Log2(FoldChange))") + xlab("Timepoint") +
  scale_colour_manual(values = up.colors) + 
  geom_errorbar(limits, width = 0.1, colour = "gray76") + 
  theme(axis.title.y = element_text(size = 25, vjust = 1.75), 
        axis.title.x = element_text(size = 25, vjust = 0.01), 
        legend.title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.justification = c(1,0), 
        legend.position = c(1,0), 
        title = element_text(size = 30, vjust = 0.99), 
        legend.key.size = unit(1.5, "cm"))
dev.off()



png(filename = "analysis_figures/downCluster_line.png", width = 1000, height = 1000)
ggplot(down.means, aes(x = Timepoint, y = aveFoldChange, group = Cluster, colour = Cluster)) + 
  geom_line(size = 2.5) +
  theme_bw() + 
  ggtitle("Downregulated Gene Cluster Expression Patterns") + 
  ylab("Mean(Log2(FoldChange))") + xlab("Timepoint") + 
  scale_colour_manual(values = down.colors) + 
  geom_errorbar(limits, width = 0.1, colour = "gray76") + 
  theme(axis.title.y = element_text(size = 25, vjust = 1.75), 
        axis.title.x = element_text(size = 25, vjust = 0.01), 
        legend.title = element_text(size = 20), 
        text = element_text(size = 20), 
        legend.justification = c(1,0), legend.position = c(1,0), 
        title = element_text(size = 30, vjust = 0.99), 
        legend.key.size = unit(1.5, "cm"))
dev.off()

