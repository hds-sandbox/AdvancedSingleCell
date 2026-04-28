#a function to nicely plot correlations of cell features and principal components of a seurat object
plotCorrelations <- function(object, measures=c('nCount_RNA'), nfeatures=1000){
    options(repr.plot.width=14, repr.plot.height=5)
    obj <- FindVariableFeatures(object, nfeatures = nfeatures)
    var <- obj@assays$RNA@meta.features$vst.variable
    mat <- t( obj@assays$RNA@counts[var,] )
    pca <- prcomp(mat,  center = F, scale. = F)
    meta <- object@meta.data
    r2 <- matrix(0, nrow=length(measures), ncol=10)
    r2max <- c()
    for(ROW in 1:length(measures)){
        for(COL in 1:10){
            r2[ROW, COL] <- summary( lm( meta[, measures[ROW]] ~ pca$x[,COL]) )$r.squared
            }
        r2max <- c(r2max, which.max(r2[ROW,]))
        }
    
    plot_list = list()
    for(ROW in 1:length(measures)){
        p <- ggplot( meta, aes(x=pca$x[,r2max[ROW]], y=meta[,measures[ROW]])) + 
        geom_point(alpha=0.1, size=5)+
        geom_smooth(se=TRUE, method="lm") +
        ggtitle( paste("Regression of",measures[ROW],"~ PCA",r2max[ROW], "with maximum R^2",r2[ROW,r2max[ROW]]) ) + 
        xlab(paste("PCA",r2max[ROW])) + ylab(measures[ROW])

        print(p)
        }
    }

#a function to plot scores from a feature list on a UMAP, where the seurat object already contains marker scores for markers_list
plotScoresUMAP <- function(markers_list, 
                           seurat_data, 
                           repr.plot.width=20,
                           repr.plot.height=18){
    L <- length(names(seurat_data@meta.data))
    F <- length(names(markers_list))-1
    newnames <- names(seurat_data@meta.data)
    newnames[(L-F):L] <- names(markers_list)
    names(seurat_data@meta.data) <- newnames
    for(N in names(seurat_data@meta.data)[(L-F):L])
        seurat_data@meta.data[N] <- seurat_data@meta.data[N] / 
                                            max(seurat_data@meta.data[N])
    
    
    options(repr.plot.width=repr.plot.width, repr.plot.height=repr.plot.height)

    p <- FeaturePlot(seurat_data, 
            reduction = "umap", 
            features = names(markers_list), 
            order = TRUE,
            min.cutoff = 0,
            label = TRUE, 
            label.size = 7) + theme(legend.position = "right")
    
    return(p)
    }

#a function to correctly name the scores as in the function above, where the names are according to the ones in the markers list
renameScores <- function(markers_list, 
                           seurat_data, 
                           repr.plot.width=20,
                           repr.plot.height=18){
    L <- length(names(seurat_data@meta.data))
    F <- length(names(markers_list))-1
    newnames <- names(seurat_data@meta.data)
       message("Scores renamed FROM\n")
       message(cat(newnames[(L-F):L], sep="\n"))
       message("TO\n")
       message(cat(names(markers_list), sep="\n"))
    newnames[(L-F):L] <- names(markers_list)
    names(seurat_data@meta.data) <- newnames
    for(N in names(seurat_data@meta.data)[(L-F):L])
        seurat_data@meta.data[N] <- seurat_data@meta.data[N] / 
                                            max(seurat_data@meta.data[N])

    return(seurat_data)
    }



## Function that check for the presence of data and
## downloads it if necessary

downloadData <- function(){
# Check if the folder exists
if (!dir.exists("../Data")) {
  message("Data folder does not exists! Create one")
  dir.create("../Data")
  dir.create("../Data/faba1")
  dir.create("../Data/faba2")
  dir.create("../Data/faba3")
} else {
  message("Data folder exists. Check for files and eventually downloads them. Please wait.")
}

if (!file.exists("../Data/faba1.tar") ){
    message("Download ../Data/faba1.tar and unzip")
    system("wget https://github.com/hds-sandbox/AdvancedSingleCell/raw/refs/heads/main/Modules/Wells_data_analysis/Datasets/faba1.tar -O ../Data/faba1.tar")
    system("tar -xvf ../Data/faba1.tar -C ../Data/faba1 --strip-components=1")}

if (!file.exists("../Data/faba2.tar") ){
    message("Download ../Data/faba2.tar and unzip")
    system("wget https://github.com/hds-sandbox/AdvancedSingleCell/raw/refs/heads/main/Modules/Wells_data_analysis/Datasets/faba2.tar -O ../Data/faba2.tar")
    system("tar -xvf ../Data/faba2.tar -C ../Data/faba2 --strip-components=6")}

if (!file.exists("../Data/faba3.tar") ){
    message("Download ../Data/faba3.tar and unzip")
    system("wget https://github.com/hds-sandbox/AdvancedSingleCell/raw/refs/heads/main/Modules/Wells_data_analysis/Datasets/faba3.tar -O ../Data/faba3.tar")
    system("tar -xvf ../Data/faba3.tar -C ../Data/faba3 --strip-components=6")}

message("Done!")
}

## Converts the genes of the faba bean dataset to the corresponding orthologs
##
geneConvert <- function(seuratObject, genes, from, to){
    message("Converting gene names to orthologs. Please wait.")
    ortho <- readr::read_delim("https://github.com/hds-sandbox/AdvancedSingleCell/raw/refs/heads/main/Modules/Wells_data_analysis/Datasets/orthologTable.csv", ";")
    #error if from and to do not match one of the columns of the ortholog table
    if( !(from %in% colnames(ortho)) | !(to %in% colnames(ortho)) ){
        stop("Error: from and to must be one of the columns of the ortholog table. Names are: ", paste(colnames(ortho), collapse=", "))
    }   

    oldNames <- rownames(seuratObject)
    oldNames <- unlist( sapply( oldNames, function(x) gsub("GeneNAME", "", x) ) )
    newNames <- ortho %>%
            filter(str_detect(from, genes)) %>%
            pull(to)
    
    return(newNames)
}
    
## Function that assigns cluster names
## based on the highest score for the markers
## unnamed clusters must be in the Ident() of the object
## markersList is a list of markers where each element has name "cellType_scoring"
## (or any other name than "scoring", but it is important that the cellType comes
## first and that there is an underscore between the two words)
                                
clusterNames <- function(seuratObject, markersList){

    message("Cluster assignment started")
    
    clusterAssign <- as.vector(Idents(seuratObject))
    newClusterAssign <- clusterAssign
    for(CLST in unique(clusterAssign)){
        totalScore <- list()
        totalScoreNames <- as.vector( sapply(names(features_list), function(x) unlist(strsplit(x,"_"))[1]) )
        for(NAME in totalScoreNames){
            longNAME <- paste(NAME, "scoring", sep="_")
            totalScore[NAME] = sum(as.vector(seuratObject@meta.data[clusterAssign == CLST,longNAME]))
            }
        assignment <- totalScoreNames[ which.max(as.vector(totalScore)) ]
        newClusterAssign[ newClusterAssign==CLST ] <- assignment
        message(paste0("--- ", assignment, " assigned to ", CLST))
        }

    message("Cluster assignment finished")
    newClusterAssign
}