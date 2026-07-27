


# ImportH5ad
# import data to Seurat format
# Input: folder. Folder containing the processed files from h5ad. Should include a .mtx sparce matrix count table, cellID, geneID, metadata. 
# org: organism code.
# Output: a seurat object.
# note: in the output file, the .mtx is loaded to RNA assay, count slot. Optimally it is umi count or raw counts.

ImportH5ad <- function (folder, org) {
  # load cell barcode and gene ID
  cellID <- read.csv(file = paste0(folder, "cellID_", org, ".csv"), header=F)
  geneID <- read.csv(file = paste0(folder, "geneID_", org, ".csv"), header=F)
  #load metadata
  metadata <- read.csv(file = paste0(folder, "metadata_",org, ".csv"), header=T, row.names = 1)
  #load count table
  counts <- readMM(file = paste0(folder, "scdata_", org, ".mtx"))
  rownames(counts) <- cellID[,1]
  colnames(counts) <- geneID[,1]
  counts<-t(counts)
  scdata <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = metadata, project = paste(org))
  scdata <- DietSeurat(scdata)
  return(scdata)
}

# ReadSTAR
# import gene.out.tab from a batch run of star to a count matrix.
# input: folder. the folder containing a batch run result of star, e.g., a foler with many subfolders of well names , "A1", "A2", ..., "P24".
# input: plate.name. a string of plate name. E.g.: "PN_plate7".
# output: data.raw. a count matrix. row is gene count, column is cell/well names.
# example output
#                     plate7_P5 |  plate7_P6 | plate7_P7 | plate7_P8 | plate7_P9
# Xkr4                    0            0           0           0          0
# Rp1                     0            0           0           0          0
# Sox17                   0            0           0           0          0
# Mrpl15                  0            0           0           0          0
# Lypla1                  0            0           0           0          0
# ERCC-00170              0            0           0           0          0
# ERCC-00171              0            0           0           0          0
# Gfp_transgene           1            0           2           2          2
# Cre_transgene           0            0           0           0          0
# Tdtom_transgene        20            1          89           0          0
# note: this function assumes each cell is placed in a well in 384 plate format. There is no duplation of well names.


ReadSTAR <- function(plate.name, folder) {
  files <- list.files(path = folder, pattern = "Gene.out.tab",recursive = T,full.names = T)
  data.raw<-data.frame()
  for (file in files) {
    cell.data<-read.table(file=file)
    cell.data<-cell.data[5:23436,]
    cell.data.df<-data.frame(counts=cell.data[,2])
    rownames(cell.data.df)<-cell.data$V1
    well.name<-unlist(strsplit(file,split = "/"))
    colnames(cell.data.df)<-paste(plate.name, "_",  well.name[length(well.name)-1], sep="")
    if(file==files[1]){
      data.raw <- cell.data.df
    }else{
      data.raw<-cbind(data.raw, cell.data.df)
    }
  }
  return(data.raw)
  
}

# LogCPM 
# calculate log(cpm+1) from umi counts for a seurat object
# input: scdata. seurat object. must have "RNA" assay and slot "counts"
# output: the same seurat object with log(cpm+1) added as "logCPM" slot under assay "RNA".

LogCPM <- function(scdata) {
  data.count <- as(scdata[["RNA"]]$counts, "sparseMatrix")
  data.cpm <- as(apply(data.count,2, function(x) (x/sum(x))*1000000), "sparseMatrix")
  scdata.logCPM <- as(log(data.cpm + 1), "sparseMatrix")
  scdata@assays$RNA$logCPM <- scdata.logCPM
  return(scdata)
  
}

# SeuratProcess
# run seurat data normalization with scTransform, vst.flavor = v2, run PCA with 30 pcs. run harmony across batch with 30 dims, run umap and run tsne with 30 dims, find neighbors and clusters.
# input: scdata. A seurat object
# input: batch.name. The metadata column name where batch label is stored.
# output: a seurat object with the described process.
SeuratProcess <- function(scdata, batch.name="key", npc=30, clustering.resolution=1) {
  
  # Split the object
  DefaultAssay(scdata) <- "RNA"
  object_list <- SplitObject(scdata, split.by = batch.name)
  object_list2 <- list()
  # Normalize counts data
  for (i in 1:length(object_list)) {
    object_list2[[i]] <- NormalizeData(object_list[[i]])
  }
  # Apply SCTransform to each object in the list
  object_list2 <- lapply(object_list2, function(x) {
    x <- SCTransform(x, verbose = TRUE) #Set verbose to true to see progress.
    return(x)
  })  

  # Merge the transformed objects
  scdata <- merge(object_list2[[1]], y = object_list2[-1])
  var.features <- lapply(object_list2, function(item) {
    item@assays$SCT@var.features
  })
  var.features <- unique(unlist(var.features))
  VariableFeatures(scdata) <- var.features
  # join RNA assay
  scdata[["RNA"]] <- JoinLayers(scdata[["RNA"]])
  # Scale data
  DefaultAssay(scdata) <- "RNA"
  scdata <- ScaleData(scdata)
  #PCA, harmony, clustering
  DefaultAssay(scdata) <- "SCT"
  scdata <-     RunPCA(scdata, npcs = npc, verbose = FALSE) %>%
    RunHarmony(dims=1:npc, verbose= FALSE, group.by.vars = batch.name) %>%
    RunUMAP(reduction = "harmony", dims = 1:npc, verbose = FALSE) %>%
    RunTSNE(dims = 1:npc, reduction = "harmony") %>%
    FindNeighbors(reduction = "harmony", dims = 1:npc, verbose = FALSE) %>%
    FindClusters(resolution = clustering.resolution)
  return(scdata)
}


# DotPlotMM
# plot dotplot using ensembl official gene symbol in allen insitute mouse single cell data.
# the gene symbol is looked up from ensemble look up table first.
# example ensembl table
# head(ensembl.table)
# ensembl_gene_id external_gene_name
# 1 ENSMUSG00000064336              mt-Tf
# 2 ENSMUSG00000064337            mt-Rnr1
# 3 ENSMUSG00000064338              mt-Tv
# 4 ENSMUSG00000064339            mt-Rnr2
# 5 ENSMUSG00000064340             mt-Tl1
# 6 ENSMUSG00000064341             mt-Nd1
# column names are ensembl_gene_id and external_gene_name. Order matters.
# scdata.MM is the allen mouse data converted to seurat. 

DotPlotMM <- function (scdata.MM, feature, ensembl.table, group.by = "subclass") {
  if (sum(feature %in% ensembl.table[,2])==length(feature)) {
    geneID <- c(ensembl.table[ensembl.table[,2] %in% feature, 1])
    dotplot <- DotPlot(scdata.MM, features = geneID, cols = c("lightgrey","darkred"), group.by=group.by)
    print(paste(feature, " = ", geneID))
    return(dotplot)
  } else {
    print("Some features not found in ensembl.table")
  }
  
}



#Function MergeCluster
#merge clusters based on DE score. Only clusters with DE score passing the set threshold will be mainteined. smaller clusters are merged to its nearest neighbour. 
#The function will use the "seurat_cluster" slot in the meta data to merge.
#Input: Object, padj.th, lfc.th, q.diff.th, de.score.th, min.cells
#Object, the seurat object to use. Must have "seurat_cluster" slot.


MergeCluster <- function(object, padj.th, lfc.th, q.diff.th, de.score.th, min.cells, assay="RNA") {
  #create a new column in metadata called cluster_merge for making merged cluster labels.
  cluster_merge<-object@meta.data$seurat_clusters
  #if cluster_merge already exists, remove the existing column.
  if (is.null(object@meta.data$cluster_merge)==FALSE) {
    print("Removing current cluster_merge column")
    object@meta.data<-object@meta.data[,c(names(object@meta.data)[names(object@meta.data)!="cluster_merge"])]
  }
  object@meta.data<-cbind(object@meta.data,cluster_merge)
  nclusters = length(unique(cluster_merge))
  
  # if number of clusters = 1, output only 1 cluster no merging
  if (nclusters == 1) {
    print("only 1 cluster, not merging")
  } else if (nclusters == 2)   { # if number of clusters = 2, test between the 2. 
    names = unique(cluster_merge)
    # if the smaller cluster has less than min.cells, merge the two.
    if (min(table(cluster_merge))<min.cells) {
      object@meta.data$cluster_merge[object@meta.data$cluster_merge == names[1]] <- names[2]
      print(paste("1 of the 2 clusters smaller than min.cells, Merged ", names[1], " with ", names[2], sep=""))
      object@meta.data$seurat_clusters<-object@meta.data$cluster_merge
    } else { # if both clusters are big enough, test de.score.
      marker.1 <- FindMarkers(object, ident.1 = names[1], ident.2 = names[2], logfc.threshold = lfc.th, verbose = FALSE, recorrect_umi=FALSE)
      marker.1<-marker.1[marker.1$p_val_adj<padj.th,]
      marker.1<-marker.1[abs(marker.1$pct.1-marker.1$pct.2)/apply(marker.1[,c("pct.1","pct.2")], 1, max)>q.diff.th,]
      norm.pvalue <- -log10(marker.1$p_val_adj)
      norm.pvalue[norm.pvalue>20]<-20
      DE.score<-sum(norm.pvalue)
      #if DE gene thresholds not met, merge to its nearest neighbour
      if (DE.score < de.score.th) {
        object@meta.data$cluster_merge[object@meta.data$cluster_merge == names[1]] <- names[2]
        print(paste("DE.score is smaller than DE.score.th, Merged ", names[1], " with ", names[2], sep=""))
        object@meta.data$seurat_clusters<-object@meta.data$cluster_merge
        #if DE gene thresholds met, keep both
      } 
    }
  
  } else {  # if number of clusters >= 3, find 2 nearest neighbours  to test.
    #loop through all cluster names
    merged_this_round <- TRUE
    
    # Start a while loop that continues as long as a merge occurs
    while (merged_this_round) {
      merged_this_round <- FALSE # Assume no merges will happen this round
      
      # Recalculate everything inside the loop to get the current state
      object <- FindVariableFeatures(object, verbose = FALSE)
      matrix.aggr <- AggregateExpression(object = object, features = object@assays[[assay]]@var.features, group.by = "cluster_merge")
      matrix.aggr <- as.matrix(matrix.aggr[[assay]])
      dist.mat <- as.matrix(dist(t(matrix.aggr), diag = T, upper = T))
      
      # If there's only one cluster left, stop
      if (nrow(dist.mat) <= 1) {
        break
      }
      
      # Find the closest pair of clusters based on the distance matrix
      # Ignore the diagonal (distance to self)
      min_dist_val <- min(dist.mat[dist.mat != 0])
      closest_pair_indices <- which(dist.mat == min_dist_val, arr.ind = TRUE)
      
      # In case of ties, just take the first one
      pair_to_merge <- rownames(dist.mat)[closest_pair_indices[1, ]]
      cluster_name_1 <- unlist(strsplit(pair_to_merge[1], split = "g"))[2]
      cluster_name_2 <- unlist(strsplit(pair_to_merge[2], split = "g"))[2]
      
      # Ensure they are not the same cluster
      if (cluster_name_1 == cluster_name_2) {
        dist.mat[closest_pair_indices[1,1], closest_pair_indices[1,2]] = 99999
        next
      }
      
      print(paste("Testing potential merge between clusters", cluster_name_1, "and", cluster_name_2))
      
      # Check if a merge is warranted based on your criteria (cell count, DE score, etc.)
      # This logic from your original function needs to be adapted here.
      
      # Example logic for a merge based on DE score:
      marker <- FindMarkers(object, ident.1 = cluster_name_1, ident.2 = cluster_name_2, logfc.threshold = lfc.th, verbose = FALSE, recorrect_umi = FALSE)
      marker <- marker[marker$p_val_adj < padj.th, ]
      marker <- marker[abs(marker$pct.1 - marker$pct.2) / apply(marker[, c("pct.1", "pct.2")], 1, max) > q.diff.th, ]
      norm.pvalue <- -log10(marker$p_val_adj)
      norm.pvalue[norm.pvalue > 20] <- 20
      DE.score <- sum(norm.pvalue)
      
      # If criteria met, merge them and set the flag to continue the loop
      if (DE.score < de.score.th) {
        object@meta.data$cluster_merge[object@meta.data$cluster_merge == cluster_name_2] <- cluster_name_1
        print(paste("Merged", cluster_name_2, "with", cluster_name_1))
        object@meta.data$seurat_clusters <- object@meta.data$cluster_merge
        merged_this_round <- TRUE
      }
    }
    }
  return(object)
}





#Function iter.clust will perform iterative clustering
#Input
#object must have meta data column "cluster_merge", and "seurat_cluster".
iter.clust<-function(object, de.param, var.features, res=1, assay="RNA") {
  DefaultAssay(object) <- assay
  object.tokeep<-list()
  is.increased=T
  object.tosubcluster <- SplitObject(object, split.by = "cluster_merge")
  k=0
  while (length(object.tosubcluster)>1) {
    print(paste("started iteration",k))
    object.tosubcluster.k<-object.tosubcluster
    print(object.tosubcluster.k)
    for (i in 1:length(object.tosubcluster)) {
          scdata <- object.tosubcluster.k[[i]]
          scdata <- ScaleData(scdata, verbose = FALSE)
          print("")
          print("*************************************")
          name.i <- names(object.tosubcluster.k)[i]
          print(paste(i+k, "th subclustering, cluster name:", name.i))
          print("*************************************")
          print("")
          scdata<-FindVariableFeatures(scdata, assay=assay)
          tmp.var.features <- VariableFeatures(scdata, assay = assay)
          

      if (length(tmp.var.features)>30) {
        if (ncol(scdata)<100) {
          res=1
        }
        print("PCA analysis")
        scdata <- RunPCA(scdata, npcs = 29, verbose=FALSE)
        #debug Error in irlba(A = t(x = object), nv = npcs, ...) : 
        #   max(nu, nv) must be strictly less than min(nrow(A), ncol(A))
        scdata <- FindNeighbors(scdata, dims = 1:29,reduction = "pca")
        print("clustering")
        scdata <- FindClusters(scdata, resolution = res)
        cluster_merge<-paste(names(object.tosubcluster)[i], ".", scdata@meta.data$seurat_clusters,sep="")
        scdata@meta.data<-scdata@meta.data[,1:(ncol(scdata@meta.data)-1)]
        scdata@meta.data$seurat_clusters<-cluster_merge
        if (length(unique(scdata@meta.data$seurat_clusters))==1) {
          print("Only 1 cluster, stop from further subclustering.")
          scdata@meta.data<-cbind(scdata@meta.data, cluster_merge=scdata@meta.data$seurat_clusters)
          object.tokeep[[name.i]]<-scdata
          object.tosubcluster[[name.i]] <- NULL
        } else {
          print("More than 1 clusters, Merging clusters based on DE gene threshold")
          scdata<-merge_cl(object=scdata, de.param, assay=assay)
          if (length(unique(scdata@meta.data$cluster_merge))==1) {
            print("Only 1 cluster after merging, stop from further subclustering.")
            object.tokeep[[name.i]]<-scdata
            object.tosubcluster[[name.i]] <- NULL
          } else {
            object.tosubcluster[[name.i]] <- scdata
          }
        }
      } else {
        print("Less features than 30, stop further subclustering.")
        object.tokeep[[name.i]]<-scdata
        object.tosubcluster[[name.i]] <- NULL
      }
      
      
    }
    object.tosubcluster.list <- map(object.tosubcluster, ~SplitObject(., split.by = "cluster_merge"))
    object.tosubcluster <- flatten(object.tosubcluster.list)
    k=k+1
  }
  # merge the last cluster to tokeep
  if (length(object.tosubcluster)==1){
    object.tokeep[[names(object.tosubcluster[1])]]<-object.tosubcluster[1]
  }
  # merge all clusters that can't be subdivided.
  object.final <- Merge_Seurat_List(list_seurat = object.tokeep)
  object.final@meta.data$seurat_clusters <- object.final@meta.data$cluster_merge
  return(object.final)
}

  

#Bootstraping
bootstrap_iter_clust <- function (seurat.obj, n.iter=100, de.param, var.features = var.features, res=1, assay="RNA") {
  #create a list containing coclustering matrix of each round
  n.cluster <- c()
  cocluster.matrix<-data.frame()
  start.time<-Sys.time()
  for (i in 1:n.iter) {
    print("")
    print("*************************************")
    print("*************************************")
    print(paste("Start iteration",i))
    print(paste("Progress:",i/n.iter*100,"%"))
    print(paste("Elapsed time:", difftime(Sys.time(), start.time, units = "mins"), " mins"))
    print("*************************************")
    print("*************************************")
    print("")
    #subsample by 80%
    print("subsampling 80%")
    names<-colnames(seurat.obj)
    pred=sample(names, ncol(seurat.obj)*0.8)
    print(head(pred))
    scdata.sample<-seurat.obj[,pred]
    #iterative clustering the subsampled data.
    scdata.sample <- DietSeurat(scdata.sample)
    print("Start initial clustering")
    all.genes<-rownames(rownames(scdata.sample))
    scdata.sample <- ScaleData(scdata.sample, features = all.genes)
    scdata.sample <- RunPCA(scdata.sample, assay = assay, features = var.features, verbose = FALSE, npcs = 29)
    scdata.sample <- RunTSNE(scdata.sample, dims = 1:29, reduction = "pca" ,check_duplicates = FALSE, verbose = FALSE)
    scdata.sample <- RunUMAP(scdata.sample, dims = 1:29, reduction = "pca" ,verbose = FALSE)
    scdata.sample <- FindNeighbors(scdata.sample, dims = 1:29, verbose = FALSE,reduction = "pca")
    scdata.sample <- FindClusters(scdata.sample, verbose = FALSE, res=res)
    print("test initial clustering stability")
    scdata.sample <- merge_cl(object=scdata.sample, de.param, assay=assay)
    print("Start iterative clustering.")
    scdata.sample <- iter.clust(object = scdata.sample, de.param, var.features = var.features, res=res, assay=assay)
    #save cocluster matrix
    print("Generat co-clustering matrix.")
    cluster.i<-data.frame(name=rownames(scdata.sample@meta.data),cluster=paste(scdata.sample@meta.data$cluster_merge,".iter",i, sep=""))
    cocluster.matrix<-rbind(cocluster.matrix, cluster.i)
    print(paste("iteration",i,"done, merged data with the last iteration"))
    n.cluster.i <- length(unique(scdata.sample@meta.data$cluster_merge))
    n.cluster <- c(n.cluster, n.cluster.i)
  }
  print(paste("finished at",Sys.time()))
  print(paste("Total elapsed time:",difftime(Sys.time(), start.time, units = "mins"), " mins"))
  cocluster.matrix.average <- crossprod(table(cocluster.matrix[c(2,1)]))
  boot.list<-list(cocluster.matrix.average, n.cluster)
  return(boot.list)
}

# parallele version
bootstrap_iter_clust_parallel <- function (seurat.obj, n.iter = 100, de.param, var.features, res = 1, ncore=4, assay = "RNA") {

  # Use all available cores except one
  # Adjust the number of cores as needed
  num_cores <- ncore
  registerDoParallel(num_cores)

  # Replace the for loop with foreach
  results <- foreach(i = 1:n.iter, .combine = 'rbind', .packages = c("Seurat")) %dopar% {
    # Subsample by 80%
    names <- colnames(seurat.obj)
    pred <- sample(names, ncol(seurat.obj) * 0.8)
    scdata.sample <- seurat.obj[, pred]
    
    # Process the subsampled data
    scdata.sample <- DietSeurat(scdata.sample)
    all.genes <- rownames(scdata.sample)
    scdata.sample <- ScaleData(scdata.sample, features = all.genes)
    scdata.sample <- RunPCA(scdata.sample, assay = assay, features = var.features, verbose = FALSE, npcs = 29)
    scdata.sample <- RunTSNE(scdata.sample, dims = 1:29, reduction = "pca", check_duplicates = FALSE, verbose = FALSE)
    scdata.sample <- RunUMAP(scdata.sample, dims = 1:29, reduction = "pca", verbose = FALSE)
    scdata.sample <- FindNeighbors(scdata.sample, dims = 1:29, verbose = FALSE, reduction = "pca")
    scdata.sample <- FindClusters(scdata.sample, verbose = FALSE, res = res)
    
    # Assume merge_cl() and iter.clust() are defined elsewhere
    scdata.sample <- merge_cl(object = scdata.sample, de.param, assay = assay)
    scdata.sample <- iter.clust(object = scdata.sample, de.param, var.features = var.features, res = res, assay = assay)
    
    # Generate the cocluster matrix for this iteration
    cluster.i <- data.frame(name = rownames(scdata.sample@meta.data),
                            cluster = paste(scdata.sample@meta.data$cluster_merge, ".iter", i, sep = ""))
    
    # Return the data frame and number of clusters for this iteration
    list(cluster.i, length(unique(scdata.sample@meta.data$cluster_merge)))
  }
  
  # Stop the parallel cluster
  stopImplicitCluster()
  
  # Post-processing outside the loop
  cocluster.matrix <- do.call(rbind, lapply(results, `[[`, 1))
  n.cluster <- unlist(lapply(results, `[[`, 2))
  
  cocluster.matrix.average <- crossprod(table(cocluster.matrix[c(2,1)]))
  boot.list <- list(cocluster.matrix.average, n.cluster)
  
  return(boot.list)
}


# plot co-clustering matrix

plot.cocluster <- function(co.result) {
  MO.comatrix <- co.result[[1]]
  #remove boundary cells.
  MO.comatrix<- MO.comatrix[rowMaxs(MO.comatrix)>50,]
  MO.comatrix<- MO.comatrix[,colMaxs(MO.comatrix)>50]
  cocluster.heatmap <- pheatmap(MO.comatrix, show_colnames = F, show_rownames = F, cutree_rows = median(co.result[[2]]), cutree_cols = median(co.result[[2]]))
  return(cocluster.heatmap)
}


# add bootstrapped cluster id to meta data seurat.

assign.bootedclusterID <- function(co.result, seurat.obj) {
  MO.comatrix <- co.result[[1]]
  #remove boundary cells.
  MO.comatrix<- MO.comatrix[rowMaxs(MO.comatrix)>50,]
  MO.comatrix<- MO.comatrix[,colMaxs(MO.comatrix)>50]
  #cluster the co-cluster matrix to the median number of clusters.
  cocluster.heatmap <- pheatmap(MO.comatrix, show_colnames = F, show_rownames = F, cutree_rows = median(co.result[[2]]), cutree_cols = median(co.result[[2]]))
  concensus.cluster<- cutree(cocluster.heatmap$tree_row, k = median(co.result[[2]]))
  concensus.cluster.df<-data.frame(concensus.cluster)
  rownames(concensus.cluster.df) <- sub(pattern = "_", replacement = "", x = rownames(concensus.cluster.df))
  #subset the cells and add cluster name to metadata
  scdata.labeled <- subset(seurat.obj, cells = rownames(concensus.cluster.df))
  concensus.cluster.df<- concensus.cluster.df[rownames(scdata.labeled@meta.data),]
  scdata.labeled@meta.data<-cbind(scdata.labeled@meta.data, cluster_merge_booted=concensus.cluster.df)
  return(scdata.labeled)
}





#save seurat to anndata

save.h5ad <- function(seurat.object, assay="RNA", layer="counts", filename){
  counts_matrix <- GetAssayData(seurat.object, assay = "RNA", layer = "counts")
  obs_metadata <- seurat.object@meta.data
  gene_names <- rownames(counts_matrix)
  var_metadata <- data.frame(
    gene_names = gene_names,
    row.names = gene_names
  )
  # Create AnnData object with the matrix and metadata
  adata <- AnnData(
    X = t(counts_matrix),
    obs = obs_metadata,
    var = var_metadata
  )
  write_h5ad(adata, filename)
}


# helper function to test if a pair of clusters needs to be merged based on de.param.

test_merge <- function(de.pair, de.param, merge.type="undirectional")
{
  if(length(de.pair)==0){
    return(TRUE)
  }
  to.merge = FALSE
  if(merge.type=="undirectional"){
    if(!is.null(de.param$de.score.th)){
      to.merge=de.pair$score < de.param$de.score.th
    }
    if(!to.merge & !is.null(de.param$min.genes)){
      to.merge=de.pair$num < de.param$min.genes
    }
  }
  else{
    if(!is.null(de.param$de.score.th)){
      to.merge=de.pair$up.score < de.param$de.score.th | de.pair$down.score < de.param$de.score.th
    }
    if(!to.merge & !is.null(de.param$min.genes)){
      to.merge=de.pair$up.num < de.param$min.genes | de.pair$down.num < de.param$min.genes
    }
  }
  return(to.merge)
}

# From hicat merge clusters, convert to Seurat

merge_cl<- function(#norm.dat,
                    object,
                    #cl, 
                    #rd.dat = NULL,
                    #rd.dat.t = NULL, 
                    de.param = de_param(), 
                    merge.type = c("undirectional","directional"), 
                    #max.cl.size = 300,
                    #de.method = "limma",
                    #de.genes = NULL, 
                    #return.markers = FALSE,
                    #pairBatch =40,
                    #verbose = 0,
                    assay="RNA"
                    )
{
  print("start merging")
  #create a new column in metadata called cluster_merge for making merged cluster labels.
  cl<-droplevels(as.factor(object@meta.data$seurat_clusters))
  init_cluster<-length(unique(cl))
  #if cluster_merge already exists, remove the existing column.
  if (is.null(object@meta.data$cluster_merge)==FALSE) {
    print("Removing current cluster_merge column")
    object@meta.data<-object@meta.data[,c(names(object@meta.data)[names(object@meta.data)!="cluster_merge"])]
  }
  object@meta.data<-cbind(object@meta.data, cluster_merge=cl)
  merge.type=merge.type[1]
  ###Merge small clusters with the closest neighbors first.
  print("merge small clusters")
  while(TRUE){
    cl<-droplevels(as.factor(object@meta.data$cluster_merge))
    cl.size = table(cl)
    #if only 1 cluster, break
    if(length(cl.size)==1){
      break
      print("only 1 cluster, no merging")
    }
    # if no small cluster, break
    cl.small =  names(cl.size)[cl.size < de.param$min.cells]
    if(length(cl.small)==0){
      break
      print("no small clusters")
    }
    object <- FindVariableFeatures(object, verbose = FALSE)
    matrix.aggr <- AggregateExpression(object = object, features = object@assays[[assay]]@var.features, group.by = "cluster_merge")
    matrix.aggr <- as.matrix(matrix.aggr[[assay]])
    dist.mat <- as.matrix(dist(t(matrix.aggr), diag = F, upper = T))
    rownames(dist.mat)<- gsub("g", "", rownames(dist.mat))
    colnames(dist.mat)<- gsub("g", "", colnames(dist.mat))
    cl.sim = dist.mat
    smallest <- names(cl.size)[cl.size==min(cl.size)][1]
    tmp <- cl.sim[smallest,]
    nn.tmp <- names(tmp)[tmp==min(tmp[tmp>0])]
    print(paste("merging ",nn.tmp," and ", smallest))
    object@meta.data$cluster_merge[object@meta.data$cluster_merge == smallest] <- nn.tmp
    cl<-droplevels(as.factor(object@meta.data$cluster_merge))
  }		
  print("Merging DE.score < DE.score.th clusters")
  while(length(unique(cl)) > 1) { # only continue the following step if clusters > 1 after the above step.
    cl<-droplevels(as.factor(object@meta.data$cluster_merge))
    if(length(unique(cl)) == 2) { # if 2 cluster from the previous step, compare the 2 with de.score.
      print("only 2 clusters")
      names = unique(cl)
      marker.1 = FindMarkers(object, ident.1 = names[1], ident.2 = names[2], logfc.threshold = de.param$lfc.th, verbose = FALSE, recorrect_umi=FALSE, group.by = "cluster_merge")
      marker.1 = marker.1[marker.1$p_val_adj<de.param$padj.th,]
      marker.1 = marker.1[abs(marker.1$pct.1-marker.1$pct.2)/apply(marker.1[,c("pct.1","pct.2")], 1, max)>de.param$q.diff.th,]
      norm.pvalue = -log10(marker.1$p_val_adj)
      norm.pvalue[norm.pvalue>20]<-20
      DE.score = sum(norm.pvalue)
      #if DE gene thresholds not met, merge to its nearest neighbour
      if (DE.score < de.param$de.score.th) {
        object@meta.data$cluster_merge <- names[2]
        print(paste("DE.score is smaller than DE.score.th, Merged ", names[1], " with ", names[2], sep=""))
        object@meta.data$seurat_clusters<-object@meta.data$cluster_merge
        #if DE gene thresholds met, keep both
      } else {
        print("DE.score met, not merging")
        break
      }
    }
    
    else { # if 3 or more clusters from the previous step, compare nearest neighbours.
      print("more than 2 clusters, testing DE.score.th")
      #calculate cluster distance and determine neighbours.
      object <- FindVariableFeatures(object, verbose = FALSE)
      matrix.aggr <- AggregateExpression(object = object, features = object@assays[[assay]]@var.features, group.by = "cluster_merge")
      matrix.aggr <- as.matrix(matrix.aggr[[assay]])
      dist.mat <- as.matrix(dist(t(matrix.aggr), diag = F, upper = T))
      rownames(dist.mat)<- gsub("g", "", rownames(dist.mat))
      colnames(dist.mat)<- gsub("g", "", colnames(dist.mat))
      clusters <- rownames(dist.mat)
      # Use lapply to iterate through each cluster and create a small data frame
      results_list <- lapply(clusters, function(cluster) {
        # Sort the distances for the current cluster
        sorted_distances <- sort(dist.mat[cluster, ], index.return = TRUE)
        # Get the names of the closest and second-closest neighbours (indices 2 and 3)
        neighbours <- names(sorted_distances$x)[2:3]
        # Create a data frame for this cluster and its two neighbours
        df <- data.frame(
          Cluster = rep(cluster, 2),
          Neighbor = neighbours,
          row.names = NULL
        )
        return(df)
      })
      # Combine the list of data frames into a single data frame
      tomerge.df <- do.call(rbind, results_list)
      # remove duplicated cluster pairs
      sorted_pairs <- data.frame(
        pmin(tomerge.df$Cluster, tomerge.df$Neighbor),
        pmax(tomerge.df$Cluster, tomerge.df$Neighbor)
      )
      # Identify the rows that are duplicated based on the sorted pairs
      unique_rows <- !duplicated(sorted_pairs)
      # Filter the original data frame to keep only the unique rows
      tomerge.df <- tomerge.df[unique_rows, ]
      tomerge.df$merge <- rep("FALSE", nrow(tomerge.df))
      # determine if each pair of clusters should be merged based on DE.score.
      for (j in 1:nrow(tomerge.df)) {
        names=c(tomerge.df[j,1], tomerge.df[j,2])
        marker.1 <- FindMarkers(object, ident.1 = names[1], ident.2 = names[2], logfc.threshold = de.param$lfc.th, verbose = FALSE, recorrect_umi=FALSE, group.by = "cluster_merge")
        marker.1<-marker.1[marker.1$p_val_adj < de.param$padj.th,]
        marker.1<-marker.1[abs(marker.1$pct.1-marker.1$pct.2)/apply(marker.1[,c("pct.1","pct.2")], 1, max) > de.param$q.diff.th,]
        norm.pvalue <- -log10(marker.1$p_val_adj)
        norm.pvalue[norm.pvalue > 20] <- 20
        DE.score <- sum(norm.pvalue)
        print(paste(names[1],names[2],DE.score))
        #if DE gene thresholds not met, merge to its nearest neighbour
        if (DE.score < de.param$de.score.th) {
          tomerge.df$merge[j] <- TRUE
          print(paste("DE.score is smaller than DE.score.th, Merged ", names[1], " with ", names[2], sep=""))
          #if DE gene thresholds met, keep both
        }
      }
      # use the resulting tomerge.df to determine which clusters to merge, as below.
      # if no merging break out.
      if (sum(as.logical(tomerge.df$merge))==0) {
        print(paste("All DE.score pass threshold. No cluster needs to be merged."))
        break
      }
      # merge clusters.
      # create cluster pairs to merge
      tomerge.df.do <- tomerge.df[tomerge.df$merge==TRUE, c("Cluster", "Neighbor")]
      # Use graph to finalize which clusters to merge together.
      g <- graph_from_data_frame(tomerge.df.do, directed = FALSE)
      # Get the connected components of the graph
      components <- components(g)
      # Find the lowest ID cluster for each component
      find_community <- function(components) {
        df = data.frame()
        for (i in unique(components$membership)) {
          # Get all clusters in the current component
          clusters_in_component <- names(components$membership[components$membership == i])
          # Find the lowest numerical ID among them (the "root")
          root_name <- clusters_in_component[which.min(as.numeric(gsub("g", "", clusters_in_component)))]
          # Return a vector of the root and all other clusters in the component
          other_members <- clusters_in_component[clusters_in_component != root_name]
          # Return a two-column data frame for this component
          df.i = data.frame(Cluster = rep(root_name, length(other_members)), Neighbor = other_members)
          df = rbind(df, df.i)
        }
        return(df)
      }
      tomerge.df.do <- find_community(components)
      tomerge.df.do <- tomerge.df.do[order(tomerge.df.do$Cluster),]
      tomerge.df.do$Cluster<- gsub("g", "", tomerge.df.do$Cluster)
      tomerge.df.do$Neighbor<- gsub("g", "", tomerge.df.do$Neighbor)
      # perform the merge
      for (l in 1:nrow(tomerge.df.do)) {
        print(paste("DE.score is smaller than DE.score.th, Merged ", tomerge.df.do$Neighbor[l], " with ", tomerge.df.do$Cluster[l], sep=""))
        object@meta.data$cluster_merge[object@meta.data$cluster_merge == tomerge.df.do$Neighbor[l]] <- tomerge.df.do$Cluster[l]
      }
    }
  }
  print(paste("Merging done. Merged ", init_cluster, "clusters into", length(unique(object@meta.data$cluster_merge))))
  object@meta.data$seurat_clusters <- object@meta.data$cluster_merge
  return(object)
  
}



  
