library(R.matlab)
library(Hmisc)
library(pheatmap)
library(RColorBrewer)
library(dbscan)
rm(list=ls())


fdir = "../matlab/DataTables/bootstrap_Both_WithParfor_Selection_Progbar__maxfcnt1001k_100pat_0.01nstd_se-ard.kern/"
dir.create("../Correlation_plots/")
dirname = paste0("../Correlation_plots/CorrPlots_DBSCAN_findKeepAgain_", substr(fdir,start = 22,stop = nchar(fdir)))
dir.create(dirname)

clusternum = c(1,2,3,4,5,6,7,8,9.10)
startnum = 1 # only used if we are not doing all the tumors, 3-5 is LIHC, UCEC, BLCA
sigvalue = 5e-2
paletteLength <- 64
tickfont = 5
mingroupsize = 5;


getDataMat <- function(dataList){
  type = dataList$type[[1]][[1]][1]
  genes = unlist(lapply(dataList$genes,function(x){
    x[[1]][[1]][1]
  }))
  corr = rcorr(t(dataList$ard),type = "spearman")
  spmat = corr$r
  pmat = corr$P
  spdist = as.dist((1-spmat)/2) # make the distance between 0 and 1
  dbclust = hdbscan(x = spdist,minPts = mingroupsize)
  clusterGeneInd = which(sort(dbclust$cluster,index.return=T)$ix > 0) # do not show genes hdbscan discarded as "noise points"
  clusterGeneName = genes[clusterGeneInd]
  spdist_clust = spmat[clusterGeneInd,clusterGeneInd]
  spmat_clust = spmat[clusterGeneInd,clusterGeneInd]
  
  # sort the data for visual inspection
  sorted = apply(dataList$ard,2,
                 function(x){check = sort(x,index.return=T)
                 check$ix})
  sortedVals = apply(dataList$ard,2,
                     function(x){check = sort(x,index.return=T)
                     check$x})
  
  pmat_clust = pmat[clusterGeneInd,clusterGeneInd]
  clusters = sort(dbclust$cluster,index.return=T)$x[which(sort(dbclust$cluster,index.return=T)$x > 0)]
  
  colnames(spmat) = rownames(spmat) = colnames(pmat) = rownames(pmat) = genes
  colnames(spmat_clust) = rownames(spmat_clust) = colnames(spdist_clust) = rownames(spdist_clust) = colnames(pmat_clust) = rownames(pmat_clust) = clusterGeneName
  return(list(type=type,spcor=spmat_clust,P=pmat_clust,ard=dataList$ard,sortedArd=sorted,sortedVals=sortedVals,acc=dataList$loo,geneNames = clusterGeneName))
  
}


getnames <- function(data){
  names = c()
  for(i in 1:dim(data)[3]){
    names[i] = data[,,i]$type[[1]][[1]][1]
  }
  return(names)
}
files = list.files(fdir)
data = list()
namesvec = c()
for(fname in grep(".mat",files)){
  allData = readMat(paste0(fdir,files[fname]))
  #allData = readMat("./UMAP_Plots_GP_THICK_OneDFI_Bootstrap/meta.mat")
  
  namesvec[fname] =   getnames(allData[[1]])
  data[[fname]] = getDataMat(allData[[1]][,,1])
  
}  

names(data) = namesvec

for(i in startnum:length(data)){
  data[[i]]$spcor[which(data[[i]]$P > sigvalue)] = NA
  
  pdf(paste0(dirname, "/Hist_", names(data)[i],".pdf"))
  par(mfrow=c(3,1))    # set the plotting area into a 1*2 array
  hist(data[[i]]$sortedArd[which(data[[i]]$geneNames == "WNT2"),],breaks = 10,main = "SE-ARD Lengthscale Ranks",xlab = paste0("WNT2 Lengthscale Rank, ",names(data)[i]),ylim=c(0,1000))
  hist(data[[i]]$sortedArd[which(data[[i]]$geneNames == "WNT11"),],breaks = 10,main = "",xlab = paste0("WNT11 Lengthscale Rank, ",names(data)[i]),ylim=c(0,1000))
  hist(data[[i]]$sortedArd[which(data[[i]]$geneNames == "WNT5A"),],breaks = 10,main = "",xlab = paste0("WNT5A Lengthscale Rank, ",names(data)[i]),ylim=c(0,1000))
  dev.off()
  
  sortedNames = rep("",length(data[[i]]$geneNames))
  for(ii in 1:ncol(data[[i]]$sortedArd)){
    names = c()
    for(jj in 1:nrow(data[[i]]$sortedArd)){
      names[jj] = data[[i]]$geneNames[data[[i]]$sortedArd[jj,ii]]
    }
    sortedNames = cbind(sortedNames,names)
  }
  
  means = rowMeans(data[[i]]$sortedArd)
  annotations = data.frame(Mean_rank = means)
  rownames(annotations) = data[[i]]$geneNames
  for(j in 1:length(clusternum)){
    
    pdf(paste0(dirname, "/Dendro_", names(data)[i], "_", clusternum[j],".pdf"))
    pheatmap(data[[i]]$spcor, clustering_distance_rows = "euclidean", 
             treeheight_col = 0, 
             cutree_rows = clusternum[j],
             na_col = "black",
             fontsize = tickfont,annotation_row = annotations)
    
    dev.off()
  }
}  


