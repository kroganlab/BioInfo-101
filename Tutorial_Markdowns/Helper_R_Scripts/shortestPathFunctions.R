# shortest path type beat
# functions included here:
#
# find_shortestPaths
# gen_NodeInfo
# gen_RandomShortestPaths
# normalize_ShortestPaths

### LOAD PACKAGES / HELPER FUNCTIONS###
library(data.table)
library(igraph)
library(parallel)
library(doParallel)
library(pbapply)
source("C:\\Users\\isido\\github\\TBRU_network_analysis\\Helper_R_Scripts\\enrichmentTestFunctions.R")
pathToData <- "C:\\Users\\isido\\github\\TBRU_network_analysis\\Input_Data\\enrichment\\"

find_ShortestPaths <- function(startpoints, targets, graph, distances = NULL, trackStartEnds = FALSE){    # startpoints and targets are both character vectors of nodes in the graph, graph is an igraph object, distances is a 
                                                                                  # distance matrix of the  
  if (is.null(distances)){
    distances <- igraph::distances(graph)
  }
  
  closestTarget <- data.table(interactor = startpoints)
  closestTarget[, knockdown := lapply(closestTarget$interactor, FUN = function(x, candidates, distances){
    candidateNodeDistances <- distances[x,candidates][!distances[x,candidates] == 0]
    return (data.table(candidateNodeDistances, names(candidateNodeDistances))[candidateNodeDistances == min(candidateNodeDistances),V2] )
  }, candidates = targets, distances = distances)]
  
  
  paths <- lapply( 1:nrow(closestTarget), FUN = function(i, targets, graph){
    tempPath <- igraph::shortest_paths( graph, closestTarget[i, interactor], closestTarget[i, knockdown][[1]] )$vpath
    }, targets = closestTargets, graph = graph)
  print(length(paths))
  paths <- unlist(paths, recursive = F)
  
  nodes <- unique(names(unlist(paths)))
  
  edges <- unlist(lapply(paths, FUN = function(path){
    return(lapply(1:(length(path) -1), FUN = function(n, path){
      if (length(path) == 0){
        return(NULL)
      } else{ return( c(names(path[n]),names(path[n+1])) ) }
    }, path = path))
  }))
  
  if (!trackStartEnds){
    return(list(paths = paths, nodes = nodes, edges = edges))
  } else{
    
    startEnds <- suppressWarnings( do.call(rbind, lapply(paths, FUN = function(p){
      names <- names(p)
      len <- length(names)
      if (len > 2){
        return(data.table("inpath" = c(NA, names[2:(len-1)]), "start" = names[1], "end" = names[len] ))
      } else{
        return(data.table("inpath" = c(NA), "start" = names[1], "end" = names[len] ))
      }
    })))
    
    startEnds <- startEnds[!is.na(start)]
      
    return(list(paths = paths, nodes = nodes, edges = edges, startEnds = startEnds))
  }
}


gen_NodeInfo <- function(pathsobj){
  
  pathsobj$edgemat <- matrix(pathsobj$edges, ncol = 2, byrow = T)
  temp.graph <- igraph::simplify(igraph::graph_from_edgelist(pathsobj$edgemat) )
  
  nodeInfo <- data.table(string = pathsobj$nodes, betweenness = igraph::betweenness(temp.graph), degree = igraph::degree(temp.graph), eigen = igraph::evcent(temp.graph)$vector, reach = (ego_size(temp.graph, 2)-1)/(vcount(temp.graph)-1))
  nodeInfo[data.table(pathsobj$edges)[, .N, by = V1], npaths := i.N, on = .(string = V1)] 
  
  return(nodeInfo)
}


gen_RandomShortestPaths <- function(nstarts, ntargets, niter, distances, graph){

  # Setting up multithreaded process to speed up run time -- this is written to work for Windows and Mac ( I think Linux will work also..?)
  numCores <- round(parallel::detectCores() * .70)
  
  # Prepare the Clustering
  if (Sys.info()["sysname"] == "Windows"){
    cl <- makeCluster(numCores) 
    registerDoParallel(cl)
    clusterExport(cl, list( "find_ShortestPaths", "nstarts", "ntargets", "distances", "graph", "gen_NodeInfo" ), envir = environment())
    clusterEvalQ(cl,library("igraph"))
    clusterEvalQ(cl,library("data.table"))
    
    # for n iterations, pick x random starts and y random targets as supplied at the start of fxn, find shortest paths as above, and return centrality measures for each node in nodeInfo
    simulated <- pbapply::pblapply( 1:niter, FUN = function(n){
      paths <- suppressWarnings( find_ShortestPaths(rownames(distances)[sample.int(nrow(distances), nstarts)], 
                                                   rownames(distances)[sample.int(nrow(distances), ntargets)], 
                                                   graph, distances) ) 
      return(gen_NodeInfo(paths))
    }, cl = cl)
    stopCluster(cl)
  } else {
    simulated <- pbapply::pblapply( 1:niter, FUN = function(n){
      paths <- suppressWarnings( find_ShortestPaths(rownames(distances)[sample.int(nrow(distances), nstarts)], 
                                                    rownames(distances)[sample.int(nrow(distances), ntargets)], 
                                                    graph, distances) ) 
      return(gen_NodeInfo(paths))
    }, cl = numCores)
  }
  
  return(do.call(rbind, simulated))
}

normalizeToNetwork <- function(){
  
}

### THESE PVALUE FUNCTIONS ARE ALL MESSED UP
calcPval <- function(value, dist, cond = "greater", background = NULL){
  if (!is.null(background)){
    zeros <- background - length(dist)
    dist <- c(dist, rep(0, zeros))
  }
  
  if (cond == "greater"){
    dist <- as.numeric(dist)
    pv <- (length(dist[dist > value]) + 1 ) / (length(dist)+1)
  }
  return(pv)
}


applyPval <- function(valuecol, applyTo, background, keycol = "string"){
  #maybe incl p.adjust?
  pvals <- lapply(applyTo[[keycol]], FUN = function(key){
    return(calcPval( applyTo[ applyTo[[keycol]] == key, ][[valuecol]], as.numeric(background[ background[[keycol]] == key, ][[valuecol]]), background = 1000) )
  } )
  return(pvals)
}


imputeToMin <- function(vect, includeZeros = T){
  
  if (includeZeros == T){
    vect[is.na(vect) | vect == 0] <- min(vect[vect  != 0 ][is.finite(vect[vect  != 0 ])])
  } else {
    vect[is.na(vect)] <- min(vect[is.finite(vect)])
  }
  
  return(vect)
}


normalize_ShortestPaths <- function(simulatedInfo, testInfo, useSim = TRUE, keycol = "string"){
  
  if (useSim){
    meanSimulatedInfo <- copy(simulatedInfo)
    meanSimulatedInfo[, betweenness := imputeToMin(betweenness)]
    meanSimulatedInfo <- meanSimulatedInfo[,.(betwn = mean(betweenness), degree = mean(degree), eigen = mean(eigen), reach = mean(reach), npaths = mean(npaths)), by = string]
    meanSimulatedInfo <- meanSimulatedInfo[ meanSimulatedInfo[[keycol]] %in% testInfo[[keycol]] ]
    
    testInfo[, c("normBetwn", "normDegree", "normEigen", "normReach", "normPaths") := .(betweenness/meanSimulatedInfo$betwn, degree/meanSimulatedInfo$degree, eigen/meanSimulatedInfo$eigen, reach/meanSimulatedInfo$reach, npaths/meanSimulatedInfo$npaths)]
    
    pvals <- do.call(cbind, lapply( c("betweenness", "degree", "eigen", "reach", "npaths"), FUN = function(measureCol){
      pvalCol <- data.table( applyPval(measureCol, testInfo, simulatedInfo, keycol = keycol))
      setnames(pvalCol, "V1", paste0(measureCol, "_pval"))
      return(pvalCol)
    }))
    
    testInfo <- cbind(testInfo, pvals)
    return(testInfo)
  } else {
    
  }
}


shortenChar <- function(stringVect, numchar=60){
 return( lapply(stringVect, function(string){
    if (nchar(string) < numchar) {
      return(string)
    } else {
      return(string[1:numchar])
    }
  }) )
}


prepare_subnetworkForEnrichment <- function(geneGroups, geneColName, groupColName, geneLabel = "symbol", subnet = "Proteome", reassign = FALSE){
  
  setnames(geneGroups, c(geneColName, groupColName), c("Protein", "Group"))
  groups <- geneGroups
  # Prepare gene groups
#  humanUPIDmap <- fread("..\\Input_Data\\enrichment\\humanSymbol-uniprotQueryResults.tsv")
  
#  
#  if (sum(geneGroups[[geneColName]] %in% humanUPIDmap$To) < sum(geneGroups[[geneColName]] %in% humanUPIDmap$From)){
#    setnames(geneGroups, c(geneColName, groupColName), c("Protein", "Group"))
#    
#    uniprotGroups <- copy(geneGroups)
#    uniprotGroups[humanUPIDmap, upid := i.To, on = .(Protein = From)] [, Protein := NULL]
#    setnames(uniprotGroups, c("upid"), c("Protein"))
#    
#    groups <- list("gene" = geneGroups, "uniprot" = uniprotGroups) 
#  } else {
#    setnames(geneGroups, c(geneColName, groupColName), c("Protein", "Group"))
#    symbolGroups <- copy(geneGroups)
#   symbolGroups[humanUPIDmap, symbol := i.From, on = .(Protein = To)] [, Protein := NULL]
#    setnames(symbolGroups, c("symbol"), c("Protein"))
#    
#    groups <- list("gene" = symbolGroups, "uniprot" = geneGroups)
#  }
  
  
  # Load default gene universe
  if (!exists("unv") | reassign == T){
      stringMap <-fread( paste0(pathToData, "..\\mapping\\swissprot-STRING.csv.gz")) #[, `Gene Names` := tstrsplit(`Gene Names`, split = ";| " , keep = 1)]
      ref <- unique(net$protein1)
      unv <<- stringMap[ From %in% ref, From]
  }
  
  # Load gene sets
  if (!exists("gmt") | reassign == T){
    
    if (!exists("aliases")){
      aliases <- fread("https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz")
      setnames(aliases, "#string_protein_id", "stringID")
    }
    kegg <- fread(paste0(pathToData, "KEGGgmt.csv.gz"))[aliases, string := i.stringID, on = .(gene=alias)] [, gene := NULL]  
    go <- fread(paste0(pathToData, "GOgmt.csv.gz"))[aliases, string := i.stringID, on = .(gene=alias)] [, gene := NULL]  
    c2 <- fread(paste0(pathToData, "GSEA.C2gmt.csv.gz"))[aliases, string := i.stringID, on = .(gene=alias)] [, gene := NULL]  
    ipa <- fread(paste0(pathToData, "IPAgmt.csv.gz"))[aliases, string := i.stringID, on = .(gene=alias)] [, gene := NULL]  
    ipap <- fread(paste0(pathToData, "IPA_pathways_KEGG_MSigDBgmt.csv.gz"))[aliases, string := i.stringID, on = .(gene=alias)] [, gene := NULL]  
    
    kegg <- kegg[!is.na(string), .(ont, string)] 
    go <- go[!is.na(string), .(ont, string) ] 
    
    
    
    gmt <<- lapply( list("kegg" = kegg, "go" = go, "c2" = c2, "ipa" = ipa, "ipap" = ipap) , FUN = function(db){
      setnames(db, "string", "gene")
      return(db)
    })
  }
  
  return(groups)
}

calculate_subnetworkEnrichment <- function(groups, proteinUnv = NULL, additionalGMT = NULL){
  
  if (!is.null(proteinUnv)){
    unv <- proteinUnv
  } 
  if (!is.null(additionalGMT)){
    if (exists("gmt")){
      gmt <- c(gmt, additionalGMT)
    } else {
      gmt <- additionalGMT
    }
  }
  if (!exists("gmt") | !exists("unv")){ stop("No gmt or unv data structures found in environment, try running prepare_subnetworkForEnrichment() or supply in optional arguments") }
  
  enrichments <- lapply(names(gmt), function(db, groups, unvs){
      univ <- unvs
      gg <- groups
    return( enricherOnGroups(groupTable = gg, geneColumn = "Protein", groupColumns = c("Group"), term2gene.gmt = gmt[[db]], universe = univ))
  }, groups = groups, unvs = unv )
  names(enrichments) <- names(gmt)

  return(enrichments)
}


substrRight <- function(x,n){
  substr(x,nchar(x)-n+1, nchar(x))
}


enrichmentNicelyFormatted <- function(enrichment , topn = 2, dataSrc = "", groupSrc = "", otherAnt = "", save = FALSE, subDir = "", dims = NULL, ...){
  title <- paste(dataSrc, groupSrc)
  hminfo <- enrichHeatmapBestPerGroup(enrichment, NULL, groupColumn = "Group", topN = topn, max_pAdjust = 0.05, cluster_columns = FALSE, title = title, ...)

  if (save){
    Prefix = paste(dataSrc, groupSrc, "enrichmentHeatmap", paste0("topn=",topn), sep = "_")
    BackupAsPDF(hminfo, prefix  = Prefix, subDir = subDir, dimensions = dims )
  }
  return(hminfo)
}


plot_subnetworkEnrichments <- function(enrichments, genesets = NULL, topn = 2, save = F, name = "", subd = "", dims = NULL, ...){
  
  if (is.null(genesets)){
    genesets <- names(enrichments)
  }
  
  lapply(genesets, FUN = function(en){
    tryCatch( enrichmentNicelyFormatted(enrichments[[en]], topn = topn, name, en, save = save, subDir = subd, dims = dims, ...), error=function(e) NULL)
  })
  
}


sharedGenes_Helper <- function(genesets){
  
  sharedset <- genesets[1][[1]][ genesets[1][[1]] %in% genesets[2][[1]]]
  
  if (length(genesets) > 2){
  sharedGenes_Helper( c(list(sharedset), genesets[3:length(genesets)]) )
  } else {
    return(sharedset)
  }
  
}


sharedGenes_BetwAnnotations <- function(en, group, annotations, showall = F){
  
  genesets <- lapply(annotations, FUN = function(ant){
   return( strsplit( en[ toupper(ID) == toupper(ant) & Group == group, geneID],  split = "/")[[1]] )
  })
  
  if (showall){
    return(unique(unlist(genesets)))
  } else {
    return(sharedGenes_Helper(genesets))
  }
}
