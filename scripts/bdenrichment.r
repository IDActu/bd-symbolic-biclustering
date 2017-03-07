GENE_ANNOTATION_LIB<-"drosophila2.db" # R library with gene annotations

#index - index of the split, the splits were predefined
#bdParams- (pvalGenes=0.05,pvalLocs=0.1,thresScoreGenes=1,thresScoreLocations=1) # set the parameters of bidirectional enrichment
#dataset - gene expression binary dataset
#locOntologyFile - location ontology file
#geneOntology - mapping between the KEGG pathways and gene annotations
#splitIndexes - train test splits
#trainFile - train file prefix
#testFile - test file prefix
#pandaPath - path to panda file


bdenrichment <- function(index = 1, trainFile, testFile, splitIndexes, 
                         bdParams = c(pvalGenes=0.05,pvalLocs=0.1,thresScoreGenes=1,thresScoreLocations=1),
                         dataset, locOntologyFile, geneOntology, pandaPath = "panda")
{
  source("scripts/bc_library.r") # load supporting functions
  
  ## COMMON FOR ALL SPLITS
  
  #d <-readRDS(GE_DATASET)
  d <-read.csv(dataset, row.names = 1)
  
  fullnames = as.character(names(d))
  fullnames <- gsub("#",".",fullnames) # change the location names so that they match names read from matrix headers
  l=strsplit(names(d),"[#.]")
  names(d) = sapply(l,function(x){paste(x[1],x[3],sep="")})
  fullnames.df = data.frame(long=fullnames,short=seq(0,length(d)-1))
  
  fullflynames<-unique(rownames(d))
  fly2GO <- annotData(GENE_ANNOTATION_LIB,fullflynames)
  locOnt<-readLocOnt(locOntologyFile)
  names(locOnt)<-gsub("#",".",names(locOnt)) # change the location names so that they match names read from matrix headers
  keggOnt<-readLocOnt(geneOntology)
  keggTerms <- unique(unlist(keggOnt))
  keggOntInv<-lapply(keggTerms,function(x) names(keggOnt)[sapply(keggOnt,function(y) x %in% y)])
  names(keggOntInv)<-keggTerms
  
  # universal objects, no real bcList file
  bcList<-c(0,rep(1,length(fullflynames)-1))
  bcList<-factor(bcList)
  names(bcList)<-fullflynames
  # build topGO objects to be used in applyBicluster (general ontology graphs)
  bcBP<-new("topGOdata",ontology = "BP",allGenes = bcList, nodeSize = 10, annot = annFUN.gene2GO, gene2GO=fly2GO)
  bcMF<-new("topGOdata",ontology = "MF",allGenes = bcList, nodeSize = 10, annot = annFUN.gene2GO, gene2GO=fly2GO)
  bcCC<-new("topGOdata",ontology = "CC",allGenes = bcList, nodeSize = 10, annot = annFUN.gene2GO, gene2GO=fly2GO)
  
  # load 10 train-test splits into samples
  #samples <-readRDS(SPLIT_FILE)
  samples <-splitIndexes
  
  ## BICLUSTERING
    
  fpanda<-genPandaInput(trainFile,index)
  x<-paste(pandaPath," -d ",fpanda," -o fimi_out_",index,".csv -y 0.3 -t 0.3",sep="")
  system(x)
  
  trainRows<-unique(samples[[index]]$index_row)
  trainCols<-unique(samples[[index]]$index_col)
  testA = read.csv(paste(testFile,index,".csv",sep=""), row.names = 1)
  testRows<-!(seq(dim(testA)[1]) %in% trainRows)
  testCols<-!(seq(dim(testA)[2]) %in% trainCols)
  locInd <- fullnames.df[trainCols,]
  
  bc<-readPanda(paste("fimi_out_",index,".csv",sep=""))
  actflynames<-read.csv(paste("fimi_mat_",index,"_names.csv",sep=""),header=F) # read flynames that corerspond to the given train file
  bci<-lapply(bc,function(x){bcflynames=as.character(actflynames$V1[x[[3]]+1]);list(as.character(fullnames.df$long[fullnames.df$short %in% as.character(x[[1]])]),bcflynames[nchar(bcflynames)>0])}) # +1 as Panda counts rows from 0
  
  # call for all (non-trivial) biclusters
  ntbc<-bci[sapply(bci,function(x) length(x[[2]]))>2]
  bciSizes<-sapply(bci,function(x) paste(length(x[[1]]),length(x[[2]]),sep="x"))
  
  # analyze the gene direction
  sumRes <- sapply(ntbc,runBcGOTest,fullflynames=fullflynames,fly2GO=fly2GO)
  names(sumRes) <- paste(rep(seq(length(ntbc)),each=3),rep(c("BP","CC","MF"),length(ntbc)))
  sumResOnt<-sapply(ntbc,runBcLocTest,locOnt=locOnt,locInd=locInd,locTerms=unique(unlist(locOnt)))
  sumResKEGG<-sapply(ntbc,runBcKEGGTest,keggOnt=keggOnt,keggTerms=keggTerms)
  
  ## CLASSIFICATION
  
  # apply the biclusters to classify the test data
  # take the full dataset and replace all the TRAIN entries by NAs
  testA[(rowSums(testA,na.rm=T)==sum(testCols) | rowSums(testA,na.rm=T)==ncol(testA)),] <- NA
  fullgnames<-rownames(testA)
  #names(testA) = sapply(names(testA),function(x) substr(x,2,nchar(x)))
  
  # tests limited to specific test matrix subregions
  testG <- testA
  testG[testRows,]<-NA # keep genes, generalization in terms of locations
  testS <- testA
  testS[,testCols]<-NA # keep locations, generalization in terms of genes
  test2f <- testA
  test2f[testRows,trainCols]<-NA;test2f[trainRows,testCols]<-NA  
  
  desc<-descBiclusters(sumRes,bdParams[["pvalGenes"]])
  desc<-keggBiclusters(desc,sumResKEGG,bdParams[["pvalGenes"]]) # extend desc with the KEGG part
  iGenes<-lapply(desc,function(x){applyBicluster(fullflynames,x,bdParams[["thresScoreGenes"]],bcBP,bcMF,bcCC,keggOntInv)})
  #iGenes<-descBiclusterKegg(sumResKEGG,keggOnt,bdParams[[pvalGenes]],thresScoreGenes) # independent KEGG test
  iLocs<-descBiclusterOnt(sumResOnt,locOnt,bdParams[["pvalLocs"]],bdParams[["thresScoreLocations"]])
  # store the description size
  sizeL<-data.frame(goTerms=sapply(desc,function(x) sum(sapply(x,length))),locTerms=colSums(sumResOnt<bdParams[["pvalLocs"]]),genes=sapply(iGenes,length),locs=sapply(iLocs,length))
  pred<-testA
  pred[!is.na(testA)]<-0
  predG<-pred;predS<-pred
  # use the biclusters to predict
  for (i in seq(1,length(iGenes))){
    pred[fullgnames %in% iGenes[[i]],names(testA) %in% iLocs[[i]]] <- 1
    predG[fullgnames %in% bci[[i]][[2]],names(testA) %in% iLocs[[i]]] <- 1
    predS[fullgnames %in% iGenes[[i]],names(testA) %in% bci[[i]][[1]]] <- 1
  }
  predfac<-factor(t(pred),levels=c(0,1))
  predfacG<-factor(t(predG),levels=c(0,1));predfacS<-factor(t(predS),levels=c(0,1))
  resA<-table(predfac,as.vector(t(testA)))
  resG<-table(predfacG,as.vector(t(testG)))
  resS<-table(predfacS,as.vector(t(testS)))
  res2f<-table(predfac,as.vector(t(test2f)))
  
  prec<-resA[2,2]/(resA[2,2]+resA[2,1]) # precision
  AUCL<-c(all=myAUC4f(resA),gen2f=myAUC4f(res2f),keepgenes=myAUC4f(resG),keeplocations=myAUC4f(resS))
  resTL<-list(all=resA,gen2f=res2f,keepgenes=resG,keeplocations=resS)
  return(list(prec=prec,AUCL=AUCL,resTL=resTL,sizeL=sizeL))
}
