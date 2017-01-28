## biclustering functions
library(topGO)

## convert into FIMI sparse format -- column indices for 1s
## A helper function that tests whether an object is either NULL _or_ a list of NULLs
## in here serves to remove empty rows
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
## Recursively step down into list, removing all such objects 
rmNullObs <- function(x) {
   x <- Filter(Negate(is.NullOb), x)
   lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

## generate a single input file for PANDA
genPandaInput<-function(TRAIN_FILE,INDEX){
  d <- read.csv(paste(TRAIN_FILE,INDEX,".csv",sep=""), row.names = 1)
  names(d) <- seq(0,length(d)-1)
  m <- as.matrix(d)
  mm <- apply(m,1,function(x){if (sum(x)) colnames(m)[x==1]})
  fpanda <- paste("fimi_mat_",INDEX,".csv",sep="")
  fpanda_names <- paste("fimi_mat_",INDEX,"_names",".csv",sep="")
  #file.remove(fpanda)
  mm_notnull<-rmNullObs(mm)
  lapply(mm_notnull, write, fpanda, append=TRUE, ncolumns=1000)
  write.table(names(mm_notnull),fpanda_names,row.names=FALSE,col.names=FALSE,quote=FALSE)
  return(fpanda)
}  

## generate input files for PANDA
## work with dedicated column names in OVARY dataset ... obsolete
genPandaInputOvary<-function(fnames){
  for (f in fnames) {
    d <- read.csv(f)
    l <- strsplit(names(d),"[.]")
    names(d) <- sapply(l,FUN = "[[",2)
    m <- as.matrix(d)
    mm <- apply(m,1,function(x){if (sum(x)) colnames(m)[x==1]})
    fpanda <- paste("fimi_mat_",strsplit(f,"[_.]")[[1]][3],".csv",sep="")
    fpanda_names <- paste("fimi_mat_",strsplit(f,"[_.]")[[1]][3],"_names",".csv",sep="")
    file.remove(fpanda)
    mm_notnull<-rmNullObs(mm)
    lapply(mm_notnull, write, fpanda, append=TRUE, ncolumns=1000)
    write.table(names(mm_notnull),fpanda_names,row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  return(fpanda)
}  
## generate input files for PANDA for discs, treat names differently than in ovaries ... obsolete
genPandaInputDiscs<-function(fnames){
  for (f in fnames) {
    d <- read.csv(f)
    l <- strsplit(names(d),"[.]")
    names(d) <- sapply(l,function(x){paste(substring(x[1],2),x[3],sep="")})
    m <- as.matrix(d)
    mm <- apply(m,1,function(x){if (sum(x)) colnames(m)[x==1]})
    fpanda <- paste("../franta/disky/fimi_mat_",strsplit(f,"[_.]")[[1]][5],".csv",sep="")
    fpanda_names <- paste("../franta/disky/fimi_mat_",strsplit(f,"[_.]")[[1]][5],"_names",".csv",sep="")
    file.remove(fpanda)
    mm_notnull<-rmNullObs(mm)
    lapply(mm_notnull, write, fpanda, append=TRUE, ncolumns=1000)
    write.table(names(mm_notnull),fpanda_names,row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
}  

## generate input files for PANDA for discs, treat names differently than in ovaries
## generate 4 separate files
genPandaInputDiscs4<-function(fnames){
  for (f in fnames) {
    d <- read.csv(f)
    l <- strsplit(names(d),"[.]")
    gind <- factor(unlist(lapply(l,function(x) x[1]))) # split by disc name
    gsizes <- sapply(split(gind,gind),length) # find indeces to address the individual discs
    gsizes <- c(0,cumsum(gsizes))
    names(d) <- sapply(l,function(x){paste(substring(x[1],2),x[3],sep="")})
    m <- as.matrix(d)
    for (disc in seq(2,length(gsizes))){
      mm <- apply(m[,(gsizes[disc-1]+1):gsizes[disc]],1,function(x){if (sum(x)) colnames(m[,(gsizes[disc-1]+1):gsizes[disc]])[x==1]})
      fpanda <- paste("../franta/disky/fimi_mat_",strsplit(f,"[_.]")[[1]][5],"_X",(disc-1),".csv",sep="")
      fpanda_names <- paste("../franta/disky/fimi_mat_",strsplit(f,"[_.]")[[1]][5],"_X",(disc-1),"_names",".csv",sep="")
      file.remove(fpanda)
      mm_notnull<-rmNullObs(mm)
      lapply(mm_notnull, write, fpanda, append=TRUE, ncolumns=1000)
      write.table(names(mm_notnull),fpanda_names,row.names=FALSE,col.names=FALSE,quote=FALSE)
    }  
  }
}  


## read and deparse Panda file
readPanda<-function(f){
  fc <- file(f)
  mylist <- strsplit(readLines(fc), split=" [(]|[)] [[]|[]]")
  close(fc)
  mylist<-lapply(mylist,function(x){strsplit(x," ")})
  mylist<-lapply(mylist,function(x){lapply(x,function(x) sort(as.numeric(x)))})
  return(mylist)
}

# get annotations to GO mapping (use FlyBase Ids, 4113 out of 4490 available in drosophila2.db (92%), 3995 available in drosgenome1.db, 4125 when merged) 
# generates a named list of GO Ids (a list of Ids for every FlyBase Id in flyids)
annotData<-function(lib,flyids){
    library(lib,character.only=T)
    libshort<-unlist(strsplit(lib,"[.]"))[1]
    x <- get(paste(libshort,"FLYBASE",sep="")) # access drosophila2FLYBASE object to translate probes into FlyBase
    mapped_genes <- mappedkeys(x)
    xx <- unlist(as.list(x[mapped_genes]))
    probe2fly <- xx[xx %in% flyids]
    x <- get(paste(libshort,"GO",sep="")) # access drosophila2GO object to obtain GO annotations
    mapped_genes <- mappedkeys(x)
    yy <- as.list(x[mapped_genes])
    gene2GO <- lapply(yy,function(y) unique(unlist(lapply(y,function(x)x[[1]])))) # simplify the output, use GO ids only
    gene2GO <- gene2GO[names(gene2GO) %in% names(probe2fly)] # remove translations that do not relate to flyids
    names(gene2GO) <- probe2fly[match(names(gene2GO),names(probe2fly))]
    return(gene2GO)
}

# multiple calls to all biclusters
# first define function that evaluates single bc
runBcGOTest<-function(bc, fullflynames){
    bcList<-rep(0,length(fullflynames))
    bcList[fullflynames %in% bc[[2]]]<-1
    bcList<-factor(bcList)
    names(bcList)<-fullflynames
    bcBP<-new("topGOdata",ontology = "BP",allGenes = bcList, nodeSize = 20, annot = annFUN.gene2GO, gene2GO=fly2GO)
    bcCC<-new("topGOdata",ontology = "CC",allGenes = bcList, nodeSize = 30, annot = annFUN.gene2GO, gene2GO=fly2GO)
    bcMF<-new("topGOdata",ontology = "MF",allGenes = bcList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO=fly2GO)
    resFisherBP <- runTest(bcBP, algorithm = "classic", statistic = "fisher")
    resFisherCC <- runTest(bcCC, algorithm = "classic", statistic = "fisher")
    resFisherMF <- runTest(bcMF, algorithm = "classic", statistic = "fisher")
    return(list(BP=GenTable(bcBP,classic=resFisherBP,topNodes=100),CC=GenTable(bcCC,classic=resFisherCC,topNodes=30),MF=GenTable(bcMF,classic=resFisherMF,topNodes=30)))
}

## learn the bicluster description 
# take all terms that describe the bicluster at the given p-value -- compute their scores (the score is log pval)
# apply for all the results, merge BP, MF and CC
descBiclusters<-function(sumRes,pval){
  ibc<-seq(1,length(sumRes),3)
  sumRes<-lapply(seq(length(sumRes)),function(x) {sumRes[[x]]$classic<-sapply(sumRes[[x]]$classic,function(y) {gsub("< ","",y)});return(sumRes[[x]])}) # remove < symbol from pvals
  GOList<-sapply(ibc,function(x){sapply(seq(0,2),function(y){takeRows<-sumRes[[x+y]][as.numeric(sumRes[[x+y]]$classic)<pval,];return(list(takeRows$GO.ID,-log(as.numeric(takeRows$classic),10)))})})
  GOList<-apply(GOList,2,function(x){lapply(seq(1,5,2),function(y){lapply(seq(length(x[[y]])),function(z){c(x[[y]][z],x[[y+1]][z])})})})
  names(GOList)<-paste("bc",seq(1,length(GOList)),sep="")
  for (x in 1:length(GOList)) {names(GOList[[x]])<-c("BP","CC","MF")}
  #GOList<-unlist(unlist(GOList,recursive=F),recursive=F)
  #GOList<-GOList[sapply(GOList,function(x){is.na(x[1])})==0]
  return(GOList)
}

keggBiclusters<-function(desc,sumResKEGG,pvalGenes){
# extend desc with KEGG part
  for (i in seq(length(desc))){
    tKEGG<--log(sumResKEGG[sumResKEGG[,i]<pvalGenes,i],10)
    tKEGG<-lapply(seq(length(tKEGG)),function(x) c(names(tKEGG)[x],as.character(tKEGG[x])))
    desc[[i]]$KEGG<-tKEGG
  }
  return(desc)
}

## decompose a list of terms with their scores
decTList<-function(tSubList){
  terms<-unlist(lapply(tSubList,function(x) x[1]))
  scores<-as.numeric(unlist(lapply(tSubList,function(x) x[2])))
  return(list(terms=terms,scores=scores))
}

## aggregates scores reached by the individual genes
# modifies the existing vector of gScores, adds the scores for the genes contained in gList, the scores are available in tSubList
aggScore<-function(gScores,gTList,tSubList){
  if (length(gTList)>0){
    for (x in seq(length(gTList))){
      relGenes<-names(gScores) %in% gTList[[x]]
      tScore<-tSubList$scores[which(tSubList$terms==names(gTList)[x])]
      gScores[relGenes]<-gScores[relGenes]+tScore
    }
  }  
  return(gScores)
}

## apply the bicluster description to test data
# obtain all the genes from gList covered by the given (structured) list of terms
# covered means that their aggregated score exceeds the threshold
applyBicluster<-function(gList,tList,thres){
  # change the structure of the term list 
  tListRef<-lapply(tList,decTList)
  # initialize gScores
  gScores<-rep(0,length(gList))
  names(gScores)<-gList
  # genesInTerm is a topGO function that finds all the genes assigned to the list of terms in the given ontology
  # bcBP, bcMF and bcCC must exist as global objects
  gScores<-aggScore(gScores,genesInTerm(bcBP,tListRef$BP$terms),tListRef$BP)
  gScores<-aggScore(gScores,genesInTerm(bcCC,tListRef$CC$terms),tListRef$CC)
  gScores<-aggScore(gScores,genesInTerm(bcMF,tListRef$MF$terms),tListRef$MF)
  if (length(tList)>3) {gScores<-aggScore(gScores,keggOntInv[names(keggOntInv) %in% tListRef$KEGG$terms],tListRef$KEGG)}
  return(names(gScores)[gScores>thres])
}

# the individual terms to test defined in a data frame
runBcMyOntTest<-function(bc,myOntDF){
    apply(myOntDF,2,function(x){locs<-rep(0,nrow(myOntDF));locs[which(names(d) %in% bc[[1]])]<-1;t<-table(locs,x);fisher.test(t,alternative="greater")$p.value})
}

runFTest<-function(myTerm,locs,myOnt){
# runs Fisher test for the given location term and the given bicluster location indices
    fft<-table(data.frame(biclust=names(myOnt) %in% locs,termOcc=sapply(myOnt, function(x) myTerm %in% x)))
    if (min(dim(fft))<2) {return(1)}
    else return(fisher.test(fft,alternative="greater")$p.value)
}

runBcLocTest<-function(bc){
# gets all the enriched location terms for the given set of locations (one element of bicluster)
    locs <- locInd$long[locInd$long %in% bc[[1]]]
    resFisherLoc <- sapply(locTerms,function(x) runFTest(x,locs,locOnt))
    return(resFisherLoc)
}

# the ontology from Franta
readLocOnt<-function(fname){
    locOnt <- scan(fname, what="", sep="\n")
    locOnt <- strsplit(locOnt, "[[:space:]]+")
    # Extract the first vector element and set it as the list element name
    names(locOnt) <- sapply(locOnt, `[[`, 1)
    # Remove the first vector element from each list element
    locOnt <- lapply(locOnt, `[`, -1)
}    

runBcKEGGTest<-function(bc){
# gets all the KEGG pathways for the given set of genes (one element of bicluster)
    genes <- names(keggOnt)[names(keggOnt) %in% bc[[2]]]
    resFisherLoc <- sapply(keggTerms,function(x) runFTest(x,genes,keggOnt))
    return(resFisherLoc)
}

## aggregates scores reached by the individual genes
# modifies the existing vector of gScores, adds the scores for the genes contained in gList, the scores are available in tSubList
aggScoreOnt<-function(locLists,tSubList,thres,myOnt){
  # initialize lScores
  lScores<-rep(0,length(names(myOnt)))
  names(lScores)<-names(myOnt)
  if (length(tSubList)>0){
    for (x in seq(length(tSubList))){
      relLocs<-locLists[[which(names(locLists)==names(tSubList)[x])]]
      lScores[relLocs]<-lScores[relLocs]+tSubList[x]
    }
  }  
  return(names(lScores)[lScores>thres])
}

# obtain all the locations that correspond to the given term at the given pval
descBiclusterOnt<-function(sumResOnt,pval,thres){
  # consruct location lists for the individual terms
  locLists<-sapply(rownames(sumResOnt),function(x) c(x,names(locOnt)[sapply(locOnt, function(y) x %in% y)]))
  # for each bicluster get the terms exceeding the pval threshold
  bcdesc<-apply(sumResOnt,2,function(x) sapply(x[x<pval],function(y) -log(y,10)))
  # only one bidluster ...
  if (typeof(bcdesc)!="list") {bcdesc<-list(bcdesc[,1])}
  # collect the locations that exceed the aggregate score threshold
  lScoresList<-lapply(bcdesc,function(x) aggScoreOnt(locLists,x,thres,locOnt))
  return(lScoresList)
}

# obtain all the locations that correspond to the given term at the given pval
descBiclusterKegg<-function(sumResKegg,pval,thres){
  # consruct location lists for the individual terms
  locLists<-sapply(rownames(sumResKegg),function(x) c(x,names(keggOnt)[sapply(keggOnt, function(y) x %in% y)]))
  # for each bicluster get the terms exceeding the pval threshold
  bcdesc<-apply(sumResKegg,2,function(x) sapply(x[x<pval],function(y) -log(y,10)))
  # only one bidluster ...
  if (typeof(bcdesc)!="list") {bcdesc<-list(bcdesc[,1])}
  # collect the locations that exceed the aggregate score threshold
  lScoresList<-lapply(bcdesc,function(x) aggScoreOnt(locLists,x,thres,keggOnt))
  return(lScoresList)
}


# pounit AUC estimate, only a single pair of TPr and FPr values
myAUCpoint<-function(TPr,FPr){
  return(FPr*TPr/2+TPr*(1-FPr)+(1-FPr)*(1-TPr)/2)
}

# pounit AUC estimate, only a single pair of TPr and FPr values
myAUC4f<-function(ft){
  TN<-ft[1,1];TP<-ft[2,2];FN<-ft[1,2];FP<-ft[2,1];
  return(myAUCpoint(TP/(TP+FN),FP/(FP+TN)))
}

# uses the fact that the AUC is equal to the probability that a true positive is scored greater than a true negative
# pos.scores is a vector containing a score of the positive examples, and neg.scores is a vector containing the negative examples
myAUCstoch<-function(pscores,nscores){
  prand<-sample(pscores[!is.na(pscores)],1000,replace=T)
  nrand<-sample(nscores[!is.na(nscores)],1000,replace=T)
  return(mean(prand > nrand)+mean(prand==nrand)/2)
  # for bootstrap
  #aucs = replicate(1000,mean(sample(pos.scores,1000,replace=T) > sample(neg.scores,1000,replace=T)))
}