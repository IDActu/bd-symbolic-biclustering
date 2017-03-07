####################### Semantic biclustering #######################

#### Parameters
MEMGB <- 2  #memory allocated for JVM
NUM_FOLDS <- 10 #number of folds
MIN_COVER <- 500 #min term coverage on GO ontology

GENE_FILE <- "00-data/pathwayList"  #path to the row ontology file (kegg)
LOCATION_FILE <- "00-data/columnOntology" #path to the column ontology file
MATRIX <- "00-data/discMatrix.csv"  #original gene expression matrix
GO_OBO <- "00-data/go-basic.obo"  #GO in obo format
FBGN2GO <- "00-data/gene_association.fb"  #mapping FBgn to GO ID

OUTPUT_DIR_CV <- "01-crossvalidation/"
OUTPUT_DIR_ARFF <- "02-arff_files/"
OUTPUT_PREFIX_TRAIN_CV <- paste0(OUTPUT_DIR_CV, "matrixTRAIN")
OUTPUT_PREFIX_TEST_CV <- paste0(OUTPUT_DIR_CV, "matrixTEST")
OUTPUT_PREFIX_TRAIN_ARFF <- paste0(OUTPUT_DIR_ARFF, "train")
OUTPUT_PREFIX_TEST_ARFF <- paste0(OUTPUT_DIR_ARFF, "test")


#### Crossvalidation step
source("scripts/crossvalidation.R")
#make a directory for crossvalidation
system(paste("mkdir -p",OUTPUT_DIR_CV))
#build n training and testing datasets, safe them and return indexes
splitIndexes <- crossvalidation(folds = NUM_FOLDS, matrixPath = MATRIX,
                                prefixTrain = OUTPUT_PREFIX_TRAIN_CV, prefixTest = OUTPUT_PREFIX_TEST_CV)

#make a directory for results
system(paste("mkdir -p",OUTPUT_DIR_ARFF))
#### RULE&TREE LEARNING
for(ifold in 1:NUM_FOLDS)
{
  #### Create ontology files
  
  #### Build ARFF file
  #command for creating an ARFF file
  command <- paste("perl scripts/buildARFFdataset.pl --pathTrain", paste0(OUTPUT_PREFIX_TRAIN_CV, ifold, ".csv"),
                   "--pathTest", paste0(OUTPUT_PREFIX_TEST_CV,ifold, ".csv"), "--min_cover", MIN_COVER,
                   "--outputTrain", paste0(OUTPUT_PREFIX_TRAIN_ARFF, ifold, ".arff"),
                   "--outputTest",  paste0(OUTPUT_PREFIX_TEST_ARFF, ifold, ".arff"),
                   "--gene", GENE_FILE, "--location", LOCATION_FILE, "--rowOBO", GO_OBO,
                   "--fbgn2GO", FBGN2GO)
  
  #execute the command
  system(command)
  
  #memory restriction
  options(java.parameters = paste0("-Xmx",MEMGB,"g"))
  library(RWeka)
  
  #load train and test datasets
  trainDataset <- read.arff(paste0(OUTPUT_PREFIX_TRAIN_ARFF, ifold, ".arff"))
  testDataset <- read.arff(paste0(OUTPUT_PREFIX_TEST_ARFF, ifold, ".arff"))
  
  #J48
  trainJ48 <- J48(classification ~ ., data = trainDataset, control = Weka_control(C = 0.01, M = 20))
  testJ48 <- evaluate_Weka_classifier(trainJ48, newdata = testDataset, 
                                      class = TRUE, complexity = TRUE)
  write(testJ48$string, stdout())
  
  #JRip
  trainJRip <- JRip(classification ~ ., data = trainDataset,
                    control = Weka_control(F = 3, N = 2.0, O = 2, S = 1))
  testJRip <- evaluate_Weka_classifier(trainJRip, newdata = testDataset,
                                       class = TRUE, complexity = TRUE)
  write(testJRip$string, stdout())
}

#### BI-DIRECTIONAL ENRICHMENT
source("scripts/bdenrichment.r")
resTLall<-list()
AUCLall<-list()
precLall<-list()
sizeLall<-list()
for(ifold in 1:NUM_FOLDS)
{
  res <-  bdenrichment(index = ifold, trainFile = OUTPUT_PREFIX_TRAIN_CV, testFile = OUTPUT_PREFIX_TEST_CV,
               splitIndexes = splitIndexes, dataset = MATRIX, locOntologyFile = LOCATION_FILE,
               geneOntology = GENE_FILE, pandaPath = "scripts/panda",
                bdParams = c(pvalGenes=0.05,pvalLocs=0.1,thresScoreGenes=1,thresScoreLocations=1))
  
  resTLall[[length(resTLall)+1]] <- res[["resTL"]]
  AUCLall[[length(AUCLall)+1]] <- res[["AUCL"]]
  precLall[[length(precLall)+1]] <- res[["precL"]]
  sizeLall[[length(sizeLall)+1]] <- res[["sizeL"]]
  
  write(res[["AUCL"]][1], stdout())
}
