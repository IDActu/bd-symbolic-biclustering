crossvalidation <- function(folds = 10, matrixPath = "matrix.csv", prefixTrain = "matrix_TRAIN", prefixTest = "matrixTest")
{
  train_ration <- sqrt(0.7)
  samples <- list(list())
  myNames <- vector()
  
  data <- read.csv(matrixPath, row.names = 1)
    
  for(iCYCLES in 1:folds)
  {
    #sample rows
    row_sample_TRAIN <- sample(1:nrow(data), size = round(nrow(data)*train_ration), replace = FALSE)
    row_sample_TRAIN <- sort(row_sample_TRAIN, decreasing = FALSE)
    row_sample_TRAIN_Templ <- rep(FALSE, nrow(data))
    row_sample_TRAIN_Templ[row_sample_TRAIN] <- TRUE
    
    #sample columns
    col_sample_TRAIN <- sample(1:ncol(data), size = round((ncol(data))*train_ration), replace = FALSE) #minus flybaseID
    col_sample_TRAIN <- sort(col_sample_TRAIN, decreasing = FALSE)
    col_sample_TRAIN_Templ <- rep(FALSE, ncol(data))
    col_sample_TRAIN_Templ[col_sample_TRAIN] <- TRUE
    
    template <- row_sample_TRAIN_Templ %*% t(col_sample_TRAIN_Templ)
    #train index
    index <- which(template == 1)
    index_row <- index %% nrow(data)
    index_row[index_row == 0] <- nrow(data) #row index
    index_col <- ((index - 1) %/% nrow(data))+1 #column index
    
    #test index
    indexTEST <- which(template == 0)
    index_row_TEST <- indexTEST %% nrow(data)
    index_row_TEST[index_row_TEST == 0] <- nrow(data) #row index
    index_col_TEST <- ((indexTEST - 1) %/% nrow(data))+1 #column index
    
    #train dataset
    TRAIN_DATASET <- data[row_sample_TRAIN, (col_sample_TRAIN)]
    
    #test dataset
    TEST_DATASET <- data
    TEST_DATASET[,] <- NaN
    for(i in 1:length(index_row_TEST))
    {
      TEST_DATASET[index_row_TEST[i], index_col_TEST[i]] <- data[index_row_TEST[i], index_col_TEST[i]]
    }
    
    write.csv(TRAIN_DATASET, file = paste0(prefixTrain,iCYCLES,'.csv', sep = ''))
    write.csv(TEST_DATASET, file = paste0(prefixTest,iCYCLES,'.csv', sep = ''))
    
    #save indexes
    samples[iCYCLES] <- list(list(index_col = index_col, index_row = index_row, index_col_TEST = index_col_TEST, index_row_TEST  = index_row_TEST))
    myNames[length(myNames)+1] <- paste0('run', iCYCLES, sep = '')
  }
  
  names(samples) <- myNames
  #saveRDS(samples, "samples.rds")
  
  
  #create 3 testing sets
  samples_3tests <- samples
  for(irun in 1:length(samples))
  {
    #test columns
    col <- samples[[irun]]$index_col_TEST
    row <- samples[[irun]]$index_row_TEST
    
    #train columns
    train_col <- samples[[irun]]$index_col
    train_row <- samples[[irun]]$index_row
    
    onlyTestCOL <- setdiff(col, train_col) #not in train columns
    onlyTestROW <- setdiff(row, train_row) #not in train rows
    
    testUnique_COL <- sapply(col, function(x){if(x %in% onlyTestCOL){x <- TRUE} else {x <- FALSE} })
    testUnique_ROW <- sapply(row, function(x){if(x %in% onlyTestROW){x <- TRUE} else {x <- FALSE} })
    
    testUniqueCOL_INDEX <- xor(testUnique_COL, testUnique_ROW) & testUnique_COL
    testUniqueROW_INDEX <- xor(testUnique_COL, testUnique_ROW) & testUnique_ROW
    testCOMMON_INDEX <- testUnique_COL & testUnique_ROW
    
    testUniqueLOC_COL <- col[testUniqueCOL_INDEX]
    testUniqueLOC_ROW <- row[testUniqueCOL_INDEX]
    
    testUniqueGENE_COL <- col[testUniqueROW_INDEX]
    testUniqueGENE_ROW <- row[testUniqueROW_INDEX]
    
    testCOMMON_COL <- col[testCOMMON_INDEX]
    testCOMMON_ROW <- row[testCOMMON_INDEX]
    
    #create & save object
    samples_3tests[[irun]][["testUniqueLOC_COL_INDEX"]] <- testUniqueLOC_COL
    samples_3tests[[irun]][["testUniqueLOC_ROW_INDEX"]] <- testUniqueLOC_ROW
    
    samples_3tests[[irun]][["testUniqueGENE_COL_INDEX"]] <- testUniqueGENE_COL
    samples_3tests[[irun]][["testUniqueGENE_ROW_INDEX"]] <- testUniqueGENE_ROW
    
    samples_3tests[[irun]][["testCOMMON_COL_INDEX"]] <- testCOMMON_COL
    samples_3tests[[irun]][["testCOMMON_ROW_INDEX"]] <- testCOMMON_ROW
    
    #print final datasets
    TEST_DATASET_COL <- data
    TEST_DATASET_COL[,] <- NaN
    
    TEST_DATASET_ROW <- data
    TEST_DATASET_ROW[,] <- NaN
    
    TEST_DATASET_COMMON <- data
    TEST_DATASET_COMMON[,] <- NaN
    
    for(i in 1:length(samples_3tests[[irun]][["testUniqueLOC_ROW_INDEX"]]))
    {
      TEST_DATASET_COL[samples_3tests[[irun]][["testUniqueLOC_ROW_INDEX"]][i], samples_3tests[[irun]][["testUniqueLOC_COL_INDEX"]][i]] <- data[samples_3tests[[irun]][["testUniqueLOC_ROW_INDEX"]][i], samples_3tests[[irun]][["testUniqueLOC_COL_INDEX"]][i]]
    }
    
    for(i in 1:length(samples_3tests[[irun]][["testUniqueGENE_ROW_INDEX"]]))
    {
      TEST_DATASET_ROW[samples_3tests[[irun]][["testUniqueGENE_ROW_INDEX"]][i], samples_3tests[[irun]][["testUniqueGENE_COL_INDEX"]][i]] <- data[samples_3tests[[irun]][["testUniqueGENE_ROW_INDEX"]][i], samples_3tests[[irun]][["testUniqueGENE_COL_INDEX"]][i]]
    }
    
    for(i in 1:length(samples_3tests[[irun]][["testCOMMON_ROW_INDEX"]]))
    {
      TEST_DATASET_COMMON[samples_3tests[[irun]][["testCOMMON_ROW_INDEX"]][i], samples_3tests[[irun]][["testCOMMON_COL_INDEX"]][i]] <- data[samples_3tests[[irun]][["testCOMMON_ROW_INDEX"]][i], samples_3tests[[irun]][["testCOMMON_COL_INDEX"]][i]]
    }
    write.csv(TEST_DATASET_COL, file = paste0(prefixTest,irun,'_keepGenes','.csv', sep = ''))
    write.csv(TEST_DATASET_ROW, file = paste0(prefixTest,irun,'_keepLocations','.csv', sep = ''))
    write.csv(TEST_DATASET_COMMON, file = paste0(prefixTest,irun,'_bd','.csv', sep = ''))
  }
  
  #saveRDS(samples_3tests, file = "samples_3tests.RDS")
  return(samples_3tests)
}
