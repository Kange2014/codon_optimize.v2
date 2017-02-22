#' Linear regression using lm function in R
#'
#' @param trainingdata: This is a dataframe for model building
#' @param testingdata: This is also a dataframe for model testing
#'
#' @return A list of all variables: variables used in the regress function as an object of class regress
#'

regress <- function(trainingdata, testingdata){
   
  lm.fit <- lm(Expression ~., data=trainingdata[,2:ncol(trainingdata)])
  model <- step(lm.fit,trace = 0)
  
  #model <- lm(Expression ~., data=trainingdata[,2:ncol(trainingdata)]) %>%
   #   step(k = 2, direction = "both",trace = 0)
  
  model$trainingdata <- trainingdata
  model$testingdata <- testingdata
  
  model$training_pred <- data.frame(trainingdata,"Prediction" = predict(model,trainingdata[,2:ncol(trainingdata)]))
  
  as.list(environment()) %>% add_class(c("regress","model"))
}

# combinations
f <- function(x,y) paste(x,y,sep="")

#' dataNumeric: to convert string (DNA sequence) into numeric vectors by defining features
#'
#' @param string: original coding sequence
#' @param replaceString: coding sequence part to be optimized
#'
#' @return A dataframe containing all codon permutations of replaceString and their numeric vectors
#'
dataNumeric <- function(string = "", replaceString = ""){

  if(is_empty(string)) return("Original sequence needed. Please input the full-length coding sequence.")
  if(!is_empty(replaceString)){
	len = nchar(replaceString)
	substr(string,1,nchar(replaceString)) <- replaceString
	
	dna_string <- substr(string,1,len)
	fastaSeq <- s2c(dna_string)
	syn_codons <- syncodons(splitseq(fastaSeq))
	combn <- ""	
	for(j in 1:length(syn_codons)){
		combn <- outer(combn,syn_codons[[j]],f)
		combn <- as.vector(combn)
	}

	all_seqs <- toupper(paste(combn,substr(string,len+1,nchar(string)),sep=""))
	deltaG <- folding_cal(substr(all_seqs,1,48))

	data(caitab)
  
	fulldata <- c()
	for(i in 1:length(all_seqs)){
		a <- s2c(tolower(all_seqs[i]))
		freq1 <- seqinr::count(a,wordsize=1,freq=TRUE)
		freq2 <- seqinr::count(a,wordsize=2,freq=TRUE)
		freq3 <- seqinr::count(a,wordsize=3,freq=TRUE)
		cai_index <- cai(a,w=caitab$ec)
		fulldata<- rbind(fulldata,c(freq1,freq2,freq3,cai_index,deltaG[i,2]))
	}
	colnames(fulldata)[(ncol(fulldata)-1):ncol(fulldata)] = c("Codon Adaptation Index","Minimum free energy")
	fulldata <-data.frame(substring(all_seqs,1,len),fulldata)
	colnames(fulldata)[1] = "Sequence"
	
	fulldata
  }
  else{
	return("Please specify coding sequence needing permutation.")
  }
  
}

#' getData: to get numeric data from upload sequence dataset and split them into training data and testing data
#'
#' @param dataset: upload data file
#' @param string: original coding sequence
#' @param met: whether the uoload data file contains the initial codon Met
#'
#' @return A list of three dataframes containing numeric vectors, corresponding to full data, training data and testing data
#'

getData <- function(dataset, string = "", met = FALSE
                    ) {	
  
  if(is_empty(string)) return("Original sequence needed. Please input the full-length coding sequence.")
  
  num_col = ncol(dataset)
  
  ### coding sequence length of investigating amino acids
  ## if(met) { len = 3*(num_col-2) }
  ## else  {len = 3*(num_col-2+1) }
  
  known_seqs <- c()
  for(i in 1:nrow(dataset)){
	codons <- paste(dataset[i,2:(num_col-1)],collapse="")
	if(!met) { codons <- paste0("ATG",codons) }
	known_seqs <- c(known_seqs,codons)
  }

  fulldata <- dataNumeric(string,known_seqs[1])
  
  knowndata <- data.frame()
  for(i in 1:length(known_seqs)){
	loc = which(fulldata[,1] == toupper(known_seqs[i]))
	knowndata = rbind(knowndata,fulldata[loc,])
  }

  knowndata <- data.frame(knowndata,Expression = dataset[,num_col])
  newdata = fulldata[!(fulldata[,1] %in% known_seqs),]
  
  return(list("fulldata" = fulldata,"trainingdata" = knowndata,"testingdata" = newdata))
}

#' Summary method for the regress function
#'

summary.regress <- function(object,
                            ...) {

  if (is.character(object)) return(object)
  if (class(object$model)[1] != 'lm') return(object)

  cat("Linear regression (LM)\n")
  #cat("Data     :", object$model$model, "\n")
  
  print(summary(object$model))
  
  ##predict.regress(object,object$model$trainingdata,dec = 3)
}

#' Plot method for the regress function
#' Plot the actual values against predicted values by linear regression model
#'
plot.regress <- function(object,...){
  if (is.character(object)) return(object)
  if (class(object$model)[1] != 'lm') return(object)
  
  predictdata <- data.frame("LM" = object$model$training_pred$Prediction,"Actual"=object$model$trainingdata$Expression)
  p <- ggplot2::ggplot(predictdata,aes(LM,Actual))+ geom_point(shape=21,colour="blue")+geom_smooth(aes(LM,Actual,colour="red"),method="lm")
  
  sshhr(p)
}

#' Predict method for the regress function
#'
#'
#' @param object: Return value from regress
#' @param pred_data: a sequence numeric dataframe to used for prediction 
#' @param dec: Number of decimals to show
#' @param ...: further arguments passed to or from other methods
#'
predict.regress <- function(object,
                            pred_data = "",
                            dec = 3,
                            ...) {

 if (is.character(object)) return(object)

 pfun <- function(model, pred) {

    pred_val <-
      try(sshhr(
        predict(model, pred)
      ))

    if (!is(pred_val, 'try-error')) {
        pred_val %<>% data.frame %>% select(1)
        colnames(pred_val) <- "Prediction"
    }
    pred_val
  }

  if(is_empty(pred_data)){
	predict_model(object, pfun, "regress.predict", object$model$testingdata[,1],object$model$testingdata[,2:ncol(object$model$testingdata)], dec)
  }
  else{
	predict_model(object, pfun, "regress.predict", pred_data[,1],pred_data[,2:ncol(pred_data)], dec)
  }
  
}

#' Predict method for model functions
#'
#' @param object: Return value from regress
#' @param pfun: Function to use for prediction
#' @param mclass: Model class to attach
#' @param pred_data: a sequence numeric dataframe to used for prediction
#
#' @param dec: Number of decimals to show
#' @param ...: further arguments passed to or from other methods
#'
predict_model <- function(object, pfun, mclass,seqs = c(),
                          pred_data = "",
                          dec = 3,
                          ...) {

  if (is.character(object)) return(object)

  pred_val <- pfun(object$model,pred_data)

  
  pred <- data.frame("Seqs" = seqs,pred_data, pred_val, check.names = FALSE)
  pred <- pred[order(pred[,ncol(pred)],decreasing = TRUE),]
  vars <- colnames(pred)

  ## adding attributes used by other methods
  pred <- set_attr(pred, "dataset", "Model") %>%
	  set_attr("vars", vars) %>%
      set_attr("dec", dec) %>%
      set_attr("pred_data", pred_data)

  return(add_class(pred, c(mclass, "model.predict")))
}

#' Print method for the model prediction
#'
#' @param x: Return value from prediction method
#' @param ...: further arguments passed to or from other methods
#' @param n: Number of lines of prediction results to print. Use -1 to print all lines
#' @param header Header line
#'
print_predict_model <- function(x, ..., n = 10, header = "") {

  class(x) <- "data.frame"
  vars <- attr(x, "vars")
  vars <- vars %in% c("Seqs", "g", "c", "Codon.Adaptation.Index", "Minimum.free.energy","Prediction")
  pred_data <- attr(x, "pred_data")

  cat(header)
  cat("\n")
  
  if (!is.character(pred_data)) pred_data <- "-----"
  cat("Prediction dataset   :", pred_data,"\n")

  if (n == -1) {
    cat("\n")
    formatdf(x[,vars, drop = FALSE], attr(x, "dec")) %>% print(row.names = FALSE)
  } else {
    if (nrow(x) > n)
      cat("Top shown            :", n, "\n")
    cat("\n")
    head(x[,vars, drop = FALSE], n) %>% formatdf(attr(x, "dec")) %>% print(row.names = FALSE)
  }
}

#' Print method for predict.regress
#'
#' @param x: Return value from prediction method
#' @param ...: further arguments passed to or from other methods
#' @param n: Number of lines of prediction results to print. Use -1 to print all lines
#'
#' @export
print.regress.predict <- function(x, ..., n = 10)
  print_predict_model(x, ..., n = n, header = "Linear regression (LS)")
