#' Alias used to add an attribute
#'
#' @param x Object
#' @param which Attribute name
#' @param value Value to set
#
#' @examples
#' foo <- data.frame(price = 1:5) %>% set_attr("desc", "price set in experiment ...")
#'
#' @export
set_attr <- function(x, which, value) `attr<-`(x, which, value)

#' Convenience function to add a class
#'
#' @param x Object
#' @param cl Vector of class labels to add
#'
#' @examples
#' foo <- "some text" %>% add_class("text")
#' foo <- "some text" %>% add_class(c("text","another class"))
#'
#' @export
add_class <- function(x, cl) `class<-`(x, c(cl, class(x)))


#' Hide warnings and messages and return invisible
#'
#' @details Adapted from \url{http://www.onthelambda.com/2014/09/17/fun-with-rprofile-and-customizing-r-startup/}
#'
#' @param ... Inputs to keep quite
#'
#' @examples
#' sshh( library(dplyr) )
#'
#' @export
sshh <- function(...) {
  suppressWarnings( suppressMessages( ... ) )
  invisible()
}

#' Hide warnings and messages and return result
#'
#' @details Adapted from \url{http://www.onthelambda.com/2014/09/17/fun-with-rprofile-and-customizing-r-startup/}
#'
#' @param ... Inputs to keep quite
#'
#' @examples
#' sshhr( library(dplyr) )
#'
#' @export
sshhr <- function(...) suppressWarnings( suppressMessages( ... ) )


#' Is a character variable defined
#'
#' @details Is a variable NULL or an empty string
#'
#' @param x Character value to evaluate
#' @param empty Indicate what 'empty' means. Default is empty string (i.e., "")
#'
#' @return TRUE if empty, else FALSE
#'
#' @examples
#' is_empty("")
#' is_empty(NULL)
#' is_empty(NA)
#' is_empty(c())
#' is_empty("none", empty = "none")
#' is_empty("")
#' is_empty("   ")
#' is_empty(" something  ")
#'
#' @export

is_empty <- function(x, empty = "\\s*") if (is_not(x) || grepl(paste0("^",empty,"$"), x)) TRUE else FALSE

#' Is input a string?
#'
#' @details Is input a string
#'
#' @param x Input
#'
#' @return TRUE if string, else FALSE
#'
#' @examples
#' is_string("   ")
#' is_string("data")
#' is_string(c("data","data"))
#' is_string(NULL)
#'
#' @export
is_string <- function(x) if (length(x) == 1 && is.character(x) && !is_empty(x)) TRUE else FALSE

## check if a button was NOT pressed
not_pressed <- function(x) if (is.null(x) || x == 0) TRUE else FALSE

pressed <- function(x) if (!is.null(x) && x > 0) TRUE else FALSE

## drop elements from .._args variables obtained using formals
r_drop <- function(x, drop = c("dataset","data_filter")) x[-which(x %in% drop)]

#' Create a vector of interaction terms
#'
#' @param vars Variables lables to use
#' @param nway 2-way (2) or 3-way (3) interactions labels to create
#' @param sep Separator between variable names (default is :)
#'
#' @return Character vector of interaction term labels
#'
#' @examples
#' paste0("var", 1:3) %>% iterms(2)
#' paste0("var", 1:3) %>% iterms(3)
#' paste0("var", 1:3) %>% iterms(2, sep = ".")
#'
#' @export
iterms <- function(vars, nway, sep = ":") {
  if (!nway %in% c(2,3)) return(character(0))
  it <- c()
  for (i in 2:nway) {
    #it %<>% {c(., combn(vars, i) %>% apply(2, paste, collapse = sep))}
	it <- c(it, combn(vars, i) %>% apply(2, paste, collapse = sep))
    ## lm doesn't evaluate a:a
    # if (i == 2) it <- c(it, paste(vars, vars, sep = "*"))
    # if (i == 3) it <- c(it, paste(vars, vars, vars, sep = "*"))
  }
  it
}


#' Format a data.frame with a specified number of decimal places
#'
#' @param tbl Data.frame
#' @param dec Number of decimal places
#' @param perc Display numbers as percentages (TRUE or FALSE)
#' @param mark Thousand separator
#'
#' @return Data.frame for printing
#'
#' @examples
#' data.frame(x = c("a","b"), y = c(1L, 2L), z = c(-0.0005, 3)) %>%
#'   formatdf(dec = 3)
#' data.frame(x = c(1L, 2L), y = c(0.05, 0.8)) %>%
#'   formatdf(dec = 2, perc = TRUE)
#'
#' @export
formatdf <- function(tbl, dec = 3, perc = FALSE, mark = "") {

  frm <- function(x) {
    if (is.double(x)) {
      formatnr(x, dec = dec, perc = perc, mark = mark)
    } else if (is.integer(x)) {
      formatnr(x, dec = 0, mark = mark)
    } else {
      x
    }
  }

  mutate_each(tbl, funs(frm))
}

#' Format a number with a specified number of decimal places, thousand sep, and a symbol
#'
#' @param x Number or vector
#' @param sym Symbol to use
#' @param dec Number of decimal places
#' @param perc Display number as a percentage
#' @param mark Thousand separator
#'
#' @return Character (vector) in the desired format
#'
#' @examples
#' formatnr(2000, "$")
#' formatnr(2000, dec = 4)
#' formatnr(.05, perc = TRUE)
#' formatnr(c(.1, .99), perc = TRUE)
#' formatnr(data.frame(a = c(.1, .99)), perc = TRUE)
#' formatnr(data.frame(a = 1000), sym = "$", dec = 0)
#'
#' @export
formatnr <- function(x, sym = "", dec = 2, perc = FALSE, mark = ",") {
  if ("data.frame" %in% class(x)) x <- x[[1]]
  if (perc)
    paste0(sym, formatC(100 * x, digits = dec, big.mark = mark, format = "f"), "%")
  else
    paste0(sym, formatC(x, digits = dec, big.mark = mark, format = "f"))
}

#' Convenience function for is.null or is.na
#'
#' @param x Input
#'
#' @examples
#' is_not(NA)
#' is_not(NULL)
#' is_not(c())
#'
#' @export
is_not <- function(x) length(x) == 0 || is.na(x)


## fun_name is a string of the main function name
## rfun_name is a string of the reactive wrapper that calls the main function
## out_name is the name of the output, set to fun_name by default
register_print_output <- function(fun_name, rfun_name,
                                  out_name = fun_name) {

  ## Generate output for the summary tab
  output[[out_name]] <- renderPrint({
    ## when no analysis was conducted (e.g., no variables selected)
    get(rfun_name)() %>%
      {if (is.character(.)) cat(.,"\n") else .} %>%
      rm(.)
  })
  return(invisible())
}

# fun_name is a string of the main function name
# rfun_name is a string of the reactive wrapper that calls the main function
# out_name is the name of the output, set to fun_name by default
register_plot_output <- function(fun_name, rfun_name,
                                 out_name = fun_name,
                                 width = 650,
                                 height = 500) {

  ## Generate output for the plots tab
  output[[out_name]] <- renderPlot({

    ## when no analysis was conducted (e.g., no variables selected)
    get(rfun_name)() %>% { if (is.null(.)) " " else . } %>%
    { if (is.character(.)) {
        plot(x = 1, type = 'n', main = paste0("\n\n\n\n\n\n\n\n",.) ,
             axes = FALSE, xlab = "", ylab = "")
      } else {
        withProgress(message = 'Making plot', value = 1, print(.))
      }
    }
  }, width, height)

  return(invisible())
}

# minimum free energy of folding

## ViennaRNA RNAfold
folding_cal <- function(all_seqs = c()){
		num = length(all_seqs)
		ids <- paste(">tmp",1:num,sep="")
		new_seqs <- c()
		new_seqs[seq(1,2*num,2)] = ids
		new_seqs[seq(2,2*num,2)] = all_seqs
		write(new_seqs,"~/tmpseq.fasta")
		system("RNAfold < ~/tmpseq.fasta > ~/tmpseq.rnafold.txt");
	
		raw_fold_data <- readLines("~/tmpseq.rnafold.txt");
	
		mfe_vec <- c()
		for(i in seq(1,length(raw_fold_data),3)){
			mfes <- strsplit(raw_fold_data[i+2]," ")
			len <- length(mfes[[1]])
			mfe <- sub("[(]","",mfes[[1]][len])
			mfe <- sub("[)]","",mfe)
			mfe_vec <- c(mfe_vec,mfe)
		}
	
		deltaG <- data.frame(all_seqs,as.numeric(mfe_vec));

		system("rm ~/tmpseq.rnafold.txt");
		system("rm ~/tmpseq.fasta");
		system("rm ~/*.ps");

		return(deltaG);
}
