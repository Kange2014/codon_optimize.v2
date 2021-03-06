#' Create (partial) factorial design
#'
#'
#' @param factors: Categorical variables used as input for design
#' @param int: Vector of interaction terms to consider when generating design
#' @param trials: Number of trial to create. If NA then all feasible designs will be considered until a design with perfect D-efficiency is found
#' @param seed: Random seed to use as the starting point
#'
#' @return A list with all variables defined in the function as an object of class doe
#'
#'
#' @importFrom AlgDesign optFederov
#' @importFrom mvtnorm pmvnorm
#' @importFrom polycor hetcor
#'

doe <- function(factors, int = "", trials = NA, seed = NA) {

  df_list <-
    gsub("[ ]{2,}"," ",factors) %>%
    gsub("/","",.) %>%
    gsub("\\\\n","\n",.) %>%
    gsub("[ ]*;[ ]*",";",.) %>%
    gsub(";{2,}",";",.) %>%
    gsub("[;]+[ ]{0,}\n","\n",.) %>%
    gsub("[ ]{1,}\n","\n",.) %>%
    gsub("\n[ ]+","\n",.) %>%
    gsub("[\n]{2,}","\n",.) %>%
    gsub("[ ]+","_",.) %>%
    strsplit(.,"\n") %>%
    .[[1]] %>% strsplit(";")

  df_names <- c()
  for (i in seq_len(length(df_list))) {
    dt <- df_list[[i]] %>% gsub("^\\s+|\\s+$", "", .)
    df_names <- c(df_names, dt[1])
    df_list[[i]] <- dt[-1]
  }
  names(df_list) <- df_names
  model <- paste0("~ ", paste0(df_names, collapse = " + "))
  nInt <- 0
  if (!is_empty(int)) {
    model <- paste0(model, " + ", paste0(int, collapse = " + "))
    nInt <- length(int)
  }

  part_fac <- function(df, model = ~ ., int = 0, trials = NA, seed = 172110) {

	  full <- expand.grid(df)

	  ###############################################
	  # eliminate combinations from full
	  # by removing then from the variable _experiment_
	  # http://stackoverflow.com/questions/18459311/creating-a-fractional-factorial-design-in-r-without-prohibited-pairs?rq=1
	  ###############################################

	  levs <- sapply(df, length)
	  nr_levels <- sum(levs)
	  min_trials <- nr_levels - length(df) + 1
	  max_trials <- nrow(full)

	  if (!is.null(trials) && !is.na(trials)) max_trials <- min_trials <- trials

	  eff <-
	    data.frame(
	      Trials = min_trials:max_trials,
	      "D-efficiency" = NA,
	      # "Determinant" = NA,
	      "Balanced" = NA,
	      check.names = FALSE
	    )

	  for (i in min_trials:max_trials) {
      seed %>% gsub("[^0-9]","",.) %>% { if (!is_empty(.)) set.seed(seed) }
	    design <- try(AlgDesign::optFederov(model, data = full, nRepeats = 50,
                    nTrials = i, maxIteration=1000,
                    approximate = FALSE), silent = TRUE)

	    if (is(design, 'try-error')) next
	    # cor_mat <- cor(data.matrix(design$design))
	    # detcm <- det(cor_mat)
	    ind <- which(eff$Trials %in% i)
	    eff[ind,"D-efficiency"] <- design$Dea
	    # eff[ind,"Determinant"] <- round(detcm,3)
	    eff[ind,"Balanced"] <-  all(i %% levs == 0)

	    if (design$Dea == 1) break
	  }

    cor_mat <- sshhr(polycor::hetcor(design$design, std.err = FALSE)$correlations)
    detcm <- det(cor_mat)

	  if (exists("cor_mat")) {
	    list(df = df, cor_mat = cor_mat, detcm = detcm, Dea = design$Dea,
	         part = arrange_(design$design, .dots = names(df)),
	         full = arrange_(full, .dots = names(df)),
	         eff = na.omit(eff),
	         seed = seed)
	  } else if (!is.na(trials)) {
	    "No solution exists for the selected number of trials"
	  } else {
	    "No solution found"
	  }
	}

  part_fac(df_list, model = as.formula(model), int = nInt, trials = trials, seed = seed) %>%
    add_class("doe")
}

#' Summary method for doe function
#'
#' @param object: Return value from \code{\link{doe}}
#' @param eff: If TRUE print efficiency output
#' @param part: If TRUE print partial factorial
#' @param full: If TRUE print full factorial
#' @param ...: further arguments passed to or from other methods.
#'

summary.doe <- function(object, eff = TRUE, part = TRUE, full = TRUE, ...) {

  if (!is.list(object)) return(object)

  cat("Experimental design\n")
  cat("# trials for partial factorial:", nrow(object$part),"\n")
  cat("# trials for full factorial   :", nrow(object$full),"\n")
  if (!is.null(object$seed) && !is.na(object$seed))
    cat("Random seed                   :", object$seed,"\n")

  cat("\nAttributes and levels:\n")
  nl <- names(object$df)
  for (i in nl) {
    cat(paste0(i, ":"), paste0(object$df[[i]], collapse = ", "), "\n")
  }

  if (eff) {
    cat("\nDesign efficiency:\n")
    print(formatdf(object$eff, dec = 3), row.names = FALSE)
  }

  if (part) {
    cat("\nPartial factorial design correlations:\n")
    nrdec <- ifelse (object$detcm == 1, 0, 3)
    # print(formatdf(data.frame(object$cor_mat), dec = nrdec) , row.names = FALSE)
    print(round(object$cor_mat, nrdec) , row.names = FALSE)

    cat("\nPartial factorial design:\n")
    print(object$part)
  }

  if (full) {
    cat("\nFull factorial design:\n")
    print(object$full)
  }
}
