#' Function that takes a formula and data objects and creates a flexible `GenCovINLA_re` object, storing
#' all the information needed to run the `GenCovINLA` model, e.g. fixed and random effects.
#' @param formula A `GenCovINLA` formula object specifying a model.
#' @param data A data.frame with the data to use in the model
#' @param cov_data A named list containing the covariance matrices referenced in the formula. 
#' List element names should be match those reference in the formula (e.g. after any "__"). Note that
#' all matrices in `cov_data` must have named columns and rows corresponding to labels in `data` (e.g. before any "__")
#' @importFrom stats as.formula
#' @export
GenCov_re_construct <- function(formula, data, cov_data) {
  data <- as.data.frame(data)
  # formula <- abund ~ 1 + (1|sp__phy) + (1 + x|sp__phy@site__iid) + (0 + k|sp__iid) + (1|island__iid)
  # formula <- abund ~ 1 + (1|sp__phy) + (1 + x|sp__phy@site__phy) + (0 + k|sp__iid)
  fm <- lme4::findbars(formula)
  fm_char <- lapply(fm, as.character)
  fm_char = lapply(fm_char, function(x) gsub(pattern = "^0 ?[+] ?", replacement = "", x))
  pluses <- sapply(fm_char, function(x) any(grepl("[+]", x)))
  if(any(pluses)) {
    warning("GenCovINLA currently does not support correlated slopes in random effects (e.g. (1+x|site)), all slopes will be modeled independently (e.g. (1|site) + (x|site))...")
    for(i in 1:sum(pluses)) {
      form <- fm[[which(pluses)[i]]]
      new_form <- paste0("y ~ (",Reduce(paste, deparse(form)), ")")
      new_form <- gsub("|", "||", new_form, fixed = TRUE)
      split_fm <- lme4::findbars(as.formula(new_form))
      split_fm = lapply(split_fm, function(x) gsub(pattern = "^0 ?[+] ?", replacement = "", x))
      fm_char <- fm_char[-which(pluses)[i]]
      fm_char <- c(fm_char, split_fm)
    }
  }
  
  random_effects <- list()
  #nested <- numeric(length(fm_char))
  for(i in seq_along(fm_char)) {
    weights <- NULL
    replicates <- NULL
    covar_1 <- NULL
    covar_2 <- NULL
    ids_1 <- NULL
    id2_2 <- NULL
    ids_big <- NULL
    covar_big <- NULL
    side_2_is_covar <- FALSE
    
    fm_curr <- fm_char[[i]]
    
    if(fm_curr[2] == "1") {
      weights <- NULL
    } else {
      weights <- data[ , fm_curr[2]]
    }
    
    two_sides <- strsplit(fm_curr[3], "[@]")[[1]]
    if(length(two_sides) > 1) {
      nested <- TRUE
    } else {
      nested <- FALSE
    }
    two_sides <- strsplit(two_sides, "__")
    
    side_1 <- two_sides[[1]]
    if(length(side_1) == 1) {
      iid_1 = TRUE
    } else {
      if(side_1[2] == "iid") {
        iid_1 <- TRUE
      } else {
        iid_1 <- FALSE
      }
    }
    
    if(iid_1) {
      covar_1 <- diag(dplyr::n_distinct(data[ , side_1[1]]))
      ids_1 <- as.numeric(as.factor(data[ , side_1[1]]))
    } else {
      covar_1 <- cov_data[[side_1[2]]]
      by_it <- "var"
      names(by_it) <- side_1[1]
      new_dat <- data[ , side_1[1], drop = FALSE]
      new_dat <- dplyr::left_join(new_dat, dplyr::data_frame(var = rownames(covar_1),
                                    ids = seq_len(nrow(covar_1))), by = by_it)
      ids_1 <- new_dat$ids
    }
    
    if(nested) {
      
      side_2 <- two_sides[[2]]
      if(length(side_2) == 1) {
        iid_2 = TRUE
      } else {
        if(side_2[2] == "iid") {
          iid_2 <- TRUE
        } else {
          iid_2 <- FALSE
        }
      }
      
      if(iid_2) {
        replicates <- as.numeric(as.factor(data[ , side_2[1]]))
        covar_2 <- NULL
        side_2_is_covar <- FALSE
      } else {
        
        if(dplyr::n_distinct(data[ , side_2[1]]) * dplyr::n_distinct(data[ , side_1[1]]) != nrow(data)) {
          stop(paste("nesting with right-side covariance currently 
                     doesn't support hierarchical factors. Try removing", side_2[1], " | ", side_2[2]))
        } else {
          covar_2 <- cov_data[[side_2[2]]]
          by_it <- "var"
          names(by_it) <- side_2[1]
          new_dat <- data[ , side_2[1], drop = FALSE]
          new_dat <- dplyr::left_join(new_dat, dplyr::data_frame(var = rownames(covar_2),
                                                                 ids = seq_len(nrow(covar_2))), by = by_it)
          ids_2 <- new_dat$ids
          side_2_is_covar <- TRUE
          replicates <- NULL
        }
      }
    }
    
    if(!side_2_is_covar) {
      random_effects[[i]] <- list(ids = ids_1, weights = weights, covar = covar_1, replicates = replicates)
    } else {
      new_dat <- data
      new_dat$id <- 1:nrow(new_dat)
      new_dat <- dplyr::arrange(new_dat, ids_1, ids_2)
      new_dat$ids <- 1:nrow(new_dat)
      new_dat <- dplyr::arrange(new_dat, id)
      covar_big <- kronecker(covar_1, covar_2)
      ids_big <- new_dat$ids
      random_effects[[i]] <- list(ids = ids_big, weights = weights, covar = covar_big, replicates = NULL)
    }
  }
  
  names(random_effects) <- sapply(fm_char, function(x) paste(x[2], x[3], sep = " | "))
  
  list(formula = lme4::nobars(formula), random_effects = random_effects)
  
}