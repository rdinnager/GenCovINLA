#' Function that takes a formula and data objects and creates a flexible `GenCovINLA_re` object, storing
#' all the information needed to run the `GenCovINLA` model, e.g. fixed and random effects.
#' @param formula A `GenCovINLA` formula object specifying a model.
#' @param data A data.frame with the data to use in the model
#' @param cov_data A named list containing the covariance matrices referenced in the formula. 
#' List element names should be match those reference in the formula (e.g. after any "__")
#' @export
GenCov_re_construct <- function(formula, data, cov_data) {
  # formula <- y ~ 1 + (1|sp__phy) + (1 + x|sp__phy2@site__iid) + (0 + k|sp__iid)
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
    fm_curr <- fm_char[[i]]
    two_sides <- strsplit(fm_curr[3], "[@]")
    if(length(two_sides) > 1) {
      nested <- TRUE
    } else {
      nested <- FALSE
    }
    two_sides <- lapply(two_sides, function(x) unlist(strsplit(x, "__"))) 
    
    side_1 <- two_sides[[1]]
    if(length(side_1) == 1) {
      iid_1 = TRUE
    } else {
      if(side_1[2] == "iid") {
        iid_1 <- TRUE
      }
    }
    
    if(iid) {
      covar_1 <- diag(dplyr::n_distinct(data[ , side_1[1]))
    } else {
      covar_1 <- cov_data[[side_1[2]]]
    }
    
    if(nested) {
      
      
    } 
  }
}