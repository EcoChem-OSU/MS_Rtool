
### Modified version of NIST function
## https://github.com/usnistgov/NISTPFAS/blob/main/suspectlist/fn/calculate_residual.R

### 2024-03-10 - B Droz @ Oregon State University
################################################################################
####  **** Modified NIST version 
require(MetaboCoreUtils)

extract_elements <- function(composition.str) {
  single.elem <- gregexpr("[A-Z]", composition.str)[[1]]
  double.elem <- gregexpr("[A-Z][a-z]", composition.str)[[1]]
  single.elem <- single.elem[! single.elem %in% double.elem]
  if (length(single.elem) > 0) {
    elements <- substring(composition.str, single.elem, single.elem)
  }
  if (double.elem[1] != -1) {
    elements <- c(elements, substring(composition.str, double.elem, double.elem + 1))
  }
  ecounts <- rep(1, length(elements))
  nums <- gregexpr("[0-9]+", composition.str)[[1]]
  nums.count <- attr(nums, "match.length")
  for (i in 1:length(nums)) {
    if (substr(composition.str, nums[i] - 1, nums[i] - 1) %in% elements) {
      ecounts[which(elements == substr(composition.str, nums[i] - 1, nums[i] - 1))] <- substr(composition.str, nums[i], nums[i] + nums.count[i] - 1)
    }
    if (substr(composition.str, nums[i] - 2, nums[i] - 1) %in% elements) {
      ecounts[which(elements == substr(composition.str, nums[i] - 2, nums[i] - 1))] <- substr(composition.str, nums[i], nums[i] + nums.count[i] - 1)
    }
  }
  results <- c()
  results$elements <- elements
  results$counts <- as.integer(ecounts)
  results
}

collapse_ecomp <- function(ecomp) {
  if (length(ecomp$elements) == 0) {return(NA)}
  paste(sapply(1:length(ecomp$elements), function(x) {if (ecomp$counts[x] > 0)
    {o <- paste(ecomp$elements[x], ecomp$counts[x], sep = "")}; if (ecomp$counts[x] <= 0) {o <- ""};
    o}), collapse = "")
}


# modified version work only for neutral
calculate_residual <- function(ecomp, rep_unit = "CF2", exactmasses, adduct = "neutral") {
  
      parent_ecomp <- extract_elements(ecomp)
      rep_ecomp <- extract_elements(rep_unit)
      if (length(which(rep_ecomp$elements %in% parent_ecomp$elements)) != length(rep_ecomp$elements)) 
        {return("You do not have all of the elements in the repeating units in the elemental formula")}
      parent_repeating <- list(elements = parent_ecomp$elements[which(parent_ecomp$elements %in% rep_ecomp$elements)], 
                               counts = parent_ecomp$counts[which(parent_ecomp$elements %in% rep_ecomp$elements)])
      
      max_rep <- sapply(parent_repeating$elements, 
                        function(x) floor(parent_repeating$counts[which(parent_repeating$elements == x)]/rep_ecomp$counts[which(rep_ecomp$elements == x)]))
      rep_num <- min(max_rep)
      parent_remove <- list(elements = rep_ecomp$elements, 
                            counts = rep_num*rep_ecomp$counts)
      residual_ecomp <- parent_ecomp
      for (i in parent_remove$elements) {
        residual_ecomp$counts[which(residual_ecomp$elements == i)] <-  
          residual_ecomp$counts[which(residual_ecomp$elements == i)] - 
          parent_remove$counts[which(parent_remove$elements == i)]
      }
      if (length(which(residual_ecomp$counts == 0)) > 0) {
        residual_ecomp$elements <- residual_ecomp$elements[-which(residual_ecomp$counts == 0)]
        residual_ecomp$counts <- residual_ecomp$counts[-which(residual_ecomp$counts == 0)]
      }
      data.frame(parent = ecomp, rep_unit = rep_unit, 
                 rep_mass = calculateMass(collapse_ecomp(rep_ecomp)),
                   # calculate_monoisotope(rep_ecomp, exactmasses, adduct = "neutral"), 
                  rep_num = rep_num,
                 residual = collapse_ecomp(residual_ecomp), residual_charge = adduct, 
                 residual_mass = calculateMass(collapse_ecomp(residual_ecomp)) )
                 # calculate_monoisotope(residual_ecomp, exactmasses, adduct))
    }