#' Calculate CSR Scores and Assign a CSR Strategy Using the Hodgson et al. (1999) Model
#'
#' This function calculates and assigns CSR strategy types to plant species
#' based on the Hodgson et al. (1999) model, as well as performing
#' the required trait transformations for the model.
#'
#' @param data A data frame containing species names and trait data. Must include (LFW and LA optional instead of LDMC and SLA respectively):
#' \describe{
#'   \item{CH}{Canopy height (mm)}
#'   \item{LDMC}{Leaf dry matter content (\%)}
#'   \item{FP}{Flowering period (months)}
#'   \item{LS}{Lateral spread (categorical scale: 1–6, based on Hodgson et al. 1999)}
#'   \item{LDW}{Leaf dry weight (mg)}
#'   \item{SLA}{Specific leaf area (mm²/mg)}
#'   \item{LFW}{Leaf fresh weight (mg)}
#'   \item{LA}{Leaf area (mm²)}
#' }
#' @param calcLDMC Logical. Whether to calculate LDMC from LDW and LFW. Default is FALSE.
#' @param calcSLA Logical. Whether to calculate SLA from LA and LDW. Default is FALSE.
#'
#' @return A data frame with CSR scores (-2.5 to 2.5), CSR percentages (C%, S%, R%), and assigned strategy classification.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' # hodgson(exampleData)

hodgson <- function(data, calcLDMC = FALSE, calcSLA = FALSE) {

  # column renaming helper function
  renameColumnsHodgson <- function(data) {
    renameDict <- list(
      "species" = c("^species$", "^sp\\.$", "^sp$", "^taxon$", "^species[._ ]?name$", "^taxon[._ ]?name$"),
      "CH" = c("^canopy[._ ]?height$", "^CH$"),
      "LDMC" = c("^leaf[._ ]?dry[._ ]?matter[._ ]?content$", "^LDMC$"),
      "FP" = c("^flowering[._ ]?period$", "^FP$"),
      "LS" = c("^lateral[._ ]?spread$", "^LS$"),
      "LDW" = c("^leaf[._ ]?dry[._ ]?weight$", "^LDW$"),
      "SLA" = c("^specific[._ ]?leaf[._ ]?area$", "^SLA$"),
      "LFW" = c("^leaf[._ ]?fresh[._ ]?weight$", "^LFW$"),
      "LA"  = c("^leaf[._ ]?area$", "^LA$")
    )

    for (colName in names(renameDict)) {
      matchedCols <- grep(paste(renameDict[[colName]], collapse = "|"), names(data), ignore.case = TRUE, value = TRUE)
      if (length(matchedCols) > 0) {
        names(data)[names(data) %in% matchedCols] <- colName
      }
    }

    return(data)
  }

  # rename columns
  data <- renameColumnsHodgson(data)

  # check species column
  if (!"species" %in% names(data)) {
    stop("Data must contain a 'species' column.")
  }

  # required columns
  if (calcLDMC && calcSLA) {
    requiredCols <- c("species", "CH", "FP", "LS", "LDW", "LFW", "LA")
  } else if (calcLDMC) {
    requiredCols <- c("species", "CH", "FP", "LS", "LDW", "LFW", "SLA")
  } else if (calcSLA) {
    requiredCols <- c("species", "CH", "FP", "LS", "LDW", "LDMC", "LA")
  } else {
    requiredCols <- c("species", "CH", "FP", "LS", "LDW", "LDMC", "SLA")
  }

  # check column presence
  if (!all(requiredCols %in% names(data))) {
    stop(paste("Data must contain:", paste(requiredCols, collapse = ", ")))
  }

  # remove rows with missing values
  missingRows <- data %>% filter(if_any(all_of(requiredCols), ~ is.na(.)))
  if (nrow(missingRows) > 0) {
    removedSpecies <- paste(unique(missingRows$species), collapse = ", ")
    warning(sprintf(
      "%d row(s) removed due to missing trait values. affected species: %s",
      nrow(missingRows),
      removedSpecies
    ))
  }
  data <- anti_join(data, missingRows, by = "species")

  # ldmc calculation
  if (calcLDMC) {
    data <- data %>% mutate(LDMC = (LDW * 100) / LFW)
  }

  # sla calculation
  if (calcSLA) {
    data <- data %>% mutate(SLA = LA / LDW)
  }

  # trait transformations
  calculations <- data %>%
    mutate(
      chPr  = case_when(
        CH > 999 ~ 6,
        CH > 599 ~ 5,
        CH > 299 ~ 4,
        CH > 99  ~ 3,
        CH > 49  ~ 2,
        TRUE     ~ 1
      ),
      ldmcPr = sqrt(LDMC),
      fpPr   = FP,
      lsPr   = LS,
      ldwPr  = log(LDW) + 3,
      slaPr  = sqrt(SLA)
    )

  # raw csr dimension scores
  calculations <- calculations %>%
    mutate(
      rawC = (0.141 * chPr^2) + (0.09061 * lsPr^2),
      rawS = 54.6 - (1.666 * chPr^2) + (1.069 * ldmcPr^2) - (2.732 * slaPr^2) + (1.722 * lsPr^2),
      rawR = (2.518 * fpPr) - (2.748 * ldwPr) + (5.37 * slaPr)
    )

  # csr score calculation
  calculations <- calculations %>%
    mutate(
      cScore = pmin(pmax(-2.5 + 0.839 * rawC, -2.5), 2.5),
      sScore = pmin(pmax(-1.103 + 0.0474 * rawS, -2.5), 2.5),
      rScore = pmin(pmax(-2.5 + 0.119 * rawR, -2.5), 2.5)
    )

  # convert scores to percentages
  calculations <- calculations %>%
    mutate(
      cPercent = 100 * (cScore + 2.5) / (cScore + sScore + rScore + 7.5),
      sPercent = 100 * (sScore + 2.5) / (cScore + sScore + rScore + 7.5),
      rPercent = 100 * (rScore + 2.5) / (cScore + sScore + rScore + 7.5)
    )

  # strategy classification
  csrReference <- data.frame(
    strategy = c("C", "C/CR", "C/SC", "CR", "C/CSR", "SC", "CR/CSR", "SC/CSR", "R/CR",
                 "CSR", "S/SC", "R/CSR", "S/CSR", "R", "SR/CSR", "S", "R/SR", "S/SR", "SR"),
    cRef = c(2, 1, 1, 0, 1, 0, 0, 0, -1, 0, -1, -1, -1, -2, -1, -2, -2, -2, -2),
    sRef = c(-2, -2, -1, -2, -1, 0, -1, 0, -2, 0, 1, -1, 1, -2, 0, 2, -1, 1, 0),
    rRef = c(-2, -1, -2, 0, -1, -2, 0, -1, 1, 0, -2, 1, -1, 2, 0, -2, 1, -1, 0)
  )

  classifySpecies <- function(c, s, r) {
    distances <- sqrt((csrReference$cRef - c)^2 + (csrReference$sRef - s)^2 + (csrReference$rRef - r)^2)
    return(csrReference$strategy[which.min(distances)])
  }

  calculations <- calculations %>%
    rowwise() %>%
    mutate(strategyClass = classifySpecies(cScore, sScore, rScore)) %>%
    ungroup()

  # output columns
  outputCols <- c("cScore", "sScore", "rScore", "cPercent", "sPercent", "rPercent", "strategyClass")
  if (calcLDMC) outputCols <- c("LDMC", outputCols)
  if (calcSLA) outputCols <- c("SLA", outputCols)

  # combine with original data
  finalResult <- cbind(
    data,
    calculations %>% dplyr::select(dplyr::any_of(outputCols))
  )

  return(finalResult)
}
