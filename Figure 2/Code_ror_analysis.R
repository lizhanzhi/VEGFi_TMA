library(data.table)
library(openxlsx)
library(stringr)
library(fst)
library(dplyr)
load("Demo_maligant_cancer.Rdata")
blood_pt <- read.xlsx("PT_TMA.xlsx")
blood_pt <- blood_pt %>% distinct(PT, .keep_all = TRUE)
target_pt <- unique(blood_pt$PT)
target_pt <- tolower(target_pt)
setDT(Demo)
Demo[, `:=`(DRUGNAME = tolower(DRUGNAME))]
library(stringr)
target_vegf_drugs <- list()
target_vegf_drugs$bevacizumab <- c("bevacizumab", "avastin", "mvasi", "zirabev", "blinded bevacizumab")
target_vegf_drugs$ranibizumab <- c("ranibizumab")
target_vegf_drugs$brolucizumab <- c("brolucizumab")
target_vegf_drugs$aflibercept <- c("aflibercept", "zaltrap", "eylea")
target_vegf_drugs$conbercept <- c("conbercept")
target_vegf_drugs$pegaptanib <- c("pegaptanib")
target_vegf_drugs$ramucirumab <- c("ramucirumab", "cyramza")
target_vegf_drugs$nintedanib <- c("nintedanib", "vargatef", "ofev")
target_vegf_drugs$apatinib <- c("apatinib", "aitan")
target_vegf_drugs$axitinib <- c("axitinib", "inlyta")
target_vegf_drugs$sunitinib <- c("sunitinib", "sutent", "sunitinib malate")
target_vegf_drugs$sorafenib <- c("sorafenib", "nexavar")
target_vegf_drugs$regorafenib <- c("regorafenib", "stivarga")
target_vegf_drugs$cabozantinib <- c("cabozantinib", "cabometyx", "xl184", "cometriq")
target_vegf_drugs$pazopanib <- c("pazopanib", "votrient 200mg", "votrient")
target_vegf_drugs$lenvatinib <- c("lenvatinib", "lenvima", "lenvatinib mesylate")
target_vegf_drugs$anlotinib <- c("anlotinib", "focus v")
target_vegf_drugs$fruquintinib <- c("fruquintinib", "elunate")
target_vegf_drugs$tivozanib <- c("tivozanib")
target_vegf_drugs$cediranib <- c("cediranib")
target_vegf_drugs$brivanib <- c("brivanib")
target_vegf_drugs$vandetanib <- c("vandetanib", "caprelsa")
target_vegf_drugs$vegf <- c("bevacizumab", "avastin", "mvasi", "zirabev", "blinded bevacizumab",
                            "ranibizumab", "brolucizumab", "aflibercept", "zaltrap", "eylea",
                            "conbercept", "pegaptanib")
target_vegf_drugs$vegfr <- c("ramucirumab", "cyramza", "nintedanib", "vargatef", "ofev",
                             "apatinib", "aitan", "axitinib", "inlyta", "sunitinib", "sutent",
                             "sunitinib malate", "sorafenib", "nexavar", "regorafenib", "stivarga",
                             "vandetanib", "caprelsa", "cabozantinib", "cabometyx", "xl184", "cometriq",
                             "pazopanib", "votrient 200mg", "votrient", "lenvatinib", "lenvima",
                             "lenvatinib mesylate", "anlotinib", "focus v", "fruquintinib", "elunate",
                             "tivozanib", "cediranib", "brivanib")
all_drug_names <- c(target_vegf_drugs$vegf, target_vegf_drugs$vegfr)
target_vegf_drugs$all_drug <- all_drug_names
for (drug_class in names(target_vegf_drugs)) {
  Demo[[paste0(drug_class, ".all")]] <- 0
  for (drug_name in target_vegf_drugs[[drug_class]]) {
    escaped_drug_name <- str_replace_all(drug_name, "([\\\\\\.\\+\\*\\?\\[\\^\\]\\$\\(\\)\\{\\}\\=\\!\\<\\>\\|\\:\\-])", "\\\\\\1")
    Demo[[paste0(drug_class, ".all")]] <- pmax(Demo[[paste0(drug_class, ".all")]], as.integer(grepl(escaped_drug_name, Demo$DRUGNAME, ignore.case = TRUE)))
  }
}
vegf_columns <- grep(".all$", names(Demo), value = TRUE)
escaped_target_vegf_drugs <- sapply(target_vegf_drugs, function(drugs) {
  sapply(drugs, function(drug) {
    str_replace_all(drug, "([\\\\\\.\\+\\*\\?\\[\\^\\]\\$\\(\\)\\{\\}\\=\\!\\<\\>\\|\\:\\-])", "\\\\\\1")
  })
})
all_drug_pattern <- paste(unlist(escaped_target_vegf_drugs), collapse = "|")
Demo[, VEGF := as.integer(grepl(all_drug_pattern, DRUGNAME, ignore.case = TRUE))]
Demo$PT_split <- strsplit(tolower(Demo$PT), "\\,")
Demo$TargetPT.all <- as.integer(sapply(Demo$PT_split, function(pt_array) {
  any(target_pt %in% pt_array)
}))
table(Demo$TargetPT.all)
head(Demo)
library(data.table)
library(parallel)
result_df <- data.frame(
  DrugClass = character(),
  PT = character(),
  ROR = numeric(),
  CI_up_ROR = numeric(),
  CI_down_ROR = numeric(),
  PRR = numeric(),
  CI_up_PRR = numeric(),
  CI_down_PRR = numeric(),
  IC = numeric(),
  IC025 = numeric(),
  P = numeric(),
  E = numeric(),
  a = integer(),
  b = integer(),
  c = integer(),
  d = integer(),
  EBGM = numeric(),
  CI_up_EBGM = numeric(),
  CI_down_EBGM = numeric(),
  EBGM05 = numeric(),
  Signal = logical(),
  stringsAsFactors = FALSE
)
calculate_metrics <- function(drug_class, pt, Demo) {
  pt_lower <- tolower(pt)
  Demo[, TargetPT.all := sapply(PT_split, function(x) pt_lower %in% x)]
  num_a <- as.numeric(sum(Demo$chemo.all == 1 & Demo$TargetPT.all == 1, na.rm = TRUE))
  num_b <- as.numeric(sum(Demo$chemo.all == 1 & Demo$TargetPT.all == 0, na.rm = TRUE))
  num_c <- as.numeric(sum(Demo$chemo.all == 0 & Demo$TargetPT.all == 1, na.rm = TRUE))
  num_d <- as.numeric(sum(Demo$chemo.all == 0 & Demo$TargetPT.all == 0, na.rm = TRUE))
  num_total_a <- num_a + num_c
  total_reports <- nrow(Demo)
  ror <- NA
  CI_up_ror <- NA
  CI_down_ror <- NA
  prr <- NA
  CI_up_prr <- NA
  CI_down_prr <- NA
  ic <- NA
  P <- NA
  E <- NA
  ic025 <- NA
  if (!is.na(num_a) && !is.na(num_b) && !is.na(num_c) && !is.na(num_d) &&
      num_a > 0 & num_b > 0 & num_c > 0 & num_d > 0) {
    ror <- (num_a * num_d) / (num_b * num_c)
    ln_ror <- log(ror)
    se_ln_ror <- sqrt(1/num_a + 1/num_b + 1/num_c + 1/num_d)
    CI_up_ror <- exp(ln_ror + 1.96 * se_ln_ror)
    CI_down_ror <- exp(ln_ror - 1.96 * se_ln_ror)
  }
  if (!is.na(num_a) && !is.na(num_b) && !is.na(num_c) && !is.na(num_d) &&
      (num_a + num_c) > 0 && (num_b + num_d) > 0) {
    prr <- (num_a / (num_a + num_b)) / (num_c / (num_c + num_d))
    P <- num_a / (num_a + num_b)
    E <- num_c / (num_c + num_d)
    se_ln_prr <- sqrt(1/num_a - 1/(num_a + num_b) + 1/num_c - 1/(num_c + num_d))
    ln_prr <- log(prr)
    CI_up_prr <- exp(ln_prr + 1.96 * se_ln_prr)
    CI_down_prr <- exp(ln_prr - 1.96 * se_ln_prr)
  }
  if (!is.na(num_total_a) && !is.na(total_reports) && 
      num_total_a > 0 && total_reports > 0) {
    expected <- (num_total_a * (num_a + num_b)) / total_reports
    N_obs <- num_a
    ic <- log2((N_obs + 0.5) / (expected + 0.5))
    ic025 <- ic - 3.3 * (N_obs + 0.5)^(-1/2) - 2 * (N_obs + 0.5)^(-3/2)
  }
  ebgm <- ifelse((num_a + num_c) > 0 & (num_a + num_b) > 0, 
                 (num_a * (num_a + num_b + num_c + num_d)) / ((num_a + num_c) * (num_a + num_b)), NA)
  CI_up_ebgm <- ifelse(!is.na(ebgm), exp(log(ebgm) + 1.96 * sqrt(1/num_a + 1/num_b + 1/num_c + 1/num_d)), NA)
  CI_down_ebgm <- ifelse(!is.na(ebgm), exp(log(ebgm) - 1.96 * sqrt(1/num_a + 1/num_b + 1/num_c + 1/num_d)), NA)
  ebgm05 <- ifelse(!is.na(ebgm), ebgm / exp(1.96 * sqrt(1/num_a + 1/num_b + 1/num_c + 1/num_d)), NA)
  signal <- ifelse(!is.na(ebgm05) & ebgm05 > 2, TRUE, FALSE)
  result_list <- list(
    DrugClass = drug_class,
    PT = pt,
    ROR = ror,
    CI_up_ROR = CI_up_ror,
    CI_down_ROR = CI_down_ror,
    PRR = prr,
    CI_up_PRR = CI_up_prr,
    CI_down_PRR = CI_down_prr,
    IC = ic,
    IC025 = ic025,
    P = P,
    E = E,
    a = num_a,
    b = num_b,
    c = num_c,
    d = num_d,
    EBGM = ebgm,
    CI_up_EBGM = CI_up_ebgm,
    CI_down_EBGM = CI_down_ebgm,
    EBGM05 = ebgm05,
    Signal = signal
  )
  return(result_list)
}
no_cores <- detectCores() - 16
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(data.table))
tasks <- expand.grid(drug_class = names(target_vegf_drugs), pt = blood_pt$PT)
clusterExport(cl, varlist = c("calculate_metrics", "Demo", "tasks"))
result_list <- parLapply(cl, seq_len(nrow(tasks)), function(i) {
  drug_class <- tasks$drug_class[i]
  pt <- tasks$pt[i]
  calculate_metrics(as.character(drug_class), as.character(pt), Demo)
})
stopCluster(cl)
result_df <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
