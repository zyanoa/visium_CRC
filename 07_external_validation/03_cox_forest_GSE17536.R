#!/usr/bin/env Rscript

# Project: CRC MMR spatial immune ecology
# Purpose: Univariate Cox regression forest plot for GSE17536 clinical variables and immune signatures.
# Note: This reproduces the logic of the original ezcox/show_forest block as univariate Cox analysis.

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(ggplot2)
  library(survminer)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    score_table = "results/07_external_validation/survival/GSE17536_clinical_signature_scores.csv",
    output_dir = "results/07_external_validation/cox",
    variables = "Age,Sex,Stage_group,dMMR_TLS_program,Plasma_IgG,B-cells,Plasma cells"
  )
  for (arg in args) {
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2 && kv[1] %in% names(out)) out[[kv[1]]] <- kv[2]
  }
  out
}

args <- parse_args()
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(args$score_table, check.names = FALSE, row.names = 1)
if ("TLS" %in% colnames(df)) df[["TLS"]] <- NULL
if ("Stroma4_1" %in% colnames(df) && !"dMMR_TLS_program" %in% colnames(df)) {
  df[["dMMR_TLS_program"]] <- df[["Stroma4_1"]]
  df[["Stroma4_1"]] <- NULL
}

# Normalize clinical variables when present.
df$Time <- as.numeric(df$Time)
df$Event <- as.numeric(df$Event)

if ("Age" %in% colnames(df)) {
  df$Age <- suppressWarnings(as.numeric(df$Age))
  df$Age_group <- ifelse(df$Age < 60, "<60", ">=60")
}

if ("Stage" %in% colnames(df)) {
  stage_chr <- as.character(df$Stage)
  stage_num <- suppressWarnings(as.numeric(stage_chr))
  stage_num[is.na(stage_num)] <- match(stage_chr[is.na(stage_num)], c("I", "II", "III", "IV"))
  df$Stage_group <- ifelse(stage_num %in% c(1, 2), "I,II", "III,IV")
}

make_high_low <- function(dat, variable) {
  if (!variable %in% colnames(dat)) return(dat)
  score <- suppressWarnings(as.numeric(dat[[variable]]))
  if (all(is.na(score))) return(dat)
  tmp <- data.frame(Time = dat$Time, Event = dat$Event, score = score)
  tmp <- tmp[complete.cases(tmp), ]
  cut <- tryCatch(survminer::surv_cutpoint(tmp, time = "Time", event = "Event", variables = "score"), error = function(e) NULL)
  cutoff <- if (!is.null(cut)) summary(cut)$cutpoint[1, "cutpoint"] else median(score, na.rm = TRUE)
  new_col <- paste0(make.names(variable), "_group")
  dat[[new_col]] <- ifelse(score > cutoff, "High", "Low")
  dat[[new_col]] <- factor(dat[[new_col]], levels = c("Low", "High"))
  attr(dat[[new_col]], "cutoff") <- cutoff
  dat
}

signature_vars <- intersect(c("dMMR_TLS_program", "Plasma_IgG", "B-cells", "Plasma cells", "Tregs"), colnames(df))
for (v in signature_vars) df <- make_high_low(df, v)

write.csv(df, file.path(args$output_dir, "GSE17536_cox_input_table.csv"), row.names = TRUE, quote = FALSE)

# Default covariates corresponding to the original show_forest/ezcox block.
requested <- trimws(strsplit(args$variables, ",")[[1]])
cox_vars <- c()
for (v in requested) {
  if (v %in% c("dMMR_TLS_program", "Plasma_IgG", "B-cells", "Plasma cells", "Tregs")) {
    g <- paste0(make.names(v), "_group")
    if (g %in% colnames(df)) cox_vars <- c(cox_vars, g)
  } else if (v == "Age" && "Age_group" %in% colnames(df)) {
    cox_vars <- c(cox_vars, "Age_group")
  } else if (v %in% colnames(df)) {
    cox_vars <- c(cox_vars, v)
  }
}
cox_vars <- unique(cox_vars)

run_uni_cox <- function(dat, variable) {
  f <- as.formula(paste0("Surv(Time, Event) ~ `", variable, "`"))
  fit <- tryCatch(coxph(f, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  if (nrow(s$coef) == 0) return(NULL)
  data.frame(
    variable = variable,
    term = rownames(s$coef),
    HR = s$coef[, "exp(coef)"],
    lower_95 = s$conf.int[, "lower .95"],
    upper_95 = s$conf.int[, "upper .95"],
    p.value = s$coef[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  )
}

cox_results <- bind_rows(lapply(cox_vars, function(v) run_uni_cox(df, v)))
if (nrow(cox_results) == 0) stop("No Cox models could be fitted. Check input variables.")
cox_results$`HR (95% CI)` <- sprintf("%.2f (%.2f-%.2f)", cox_results$HR, cox_results$lower_95, cox_results$upper_95)
cox_results$p_label <- ifelse(cox_results$p.value < 0.001, "<0.001", sprintf("%.3f", cox_results$p.value))
write.csv(cox_results, file.path(args$output_dir, "GSE17536_univariate_cox_results.csv"), row.names = FALSE, quote = FALSE)

# Forest plot using ggplot2 to avoid manual CSV editing and reduce dependencies.
plot_df <- cox_results %>%
  mutate(label = gsub("`", "", term), label = gsub("_group", "", label), label = gsub("High", "High vs Low", label)) %>%
  arrange(HR) %>%
  mutate(label = factor(label, levels = label))

p <- ggplot(plot_df, aes(x = HR, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), height = 0.2) +
  geom_point(size = 2.5) +
  scale_x_log10() +
  theme_classic(base_size = 12) +
  labs(title = "Univariate Cox regression", x = "Hazard ratio (log scale)", y = NULL)

ggsave(file.path(args$output_dir, "Figure5_univariate_cox_forest_GSE17536.pdf"), p, width = 7, height = max(4, 0.35 * nrow(plot_df) + 1.5))

message("Done: Cox results written to ", args$output_dir)
