#!/usr/bin/env Rscript

# ============================================================
# Figure 1G | Reference Cohort Survival Analysis by Age Class
# File: figure1g_reference_survival_by_age_class.R
#
# Description:
#   Generates Kaplan–Meier survival curves for the external
#   reference cohort stratified by age class and exports a
#   supplementary survival summary table.
#
# Inputs:
#   - data/pediatric_aya_hgg_study_data.rds
#   - data/pediatric_aya_hgg_external_data.rds
#
# Outputs:
#   - Figure1G_reference_survival_assoc.pdf
#   - Both_supplementary_survival_table.csv
#
# Author: Nicole Tignor
# Affiliation: Icahn School of Medicine at Mount Sinai
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(RColorBrewer)
})

get_green <- function(n) {
  colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[-(1:2)])(n)
}

study.data <- readRDS("data/pediatric_aya_hgg_study_data.rds")
external.data <- readRDS("data/pediatric_aya_hgg_external_data.rds")

clinical.data <- study.data$clinical
ref.data <- external.data$reference_cohort$ref.data
ref.cdata <- external.data$reference_cohort$ref.cdata

ref.data$has.complete.mut <- ref.data$id %in% ref.cdata$id

vali.surv.data <- ref.data[, c("id", "SurvObj", "age_class_name", "First.Diagnosis", "Gender")]
vali.surv.data <- data.frame(vali.surv.data, as.matrix(vali.surv.data$SurvObj))

get_data00 <- function(
    id = clinical.data$id,
    os = clinical.data$os,
    os_status = clinical.data$os_status,
    max.day = 2000
) {
  os_status <- ifelse(os > max.day, 0, os_status)
  os <- ifelse(os > max.day, max.day, os)
  
  data0 <- data.frame(
    id = id,
    days = os,
    os.status = os_status
  )
  
  data0$SurvObj <- with(data0, Surv(days, os.status))
  data0$os.status[is.na(data0$days)] <- NA
  
  data0
}

surv.data <- vali.surv.data
surv.data$SurvObj <- with(surv.data, Surv(time, status))

get_surv.cl <- function(
    cl = NULL,
    mytitle = "Age class",
    data = surv.data,
    scheme = get_green,
    mypal0 = NULL,
    cl.name = "age_class_name"
) {
  data <- data[!is.na(data$time), ]
  data <- data[data$First.Diagnosis %in% "Yes", ]
  data$cl <- data[, cl.name]
  
  fit <- survfit(SurvObj ~ cl, data = data)
  fit.table <- as.data.frame(summary(fit)$table, check.names = FALSE)
  
  surv_table <- data.frame(
    Group = row.names(fit.table),
    n = fit.table[, "records"],
    RMST = round(fit.table[, "rmean"], 1),
    se = round(fit.table[, "se(rmean)"], 2),
    Events = fit.table[, "events"],
    Censored = fit.table[, "n.start"] - fit.table[, "events"]
  )
  
  log_rank_test <- survdiff(Surv(time, status) ~ cl, data = data)
  log_p_value <- pchisq(
    log_rank_test$chisq,
    length(log_rank_test$n) - 1,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  surv_table <- surv_table %>%
    dplyr::mutate(
      log_p_value = log_p_value,
      RMST_SE = paste0(RMST, " (", se, ")")
    ) %>%
    dplyr::select(Group, n, RMST_SE, Events, Censored, log_p_value)
  
  print(surv_table)
  
  write.csv(
    surv_table,
    paste0(mytitle, "_supplementary_survival_table.csv"),
    row.names = FALSE
  )
  
  if (!is.null(scheme)) {
    mypal0 <- scheme(length(levels(data$cl)))
    names(mypal0) <- levels(data$cl)
  }
  
  mypal <- mypal0[match(levels(data$cl), names(mypal0))]
  
  splot <- ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    palette = as.character(mypal),
    conf.int = FALSE,
    risk.table = TRUE,
    size = 0.5,
    title = mytitle
  )
  
  splot
}

surv.data$Gender <- factor(surv.data$Gender, levels = c("Male", "Female"))
surv.data <- surv.data[!surv.data$id %in% clinical.data$id, ]

splot <- get_surv.cl(mytitle = "Both")

aligned <- cowplot::align_plots(
  splot$plot + theme(plot.margin = margin(r = 10)),
  splot$table + theme(plot.margin = margin(r = 10)),
  align = "v"
)

pdf("Figure1G_reference_survival_assoc.pdf", width = 6, height = 7)
cowplot::plot_grid(aligned[[1]], aligned[[2]], ncol = 1, rel_heights = c(4, 1.5))
dev.off()