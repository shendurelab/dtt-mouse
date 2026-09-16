
#########################################################################################################
### Figure S1. Characterization of epegRNA::tape-BC::circtape and PEmax integrations across E13.5 embryos

###########################################################################################
### Fig. S1C: Tape-BC read-count rank-abundance, shown separately for embryos #2, #3 and #6

library(ggplot2)
library(dplyr)

data_path = "/Users/cxqiu/GitHub/mouse_sprint/tape_pipeline/tables"
save_path = "/Volumes/f0085ts/work/tapemouse/making_figures"

dat = read.csv(paste0(data_path, "/suppfig1a_barcode_rank_abundance.csv"))
dat = dat[dat$embryo %in% paste0("embryo ", c(2,3,6)),]
dat$log_reads = log(dat$reads)
dat$copies = factor(dat$copies, levels = names(table(dat$copies)))

dat_cutoff = read.csv(paste0(data_path, "/suppfig1a_cutoffs.csv"))
dat_cutoff = dat_cutoff[dat_cutoff$embryo %in% paste0("embryo ", c(2,3,6)),]
dat_cutoff$log_reads_at_freq_cutoff = log(dat_cutoff$reads_at_freq_cutoff)

p = dat %>%
    ggplot(aes(rank, log_reads, color = copies)) +
    geom_point(size = 0.6) +
    geom_hline(data = dat_cutoff, aes(yintercept = log_reads_at_freq_cutoff), linetype = "dashed") +
    facet_wrap(~ embryo, nrow = 1) +
    labs(x = "Barcode rank (ordered high to low)", y = "Log (Reads per barcode)", fill = NULL) +
    scale_y_continuous(breaks = seq(0, 20, by = 2)) +
    theme_classic()

ggsave(paste0(save_path, "/FigS1/FigS1_TAPE-BC_read-count_rank-abundance.pdf"), p, height = 3, width = 10)


#################################################################################
### Fig. S1D: PEmax copy number in each E13.5 embryo was assessed by duplex ddPCR

# =============================================================================
# Cas9 copy number from a Bio-Rad ddPCR CSV export
# - Computes CN = (Cas9 / GAPDH) x 2
# - Adds 95% Poisson confidence intervals, propagated through the ratio
# - Plots with ggplot2
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# --- 1. Input --------------------------------------------------------------

csv_path <- "/Volumes/f0085ts/work/tapemouse/figures_data/E13_emrby_1_through_10_CN_analysis.csv"   # <-- edit path

df <- read_csv(csv_path, show_col_types = FALSE,
               locale = locale(encoding = "UTF-8"))

# The concentration column has a µ in its name; rename by position/pattern
names(df)[grepl("^Conc", names(df))] <- "Conc"
df <- df %>% rename(Sample = `Sample description 1`)

df <- df[df$Sample != "NTC",]

# "No Call" -> NA, force numeric
df <- df %>% mutate(Conc = suppressWarnings(as.numeric(Conc)))

# --- 2. Per-target Poisson stats -------------------------------------------
# p       = k/N                     fraction positive droplets
# lambda  = -ln(1 - p)              Poisson-corrected copies/droplet
# var(l)  = (exp(l) - 1) / N        variance of lambda estimator
# var(ln l) = var(l) / l^2          used for log-space CI propagation

per_target <- df %>%
  mutate(
    p       = Positives / `Accepted Droplets`,
    lambda  = ifelse(Positives == 0, 0, -log(1 - p)),
    var_lam = ifelse(Positives == 0, NA_real_,
                     (exp(lambda) - 1) / `Accepted Droplets`),
    var_ln  = ifelse(lambda == 0 | is.na(lambda), NA_real_,
                     var_lam / lambda^2)
  ) %>%
  select(Sample, Target, Positives, `Accepted Droplets`, lambda, var_ln)

# --- 3. Pivot to one row per sample and compute CN + CI --------------------

wide <- per_target %>%
  pivot_wider(
    id_cols     = Sample,
    names_from  = Target,
    values_from = c(Positives, lambda, var_ln)
  )

results <- wide %>%
  mutate(
    CN = case_when(
      is.na(lambda_GAPDH)     ~ NA_real_,           # no reference (NTC)
      Positives_Cas9 == 0     ~ 0,                  # true zero (E4)
      TRUE                    ~ 2 * lambda_Cas9 / lambda_GAPDH
    ),
    se_ln = sqrt(var_ln_Cas9 + var_ln_GAPDH),
    CN_lo = CN * exp(-1.96 * se_ln),
    CN_hi = CN * exp( 1.96 * se_ln)
  ) %>%
  mutate(Sample = factor(Sample, levels = unique(df$Sample)))

# --- 4. Categorize for coloring & make display labels ----------------------

results <- results %>%
  mutate(
    category = case_when(
      Sample == "NTC"          ~ "No template control",
      is.na(CN) | CN == 0      ~ "No call / zero",
      TRUE                     ~ "Measured"
    ),
    CN_plot    = ifelse(is.na(CN),    0, CN),
    CN_lo_plot = ifelse(is.na(CN_lo), 0, CN_lo),
    CN_hi_plot = ifelse(is.na(CN_hi), 0, CN_hi),
    label = case_when(
      Sample == "NTC"          ~ "N/A\n(no GAPDH)",
      is.na(CN) | CN == 0      ~ "No Call\n(0 pos. droplets)",
      TRUE                     ~ sprintf("%.2f", CN)
    )
  )

print(results %>% select(Sample, CN, CN_lo, CN_hi, category))

# --- 5. Plot ---------------------------------------------------------------

palette <- c(
  "Measured"            = "#4472C4",
  "No call / zero"      = "#E07A7A",
  "No template control" = "#CCCCCC"
)

p <- ggplot(results, aes(x = Sample, y = CN_plot, fill = category)) +
  geom_col(color = "black", width = 0.65, linewidth = 0.3) +
  geom_errorbar(aes(ymin = CN_lo_plot, ymax = CN_hi_plot),
                width = 0.2, linewidth = 0.5) +
  # Numeric labels for measured samples (above the upper CI)
  geom_text(data = filter(results, category == "Measured"),
            aes(y = CN_hi_plot + 0.3, label = label),
            size = 3.5, fontface = "bold") +
  # Italic labels for control / no-call samples (near zero)
  geom_text(data = filter(results, category != "Measured"),
            aes(y = 0.15, label = label),
            size = 3, fontface = "italic", color = "grey30") +
  scale_fill_manual(values = palette, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Sample",
    y = "Cas9 copy number (Cas9 / GAPDH x 2)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position   = c(0.87, 0.85),
    legend.background = element_rect(fill = NA, color = NA),
    plot.title        = element_text(hjust = 0.5),
    plot.subtitle     = element_text(hjust = 0.5, size = 10, color = "grey30")
  )

print(p)

# --- 6. Save ---------------------------------------------------------------

ggsave("/Volumes/f0085ts/work/tapemouse/making_figures/FigS1_copy_number_plot.pdf", p, width = 7, height = 5)


