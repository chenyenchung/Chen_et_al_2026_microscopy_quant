library(ggplot2)
library(data.table)
library(ggrastr)

if (!dir.exists("result/fig5")) {
  dir.create("result/fig5", recursive = TRUE)
}

p_to_sig <- function(p_value) {
  if (!is.finite(p_value)) return("n.s.")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  "n.s."
}

sig_annotation <- function(
  y, group, label, x = c(1, 2), bar_offset = 0.05, label_offset = 0.08,
  color = "black", text_size = 6
) {
  valid <- is.finite(y) & !is.na(group)
  finite_y <- y[valid]
  finite_group <- group[valid]
  if (length(finite_y) == 0) {
    stop("sig_annotation() needs at least one finite y value")
  }
  data_max <- max(finite_y)
  se_upper <- max(
    vapply(
      split(finite_y, finite_group), function(group_y) mean_se(group_y)$ymax,
      numeric(1)
    ),
    na.rm = TRUE
  )
  annotation_base <- max(data_max, se_upper, na.rm = TRUE)
  y_ptp <- diff(range(finite_y))
  if (y_ptp == 0) {
    y_ptp <- max(abs(annotation_base), 1)
  }

  list(
    annotate(
      "segment", x = x[1], xend = x[2],
      y = annotation_base + bar_offset * y_ptp,
      yend = annotation_base + bar_offset * y_ptp,
      color = color
    ),
    annotate(
      "text", x = mean(x), y = annotation_base + label_offset * y_ptp,
      label = label, color = color, size = text_size
    )
  )
}

# Load Su(H) MARCM
marcm <- fread("Su(H)/manual_quantification.csv")
setnames(
  marcm, old = c("Condition", "Vsx2+/Hth+/GFP+", "Hth+/GFP+"),
  new = c("Condition", "Vsx2", "Hth")
)

tobj <- glm(
  Vsx2 / Hth ~ Condition, data = marcm,
  family = quasibinomial, weight = Hth
)
tsobj <- summary(tobj)
beta0 <- tsobj$coefficients[1, 1]
beta1 <- tsobj$coefficients[2, 1]
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))
phi <- tsobj$dispersion
sig_label <- p_to_sig(tsobj$coefficients[2, 4])
fwrite(
  data.table(
    p_value = tsobj$coefficients[2, 4],
    p0 = p0,
    p1 = p1,
    phi = phi
  ),
  "result/fig5/Su(H)_marcm_stats.csv"
)

ssize <- table(marcm$Condition)
marcm$Condition <- factor(
  marcm$Condition,
  levels = c("WT MARCM", "Su(H) MARCM"),
  labels = c(
    paste0("WT Clone", " (n = ", ssize["WT MARCM"], ")"),
    paste0("Su(H) Clone", " (n = ", ssize["Su(H) MARCM"], ")")
  )
)

set.seed(325)
marcm |>
  ggplot(aes(x = Condition, y = Vsx2 / Hth, color = Condition)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  sig_annotation(marcm$Vsx2 / marcm$Hth, marcm$Condition, sig_label) +
  labs(y = "# Vsx2+/Hth+ in clone /\n# Hth+ in clone") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 12)
  ) +
  scale_y_continuous(limits = c(0, NA), labels = scales::number_format()) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")
ggsave("./result/fig5/Su(H)-MARCM-Vsx2-quant.pdf", width = 3, height = 4)

## Vsx1 neuronal RNAi
res_path <- list.files(
  "Vsx1-NE-LOF/result", full.names = TRUE,
  pattern = "_preprocessed\\.csv$"
)
if (length(res_path) == 0) {
  stop("No Vsx1-NE-LOF preprocessed CSV files were found")
}

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

required_cols <- c("sample", "type", "bsh_num", "x_std", "y_std")
missing_cols <- setdiff(required_cols, names(res_tbl))
if (length(missing_cols) > 0) {
  stop(
    "Vsx1-NE-LOF data are missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

expected_types <- c("mCherry RNAi", "Vsx1 NE RNAi")
unexpected_types <- setdiff(unique(res_tbl$type), expected_types)
if (length(unexpected_types) > 0) {
  stop(
    "Unexpected Vsx1-NE-LOF condition(s): ",
    paste(unexpected_types, collapse = ", ")
  )
}
if (any(!is.finite(res_tbl$x_std)) || any(!is.finite(res_tbl$y_std))) {
  stop("Vsx1-NE-LOF standardized coordinates must be finite")
}

bsh_check <- res_tbl[
  , .(
    bsh_values = uniqueN(bsh_num),
    bsh_num = first(bsh_num)
  ),
  by = sample
]
if (
  any(bsh_check$bsh_values != 1) ||
  any(!is.finite(bsh_check$bsh_num)) ||
  any(bsh_check$bsh_num <= 0)
) {
  stop("Each Vsx1-NE-LOF sample must have one positive finite bsh_num")
}

res_tbl[, type := factor(
  type,
  levels = expected_types,
  labels = c("mCherry RNAi", "Vsx1 RNAi")
)]
setorder(res_tbl, type)

p <- res_tbl |>
  ggplot(aes(x = x_std, y = y_std, color = type)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point_rast(
    size = 2, alpha = 0.1, raster.dpi = 300
  ) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  scale_color_manual(values = c("#435274", "#ba3c3c")) +
  labs(color = "Condition") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
ggsave(
  "./result/fig5/Vsx1-NE-LOF-overlay.pdf", p,
  width = 4, height = 4
)

ttbl <- res_tbl[
  , .(
    bsh_num = first(bsh_num),
    total = .N,
    anterior = sum(x_std >= 0),
    posterior = sum(x_std < 0)
  ),
  by = c("sample", "type")
]
if (any(ttbl$anterior + ttbl$posterior != ttbl$total)) {
  stop("Anterior and posterior Vsx1 counts do not recover all rows")
}

ttbl[
  , `:=`(
    anterior_ratio = anterior / bsh_num,
    posterior_ratio = posterior / bsh_num
  )
]
n <- table(ttbl$type)
ttbl[, type_n := factor(
  paste0(type, " (n = ", n[as.character(type)], ")"),
  levels = paste0(levels(type), " (n = ", n[levels(type)], ")")
)]

regional_tests <- list(
  anterior = t.test(anterior_ratio ~ type, data = ttbl),
  posterior = t.test(posterior_ratio ~ type, data = ttbl)
)

stats_tbl <- rbindlist(lapply(names(regional_tests), function(region) {
  ratio_col <- paste0(region, "_ratio")
  test_obj <- regional_tests[[region]]
  data.table(
    region = region,
    test = test_obj$method,
    p_value = test_obj$p.value,
    control_mean = ttbl[type == "mCherry RNAi", mean(get(ratio_col))],
    rnai_mean = ttbl[type == "Vsx1 RNAi", mean(get(ratio_col))],
    n_control = unname(n["mCherry RNAi"]),
    n_rnai = unname(n["Vsx1 RNAi"])
  )
}))
fwrite(stats_tbl, "./result/fig5/Vsx1-NE-LOF_stats.csv")

regional_ratio_plot <- function(tbl, ratio_col, y_label, sig_label) {
  tbl |>
    ggplot(aes(x = type_n, y = .data[[ratio_col]], color = type_n)) +
    geom_jitter(width = 0.1, height = 0, size = 1) +
    stat_summary(fun.data = "mean_se", pch = "-", size = 2) +
    sig_annotation(
      tbl[[ratio_col]], tbl$type_n, sig_label, label_offset = 0.22, bar_offset = 0.1
    ) +
    labs(y = y_label) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(angle = 60, hjust = 1, size = 12)
    ) +
    scale_y_continuous(limits = c(0, NA), labels = scales::number_format()) +
    scale_color_manual(values = c("#435274", "#ba3c3c")) +
    guides(color = "none")
}

set.seed(325)
p <- regional_ratio_plot(
  ttbl, "anterior_ratio", "# Anterior Vsx1+ /\n# Whole-sample Bsh+",
  p_to_sig(regional_tests$anterior$p.value)
)
ggsave(
  "./result/fig5/Vsx1-NE-LOF-anterior-quant.pdf", p,
  width = 3, height = 4
)

p <- regional_ratio_plot(
  ttbl, "posterior_ratio", "# Posterior Vsx1+ /\n# Whole-sample Bsh+",
  p_to_sig(regional_tests$posterior$p.value)
)
ggsave(
  "./result/fig5/Vsx1-NE-LOF-posterior-quant.pdf", p,
  width = 3, height = 4
)
