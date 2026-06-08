library(ggplot2)
library(data.table)
library(cowplot)
library(ggrastr)

if (!dir.exists("result/fig4")) dir.create("result/fig4", recursive = TRUE)

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

res_path <- list.files(
  "Ez-LOF-Vsx1/result", full.names = TRUE
)

res_tbl <- lapply(res_path, read.csv, row.names = 1) |>
  do.call(rbind, args = _)

res_tbl$theta <- atan2(res_tbl$y_std, res_tbl$x_std) * 180 / pi
res_tbl$r <- sqrt(res_tbl$x_std^2 + res_tbl$y_std^2)
res_tbl <- as.data.table(res_tbl)
res_tbl[, type := sub("(.*)_[0-9]$", "\\1", sample)]
res_tbl$type <- factor(
  res_tbl$type, levels = c("mChi", "Ezi"),
  labels = c("mCherry RNAi", "E(z) RNAi")
)
res_tbl[, ectopic := x_std < 0]

# Statistics
ttbl <- res_tbl[, .(N = .N, ect_N = sum(ectopic)), by = c("sample", "type")][
  , type_n := paste0(type, " (n = ", .N, ")"), by = type
][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
]
stat_tbl <- list()

## Posterior abundance
tobj <- glm(ect_N / N ~ type, data = ttbl, family = quasibinomial, weights = N)

tsobj <- summary(tobj)
beta0 <- tsobj$coefficients[1, 1]
beta1 <- tsobj$coefficients[2, 1]
phi <- tsobj$dispersion
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))

stat_tbl[["posterior_ectopic"]] <- data.table(
  exp = "Numeber of posterior Vsx1 (E(z) RNAi, neuronal)",
  pvalue = tsobj$coefficients[2, 4],
  p0 = p0,
  p1 = p1,
  nctrl = 3,
  nexp = 4,
  phi = phi
)

## Temporal bias
### Whole OL
tobj <- lm(r ~ type + sample, data = res_tbl)
tsobj <- summary(tobj)
stat_tbl[["all_radial"]] <- data.table(
  exp = "Normalized distance of all Vsx1 (E(z) RNAi, neuronal)",
  pvalue = tsobj$coefficients[2, 4],
  p0 = tsobj$coefficients[1, 1],
  p1 = tsobj$coefficients[2, 1] + tsobj$coefficients[1, 1],
  nctrl = 3,
  nexp = 4,
  phi = NA
)

tobj <- lm(r ~ type + sample, data = res_tbl[ectopic == TRUE])
tsobj <- summary(tobj)
stat_tbl[["ectopic_radial"]] <- data.table(
  exp = "Normalized distance of ectopic Vsx1 (E(z) RNAi, neuronal)",
  pvalue = tsobj$coefficients[2, 4],
  p0 = tsobj$coefficients[1, 1],
  p1 = tsobj$coefficients[2, 1] + tsobj$coefficients[1, 1],
  nctrl = 3,
  nexp = 4,
  phi = NA
)

## Visualization
res_tbl[type == "mCherry RNAi"] |>
  ggplot(aes(x = x_std, y = y_std)) +
  geom_point(color = "#435274", size = 0.5) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave("./result/fig4/Ez-LOF-Vsx1-Ctrl.pdf", width = 4, height = 4)

res_tbl[type == "E(z) RNAi"] |>
  ggplot(aes(x = x_std, y = y_std)) +
  geom_point(color = "#ba3c3c", size = 0.5) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave("./result/fig4/Ez-LOF-Vsx1-KD.pdf", width = 4, height = 4)

ttbl |>
  ggplot(aes(x = type_n, y = ect_N / N, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  annotate(
    "segment", x = 1, xend = 2, y = 0.16, yend = 0.16, color = "black"
  ) +
  annotate(
    "text", x = 1.5, y = 0.17, label = "*", color = "black", size = 6
  ) +
  labs(y = "# Posterior Vsx1 / Total (All Vsx1+)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")
ggsave("./result/fig4/Ez-LOF-Vsx1-quant.pdf", width = 3, height = 4)

## Ez LOF Bi anterior abundance
res_path <- list.files(
  "Ez-LOF-Bi/result", full.names = TRUE, pattern = "preprocessed.csv$"
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(
  res_tbl$type,
  levels = c("mCherry RNAi", "Ezi RNAi"),
  labels = c("mCherry RNAi", "E(z) RNAi")
)
res_tbl[, anterior := x_std > 0]

res_tbl |>
  ggplot(aes(x = x_std, y = y_std, color = type)) +
  geom_point_rast(size = 0.5, alpha = 0.12, raster.dpi = 300) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  guides(color = "none")
ggsave("./result/fig4/Ez-LOF-Bi-scatter.pdf", width = 4, height = 4)

p <- res_tbl |>
  ggplot(aes(x = x_std, y = y_std, color = type)) +
  geom_point(size = 0.5, alpha = 0.12) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(color = "Condition")

ggsave("./result/fig4/Ez-LOF-Bi-legend.pdf", get_legend(p))

ttbl <- res_tbl[
  , .(total = .N, anterior = sum(anterior)), by = c("sample", "type")
][
  , type_n := paste0(type, " (n = ", .N, ")"), by = type
][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
]

tobj <- glm(
  anterior / total ~ type, data = ttbl, family = quasibinomial,
  weights = total
)
tsobj <- summary(tobj)
beta0 <- tsobj$coefficients[1, 1]
beta1 <- tsobj$coefficients[2, 1]
phi <- tsobj$dispersion
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))
n <- table(ttbl$type)

stat_tbl[["bi_anterior"]] <- data.table(
  exp = "Number of anterior Bi cells (E(z) RNAi, neuronal)",
  pvalue = tsobj$coefficients[2, 4],
  p0 = p0,
  p1 = p1,
  nctrl = unname(n["mCherry RNAi"]),
  nexp = unname(n["E(z) RNAi"]),
  phi = phi
)

sig_label <- p_to_sig(tsobj$coefficients[2, 4])

ttbl |>
  ggplot(aes(x = type_n, y = anterior / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  sig_annotation(
    ttbl$anterior / ttbl$total, ttbl$type_n, sig_label, label_offset = 0.16
  ) +
  labs(y = "# Anterior Bi cells\n/ Total (All Bi+ cells)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")
ggsave("./result/fig4/Ez-LOF-Bi-quant.pdf", width = 3, height = 4)

## Ez LOF Vvl radial variance
res_path <- list.files(
  "Ez-LOF-Vvl/result", full.names = TRUE, pattern = "preprocessed.csv$"
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(
  res_tbl$type,
  levels = c("mCherry RNAi", "E(z) RNAi")
)

res_tbl |>
  ggplot(aes(x = x_std, y = y_std, color = type)) +
  geom_point_rast(size = 0.5, alpha = 0.12, raster.dpi = 300) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  guides(color = "none")
ggsave("./result/fig4/Ez-LOF-Vvl-scatter.pdf", width = 4, height = 4)

ttbl <- res_tbl[
  , .(var_r = var(r)), by = c("sample", "type")
][
  , type_n := paste0(type, " (n = ", .N, ")"), by = type
][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
]

tobj <- t.test(var_r ~ type, data = ttbl)
n <- table(ttbl$type)

stat_tbl[["vvl_radial_variance"]] <- data.table(
  exp = "Variance of normalized distance of Vvl (E(z) RNAi, neuronal)",
  pvalue = tobj$p.value,
  p0 = ttbl[type == "mCherry RNAi", mean(var_r)],
  p1 = ttbl[type == "E(z) RNAi", mean(var_r)],
  nctrl = unname(n["mCherry RNAi"]),
  nexp = unname(n["E(z) RNAi"]),
  phi = NA
)

sig_label <- p_to_sig(tobj$p.value)

ttbl |>
  ggplot(aes(x = type_n, y = var_r, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se", pch = "-", size = 2
  ) +
  sig_annotation(
    ttbl$var_r, ttbl$type_n, sig_label, label_offset = 0.3, bar_offset = 0.15
  ) +
  labs(y = "Variance on the radial axis") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  guides(color = "none")
ggsave("./result/fig4/Ez-LOF-Vvl-variance.pdf", width = 3, height = 4)

## Ez LOF Toy radial variance
res_path <- list.files(
  "Ez-LOF-Toy/result", full.names = TRUE, pattern = "preprocessed.csv$"
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(
  res_tbl$type,
  levels = c("mCherry RNAi", "E(z) RNAi")
)

res_tbl |>
  ggplot(aes(x = x_std, y = y_std, color = type)) +
  geom_point_rast(size = 0.5, alpha = 0.12, raster.dpi = 300) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  guides(color = "none")
ggsave("./result/fig4/Ez-LOF-Toy-scatter.pdf", width = 4, height = 4)

ttbl <- res_tbl[
  , .(var_r = var(r)), by = c("sample", "type")
][
  , type_n := paste0(type, " (n = ", .N, ")"), by = type
][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
]

tobj <- t.test(var_r ~ type, data = ttbl)
n <- table(ttbl$type)

stat_tbl[["toy_radial_variance"]] <- data.table(
  exp = "Variance of normalized distance of Toy (E(z) RNAi, neuronal)",
  pvalue = tobj$p.value,
  p0 = ttbl[type == "mCherry RNAi", mean(var_r)],
  p1 = ttbl[type == "E(z) RNAi", mean(var_r)],
  nctrl = unname(n["mCherry RNAi"]),
  nexp = unname(n["E(z) RNAi"]),
  phi = NA
)

sig_label <- p_to_sig(tobj$p.value)

ttbl |>
  ggplot(aes(x = type_n, y = var_r, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se", pch = "-", size = 2
  ) +
  sig_annotation(
    ttbl$var_r, ttbl$type_n, sig_label, label_offset = 0.3, bar_offset = 0.15
  ) +
  labs(y = "Variance on the radial axis") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  guides(color = "none")
ggsave("./result/fig4/Ez-LOF-Toy-variance.pdf", width = 3, height = 4)

write.csv(rbindlist(stat_tbl), "result/fig4/stat.csv", row.names = FALSE)
