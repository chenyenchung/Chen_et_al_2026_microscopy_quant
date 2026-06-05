library(data.table)
library(ggplot2)
library(cowplot)
library(patchwork)

if (!dir.exists("result/fig3")) dir.create("result/fig3")
stat_tbl <- list()

quasibinom_stat <- function(
  formula, ttbl, ctrl, exp,
  group_col = "type", weights_col = "total"
) {
  tobj <- do.call(glm, list(
    formula = formula, data = ttbl, family = quasibinomial,
    weights = ttbl[[weights_col]]
  ))
  tsobj <- summary(tobj)
  beta0 <- tsobj$coefficients[1, 1]
  beta1 <- tsobj$coefficients[2, 1]
  p0 <- 1 / (1 + exp(-beta0))
  p1 <- 1 / (1 + exp(-beta0 - beta1))
  n <- table(ttbl[[group_col]])
  data.frame(
    p_value = tsobj$coefficients[2, 4],
    n_ctrl  = unname(n[ctrl]),
    n_exp   = unname(n[exp]),
    fc      = p1 / p0,
    p0      = p0,
    p1      = p1,
    phi     = tsobj$dispersion
  )
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

## Vsx1 OE Pm4/Mi4
res_path <- list.files(
  "./Vsx1-GOF-Mi4/result", full.names = TRUE
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(res_tbl$type)
res_tbl$type <- relevel(res_tbl$type, "LacZ OE")
res_tbl[
  , cell := ifelse(int_1 > int_threshold, "Mi4 (Ct/Run)", "Pm4 (Run)")
]
res_tbl$cell <- factor(res_tbl$cell, levels = c("Pm4 (Run)", "Mi4 (Ct/Run)"))

res_tbl[in_plane == TRUE & type == "LacZ OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(color = "none")
ggsave("./result/fig3/Vsx1-GOF-Mi4-Ctrl.pdf", width = 4, height = 4)

res_tbl[in_plane == TRUE & type == "Vsx1 OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(color = "none")
ggsave("./result/fig3/Vsx1-GOF-Mi4-OE.pdf", width = 4, height = 4)

p <- res_tbl[in_plane == TRUE & type == "LacZ OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell type")

ggsave("./result/fig3/Vsx1-GOF-Mi4-legend.pdf", get_legend(p))

ttbl <- res_tbl[
  x_std > 0, .(total = .N, Pm4 = sum(cell == "Pm4 (Run)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"]

stat_tbl[["Vsx1-GOF-Mi4"]] <- quasibinom_stat(
  Pm4 / total ~ type, ttbl, "LacZ OE", "Vsx1 OE"
)
sig_label <- p_to_sig(stat_tbl[["Vsx1-GOF-Mi4"]]$p_value)

ttbl |>
  ggplot(aes(x = type_n, y = Pm4 / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  sig_annotation(ttbl$Pm4 / ttbl$total, ttbl$type_n, sig_label) +
  labs(y = "# Anterior Pm4 (Run)\n/ Total (All Anterior Run+)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")

ggsave("./result/fig3/Vsx1-GOF-Mi4-quant.pdf", width = 1.5, height = 3.55)

## Vsx1/2 KD Pm4/Mi4
res_path <- list.files(
  "./Vsx1-LOF-Mi4/result", full.names = TRUE
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(res_tbl$type)
res_tbl$type <- relevel(res_tbl$type, "mCherry RNAi")
res_tbl[
  , cell := ifelse(int_0 > int_threshold, "Mi4 (Ct/Run)", "Pm4 (Run)")
]
res_tbl$cell <- factor(res_tbl$cell, levels = c("Pm4 (Run)", "Mi4 (Ct/Run)"))

res_tbl[in_plane == TRUE & type == "mCherry RNAi"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Vsx1-LOF-Mi4-Ctrl.pdf", width = 4, height = 4)

res_tbl[in_plane == TRUE & type == "Vsx1/2 RNAi"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Vsx1-LOF-Mi4-KD.pdf", width = 4, height = 4)

p <- res_tbl[in_plane == TRUE & type == "mCherry RNAi"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell type")

ggsave("./result/fig3/Vsx1-LOF-Mi4-legend.pdf", get_legend(p))

ttbl <- res_tbl[
  x_std > 0, .(total = .N, Pm4 = sum(cell == "Pm4 (Run)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"]

stat_tbl[["Vsx1-LOF-Mi4"]] <- quasibinom_stat(
  Pm4 / total ~ type, ttbl, "mCherry RNAi", "Vsx1/2 RNAi"
)
sig_label <- p_to_sig(stat_tbl[["Vsx1-LOF-Mi4"]]$p_value)

ttbl |>
  ggplot(aes(x = type_n, y = Pm4 / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  sig_annotation(ttbl$Pm4 / ttbl$total, ttbl$type_n, sig_label) +
  labs(y = "# Anterior Pm4 (Run)\n/ Total (All Anterior Run+)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")

ggsave("./result/fig3/Vsx1-LOF-Mi4-quant.pdf", width = 1.5, height = 3.55)

## Vsx1 OE Tm5e/Tm29/Tm33/TmY5a
res_path <- list.files(
  "./Vsx1-GOF-Tm5e/result", full.names = TRUE
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl[type == "Vsx1 GOF", type := "Vsx1 OE"]
res_tbl$type <- factor(res_tbl$type)
res_tbl$type <- relevel(res_tbl$type, "LacZ OE")
res_tbl[
  , cell := ifelse(
    int_0 > int_threshold, "Tm29/33/TmY5a (Ey/Kn/D)", "Tm5e (Ey/Kn)"
  )
]
res_tbl$cell <- factor(
  res_tbl$cell, levels = c("Tm5e (Ey/Kn)", "Tm29/33/TmY5a (Ey/Kn/D)")
)

res_tbl[in_plane == TRUE & type == "LacZ OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Vsx1-GOF-Tm5e-Ctrl.pdf", width = 4, height = 4)

res_tbl[in_plane == TRUE & type == "Vsx1 OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Vsx1-GOF-Tm5e-Exp.pdf", width = 4, height = 4)

p <- res_tbl[in_plane == TRUE & type == "Vsx1 OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell type")

ggsave("./result/fig3/Vsx1-GOF-Tm5e-legend.pdf", get_legend(p))

ttbl <- res_tbl[
  in_plane == TRUE,
  .(total = .N, Tm5e = sum(cell == "Tm5e (Ey/Kn)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"]

stat_tbl[["Vsx1-GOF-Tm5e"]] <- quasibinom_stat(
  Tm5e / total ~ type, ttbl, "LacZ OE", "Vsx1 OE"
)
sig_label <- p_to_sig(stat_tbl[["Vsx1-GOF-Tm5e"]]$p_value)

ttbl |>
  ggplot(aes(x = type_n, y = Tm5e / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  sig_annotation(ttbl$Tm5e / ttbl$total, ttbl$type_n, sig_label) +
  labs(y = "# Tm5e (Ey/Kn) /\nTotal (All Ey/Kn+)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")

ggsave("./result/fig3/Vsx1-GOF-Tm5e-quant.pdf", width = 1.5, height = 3.55)

## Vsx1/2 KD Tm5e/Tm29/Tm33/TmY5a
res_path <- list.files(
  "./Vsx1-LOF-Tm5e/result/", full.names = TRUE
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(res_tbl$type)
res_tbl$type <- relevel(res_tbl$type, "mCherry RNAi")
res_tbl[
  , cell := ifelse(int_0 > 800, "Tm29/33/TmY5a (Ey/Kn/D)", "Tm5e (Ey/Kn)")
]
res_tbl$cell <- factor(res_tbl$cell, levels = c("Tm5e (Ey/Kn)", "Tm29/33/TmY5a (Ey/Kn/D)"))

res_tbl[in_plane == TRUE & type == "mCherry RNAi"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(color = "none")
ggsave("./result/fig3/Vsx1-LOF-Tm5e-Ctrl.pdf", width = 4, height = 4)

res_tbl[in_plane == TRUE & type == "Vsx1/2 RNAi"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(color = "none")
ggsave("./result/fig3/Vsx1-LOF-Tm5e-KD.pdf", width = 4, height = 4)

p <- res_tbl[in_plane == TRUE & type == "mCherry RNAi"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(
    values = c("#1B9E77", "#7570B3")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell type")

ggsave("./result/fig3/Vsx1-LOF-Tm5e-legend.pdf", get_legend(p))

ttbl <- res_tbl[
  in_plane == TRUE & x_std > 0,
  .(total = .N, Tm5e = sum(cell == "Tm5e (Ey/Kn)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"]

stat_tbl[["Vsx1-LOF-Tm5e"]] <- quasibinom_stat(
  Tm5e / total ~ type, ttbl, "mCherry RNAi", "Vsx1/2 RNAi"
)
sig_label <- p_to_sig(stat_tbl[["Vsx1-LOF-Tm5e"]]$p_value)

ttbl |>
  ggplot(aes(x = type_n, y = Tm5e / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  sig_annotation(ttbl$Tm5e / ttbl$total, ttbl$type_n, sig_label) +
  labs(y = "# Tm5e (Ey/Kn)\n/ Total (All Ey/Kn+)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")

ggsave("./result/fig3/Vsx1-LOF-Tm5e-quant.pdf", width = 1.5, height = 3.55)

## Bi OE TE
res_path <- list.files(
  "./Bi-GOF-TE/result", full.names = TRUE
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(res_tbl$type)
res_tbl$type <- relevel(res_tbl$type, "LacZ OE")
res_tbl[, cell := "Mi1 (Bsh+/Dimm-)"]
res_tbl[int_0 > int_threshold, cell := "TE (Bsh+/Dimm+)"]
res_tbl$cell <- factor(
  res_tbl$cell, levels = c("Mi1 (Bsh+/Dimm-)", "TE (Bsh+/Dimm+)")
)

res_tbl[in_plane == TRUE & type == "LacZ OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_color_manual(
    values = c("#7570B3", "#D95F02")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Bi-GOF-TE-Ctrl.pdf", width = 4, height = 4)

res_tbl[in_plane == TRUE & type == "Bi OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_color_manual(
    values = c("#7570B3", "#D95F02")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Bi-GOF-TE-Exp.pdf", width = 4, height = 4)

p <- res_tbl[in_plane == TRUE & type == "LacZ OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(
    values = c("#7570B3", "#D95F02")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell type")

ggsave("./result/fig3/Bi-GOF-TE-legend.pdf", get_legend(p))

ttbl <- res_tbl[
  in_plane == TRUE, .(total = .N, TE = sum(cell == "TE (Bsh+/Dimm+)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
  ]

stat_tbl[["Bi-GOF-TE (Whole OL)"]] <- quasibinom_stat(
  TE / total ~ type, ttbl, "LacZ OE", "Bi OE"
)
sig_label <- p_to_sig(stat_tbl[["Bi-GOF-TE (Whole OL)"]]$p_value)

ttbl |>
  ggplot(aes(x = type_n, y = TE / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  sig_annotation(ttbl$TE / ttbl$total, ttbl$type_n, sig_label) +
  labs(y = "# TE (Bsh+/Dimm+) /\nTotal (All in-plane Bsh+)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")

ggsave("./result/fig3/Bi-GOF-TE-quant-sup.pdf", width = 1.5, height = 3.55)

ttbl <- res_tbl[
  in_plane == TRUE & x_std > 0,
  .(total = .N, ant_TE = sum(cell == "TE (Bsh+/Dimm+)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
]

stat_tbl[["Bi-GOF-TE (Anterior)"]] <- quasibinom_stat(
  ant_TE / total ~ type, ttbl, "LacZ OE", "Bi OE"
)
sig_label <- p_to_sig(stat_tbl[["Bi-GOF-TE (Anterior)"]]$p_value)

ttbl |>
  ggplot(aes(x = type_n, y = ant_TE / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  sig_annotation(ttbl$ant_TE / ttbl$total, ttbl$type_n, sig_label) +
  labs(y = "# Ectopic TE (Bsh+/Dimm+) /\nTotal (Anterior in-plane Bsh+)") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA), labels = scales::number_format()) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")

ggsave("./result/fig3/Bi-GOF-TE-quant.pdf", width = 1.5, height = 3.55)

## Bi KO TE
data_tbl <- fread("./Bi-LOF-TE/result/all_bsh_cells.csv")

data_tbl$type <- factor(
  data_tbl$condition, levels = c("CantonS", "sgBi"),
  labels = c("Control", "Bi sKO")
)
data_tbl[, cell_class := "Mi1 (Bsh+/Dimm-)"]
data_tbl[int_2 > 500, cell_class := "TE (Bsh+/Dimm+)"]
data_tbl$cell_class <- factor(
  data_tbl$cell_class, levels = c("Mi1 (Bsh+/Dimm-)", "TE (Bsh+/Dimm+)")
)

bi_lof_te_density_tbl <- data_tbl[
  type == "Control" & cell_class == "TE (Bsh+/Dimm+)" & in_plane == TRUE
]

data_tbl[type == "Control" & int_0 > 1000 & in_plane == TRUE] |>
ggplot(aes(pca_x, pca_y, color = cell_class)) +
  geom_point(size = 0.5) +
  geom_density_2d(
    data = bi_lof_te_density_tbl,
    aes(x = pca_x, y = pca_y),
    contour_var = "ndensity",
    breaks = 0.5,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_color_manual(
    values = c("#7570B3", "#D95F02")
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(color = "none")
ggsave("./result/fig3/Bi-LOF-TE-Ctrl.pdf", width = 4, height = 4)

data_tbl[type == "Bi sKO" & int_0 > 1000 & in_plane == TRUE] |>
ggplot(aes(pca_x, pca_y, color = cell_class)) +
  geom_point(size = 0.5) +
  geom_density_2d(
    data = bi_lof_te_density_tbl,
    aes(x = pca_x, y = pca_y),
    contour_var = "ndensity",
    breaks = 0.5,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  scale_color_manual(
    values = c("#7570B3", "#D95F02")
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(color = "none")
ggsave("./result/fig3/Bi-LOF-TE-sKO.pdf", width = 4, height = 4)

p <- data_tbl[type == "Control" & int_0 > 1000 & in_plane == TRUE] |>
ggplot(aes(pca_x, pca_y, color = cell_class)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(
    values = c("#7570B3", "#D95F02")
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell type")

ggsave("./result/fig3/Bi-LOF-TE-legend.pdf", get_legend(p))

ttbl <- fread("./Bi-LOF-TE/result/te_ratios.csv")
n <- table(ttbl$condition)
ttbl$type <- factor(
  ttbl$condition, levels = c("CantonS", "sgBi"),
  labels = c(
    paste0("Control (n = ", n["CantonS"], ")"),
    paste0("Bi sKO (n = ", n["sgBi"], ")")
  )
)

ttbl[, total := False + True]
stat_tbl[["Bi-LOF-TE"]] <- quasibinom_stat(
  True / total ~ condition, ttbl, "CantonS", "sgBi", group_col = "condition"
)
sig_label <- p_to_sig(stat_tbl[["Bi-LOF-TE"]]$p_value)

ttbl |>
  ggplot(aes(x = type, y = ratio, color = type)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  sig_annotation(ttbl$ratio, ttbl$type, sig_label) +
  labs(y = "# In-clone&ROI TE (Bsh+/Dimm+)\n/ In-clone&ROI Bsh+") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")

ggsave("./result/fig3/Bi-LOF-TE-quant.pdf", width = 1.5, height = 3.55)

## This must be run at the end of the script!!
stat_df <- do.call(rbind, stat_tbl)
write.csv(stat_df, "./result/fig3/stat_summary.csv")
