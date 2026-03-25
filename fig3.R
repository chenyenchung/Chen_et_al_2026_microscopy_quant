library(data.table)
library(ggplot2)
library(cowplot)
library(patchwork)

if (!dir.exists("result/fig3")) dir.create("result/fig3")
stat_tbl <- list()

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
  x_std > 0, .(total = .N, Mi4 = sum(cell == "Mi4 (Ct/Run)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"]

ttbl |>
  ggplot(aes(x = type_n, y = Mi4 / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  annotate(
    "segment", x = 1, xend = 2, y = 0.55, yend = 0.55, color = "black"
  ) +
  annotate(
    "text", x = 1.5, y = 0.57, label = "*", color = "black", size = 6
  ) +
  labs(y = "# Anterior Mi4 (Ct/Run)\n/ Total (All Anterior Run+)") +
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

ggsave("./result/fig3/Vsx1-GOF-Mi4-quant.pdf", width = 1.25, height = 3.55)

tobj <- glm(
  Mi4 / total ~ type, data = ttbl, family = quasibinomial, weights = total
)
tsobj <- summary(tobj)
beta0 <- tsobj$coefficients[1, 1]
beta1 <- tsobj$coefficients[2, 1]
phi <- tsobj$dispersion
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))
stat_tbl[["Vsx1-GOF-Mi4"]] <- data.frame(
  p_value = tsobj$coefficients[2, 4],
  n_ctrl = 4L,
  n_exp = 4L,
  fc = p1 / p0,
  p0 = p0,
  p1 = p1,
  phi = phi
)

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
  x_std > 0, .(total = .N, Mi4 = sum(cell == "Mi4 (Ct/Run)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"]

ttbl |>
  ggplot(aes(x = type_n, y = Mi4 / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  annotate(
    "segment", x = 1, xend = 2, y = 0.7, yend = 0.7, color = "black"
  ) +
  annotate(
    "text", x = 1.5, y = 0.72, label = "*", color = "black", size = 6
  ) +
  labs(y = "# Anterior Mi4 (Ct/Run)\n/ Total (All Anterior Run+)") +
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

ggsave("./result/fig3/Vsx1-LOF-Mi4-quant.pdf", width = 1.25, height = 3.55)

tobj <- glm(
  Mi4 / total ~ type, data = ttbl, family = quasibinomial, weights = total
)
tsobj <- summary(tobj)
beta0 <- tsobj$coefficients[1, 1]
beta1 <- tsobj$coefficients[2, 1]
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))
phi <- tsobj$dispersion
stat_tbl[["Vsx1-LOF-Mi4"]] <- data.frame(
  p_value = tsobj$coefficients[2, 4],
  n_ctrl = 3L,
  n_exp = 3L,
  fc = p1 / p0,
  p0 = p0,
  p1 = p1,
  phi = phi
)

## Vsx1 OE Tm5e/Tm29/Tm33/TmY5a
res_path <- list.files(
  "./Vsx1-GOF-Tm5e/result", full.names = TRUE
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(res_tbl$type)
res_tbl$type <- relevel(res_tbl$type, "LacZ OE")
res_tbl[
  , cell := ifelse(
    int_0 > D_threshold, "Tm29/33/TmY5a (Ey/Kn/D)", "Tm5e (Ey/Kn)"
  )
]
res_tbl$cell <- factor(
  res_tbl$cell, levels = c("Tm5e (Ey/Kn)", "Tm29/33/TmY5a (Ey/Kn/D)")
)

res_tbl[with_ey_kn == TRUE & type == "LacZ OE"] |>
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

res_tbl[with_ey_kn == TRUE & type == "Vsx OE"] |>
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

p <- res_tbl[with_ey_kn == TRUE & type == "Vsx OE"] |>
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
  , .(total = .N, Tm29_et_al = sum(cell == "Tm29/33/TmY5a (Ey/Kn/D)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"]

ttbl |>
  ggplot(aes(x = type_n, y = Tm29_et_al / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  annotate(
    "segment", x = 1, xend = 2, y = 0.73, yend = 0.73, color = "black"
  ) +
  annotate(
    "text", x = 1.5, y = 0.75, label = "*", color = "black", size = 6
  ) +
  labs(y = "# Tm29/33/TmY5a (Ey/Kn/D) /\nTotal (All Ey/Kn+)") +
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

ggsave("./result/fig3/Vsx1-GOF-Tm5e-quant.pdf", width = 1.6, height = 3.55)

tobj <- glm(
  Tm29_et_al / total ~ type, data = ttbl,
  family = quasibinomial, weights = total
)
tsobj <- summary(tobj)
beta0 <- tsobj$coefficients[1, 1]
beta1 <- tsobj$coefficients[2, 1]
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))
phi <- tsobj$dispersion
stat_tbl[["Vsx1-GOF-Tm5e"]] <- data.frame(
  p_value = tsobj$coefficients[2, 4],
  n_ctrl = 3L,
  n_exp = 3L,
  fc = p1 / p0,
  p0 = p0,
  p1 = p1,
  phi = phi
)

## Bi OE DRA-Dm
res_path <- list.files(
  "./Bi-GOF-Dm/result", full.names = TRUE
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(res_tbl$type)
res_tbl$type <- relevel(res_tbl$type, "LacZ OE")
res_tbl[, dac := ifelse(int_2 > dac_threshold, "Dac+", "Dac-")]
res_tbl[, cell := "Other"]
res_tbl[dac == "Dac+" & y_std > 0, cell := "Dm-DRA (Dac+/Tj, dorsal)"]
res_tbl[dac == "Dac+" & y_std < 0, cell := "Dm8 (Dac+/Tj, ventral)"]
res_tbl$cell <- factor(res_tbl$cell, levels = c(
  "Dm-DRA (Dac+/Tj, dorsal)", "Dm8 (Dac+/Tj, ventral)", "Other"
))

res_tbl[type == "LacZ OE"] |>
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
    values = c("#D95F02", "#7570B3", "grey80")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Bi-GOF-Dm-Ctrl.pdf", width = 4, height = 4)

res_tbl[type == "Bi OE"] |>
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
    values = c("#D95F02", "#7570B3", "grey80")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Bi-GOF-Dm-Exp.pdf", width = 4, height = 4)

p <- res_tbl[type == "LacZ OE"] |>
ggplot(aes(x_std, y_std, color = cell)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "bottom"
  ) +
  scale_color_manual(
    values = c("#D95F02", "#7570B3", "grey80")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell type")

ggsave("./result/fig3/Bi-GOF-Dm-legend.pdf", get_legend(p))

ttbl_dm8 <- res_tbl[
  , .(total = .N, Dm8 = sum(cell == "Dm8 (Dac+/Tj, ventral)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
]

qdm8 <- ttbl_dm8 |>
  ggplot(aes(x = type_n, y = Dm8 / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  annotate(
    "segment", x = 1, xend = 2, y = 0.27, yend = 0.27, color = "black"
  ) +
  annotate(
    "text", x = 1.5, y = 0.29, label = "*", color = "black", size = 6
  ) +
  labs(
    title = "Dm8",
    y = "# Dm8 (Dac+/Tj, ventral) /\nTotal (All Dac+/Tj+)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")

tobj_dm8 <- glm(
  Dm8 / total ~ type, data = ttbl_dm8, family = quasibinomial, weights = total
)
tsobj_dm8 <- summary(tobj_dm8)
beta0 <- tsobj_dm8$coefficients[1, 1]
beta1 <- tsobj_dm8$coefficients[2, 1]
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))
phi <- tsobj_dm8$dispersion
stat_tbl[["Bi-GOF-Dm (Dm8)"]] <- data.frame(
  p_value = tsobj_dm8$coefficients[2, 4],
  n_ctrl = 3L,
  n_exp = 3L,
  fc = p1 / p0,
  p0 = p0,
  p1 = p1,
  phi = phi
)

ttbl_dra <- res_tbl[
  , .(total = .N, DmDRA = sum(cell == "Dm-DRA (Dac+/Tj, dorsal)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
]

qdra <- ttbl_dra |>
  ggplot(aes(x = type_n, y = DmDRA / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  annotate(
    "segment", x = 1, xend = 2, y = 0.22, yend = 0.22, color = "black"
  ) +
  annotate(
    "text", x = 1.5, y = 0.24, label = "*", color = "black", size = 6
  ) +
  labs(
    title = "DmDRA", 
    y = "# DmDRA (Dac+/Tj, dorsal) /\nTotal (All Dac+/Tj+)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_manual(
    values = c("#435274", "#ba3c3c")
  ) +
  guides(color = "none")

qdra / qdm8
ggsave("./result/fig3/Bi-GOF-Dm-quant.pdf", width = 1.8, height = 5.6)

tobj_dra <- glm(
  DmDRA / total ~ type, data = ttbl_dra, family = quasibinomial, weights = total
)
tsobj_dra <- summary(tobj_dra)
beta0 <- tsobj_dra$coefficients[1, 1]
beta1 <- tsobj_dra$coefficients[2, 1]
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))
phi <- tsobj_dra$dispersion
stat_tbl[["Bi-GOF-Dm (DmDRA)"]] <- data.frame(
  p_value = tsobj_dra$coefficients[2, 4],
  n_ctrl = 3L,
  n_exp = 3L,
  fc = p1 / p0,
  p0 = p0,
  p1 = p1,
  phi = phi
)

## Bi OE TE
res_path <- list.files(
  "./Bi-GOF-TE/result", full.names = TRUE
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()

res_tbl$type <- factor(res_tbl$type)
res_tbl$type <- relevel(res_tbl$type, "LacZ OE")
res_tbl[, cell := "Other"]
res_tbl[int_1 > dimm_threshold & int_2 > fs_threshold, cell := "TE (Dimm+/Fs+)"]
res_tbl$cell <- factor(res_tbl$cell, levels = c("TE (Dimm+/Fs+)", "Other"))

res_tbl[cell == "TE (Dimm+/Fs+)" & type == "LacZ OE"] |>
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
    values = c("#D95F02", "#7570B3")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Bi-GOF-TE-Ctrl.pdf", width = 4, height = 4)

res_tbl[cell == "TE (Dimm+/Fs+)" & type == "Bi OE"] |>
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
    values = c("#D95F02", "#7570B3")
  ) +
  guides(color = "none")
ggsave("./result/fig3/Bi-GOF-TE-Exp.pdf", width = 4, height = 4)

p <- res_tbl[cell == "TE (Dimm+/Fs+)" & type == "LacZ OE"] |>
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
    values = c("#D95F02", "#7570B3")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "Cell type")

ggsave("./result/fig3/Bi-GOF-TE-legend.pdf", get_legend(p))

ttbl <- res_tbl[
  , .(total = .N, TE = sum(cell == "TE (Dimm+/Fs+)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
  ]

ttbl |>
  ggplot(aes(x = type_n, y = TE / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  labs(y = "# TE (Dimm+/Fs+) /\nTotal (All Hoechst+)") +
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

ggsave("./result/fig3/Bi-GOF-TE-quant-sup.pdf", width = 1.6, height = 3.55)

tobj <- glm(
  TE / total ~ type, data = ttbl, family = quasibinomial, weights = total
)
tsobj <- summary(tobj)
beta0 <- tsobj$coefficients[1, 1]
beta1 <- tsobj$coefficients[2, 1]
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))
phi <- tsobj$dispersion
stat_tbl[["Bi-GOF-TE (Whole OL)"]] <- data.frame(
  p_value = tsobj$coefficients[2, 4],
  n_ctrl = 3L,
  n_exp = 3L,
  fc = p1 / p0,
  p0 = p0,
  p1 = p1,
  phi = phi
)

ttbl <- res_tbl[
  x_std > 0, .(total = .N, ant_TE = sum(cell == "TE (Dimm+/Fs+)")),
  by = c("type", "sample")
][, type_n := paste0(type, " (n = ", .N, ")", sep = ""), by = "type"][
  , type_n := factor(type_n, levels = rev(unique(type_n)))
]

ttbl |>
  ggplot(aes(x = type_n, y = ant_TE / total, color = type_n)) +
  geom_jitter(width = 0.1, height = 0, size = 1) +
  stat_summary(
    fun.data = "mean_se",  pch = "-", size = 2
  ) +
  # annotate(
  #   "segment", x = 1, xend = 2, y = 0.0033, yend = 0.0033, color = "black"
  # ) +
  # annotate(
  #   "text", x = 1.5, y = 0.0035, label = "*", color = "black", size = 6
  # ) +
  labs(y = "# Ectopic TE (Dimm+/Fs+) /\nTotal (All Hoechst+)") +
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

ggsave("./result/fig3/Bi-GOF-TE-quant.pdf", width = 1.6, height = 3.55)

tobj <- glm(
  ant_TE / total ~ type, data = ttbl, family = quasibinomial, weights = total
)
tsobj <- summary(tobj)
beta0 <- tsobj$coefficients[1, 1]
beta1 <- tsobj$coefficients[2, 1]
p0 <- 1 / (1 + exp(-beta0))
p1 <- 1 / (1 + exp(-beta0 - beta1))
phi <- tsobj$dispersion
stat_tbl[["Bi-GOF-TE (Anterior)"]] <- data.frame(
  p_value = tsobj$coefficients[2, 4],
  n_ctrl = 3L,
  n_exp = 3L,
  fc = p1 / p0,
  p0 = p0,
  p1 = p1,
  phi = phi
)

stat_df <- do.call(rbind, stat_tbl)
write.csv(stat_df, "./result/fig3/stat_summary.csv")
