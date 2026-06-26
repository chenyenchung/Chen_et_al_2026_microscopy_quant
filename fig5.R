library(ggplot2)
library(data.table)

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
