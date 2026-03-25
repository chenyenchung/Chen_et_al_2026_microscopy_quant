library(ggplot2)
library(data.table)

if (!dir.exists("result/fig5")) {
  dir.create("result/fig5", recursive = TRUE)
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
  annotate(
    "segment", x = 1, xend = 2, y = 1.1, yend = 1.1, color = "black"
  ) +
  annotate(
    "text", x = 1.5, y = 1.2, label = "*", color = "black", size = 6
  ) +
  labs(y = "# Vsx2+/Hth+ in clone /\n# Hth+ in clone") +
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
ggsave("./result/fig5/Su(H)-MARCM-Vsx2-quant.pdf", width = 3, height = 4)
