library(ggplot2)
library(data.table)

if (!dir.exists("result/fig4")) dir.create("result/fig4", recursive = TRUE)

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

write.csv(rbindlist(stat_tbl), "result/fig4/stat.csv", row.names = FALSE)

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
