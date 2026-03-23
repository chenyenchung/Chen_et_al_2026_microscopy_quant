library(data.table)
library(ggplot2)
library(cowplot)

if (!dir.exists("result/fig4")) dir.create("result/fig4", recursive = TRUE)

# Runt
res_path <- list.files(
  "./Ez-LOF-CtRun/result/", full.names = TRUE, pattern = "run"
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()
res_tbl$type <- factor(
  res_tbl$type, levels = c("mCherry RNAi", "E(z) RNAi"),
)

res_tbl |>
  ggplot(aes(x = x_std, y = y_std, color = type)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  guides(color = "none") +
  scale_color_manual(
    values = c("#435274","#ba3c3c")
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  )
ggsave("./result/fig4/Ez-LOF-Run_compare.pdf", width = 4, height = 4)

# Cut
res_path <- list.files(
  "./Ez-LOF-CtRun/result/", full.names = TRUE, pattern = "ct"
)

res_tbl <- lapply(res_path, fread) |>
  rbindlist()
res_tbl$type <- factor(
  res_tbl$type, levels = c("mCherry RNAi", "E(z) RNAi"),
)
res_tbl |>
  ggplot(aes(x = x_std, y = y_std, color = type)) +
  geom_point(size = 0.1, alpha = 0.5) +
  theme_minimal() +
  guides(color = "none") +
  scale_color_manual(
    values = c("#435274","#ba3c3c")
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  )

ggsave("./result/fig4/Ez-LOF-Ct_compare.pdf", width = 4, height = 4)

p <- res_tbl |>
  ggplot(aes(x = x_std, y = y_std, color = type)) +
  geom_point(size = 0.1, alpha = 0.5) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  scale_color_manual(
    values = c("#435274","#ba3c3c")
  ) +
  labs(color = "Condition") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  )
ggsave("./result/fig4/Ez-LOF-compare_legend.pdf",get_legend(p))
