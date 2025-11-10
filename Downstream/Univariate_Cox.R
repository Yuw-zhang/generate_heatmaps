library(survival)
library(broom)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)


colon <-  read.csv('/Users/yuw/Documents/Development/Colon_Data/data/read.csv')

###########Select effective rows#############
colon$til_status <- factor(colon$til_status, levels=c('low', 'high'))
colon$AJCC_stage <- factor(colon$AJCC_stage, levels=c('Stages I-II', 'Stages III-IV'), labels=c('I-II', 'III-IV'))
colon$tumorSize_status <-factor(colon$tumorSize_status, levels=c('low', 'high'),
                                labels = c("≤4.8cm", ">4.8cm"))
colon$MSI <- factor(colon$MSI, levels=c('MSS/MSI-L', 'MSI-H'))
colon$histoGrade <- factor(colon$histoGrade, levels=c('Well-Mod', 'Poorly'), labels=c('Well-Mod Dif', 'Poorly Dif'))
colon$immuneSubtype <- factor(colon$immuneSubtype, levels=c('C1', 'C2'), labels=c('C1*', 'C2*'))
#colon$buddingScore_status <- factor(colon$buddingScore_status, levels=c('Bd1/Bd2', 'Bd3'),
#                                    labels=c("G1/G2", "G3"))


###########Modify the column name#############
names(colon)[names(colon)=='til_status'] <- 'TILs%'
names(colon)[names(colon)=='AJCC_stage'] <- 'AJCC stages'
names(colon)[names(colon)=='tumorSize_status'] <- 'Tumor size'
names(colon)[names(colon)=='histoGrade'] <- 'Histological grades'
names(colon)[names(colon)=='immuneSubtype'] <- 'Immune subtypes'
names(colon)[names(colon)=='MSI'] <- 'MSI status'
#names(colon)[names(colon)=='buddingScore_status'] <- 'Tumor budding score'


# 定义你要分析的变量
vars <- c("`TILs%`", "`AJCC stages`", "`Tumor size`", "`Histological grades`", "`Immune subtypes`", "MSI status")
#vars <- c("`TILs%`", "`AJCC stages`", "`Tumor size`", "`Histological grades`", "`Immune subtypes`", "`Risk score`")
#vars <- c("`TILs%`", '`AJCC stages`', "`Tumor size`", "`MSI status`")

# 跑 univariate Cox 回归
uni_results <- map_dfr(vars, function(var) {
  formula <- as.formula(paste("Surv(OS.month, OS) ~", var))
  model <- coxph(formula, data = colon)
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(variable = gsub("`", "", var))
})

# 提取 term 和 label
uni_results <- uni_results %>%
  mutate(
    row_id = row_number(),
    term_clean = gsub(".*(low|high|I-II|III-IV|>4.9cm|≤4.9cm|G1/G2|G3|C2\\*|Well-Mod Dif|Poorly Dif|MSS/MSI-L|MSI-H|≤0|>0).*", "\\1", term),
    label = term_clean,
    p_label = ifelse(p.value < 0.001, "p < 0.001", paste0("p = ", round(p.value, 3)))
  )

# 获取参考组水平（每个变量的第一个 factor level）
ref_terms <- map_chr(vars, function(var) {
  var_clean <- gsub("`", "", var)
  levels(colon[[var_clean]])[1]
})
names(ref_terms) <- gsub("`", "", vars)

# 构建参考组行
ref_rows <- map_dfr(names(ref_terms), function(var) {
  data.frame(
    variable = var,
    term = NA,
    term_clean = ref_terms[[var]],
    estimate = 1,
    conf.low = 1,
    conf.high = 1,
    p.value = NA,
    p_label = "",
    label = paste0(var, ": ", ref_terms[[var]])
  )
})

# 给 uni_results 添加变量顺序编号（手动定义顺序）
var_order <- c(
  "TILs%", 
  "AJCC stages",
  "Tumor size", 
  "Histological grades",
  "Immune subtypes",
  "MSI status"
)

# 给每条记录加上变量顺序
uni_results <- uni_results %>%
  mutate(
    row_id = row_number(),  # 保留原顺序（如图中顺序）
    var_factor = factor(variable, levels = var_order)
  )

# 创建参考组（ref_rows）
ref_terms <- map_chr(var_order, function(var) {
  var_clean <- gsub("`", "", var)
  levels(colon[[var_clean]])[1]
})
names(ref_terms) <- var_order

ref_rows <- map_dfr(names(ref_terms), function(var) {
  data.frame(
    variable = var,
    term = NA,
    term_clean = ref_terms[[var]],
    estimate = 1,
    conf.low = 1,
    conf.high = 1,
    p.value = NA,
    p_label = "",
    label = paste0(var, ": ", ref_terms[[var]]),
    var_factor = factor(var, levels = var_order),  # 保持同一列
    row_id = -1  # 参考组设为 -1，以便排序在前
  )
})

# 合并并排序
uni_all <- bind_rows(ref_rows, uni_results) %>%
  arrange(var_factor, row_id) %>%
  mutate(
    label = factor(label, levels = rev(unique(label)))
  )
# 图形参数
max_hr <- max(uni_all$conf.high, na.rm = TRUE)

# 颜色设置
uni_all <- uni_all %>%
  mutate(color = ifelse(grepl(":", label), "lightblue", "orange"))

# 绘图
library(ggplot2)

png(filename = "/Users/yuw/Desktop/os-unCox-CRC.png",
    width = 12, height = 5, units = "in", res = 300)

ggplot(uni_all, aes(y = label, x = estimate)) +
  geom_point(aes(color = color), size = 4) +
  geom_segment(aes(x = conf.low, xend = conf.high, y = label, yend = label, color = color), size = 1.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_text(aes(x = max_hr + 0.4, label = p_label), hjust = 0, size = 4.2) +
  scale_color_identity() +
  labs(x = "Hazard Ratio (95% CI)", y = NULL,
       title = "Univariate Cox Regression\nOverall Survival: TCGA READ") +
  theme_minimal(base_size = 14) +
  xlim(0, max_hr + 1) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_blank(),
    axis.text.y = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_blank()
  )

dev.off()