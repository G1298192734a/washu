suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(data.table))
library(ggsci)

#fread("/work/frasergen/backup/3d/project/WGBS/241891_zhenjun/analysis_20241029/03.single/Tr01-14_B/Tr01-14_B.deduplicated_splitting_report.txt") -> mc 
read_table("ctype_classify.txt",show_col_types = 0) -> mc
round(mc$ctype_rate*100 ,2) -> mc$ctype_rate
str_c("m",mc$ctype," (", mc$mC,", ",mc$ctype_rate, "%)") -> mc$label

ggplot(mc) + # factor 用于调节顺序
  geom_col(aes(x="",y=ctype_rate,fill=factor(label,levels=label)),color="black") +
  coord_polar(theta = "y",start=90) +
  labs(title="Classification of mC") +
  scale_fill_futurama() +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5,face="bold",size=14),
    legend.title = element_blank(),
    legend.key.size = unit(0.6, "line"),      # legend 大小和间距
    legend.key.spacing.y = unit(0.3, "line"),
    legend.text = element_text(size=12)
  ) -> pl

ggsave("par_mc_percent.pdf",pl)


