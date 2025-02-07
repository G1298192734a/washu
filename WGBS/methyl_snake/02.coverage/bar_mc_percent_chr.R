suppressPackageStartupMessages(library(tidyverse))
library(ggsci)

read_table("ctype_classify_chr.txt",show_col_types = 0) -> mc
round(mc$mC_rate*100,2) -> mc$mC_rate

ggplot(mc,aes(x=chrom,y=mC_rate,fill=factor(ctype,levels=c("CG","CHG","CHH")))) + # factor 用于调节顺序
  geom_bar(color="black",position="stack",width=0.6,stat="identity") +
  labs(title="Classification of chrom mC",x="chromsome",y="percent(%)") +
  scale_fill_futurama() +
  theme_classic(base_family = "serif",base_size = 20) +
  theme(
    panel.border = element_rect(color="black",fill=NA),
    plot.margin=margin(10,30,10,30),
    plot.title = element_text(hjust = 0.5,face="bold",size=14),
    legend.title = element_blank(),
    legend.position="bottom",
    #legend.justification=c(.5,0),
    axis.line = element_blank(),
    axis.ticks = element_line(color="black"),
    axis.text.x = element_text(hjust = 1,angle=45),
  ) -> pl

ggsave("bar_mc_percent_chr.pdf",pl)


