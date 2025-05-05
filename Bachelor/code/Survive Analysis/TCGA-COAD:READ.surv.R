library(survival) # install.packages("survival")
library(survminer) # install.packages("survminer")

# 导入TCGA-COAD生存信息
coad.surv <- read.csv("/Users/mymacstudio/GitDB/CPU/Bachelor/data/TCGA/TCGA-COAD/analysis/survival_info.csv")
# 导入TCGA-READ生存信息
read.surv <- read.csv("/Users/mymacstudio/GitDB/CPU/Bachelor/data/TCGA/TCGA-READ/analysis/survival_info.csv")

## 对TCGA-COAD生存信息进行分析
df <- data.frame(
  time=coad.surv$Time/30,
  status=coad.surv$Event,
  CriticalPoint=factor(coad.surv$Risk,levels=c("before","after"))
)

# 创建生存对象
surv_object <- Surv(time=df$time, event=df$status)
# 构建生存模型(按在临界点前后分层)
fit <- survfit(surv_object~CriticalPoint, data=df)
# 生成可视化图形
p1 <- ggsurvplot(fit, 
                 pval=TRUE, # 显示Log-rank检验p值
                 pval.size=5,
                 pval.coord=c(30,0.1),
                 risk.table=FALSE, # 显示风险人数表
                 conf.int=TRUE, # 显示置信区间
                 conf.int.style="step",
                 censor.shape=124,
                 censor.size=2.5,
                 legend=c(0.9,0.1),
                 break.x.by=12,
                 break.y.by=0.1,
                 surv.scale="percent",
                 pval.method=TRUE,
                 pval.method.size=5,
                 pval.method.coord=c(10,0.1),
                 ggtheme=theme_minimal()+
                   theme(
                     panel.border=element_rect(colour="black",fill=NA,linewidth=1.5),
                     axis.text.x=element_text(family="Arial", size=12),
                     axis.text.y=element_text(family="Arial", size=12),
                     axis.title.x=element_text(family="Arial", size=16),
                     axis.title.y=element_text(family="Arial", size=16),
                     legend.justification=c(1,0),
                     legend.direction="vertical",
                     legend.background=element_rect(fill=NA,color="black",
                       linewidth=0.5
                     ),
                     legend.text=element_text(family="Arial", size=8),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()
                     ), 
                 palette="lancet" # 自定义配色
                 )+
  labs(x="Time(month)",
       y="Surviavl Probability")
p1

## 对TCGA-READ生存信息进行分析
df <- data.frame(
  time=read.surv$Time/30,
  status=read.surv$Event,
  CriticalPoint=factor(read.surv$Risk,levels=c("before","after"))
)

# 创建生存对象
surv_object <- Surv(time=df$time, event=df$status)
# 构建生存模型(按在临界点前后分层)
fit <- survfit(surv_object~CriticalPoint, data=df)
# 生成可视化图形
p2 <- ggsurvplot(fit, 
                 pval=TRUE, # 显示Log-rank检验p值
                 pval.size=5,
                 pval.coord=c(30,0.1),
                 risk.table=FALSE, # 显示风险人数表
                 conf.int=TRUE, # 显示置信区间
                 conf.int.style="step",
                 censor.shape=124,
                 censor.size=2.5,
                 legend=c(0.9,0.1),
                 break.x.by=12,
                 break.y.by=0.1,
                 surv.scale="percent",
                 pval.method=TRUE,
                 pval.method.size=5,
                 pval.method.coord=c(10,0.1),
                 ggtheme=theme_minimal()+
                   theme(
                     panel.border=element_rect(colour="black",fill=NA,linewidth=1.5),
                     axis.text.x=element_text(family="Arial", size=12),
                     axis.text.y=element_text(family="Arial", size=12),
                     axis.title.x=element_text(family="Arial", size=16),
                     axis.title.y=element_text(family="Arial", size=16),
                     legend.justification=c(1,0),
                     legend.direction="vertical",
                     legend.background=element_rect(fill=NA,color="black",
                                                    linewidth=0.5
                     ),
                     legend.text=element_text(family="Arial", size=8),
                     panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()
                   ), 
                 palette="lancet" # 自定义配色
)+
  labs(x="Time(month)",
       y="Surviavl Probability")
p2
