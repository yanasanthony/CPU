########## TCGA-COAD个体特异性网络熵分析
## 步骤1: 以正常样本为参考样本, 计算均值与标准差
## 步骤2: 对之后的时间进行分析, 将各时间点与参考样本进行合并, 计算每个基因的标准差与基因之间的pearson相关系数
## 步骤3: 按照PPI关系, 计算基因之间的因果强度指数w判断是否有指向关系，判断阈值为1
## 步骤4: 计算每个基因的波动度FT
## 步骤5：对每个中心基因计算局部因果熵
## 步骤6: 取局部因果熵排序前5%的基因为DNBs, 计算全局因果熵
## 步骤7: 对DNBs进行功能注释

# 设定工作目录
setwd("/Users/mymacstudio/GitDB/CPU/Bachelor/data/TCGA/TCGA-COAD")
# 必须数据集, 表达矩阵, PPI数据集以及临床数据集
# 参考样本表达矩阵并将表达值对数化
ref <- read.csv("ref/ref.csv",row.names=1)
ref[sapply(ref, is.numeric)] <- lapply(ref[sapply(ref, is.numeric)], function(x) log2(x + 1))
# 病例样本表达矩阵并将表达值对数化
case <- read.csv("case/case.csv",row.names=1)
case[sapply(case, is.numeric)] <- lapply(case[sapply(case, is.numeric)], function(x) log2(x + 1))
# PPI信息
PPI <- read.csv("PPI.csv")
# 临床信息
clinical_info <- read.csv("clinical/clinical_info.csv",row.names=1)

# 根据临床分期信息, 将样本归类为各时间点, 并以平均表达量作为时间点表达量
# 该数据集按照各分期数据量, 将其划分为4个时间点Stage I, Stage II, Stage III, Stage IV
stage <- clinical_info$stage
new_stage <- c()
for (i in 1:length(stage)){
  if (is.na(stage[i])){
    new_stage <- c(new_stage,NA)
  }else if (stage[i]=="Stage I" | stage[i]=="Stage IA"){
    new_stage <- c(new_stage,"Stage I")
  }else if (stage[i]=="Stage II" | stage[i]=="Stage IIA" | stage[i]=="Stage IIB" | stage[i]=="Stage IIC"){
    new_stage <- c(new_stage,"Stage II")
  }else if (stage[i]=="Stage III" | stage[i]=="Stage IIIA" | stage[i]=="Stage IIIB" | stage[i]=="Stage IIIC"){
    new_stage <- c(new_stage,"Stage III")
  }else{
    new_stage <- c(new_stage,"Stage IV")
  }
}
clinical_info$stage <- new_stage
stage.1.barcode <- rownames(clinical_info[clinical_info$stage=="Stage I" & !is.na(clinical_info$stage),])
stage.2.barcode <- rownames(clinical_info[clinical_info$stage=="Stage II" & !is.na(clinical_info$stage),])
stage.3.barcode <- rownames(clinical_info[clinical_info$stage=="Stage III" & !is.na(clinical_info$stage),])
stage.4.barcode <- rownames(clinical_info[clinical_info$stage=="Stage IV" & !is.na(clinical_info$stage),])

samples <- colnames(case)
samples <- gsub("\\.", "-", samples)
colnames(case) <- samples
stage.1.samples <- c()
stage.2.samples <- c()
stage.3.samples <- c()
stage.4.samples <- c()

for (i in 1:length(samples)){
  barcode <- substr(samples[i],1,12)
  if (barcode %in% stage.1.barcode){
    stage.1.samples <- c(stage.1.samples,samples[i])
  }else if (barcode %in% stage.2.barcode){
    stage.2.samples <- c(stage.2.samples,samples[i])
  }else if (barcode %in% stage.3.barcode){
    stage.3.samples <- c(stage.3.samples,samples[i])
  }else{
    stage.4.samples <- c(stage.4.samples,samples[i])
  }
}

# Stage I
stage.1.mean.expression <- apply(case[,stage.1.samples],1,mean)
stage.1.std.expression <- apply(case[,stage.1.samples],1,sd)
stage.1.PCC <- cor(t(case[,stage.1.samples]))
# Stage II
stage.2.mean.expression <- apply(case[,stage.2.samples],1,mean)
stage.2.std.expression <- apply(case[,stage.2.samples],1,sd)
stage.2.PCC <- cor(t(case[,stage.2.samples]))
# Stage III
stage.3.mean.expression <- apply(case[,stage.3.samples],1,mean)
stage.3.std.expression <- apply(case[,stage.3.samples],1,sd)
stage.3.PCC <- cor(t(case[,stage.3.samples]))
# Stage IV
stage.4.mean.expression <- apply(case[,stage.4.samples],1,mean)
stage.4.std.expression <- apply(case[,stage.4.samples],1,sd)
stage.4.PCC <- cor(t(case[,stage.4.samples]))

case.mean.expression <- data.frame(
  StageI=stage.1.mean.expression,
  StageII=stage.2.mean.expression,
  StageIII=stage.3.mean.expression,
  StageIV=stage.4.mean.expression
)
case.std.expression <- data.frame(
  StageI=stage.1.std.expression,
  StageII=stage.2.std.expression,
  StageIII=stage.3.std.expression,
  StageIV=stage.4.std.expression
)


# 根据PPI创建基因关系字典
gene.relations <- list()
from_symbol <- unique(PPI$from_symbol)
for (gene in from_symbol){
  to_symbol <- PPI[(PPI$from_symbol==gene),]$to_symbol
  gene.relations[[gene]] <- to_symbol
}
# 根据PPI创建基因互信息字典
gene.I <- list()
for (gene in from_symbol){
  I <- PPI[(PPI$from_symbol==gene),]$combined_score
  names(I) <- gene.relations[[gene]]
  gene.I[[gene]] <- I
}

# 创建存储各样本(可以看成是时间)SI.ICNE的容器
SI.ICNE.vec <- c()

##### 步骤1
# 以正常样本为参考样本, 计算均值和标准差
ref_means <- apply(ref,1,mean)
ref_std <- apply(ref,1,sd)
# 以参考样本的均值和标准差作为总体均值和标准差, 计算波动度作为检索
FT.df <- case.mean.expression
for (i in 1:ncol(FT.df)){
  FT.df[[i]] <- abs(FT.df[[i]]-ref_means)/ref_std
}

##### 步骤2
# 取各时间点的样本并与参考样本合并, 计算每个基因的标准差与基因之间的pearson相关系数
point <- colnames(case.mean.expression)[3]
case_std <- case.std.expression[[point]]
names(case_std) <- rownames(case.mean.expression)
if (point=="StageI"){
  case_PCC <- stage.1.PCC
}else if (point=="StageII"){
  case_PCC <- stage.2.PCC
}else if (point=="StageIII"){
  case_PCC <- stage.3.PCC
}else{
  case_PCC <- stage.4.PCC
}

##### 步骤3
# 按照PPI关系，计算基因之间的因果强度指数w
gene.w <- list()
for (gene in names(gene.relations)){
  w <- c()
  to_symbol <- gene.relations[[gene]]
  for (i in 1:length(to_symbol)){
    SD <- case_std[[to_symbol[i]]]
    PCC <- abs(case_PCC[gene,to_symbol[i]])
    w <- c(w,SD*PCC)
  }
  gene.w[[gene]] <- w
}
# 判断是否有指向关系
gene.judge.w <- list()
for (gene in from_symbol){
  judge <- c()
  values <- gene.w[[gene]]
  for (i in 1:length(values)){
    if (values[i]>1){
      judge <- c(judge,1)
    }else{judge <- c(judge,0)}
  }
  gene.judge.w[[gene]] <- judge
}
# 计算入度因果强度和出度因果强度
# 入度因果强度
gene.w_in <- list()
for (gene in names(gene.relations)){
  to_symbol <- gene.relations[[gene]]
  w <- gene.w[[gene]]
  judge <- gene.judge.w[[gene]]
  w_in <- c()
  w_in_genes <- c()
  for (i in 1:length(w)){
    if (judge[i]==1){
      w_in <- c(w_in,w[i])
      w_in_genes <- c(w_in_genes,to_symbol[i])
    }else{
      w_in <- c(w_in,0)
      w_in_genes <- c(w_in_genes,to_symbol[i])
    }
  }
  names(w_in) <- w_in_genes
  gene.w_in[[gene]] <- w_in
}
# 出度因果强度
gene.w_out <- list()
for (gene in names(gene.relations)){
  from_from_symbol <- PPI[PPI$to_symbol==gene,]$from_symbol
  w_out <- c()
  w_out_genes <- c()
  if (length(from_from_symbol)==0){
    gene.w_out[[gene]] <- c(0)
  }else{
    for (i in 1:length(from_from_symbol)){
      index <- which(gene.relations[[from_from_symbol[i]]]==gene)
      w <- gene.w[[from_from_symbol[i]]][index]
      judge <- gene.judge.w[[from_from_symbol[i]]][index]
      if (judge==1){
        w_out <- c(w_out,w)
        w_out_genes <- c(w_out_genes,from_from_symbol[i])
      }else{
        w_out <- c(w_out,0)
        w_out_genes <- c(w_out_genes,from_from_symbol[i])
      }
    }
    names(w_out) <- w_out_genes
    gene.w_out[[gene]] <- w_out
  }
}

##### 步骤4
# 计算每个基因的波动度FT
FT <- FT.df[[point]]
names(FT) <- rownames(FT.df)

##### 步骤5
# 计算每个中心基因的局部因果网络熵
# 局部入度熵的计算
gene.E.in <- c()
for (hub in names(gene.relations)){
  w_in <- gene.w_in[[hub]]
  E.in <- 0
  n <- 0
  if (sum(w_in)!=0){
    for (i in 1:length(w_in)){
      if (w_in[i]!=0){
        p.in <- w_in[i]*gene.I[[hub]][i]/sum(w_in*gene.I[[hub]])
        E.in <- E.in + (-p.in*log2(p.in)*exp(w_in[i])/(1+exp(-FT[[hub]])))
        n <- n+1
      }else{
        E.in <- E.in
      }
    }
  }else{
    E.in <- E.in
  }
  if (n == 0){
    gene.E.in <- c(gene.E.in,0)
  }else{
    gene.E.in <- c(gene.E.in,E.in/n)
  }
}
names(gene.E.in) <- names(gene.relations)

# 局部出度熵的计算
gene.E.out <- c()
for (hub in names(gene.relations)){
  w_out <- gene.w_out[[hub]]
  E.out <- 0
  n <- 0
  if (sum(w_out)!=0){
    w_out_ <- w_out[w_out!=0]
    FT_bar <- mean(FT[c(names(w_out_))])
    bottom <- 0
    for (i in 1:length(w_out_)){
      bottom <- bottom+w_out_[i]*gene.I[[names(w_out_[i])]][[hub]]
    }
    for (i in 1:length(w_out_)){
      to <- names(w_out_[i])
      p.out <- w_out_[i]*gene.I[[names(w_out_[i])]][[hub]]/bottom
      E.out <- E.out+(-p.out*log2(p.out)*exp(w_out_[i])/(1+exp(-FT_bar)))
      n <- n+1
    }
  }else{
    E.out <- E.out
  }
  if (n == 0){
    gene.E.out <- c(gene.E.out,0)
  }else{
    gene.E.out <- c(gene.E.out,E.out/n)
  }
}
names(gene.E.out) <- names(gene.relations)

#####步骤6 
# 全局熵的计算
global_E <- c()
for (hub in names(gene.relations)){
  w_in <- sum(gene.w_in[[hub]])
  w_out <- sum(gene.w_out[[hub]])
  if (w_in == 0 & w_out == 0){
    global_E <- c(global_E,0)
  }else{
    alpha <- w_in/(w_in+w_out)
    global_E <- c(global_E,log10(alpha*gene.E.in[[hub]]+(1-alpha)*gene.E.out[[hub]]))
  }
}
names(global_E) <- names(gene.relations)
global_E <- sort(global_E,decreasing=TRUE)
top <- floor(length(global_E)*0.05)
top_E <- global_E[1:top] 
SI.ICNE <- sum(top_E)


# Stage I:22.735340057943
# Stage II:21.0031581364828
# Stage III:24.535674805744
# Stage IV:22.1400162312294

# 保存基因
hubs <- names(top_E)
DNBs <- PPI[PPI$from_symbol %in% hubs,]
to <- unique(DNBs[["to_symbol"]])
DNBs <- unique(hubs,to)
DNBs <- data.frame(
  SYMBOL=DNBs
)
write.csv(DNBs,"analysis/DNBs.csv",row.names=FALSE)

# 保存临床信息以及生存信息
write.csv(clinical_info[!is.na(clinical_info$stage),!(colnames(clinical_info)%in%c("event","time"))],
          "analysis/clinical.csv",row.names=TRUE)
surv <- clinical_info[!is.na(clinical_info$event) & !is.na(clinical_info$time),]
surv <- surv[!is.na(surv$stage),]
Time <- surv$time
Event <- surv$event
Risk <- c()
for (i in 1:nrow(surv)){
  if (surv$stage[i]=="Stage I" | surv$stage[i]=="Stage II"){
    Risk <- c(Risk,"before")
  }else{
    Risk <- c(Risk,"after")
  }
}
surv.info <- data.frame(
  Time=Time,
  Event=Event,
  Risk=Risk
)
write.csv(surv.info,"analysis/survival_info.csv",row.names=FALSE)

