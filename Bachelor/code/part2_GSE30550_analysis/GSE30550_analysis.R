########## GSE30550个体特异性网络熵分析
### 时间点: Baseline Hour00 Hour005 Hour012 Hour021 Hour029 Hour036 Hour045 Hour053 Hour060 Hour077 Hour084 Hour093 Hour101 Hour108
## 步骤1: 以baseline~Hour021为参考样本, 计算均值与标准差
## 步骤2: 对之后的时间进行分析, 将各时间点与参考样本进行合并, 计算每个基因的标准差与基因之间的pearson相关系数
## 步骤3: 按照PPI关系, 计算基因之间的因果强度指数w判断是否有指向关系，判断阈值为1
## 步骤4: 计算每个基因的波动度FT
## 步骤5：对每个中心基因计算局部因果熵
## 步骤6: 取局部因果熵排序前5%的基因为DNBs, 计算全局因果熵
## 步骤7: 对DNBs进行功能注释

# 必须数据集, PPI数据集以及功能数据集
PPI <- read.csv("D:\\Github\\GitDB\\CPU\\Bachelor\\data\\GEO\\GSE30550\\PPI.csv")
func <- read.csv("D:\\Github\\GitDB\\CPU\\Bachelor\\data\\GEO\\GSE30550\\function.csv")
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
##### 步骤1
# 导入受试对象时点基因表达矩阵(行为基因, 列为时间点)
subject_matrix_path <- "D:\\Github\\GitDB\\CPU\\Bachelor\\data\\GEO\\GSE30550\\Subject 01.csv"
subject_matrix <- read.csv(subject_matrix_path,row.names=1)
colnames(subject_matrix) <- c("point0","point1","point2","point3","point4","point5","point6","point7","point8","point9",
                              "point10","point11","point12","point13","point14","point15")
# 以baseline~Hour021为参考样本, 计算均值与标准差
ref_matrix <- subject_matrix[,c("point0","point1","point2","point3","point4")]
ref_means <- apply(ref_matrix,1,mean)
ref_std <- apply(ref_matrix,1,sd)

##### 步骤2
# 取参考样本之后的时间点样本并与参考样本合并, 计算每个基因的标准差与基因之间的pearson相关系数
case_matrix <- subject_matrix[,!colnames(subject_matrix) %in% c("point0","point1","point2","point3","point4")]
point5 <- case_matrix["point15"]
point5_with_ref <- cbind(ref_matrix,point5)
point5_std <- apply(point5_with_ref,1,sd)
point5_PCC <- cor(t(point5_with_ref))

##### 步骤3
# 按照PPI关系，计算基因之间的因果强度指数w
gene.w <- list()
for (gene in from_symbol){
  w <- c()
  to_symbol <- gene.relations[[gene]]
  for (i in 1:length(to_symbol)){
    SD <- point5_std[[to_symbol[i]]]
    PCC <- point5_PCC[to_symbol[i],to_symbol[i]]
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
for (gene in from_symbol){
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
for (gene in from_symbol){
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
      }
    }
    names(w_out) <- w_out_genes
    gene.w_out[[gene]] <- w_out
  }
}

##### 步骤4
# 计算每个基因的波动度FT
case_means <- point5$point15
names(case_means) <- rownames(point5)
FT <- abs(case_means-ref_means)/ref_std

##### 计算每个中心基因的局部因果网络熵
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
  gene.E.in <- c(gene.E.in,E.in/n)
}
names(gene.E.in) <- from_symbol

# 局部出度熵的计算
gene.E.out <- c()
for (hub in names(gene.relations)){
  w_out <- gene.w_out[[hub]]
  E.out <- 0
  n <- 0
  if (sum(w_out)!=0){
    FT_bar <- mean(FT[c(names(w_out))])
    bottom <- 0
    for (i in 1:length(w_out)){
      bottom <- bottom+w_out[i]*gene.I[[names(w_out[i])]][[hub]]
    }
    for (i in 1:length(w_out)){
      to <- names(w_out[i])
      p.out <- w_out[i]*gene.I[[names(w_out[i])]][[hub]]/bottom
      E.out <- E.out+(-p.out*log2(p.out)*exp(w_out[i])/(1+exp(-FT_bar)))
      n <- n+1
    }
  }else{
    E.out <- E.out
  }
  gene.E.out <- c(gene.E.out,E.out/n)
}
names(gene.E.out) <- from_symbol

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
names(global_E) <- from_symbol
global_E <- sort(global_E,decreasing=TRUE)
top <- floor(length(global_E)*0.05)
top_E <- global_E[1:top] 
SI.ICNE <- sum(top_E)

