########## GSE30550个体特异性网络熵分析
### 时间点: Baseline Hour00 Hour005 Hour012 Hour021 Hour029 Hour036 Hour045 Hour053 Hour060 Hour077 Hour084 Hour093 Hour101 Hour108
## 步骤1: 以baseline~Hour021为参考样本, 计算均值与标准差
## 步骤2: 对之后的时间进行分析, 将各时间点与参考样本进行合并, 计算每个基因的标准差与基因之间的pearson相关系数
## 步骤3: 按照PPI关系, 计算基因之间的因果强度指数w判断是否有指向关系，判断阈值为1
## 步骤4: 计算每个基因的波动度FT
## 步骤5：对每个中心基因计算局部因果熵
## 步骤6: 取局部因果熵排序前5%的基因为DNBs, 计算全局因果熵
## 步骤7: 对DNBs进行功能注释

# 设定工作目录
setwd("/Users/mymacstudio/GitDB/CPU/Bachelor/data/GEO/GSE30550")
# 必须数据集, PPI数据集以及功能数据集
PPI <- read.csv("PPI.csv")
func <- read.csv("function.csv")
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

# 受试对象基因表达矩阵文件
subject_matrix_path <- c("Subject 01.csv","Subject 02.csv","Subject 03.csv","Subject 04.csv",
                         "Subject 05.csv","Subject 06.csv","Subject 07.csv","Subject 08.csv",
                         "Subject 09.csv","Subject 10.csv","Subject 11.csv","Subject 12.csv",
                         "Subject 13.csv","Subject 14.csv","Subject 15.csv","Subject 16.csv",
                         "Subject 17.csv")
# 时间点
timepoints <- c("point0","point1","point2","point3","point4","point5",
                "point6","point7","point8","point9","point10",
                "point11","point12","point13","point14","point15")
# 创建存储各受试对象SI.ICNE的容器
SI.ICNE.matrix <- list()
# 创建存储各受试对象中心基因的容器
hub.genes.matrix <- list()
# 创建存储各受试对象p值的容器
p.values.matrix <- list()

for (csv in 1:length(subject_matrix_path)){
  ##### 步骤1
  # 导入受试对象时点基因表达矩阵(行为基因, 列为时间点), 时间点重命名
  subject <- subject_matrix_path[csv]
  # 开始提示
  message(paste0("Starting analyzing ",substr(subject,1,10)))
  subject_matrix <- read.csv(subject,row.names=1)
  colnames(subject_matrix) <- timepoints[1:ncol(subject_matrix)]
  
  # 以baseline~Hour021为参考样本, 计算均值与标准差
  ref_matrix <- subject_matrix[,timepoints[1:4]]
  ref_means <- apply(ref_matrix,1,mean)
  ref_std <- apply(ref_matrix,1,sd)
  
  # 以参考样本的均值和标准差作为总体均值和标准差, 计算波动度作为检索
  FT_matrix <- subject_matrix
  for (p in 1:ncol(FT_matrix)){
    FT_matrix[[timepoints[p]]] <- abs(FT_matrix[[timepoints[p]]]-ref_means)/ref_std
  }
  
  ##### 创建存储各时间点SI.ICNE的容器
  SI.ICNE.subject <- c()
  ##### 创建存储各时间点中心基因的容器
  hub.genes.subject <- list()
  ##### 创建存储配对样本正态检验p值的容器
  p.values.subject <- c()
  
  ##### 步骤2
  # 取各时间点的样本并与参考样本合并, 计算每个基因的标准差与基因之间的pearson相关系数
  for (p in 1:ncol(subject_matrix)){
    point <- timepoints[p]
    point_expression <- subject_matrix[point]
    if (point %in% timepoints[1:5]){
      with_ref <- ref_matrix
    }else{
      with_ref <- cbind(ref_matrix,point_expression)
    }
    case_std <- apply(with_ref,1,sd)
    case_PCC <- cor(t(with_ref))
  
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
    FT <- FT_matrix[[point]]
    names(FT) <- rownames(FT_matrix)
  
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
    names(gene.E.in) <- from_symbol
  
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
    names(gene.E.out) <- from_symbol
  
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
    names(global_E) <- from_symbol
    global_E <- sort(global_E,decreasing=TRUE)
    top <- floor(length(global_E)*0.05)
    top_E <- global_E[1:top] 
    SI.ICNE <- sum(top_E)
    SI.ICNE.subject <- c(SI.ICNE.subject,SI.ICNE)
    hub.genes.subject[[point]] <- names(top_E)
  }
  # 获取最大SI.ICNE的时间点
  max_point <- c(timepoints[which(SI.ICNE.subject==max(SI.ICNE.subject))])[1]
  # 存储该受试对象的最大SI.ICNE时间点处的中心基因以及各时间点SI.ICNE
  hub.genes.matrix[[substr(subject,1,10)]] <- hub.genes.subject[[max_point]]
  SI.ICNE.matrix[[substr(subject,1,10)]] <- SI.ICNE.subject
  p.values.matrix[[substr(subject,1,10)]] <- p.values.subject[max_point]
  
  # 结束提示
  message("Analysis has been done successfully!")
}


# 保存SI.ICNE结果文件, 因为Subject 08, Subject 13, Subject 17缺少时间点数据, 故需要特殊处理
SI.ICNE.df <- data.frame(matrix(nrow = 16, ncol = 0))
rownames(SI.ICNE.df) <- timepoints
for (csv in subject_matrix_path){
  subject <- substr(csv,1,10)
  if (subject=="Subject 08"){
    SI.ICNE <- SI.ICNE.matrix[[subject]]
    SI.ICNE <- append(SI.ICNE,0,after=4)
    SI.ICNE.df[subject] <- SI.ICNE
  }else if (subject=="Subject 13"){
    SI.ICNE <- SI.ICNE.matrix[[subject]]
    SI.ICNE <- c(0,SI.ICNE)
    SI.ICNE <- append(SI.ICNE,0,after=6)
    SI.ICNE.df[subject] <- SI.ICNE
  }else if (subject=="Subject 17"){
    SI.ICNE <- SI.ICNE.matrix[[subject]]
    SI.ICNE <- append(SI.ICNE,0,after=6)
    SI.ICNE.df[subject] <- SI.ICNE
  }else{
    SI.ICNE.df[subject] <- SI.ICNE.matrix[[subject]]
  }
}
write.csv(SI.ICNE.df,
          "/Users/mymacstudio/GitDB/CPU/Bachelor/data/GEO/GSE30550/analysis/entropy.csv",
          row.names=TRUE)
# 保存DNB模块的中心基因
DNBs.hub <- c()
for (csv in subject_matrix_path){
  subject <- substr(csv,1,10)
  hub.genes <- hub.genes.matrix[[subject]]
  DNBs.hub <- c(DNBs.hub,hub.genes)
}
DNBs.hub <- unique(DNBs.hub)
# 保存DNB模块中心基因的功能文件
hub.function <- func[func$SYMBOL %in% DNBs.hub,]
write.csv(hub.function[,c("SYMBOL","Description")],
          "/Users/mymacstudio/GitDB/CPU/Bachelor/data/GEO/GSE30550/analysis/hubs_function.csv",
          row.names=FALSE)
# 保存DNB模块中心基因的连接文件
hub.PPI <- PPI[PPI$from_symbol %in% DNBs.hub,]
hub.PPI <- PPI[PPI$to_symbol %in% DNBs.hub,]
write.csv(hub.PPI,
          "/Users/mymacstudio/GitDB/CPU/Bachelor/data/GEO/GSE30550/analysis/hubs_PPI.csv",
          row.names=FALSE)
