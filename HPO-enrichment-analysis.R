# �����
library (org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(fdrtool)

# �����������������
setwd ( "E:/R�γ�" )
MyGeneSet<-read.table('�µ�.txt',header=FALSE)
MyGeneSet$V1<-as.character(MyGeneSet$V1)
MyGeneIDSet2=bitr(MyGeneSet$V1,fromType='SYMBOL',toType='ENTREZID',OrgDb='org.Hs.eg.db')
head(MyGeneIDSet2)

# ����n
n <- length(MyGeneIDSet2$ENTREZID)
n
# ����HPO����
BackgroundData <- read.table('genes_to_phenotypes_to_HPO_diseases.txt',header=TRUE,sep = '\t',
                             nrows = 300000,quote = "")

# ����Mֵ
HPO_ID <- unique(BackgroundData$HPO.Term.ID)
M_value <- data.frame(M_Value = character())
for(value_id in HPO_ID)
{
  allgeneInCategory <- unique(BackgroundData[BackgroundData$HPO.Term.ID == value_id,])
  M <- nrow(allgeneInCategory)
  M_Value = c(M)
  M_Frame <- data.frame(M_Value)
  M_value <- rbind(M_value,M_Frame)
}
sum_data = cbind(HPO_ID,M_value)

# ����Nֵ
BackgroundData_total_gene <- unique(BackgroundData$entrez.gene.id)
N = length(BackgroundData_total_gene)

# ����kֵ
sample_entrez <- unique(MyGeneIDSet2$ENTREZID)
genes <- data.frame(genes = character())
k_value <- data.frame(k_Value = character())
for(value_id in HPO_ID)
{
  allgeneInCategory <- unique(BackgroundData[BackgroundData$HPO.Term.ID== value_id,])
  total_hpo_gene <- allgeneInCategory$entrez.gene.id
  intersect_gene_entrezs <-intersect(sample_entrez,total_hpo_gene)
  rich_gene_symbol = NULL
  if (length(intersect_gene_entrezs) != 0)
  {
    for (intersect_gene_entrez in intersect_gene_entrezs)
      {
      intersect_genes <- MyGeneIDSet2[MyGeneIDSet2$ENTREZID == intersect_gene_entrez,]
      intersect_gene <- intersect_genes$SYMBOL
      rich_gene_symbol<- paste(rich_gene_symbol,intersect_gene,seq = ",")}
  }else
    rich_gene_symbol = " "
  gene = c(rich_gene_symbol)
  gene <- data.frame(gene)
  genes <- rbind(genes,gene)
  k = sum(sample_entrez %in% total_hpo_gene)
  k_Value = c(k)
  k_Frame <- data.frame(k_Value)
  k_value <- rbind(k_value,k_Frame)
}
sum_data <- cbind(sum_data,k_value)
sum_data <- cbind(sum_data,genes)

# ����pֵ
p_value <- data.frame(p_Value = character())
for(value_id in HPO_ID)
{
  pgeneInCategory <- sum_data[sum_data$HPO_ID== value_id,]
  k = pgeneInCategory$k_Value
  M = pgeneInCategory$M_Value
  p = phyper(k-1,M, N-M, n, lower.tail=FALSE)
  p_Value = c(p)
  p_Frame <- data.frame(p_Value)
  p_value <- rbind(p_value,p_Frame)
}
HPO_ID = data.frame(HPO_ID)
HPO_Data = cbind(HPO_ID,p_value)

# ɸѡp<0.05��HPO�Ӽ�
HPO_Data <- unique(HPO_Data[HPO_Data$p_Value<0.05,])

# ����p_adjust
HPO_Data <- HPO_Data[order(as.numeric(as.character(HPO_Data$p_Value))),]
fdr <- fdrtool(HPO_Data$p_Value,statistic = "pvalue")
q_value = fdr$qval
q_value = data.frame(q_value)
HPO_Data =cbind(HPO_Data,q_value)

# �������ݿ򱣴��ļ�,��������
P_HPO_ID <- unique(HPO_Data$HPO_ID)
SUM_HPO <- data.frame(Description = character())
for(value_id in P_HPO_ID)
{
  PgeneInCategory <- unique(BackgroundData[BackgroundData$HPO.Term.ID== value_id,])
  Description <- unique(PgeneInCategory$HPO.Term.Name)
  Description <- data.frame(Description)
  SUM_HPO <- rbind(SUM_HPO,Description)
}

# �������ݿ򱣴��ļ��������ܻ�����Ŀ��pvalue��qvalue��HPO_ID����������
gene_number <- data.frame(Count = character())
gene.name <- data.frame(gene.name = character())
for(value_id in P_HPO_ID)
{
  NUMgeneInCategory <- sum_data[sum_data$HPO_ID== value_id,]
  k = NUMgeneInCategory$k_Value
  Count = c(k)
  Count <- data.frame(Count)
  Gene_names <- NUMgeneInCategory$gene
  Gene_name <- data.frame(Gene_names)
  gene_number <- rbind(gene_number,Count)
  gene.name <- rbind(gene.name,Gene_name)
}
SUM_HPO = cbind(SUM_HPO,gene_number,HPO_Data,gene.name)
SUM_HPO <- SUM_HPO[order(as.numeric(as.character(SUM_HPO$p_Value))),]


# ��������
write.csv(SUM_HPO,"�µ�1.csv",row.names = FALSE)

# ����ͼ����
data <- read.csv("�µ�.csv",header = TRUE)
head(data)
ggplot(data = data,aes(x =Description ,y = Count,fill = -log10(p_Value))) + geom_bar(stat="identity") + 
  scale_x_discrete(limits=data$Description) + coord_flip() + labs(title = "EnrichmentHPO") + 
  theme(plot.title = element_text(size = 20,face = "bold"),axis.text = element_text(size = 12,face = "bold"),
        axis.title.x =element_text(size=14), axis.title.y=element_text(size=16),panel.background 
        = element_rect(fill="white", colour='gray')) + scale_fill_gradient(low = 'blue', high = 'red')