library(ggpubr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
setwd("E:/workfir/KEGG5")
df<-read.table("452_kegg.out_01",sep = "\t",header = T,row.names=1)
info<-read.table("all_452_mag_info.txt",sep = "\t",header = T)
ko<-read.table("ko_info",sep="\t")
SGB_info<-read.table("SGB_info",sep="\t",header = T)

df<-data.frame(t(df))
sub_230_df<-df[row.names(df) %in% SGB_info[,1],]
sub_230_df<-sub_230_df[,colnames(df) %in% ko[,1]]
sub_230_df_plot<-sub_230_df
sub_230_df_plot_ant<-data.frame(Bin=row.names(sub_230_df_plot),rawrank=c(1:nrow(sub_230_df_plot)))
sub_230_df_plot_ant_info<-merge(sub_230_df_plot_ant,SGB_info)
#here rank by greup
sub_230_df_plot_ant_info<-sub_230_df_plot_ant_info[order(sub_230_df_plot_ant_info$Grp),]

sub_230_df_plot_ant_info$grprank<-1:nrow(sub_230_df_plot_ant_info)
sub_230_df_plot_ant_info<-sub_230_df_plot_ant_info[order(sub_230_df_plot_ant_info$rawrank),]


#order for plot 
sub_230_df_plot<-sub_230_df_plot[sub_230_df_plot_ant_info$grprank,]
sub_230_df_plot_ant_info_order<-sub_230_df_plot_ant_info[order(sub_230_df_plot_ant_info$grprank),]
#plot
ha_row = rowAnnotation(df = data.frame(Abundance =sub_230_df_plot_ant_info_order$Grp),
    col = list(Abundance = c("Core" =  "#edf8b1", "M" = "#7fcdbb","Rare"="#2c7fb8")), width = unit(0.5, "cm"))

#`````````````````````````````````````````````````````````````````top 注释
colant<-data.frame(ko=colnames(sub_230_df_plot))
colant$rank<-1:nrow(colant)
colant<-merge(colant,ko,by.x="ko",by.y="V1")
colant<-colant[order(colant$V2),]
colant$frank<-1:nrow(colant)
# colant<-colant[order(colant$rank),]

sub_230_df_plot<-sub_230_df_plot[,as.numeric(colant$rank)]


 top_annotation = HeatmapAnnotation(Pathway = colant$V2)

h1=Heatmap(sub_230_df_plot, name = "Gene number",col = colorRamp2(c(0,1),c( "white","orange")),cluster_rows = F,show_column_names=T,cluster_columns=F,rect_gp = gpar(col = "gray"), row_names_gp = gpar(fontsize = 1) ,column_names_gp = gpar(fontsize = 2),column_names_side = "top",top_annotation = HeatmapAnnotation(Pathway = colant$V2)) 

Heatmap(sub_230_df_plot, name = "Gene number",col = colorRamp2(c(0,1),c( "white","orange")),cluster_rows = F,show_column_names=T,cluster_columns=F,rect_gp = gpar(col = "gray"),row_names_gp = gpar(fontsize = 0.5)) 

h1+ha_row


#不聚类
sub_230_df_plot.s<-sub_230_df_plot[order(sub_230_df_plot_rawant$Grp),]

sub_230_df_plot_rawant.s<-sub_230_df_plot_rawant[order(sub_230_df_plot_rawant$Grp),]
sub_230_df_plot_rawant.s$grprank<-1:nrow(sub_230_df_plot_rawant.s)


ha_row = rowAnnotation(df = data.frame(Abundance =sub_230_df_plot_rawant.s$Grp),
    col = list(Abundance = c("Core" =  "#edf8b1", "M" = "#7fcdbb","Rare"="#2c7fb8")), width = unit(0.5, "cm"))




h1=Heatmap(sub_230_df_plot.s, name = "Gene number",col = colorRamp2(c(0,1),c( "white","orange")),cluster_rows = F,show_column_names=F,cluster_columns=F,rect_gp = gpar(col = "gray"),row_order=)
h1+ha_row

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``new 
library(ggpubr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
setwd("E:/workfir/KEGG5")
df<-read.table("230_SGB_kegg.out",sep = "\t",header = T,row.names=1)
#info<-read.table("all_452_mag_info.txt",sep = "\t",header = T)
ko<-read.table("ko_sugar",sep="\t") #ko_others
SGB_info<-read.table("SGB_info",sep="\t",header = T)

# df<-data.frame(t(df))
sub_230_df<-df[row.names(df) %in% SGB_info[,1],]
sub_230_df<-sub_230_df[,colnames(df) %in% ko[,1]]
sub_230_df_plot<-sub_230_df
sub_230_df_plot_ant<-data.frame(Bin=row.names(sub_230_df_plot),rawrank=c(1:nrow(sub_230_df_plot)))
sub_230_df_plot_ant_info<-merge(sub_230_df_plot_ant,SGB_info,by.x="Bin",by.y="Cluster")
#here rank by greup
sub_230_df_plot_ant_info<-sub_230_df_plot_ant_info[order(sub_230_df_plot_ant_info$Grp),]

sub_230_df_plot_ant_info$grprank<-1:nrow(sub_230_df_plot_ant_info)
sub_230_df_plot_ant_info<-sub_230_df_plot_ant_info[order(sub_230_df_plot_ant_info$rawrank),]


#order for plot 
sub_230_df_plot<-sub_230_df
sub_230_df_plot<-sub_230_df_plot[order(sub_230_df_plot_ant_info$grprank),]
sub_230_df_plot_ant_info_order<-sub_230_df_plot_ant_info[order(sub_230_df_plot_ant_info$grprank),]
#plot
ha_row = rowAnnotation(df = data.frame(Abundance =sub_230_df_plot_ant_info_order$Grp),
    col = list(Abundance = c("H" =  "#edf8b1", "M" = "#7fcdbb","Rare"="#2c7fb8")), width = unit(0.5, "cm"))

#`````````````````````````````````````````````````````````````````top 注释
colant<-data.frame(ko=colnames(sub_230_df_plot))
colant$rank<-1:nrow(colant)
colant<-merge(colant,ko,by.x="ko",by.y="V1")
colant<-colant[order(colant$V2),]
colant$frank<-1:nrow(colant)
# colant<-colant[order(colant$rank),]

sub_230_df_plot<-sub_230_df_plot[,as.numeric(colant$rank)]


 # top_annotation = HeatmapAnnotation(Pathway = colant$V2)

h1=Heatmap(sub_230_df_plot, name = "Gene number",col = colorRamp2(c(0,1,10,100),c( "white","orange","#d73027","#762a83")),cluster_rows = T,show_column_names=T,cluster_columns=F,rect_gp = gpar(col = "gray",lwd = 0.01), row_names_gp = gpar(fontsize = 1) ,column_names_gp = gpar(fontsize = 2),column_names_side = "top",top_annotation = HeatmapAnnotation(Pathway = colant$V2),column_split =colant$V2 ,show_row_names=T)

# Heatmap(sub_230_df_plot, name = "Gene number",col = colorRamp2(c(0,1),c( "white","orange")),cluster_rows = F,show_column_names=T,cluster_columns=F,rect_gp = gpar(col = "gray"),row_names_gp = gpar(fontsize = 0.5)) 

h1+ha_row


#不聚类
sub_230_df_plot.s<-sub_230_df_plot[order(sub_230_df_plot_rawant$Grp),]

sub_230_df_plot_rawant.s<-sub_230_df_plot_rawant[order(sub_230_df_plot_rawant$Grp),]
sub_230_df_plot_rawant.s$grprank<-1:nrow(sub_230_df_plot_rawant.s)


ha_row = rowAnnotation(df = data.frame(Abundance =sub_230_df_plot_rawant.s$Grp),
    col = list(Abundance = c("Core" =  "#edf8b1", "M" = "#7fcdbb","Rare"="#2c7fb8")), width = unit(0.5, "cm"))




h1=Heatmap(sub_230_df_plot.s, name = "Gene number",col = colorRamp2(c(0,1),c( "white","orange")),cluster_rows = F,show_column_names=F,cluster_columns=F,rect_gp = gpar(col = "gray"),row_order=)
h1+ha_row




#---------------------------------------------------------------------------------总
df<-read.table("clipboard",sep = "\t",header = T,row.names=1)
ko<-read.table("clipboard",sep="\t")
View(ko)
SGB_info<-read.table("clipboard",sep="\t",header = T)
View(SGB_info)
sub_230_df<-df[row.names(df) %in% SGB_info[,1],]
sub_230_df<-sub_230_df[,colnames(df) %in% ko[,1]]
View(sub_230_df)
sub_230_df_plot<-sub_230_df
sub_230_df_plot_ant<-data.frame(Bin=row.names(sub_230_df_plot),rawrank=c(1:nrow(sub_230_df_plot)))
View(sub_230_df_plot_ant)
sub_230_df_plot_ant_info<-merge(sub_230_df_plot_ant,SGB_info,by.x="Bin",by.y="sample")
View(sub_230_df_plot_ant_info)
sub_230_df_plot_ant_info$grprank<-1:nrow(sub_230_df_plot_ant_info)
sub_230_df_plot<-sub_230_df
View(sub_230_df_plot)
colant<-data.frame(ko=colnames(sub_230_df_plot))
colant$rank<-1:nrow(colant)
colant<-merge(colant,ko,by.x="ko",by.y="V1")
colant<-colant[order(colant$V2),]
colant$frank<-1:nrow(colant)
View(colant)
sub_230_df_plot<-sub_230_df_plot[,as.numeric(colant$rank)]
h1=Heatmap(sub_230_df_plot, name = "Gene number",col = colorRamp2(c(0,1,10,100),c( "white","orange","#d73027","#762a83")),cluster_rows = T,show_column_names=T,cluster_columns=F,rect_gp = gpar(col = "gray",lwd = 0.01), row_names_gp = gpar(fontsize = 1) ,column_names_gp = gpar(fontsize = 2),column_names_side = "top",top_annotation = HeatmapAnnotation(Pathway = colant$V2),column_split =colant$V2,row_split = sub_230_df_plot_ant_info$Grp2,show_row_names=T)

ha_row = rowAnnotation(df = data.frame(Abundance =sub_230_df_plot_ant_info$Grp),
+col = list(Abundance = c("H" =  "#edf8b1", "Mid" = "#7fcdbb","Low"="#2c7fb8")), width = unit(0.5, "cm"))
h1+ha_row