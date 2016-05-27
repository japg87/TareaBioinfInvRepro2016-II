###Instalar Bioconductor
##Bajar el script para instalar bioconductor
source("https://bioconductor.org/biocLite.R")
##Correr el script
biocLite()
##Inscalar ComplesHeatmap
biocLite("ComplexHeatmap")

### Ejemplo intento HeatMap complejo
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

lt = readRDS(paste0(system.file(package = "ComplexHeatmap"), "/extdata/meth.rds"))
list2env(lt, envir = environment())

ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", 10), rep("Control", 10))), col = list(type = c("Tumor" = "red", "Control" = "blue")))
ha2 = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", 10), rep("Control", 10))), col = list(type = c("Tumor" = "red", "Control" = "blue")), show_legend = FALSE)

# column order of the methylation matrix which will be assigned to the expressio matrix
column_tree = hclust(dist(t(meth)))

ht_list = 
  Heatmap(meth, name = "methylation", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), cluster_columns = column_tree, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 5, column_title = "Methylation", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10)) +
  Heatmap(direction, name = "direction", col = c("hyper" = "red", "hypo" = "blue"), column_names_gp = gpar(fontsize = 8)) +
  Heatmap(expr[, column_tree$order], name = "expression", col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),cluster_columns = FALSE, top_annotation = ha2, column_names_gp = gpar(fontsize = 8), column_title = "Expression", column_title_gp = gpar(fontsize = 10)) +
  Heatmap(cor_pvalue, name = "-log10(cor_p)", col = colorRamp2(c(0, 2, 4), c("white", "white", "red")), column_names_gp = gpar(fontsize = 8)) +
  Heatmap(gene_type, name = "gene type", col = brewer.pal(length(unique(gene_type)), "Set1"), column_names_gp = gpar(fontsize = 8)) +
  Heatmap(anno, name = "anno_gene", col = brewer.pal(length(unique(anno)), "Set2"), column_names_gp = gpar(fontsize = 8)) +
  Heatmap(dist, name = "dist_tss", col = colorRamp2(c(0, 10000), c("black", "white")), column_names_gp = gpar(fontsize = 8)) +
  Heatmap(enhancer, name = "anno_enhancer", col = colorRamp2(c(0, 1), c("white", "orange")), cluster_columns = FALSE, column_names_gp = gpar(fontsize = 8), column_title = "Enhancer", column_title_gp = gpar(fontsize = 10))

ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 8, fontface = "bold"), heatmap_legend_labels_gp = gpar(fontsize = 8))

draw(ht_list, newpage = FALSE, column_title = "Correspondence between methylation, expression and other genomic features", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
invisible(ht_global_opt(RESET = TRUE))

##########################################################################################
##ahora empezar a modificarlo eliminando las columnas que no necestio para mi trabajo

###################################################################################
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

lt = readRDS(paste0(system.file(package = "ComplexHeatmap"), "/extdata/meth.rds"))
list2env(lt, envir = environment())

ha = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", 10), rep("Control", 10))), col = list(type = c("Tumor" = "red", "Control" = "blue")))
ha2 = HeatmapAnnotation(df = data.frame(type = c(rep("Tumor", 10), rep("Control", 10))), col = list(type = c("Tumor" = "red", "Control" = "blue")), show_legend = FALSE)

# column order of the methylation matrix which will be assigned to the expressio matrix
column_tree = hclust(dist(t(meth)))

ht_list = 
  Heatmap(meth, name = "methylation", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), cluster_columns = column_tree, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 5, column_title = "Methylation", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10)) +
  Heatmap(direction, name = "direction", col = c("hyper" = "red", "hypo" = "blue"), column_names_gp = gpar(fontsize = 8)) +
  Heatmap(expr[, column_tree$order], name = "expression", col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),cluster_columns = FALSE, top_annotation = ha2, column_names_gp = gpar(fontsize = 8), column_title = "Expression", column_title_gp = gpar(fontsize = 10)) +
  Heatmap(cor_pvalue, name = "-log10(cor_p)", col = colorRamp2(c(0, 2, 4), c("white", "white", "red")), column_names_gp = gpar(fontsize = 8)) 
 
ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 8, fontface = "bold"), heatmap_legend_labels_gp = gpar(fontsize = 8))

draw(ht_list, newpage = FALSE, column_title = "Correspondence between methylation, expression", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
invisible(ht_global_opt(RESET = TRUE))

################################################################################################
###Primer intento de mapa sencillo con Expresion 
####################################################################
##primero correr las librerias
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
#hacer la matriz
Expresion <- read.csv("Expresion140516.csv", row.names = 1)
#primer intento
Heatmap(Expresion)
##intente usar una paleta pero son tan pocos los colores que queda casi en blanco todo el mapa
Heatmap(Expresion, col = rev(heat.colors(9)))
###En este arme yo la paleta, aun no se ve del todo bien se ven muy cercanos los colores de los valores mas lejanos
Heatmap(Expresion, col = colorRamp2(c(0, 0.5, 1, 2, 3, 5, 10, 15, 30, 50, 100, 300, 500, 1000, 3000, 6000), c(  "khaki1", "khaki2", "gold", "orange", "orange1", "orange2", "orangered", "orangered1", "orangered2", "orangered3", "red", "red1", "red2", "red3", "red4", "firebrick4")))
### Eliminare al agrupamiento por columna para agrupar solo por patron de metilacion y no por sitio.
Heatmap(Expresion, col = colorRamp2(c(0, 0.5, 1, 2, 3, 5, 10, 15, 30, 50, 100, 300, 500, 1000, 3000, 6000), c(  "khaki1", "khaki2", "gold", "orange", "orange1", "orange2", "orangered", "orangered1", "orangered2", "orangered3", "red", "red1", "red2", "red3", "red4", "firebrick4")), cluster_columns = FALSE)


###############################
#Ahora con el mapa de metilacion
################
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
#hacer la matriz
Methylation <- read.csv("Metilacion040516.csv", row.names =1)
### primer intento
Heatmap(Methylation)
## intento con una paleta e colores
Heatmap(Methylation, col = rev(heat.colors(9)))
## paleta de colores propia
Heatmap(Methylation, col = colorRamp2(c(0, 0.5, 1, 2, 3, 5, 10, 15, 30, 50, 100, 300, 500, 1000, 3000, 6000), c(  "khaki1", "khaki2", "gold", "orange", "orange1", "orange2", "orangered", "orangered1", "orangered2", "orangered3", "red", "red1", "red2", "red3", "red4", "firebrick4")))
### corrigiendo la escala aqui solo es hasta 100
Heatmap(Expresion, col = colorRamp2(c(0, 1, 3, 5, 10, 15, 30, 40, 50, 60, 70, 80, 90, 100), c(  "lightblue", "darkslategray1", "aquamarine", "aquamarine3", "cyan", "cyan3", "deepskyblue1", "deepskyblue4", "dodgerblue1", "dodgerblue3", "blue", "blue3", "navy", "navyblue")))
### Eliminando la agrupacion por columnas
Heatmap(Expresion, col = colorRamp2(c(0, 1, 3, 5, 10, 15, 30, 40, 50, 60, 70, 80, 90, 100), c(  "lightblue", "darkslategray1", "aquamarine", "aquamarine3", "cyan", "cyan3", "deepskyblue1", "deepskyblue4", "dodgerblue1", "dodgerblue3", "blue", "blue3", "navy", "navyblue")), cluster_columns = FALSE)

###################################################################
### Ahora un mapa doble
#####################################
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
#hacer las matrices de ambos mapas
Expresion <- read.csv("Expresion140516.csv", row.names = 1)
Methylation <- read.csv("Metilacion040516.csv", row.names =1)
## hacer el maldito mapa doble o mas bien un heatmap list
f1 = colorRamp2(c(0, 1, 3, 5, 10, 15, 30, 40, 50, 60, 70, 80, 90, 100), c(  "lightblue", "darkslategray1", "aquamarine", "aquamarine3", "cyan", "cyan3", "deepskyblue1", "deepskyblue4", "dodgerblue1", "dodgerblue3", "blue", "blue3", "navy", "navyblue"))
f2 = colorRamp2(c(0, 0.5, 1, 2, 3, 5, 10, 15, 30, 50, 100, 300, 500, 1000, 3000, 6000), c(  "khaki1", "khaki2", "gold", "orange", "orange1", "orange2", "orangered", "orangered1", "orangered2", "orangered3", "red", "red1", "red2", "red3", "red4", "firebrick4"))
Heatmap(Methylation, col = f1, column_title = "Methylation", cluster_columns = FALSE) +
Heatmap(Expresion, col = f2, column_title = "Expresion", cluster_columns = FALSE)


##########################################################################################################
#
# Mapa doble con los datos invertidos columnas por filas y filas por columnas
#
###########################################################################################################

library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
#hacer las matrices de ambos mapas
Expresion <- read.csv("Expresion_160516_v2.csv", row.names = 1)
Methylation <- read.csv("Metilacion_160516_v2.csv", row.names =1)
## hacer el maldito mapa doble o mas bien un heatmap list
f1 = colorRamp2(c(0, 1, 3, 5, 10, 15, 30, 40, 50, 60, 70, 80, 90, 100), c(  "lightblue", "darkslategray1", "aquamarine", "aquamarine3", "cyan", "cyan3", "deepskyblue1", "deepskyblue4", "dodgerblue1", "dodgerblue3", "blue", "blue3", "navy", "navyblue"))
f2 = colorRamp2(c(0, 0.5, 1, 2, 3, 5, 10, 15, 30, 50, 100, 300, 500, 1000, 3000, 6000), c(  "khaki1", "khaki2", "gold", "orange", "orange1", "orange2", "orangered", "orangered1", "orangered2", "orangered3", "red", "red1", "red2", "red3", "red4", "firebrick4"))
Heatmap(Methylation, col = f1, column_title = "Methylation", cluster_rows = FALSE, row_names_gp = gpar(fontsize = 4)) +
  Heatmap(Expresion, col = f2, column_title = "Expresion", cluster_rows = FALSE, row_names_gp = gpar(fontsize = 4))

###############################################################################################
#
#Mapa doble con el primer intento de corelacion ¡¡¡¡Que pinche miedo!!!!!!
#
###############################################################################################
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
#hacer las matrices de ambos mapas
Expresion <- read.csv("Expresion_160516_v2.csv", row.names = 1)
Methylation <- read.csv("Metilacion_160516_v2.csv", row.names =1)
Correlation <- read.csv("Correlacion_160515.csv", row.names =1)
## hacer el maldito mapa doble o mas bien un heatmap list
f1 = colorRamp2(c(0, 1, 3, 5, 10, 15, 30, 40, 50, 60, 70, 80, 90, 100), c(  "lightblue", "darkslategray1", "aquamarine", "aquamarine3", "cyan", "cyan3", "deepskyblue1", "deepskyblue4", "dodgerblue1", "dodgerblue3", "blue", "blue3", "navy", "navyblue"))
f2 = colorRamp2(c(0, 0.5, 1, 2, 3, 5, 10, 15, 30, 50, 100, 300, 500, 1000, 3000, 6000), c(  "khaki1", "khaki2", "gold", "orange", "orange1", "orange2", "orangered", "orangered1", "orangered2", "orangered3", "red", "red1", "red2", "red3", "red4", "firebrick4"))
f3 = colorRamp2(c(-0.5, -0.4, -0.3, -0.2, -0.1, -0.01, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5), c("firebrick4", "red", "salmon1", "pink", "mistyrose1", "mintcream", "mintcream", "lightcyan","lightblue","lightskyblue", "blue", "navyblue" ))
Heatmap(Methylation, name = "Methylation (%)", col = f1, column_title = "Methylation", cluster_rows = FALSE, row_names_gp = gpar(fontsize = 4)) +
  Heatmap(Expresion, name = "Expresion", col = f2, column_title = "Expresion", cluster_rows = FALSE, row_names_gp = gpar(fontsize = 4)) +
    Heatmap(Correlation, name = "Correlation", col = f3, row_names_gp = gpar(fontsize = 0))

#####################
#Una version con los ejes de las x mas pequeños
##################
library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
#hacer las matrices de ambos mapas
Expresion <- read.csv("Expresion_160516_v2.csv", row.names = 1)
Methylation <- read.csv("Metilacion_160516_v2.csv", row.names =1)
Correlation <- read.csv("Correlacion_160515.csv", row.names =1)
## hacer el maldito mapa doble o mas bien un heatmap list
f1 = colorRamp2(c(0, 1, 3, 5, 10, 15, 30, 40, 50, 60, 70, 80, 90, 100), c(  "lightblue", "darkslategray1", "aquamarine", "aquamarine3", "cyan", "cyan3", "deepskyblue1", "deepskyblue4", "dodgerblue1", "dodgerblue3", "blue", "blue3", "navy", "navyblue"))
f2 = colorRamp2(c(0, 0.5, 1, 2, 3, 5, 10, 15, 30, 50, 100, 300, 500, 1000, 3000, 6000), c(  "khaki1", "khaki2", "gold", "orange", "orange1", "orange2", "orangered", "orangered1", "orangered2", "orangered3", "red", "red1", "red2", "red3", "red4", "firebrick4"))
f3 = colorRamp2(c(-0.5, -0.4, -0.3, -0.2, -0.1, -0.01, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5), c("firebrick4", "red", "pink", "mistyrose1", "mintcream", "mintcream", "mintcream", "mintcream", "lightcyan","lightblue", "blue", "navyblue" ))
Heatmap(Methylation, name = "Methylation (%)", col = f1, column_title = "Methylation", cluster_rows = FALSE, row_names_gp = gpar(fontsize = 4), column_names_gp = gpar(fontsize = 7)) +
  Heatmap(Expresion, name = "Expresion", col = f2, column_title = "Expresion", cluster_rows = FALSE, row_names_gp = gpar(fontsize = 4), column_names_gp = gpar(fontsize = 7)) +
    Heatmap(Correlation, name = "Correlation", col = f3, row_names_gp = gpar(fontsize = 0), column_names_gp = gpar(fontsize = 7))


