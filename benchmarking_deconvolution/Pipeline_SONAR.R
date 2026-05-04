library(anndata)
library(SONAR)
library(DescTools)
library(scater)
library(Matrix)
library(data.table)
library(Seurat)
library(matlabr)
library(R.matlab)

code_path <- "./core-code/"

#scrna data
SC <- read_h5ad("/data/ajita/Spatial/Datasets/Spatial_Deconvolution/DLPFC/SC.h5ad")
ref <- t(SC$X)

#annotations
cell_type <- SC$obs$cell_type
cell_type <- factor(cell_type)
CellID <- factor(colnames(ref))
cluster<- data.frame(cellname = CellID, celltype= cell_type)
typeanno <- cluster$celltype
names(typeanno) <- cluster$cellname
typeanno <- as.factor(typeanno)

st_path <- "./simulated_gene_expression_average.h5ad"
ST <- read_h5ad(st_path)

spots <-as.sparse(t(ST$X))
class(spots) <- StripAttr(class(spots))
spots <- round(spots)

row_ind <- which(rowSums(spots) > 0)
col_ind <- which(colSums(spots) > 0)
spots <- spots[row_ind,]
spots <- spots[, col_ind]

pos <- ST$obsm$spatial[col_ind,]

# row.names(pos)<-colnames(st_adata)
colnames(pos) <- c('x','y')
coords <- as.data.frame(pos)
row.names(coords)<-colnames(spots)
#get the overlap genes
overlap_gene <- intersect(rownames(spots), rownames(ref))
ref <- ref[overlap_gene,]

spots <- spots[overlap_gene,]

#calculate the nUMI and nUMI_spot
nUMI <- colSums(ref)
names(nUMI) <- colnames(ref)
nUMI_spot <- colSums(spots)
names(nUMI_spot) <- colnames(spots)

#preprocess the input data
class(ref) <- StripAttr(class(ref))
class(spots) <- StripAttr(class(spots))
processed_input<-SONAR.preprocess(sc_count=ref,sc_cell_type=typeanno,sc_nUMI=nUMI,sp_coords=coords,sp_count=spots,sp_nUMI=nUMI_spot,cores=8,type_min_cell = 0 , spot_min_UMI = 0)

#deliver the preprocessed data to SONAR
trans_data<-SONAR.deliver(processed_data=processed_input,path=code_path)

#define the bandwidth, default is 1.2 times minimal distance
temp<-dist(coords)
temp<-Matrix::Matrix(temp)
temp[temp==0] <- NA
mindist <- min(temp,na.rm = T)
h <- 1.2*mindist

#start deconvolution
SONAR.deconvolute(fname = paste0(code_path,"SONAR_main.m"),path=code_path,h,wait = TRUE)

SONAR.results <- read.csv("final.csv")
u <- fread(paste0(code_path,"u.txt"))
u[,1] <- NULL
colnames(SONAR.results) <- colnames(u)
spot_name <- read.table(file=paste0(code_path,"coord.txt"),sep=",")
rownames(SONAR.results) <- rownames(spot_name)

path <- "output.csv"
final <- SONAR.results[ , sort(colnames(SONAR.results))]
write.csv(final,file=path)


