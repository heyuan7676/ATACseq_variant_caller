


GEO_dat = read_csv('data/SraRunTable_ATAC_paired.csv')
GEO_dat = GEO_dat %>% select(c(Run, BioProject, BioSample, source_name, Cell_type, TISSUE, Cell_Line))
dim(GEO_dat)

#GEO_cell_tissue = GEO_dat[, c("Cell_type", "TISSUE", "Cell_Line")]
#GEO_cell_tissue = apply(as.data.frame(GEO_cell_tissue), 1, function(x) paste(x[complete.cases(x)], sep=','))
#GEO_cell_tissue = sapply(GEO_cell_tissue, function(x) paste(x[1], sep=','))
#table(GEO_cell_tissue)[table(GEO_cell_tissue) > 10]

GEO_sample_source = GEO_dat$source_name
length(unique(GEO_sample_source))
length(table(GEO_sample_source)[table(GEO_sample_source) >= 10])
sum(table(GEO_sample_source)[table(GEO_sample_source) >= 10])

table(GEO_sample_source)[table(GEO_sample_source) >= 10]

GEO_dat %>% filter(source_name == 'lung') %>% select(c(Run, source_name, Cell_type, TISSUE, Cell_Line))




### Haematopoiesis

GEO_dat = read.table('data/SraRunTable_Haematopoiesis.csv', sep='\t', header = T)
GEO_dat = as_tibble(GEO_dat)
GEO_dat = GEO_dat %>% select(c(Run, BioProject, BioSample, source_name, Cell_type, TISSUE, Cell_Line, differentiation_stage))
dim(GEO_dat)

length(unique(GEO_dat$BioProject))
length(unique(GEO_dat$BioSample))
length(unique(GEO_dat$Run))
length(unique(GEO_dat$Cell_type))

cell_type_count = GEO_dat %>% group_by(Cell_type) %>% count
g = ggplot(data = cell_type_count) + 
  geom_bar(aes(x = reorder(Cell_type,n), y = n), stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) + 
  xlab("Cell type") + 
  ylab("Number of samples")
ggsave(paste0('results/SraRunTable_Haematopoiesis.png'), height = 4, width = 10)

GEO_dat$Cell_type_category = "NA"


# leukemia
token = GEO_dat[grep("leukemia", GEO_dat$Cell_type, ignore.case = T), ]
GEO_dat[GEO_dat$Cell_type == "TALL-1,T cell leukemia cell line", "Cell_type_category"] = "leukemia cell line - TALL"
GEO_dat[GEO_dat$Cell_type == "acute myeloid leukemia cell line", "Cell_type_category"] = "leukemia cell line - AML2"

unique(GEO_dat$Cell_type[grep("CTV", GEO_dat$Cell_type, ignore.case = T)])
GEO_dat[GEO_dat$Cell_type == "CTV-1", "Cell_type_category"] = "leukemia cell line - CTV"
GEO_dat[GEO_dat$Cell_type == "CTV-1,T cell leukemia cell line", "Cell_type_category"] = "leukemia cell line - CTV"

unique(GEO_dat$Cell_type[grep("CLL", GEO_dat$Cell_type, ignore.case = T)])
GEO_dat[GEO_dat$Cell_type == "CLL", "Cell_type_category"] = "leukemia cell line - CLL"

unique(GEO_dat$source_name[grep("Jurkat", GEO_dat$source_name, ignore.case = T)])
GEO_dat[GEO_dat$source_name == 'Jurkat T cells', "Cell_type_category"] = "leukemia cell line - Jurkat"
GEO_dat[GEO_dat$source_name == 'Jurkat', "Cell_type_category"] = "leukemia cell line - Jurkat"
GEO_dat[GEO_dat$source_name == 'Jurkat cells', "Cell_type_category"] = "leukemia cell line - Jurkat"
GEO_dat[GEO_dat$source_name == 'K562', "Cell_type_category"] = "leukemia cell line - K562"




# PBMC
GEO_dat[GEO_dat$source_name == "Peripheral blood mononuclear cells", "Cell_type_category"] = "PBMC"



# HSPCs
unique(GEO_dat$Cell_type[grep("HSPC", GEO_dat$Cell_type, ignore.case = T)])
GEO_dat[GEO_dat$Cell_type == "HSPCs", "Cell_type_category"] = "primary hematopoietic cells"

# hematopoietic
unique(GEO_dat$Cell_type[grep("hematopoietic", GEO_dat$Cell_type, ignore.case = T)])
GEO_dat[GEO_dat$Cell_type == "hematopoietic cell line", "Cell_type_category"] = "hematopoietic cell line"

unique(GEO_dat$source_name[grep("hematopoietic", GEO_dat$source_name, ignore.case = T)])
for(alias in unique(GEO_dat$source_name[grep("hematopoietic", GEO_dat$source_name)])){
  GEO_dat[GEO_dat$source_name == alias, "Cell_type_category"] = "primary hematopoietic cells"
}

GEO_dat[GEO_dat$differentiation_stage == "Day 18 of hematopoietic differentiation", "Cell_type_category"] = "iPSC differentiated hematopoietic"



# LCL cell line
unique(GEO_dat$Cell_type[grep("lymphoblastoid", GEO_dat$Cell_type, ignore.case = T)])
unique(GEO_dat$TISSUE[grep("Lymphoblastoid", GEO_dat$TISSUE, ignore.case = T)])

for(alias in unique(GEO_dat$TISSUE[grep("Lymphoblastoid", GEO_dat$TISSUE)])){
  GEO_dat[GEO_dat$TISSUE == alias, "Cell_type_category"] = "LCL"
}
for(alias in unique(GEO_dat$Cell_type[grep("lymphoblastoid", GEO_dat$Cell_type)])){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "LCL"
}

# other cell lines
GEO_dat[GEO_dat$Cell_Line == "Ramos B", "Cell_type_category"] = "Ramos B"

# NK cells
for(alias in unique(GEO_dat$Cell_type[grep("natural killer", GEO_dat$Cell_type, ignore.case = T)])){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "NK cell"
}
for(alias in unique(GEO_dat$Cell_type[grep("NK", GEO_dat$Cell_type, ignore.case = T)])){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "NK cell"
}

for(alias in unique(GEO_dat$TISSUE[grep("natural killer", GEO_dat$TISSUE, ignore.case = T)])){
  GEO_dat[GEO_dat$TISSUE == alias, "Cell_type_category"] = "NK cell"
}


# B cell
unique(GEO_dat$Cell_type[grep("B", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("B cells", "B lymphocyte", "B-Lymphocyte and Promyeloblast", 
               "Bulk_B", "Naive_B", "Circulating B cells", "Mem_B")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "B cell"
}



# T cell
unique(GEO_dat$Cell_type[grep("T cell", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("Sorted primary Naive T cells", "Naive CD4 T cells")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "T cell"
}


# CD8 T cell
for(alias in c("sorted CD45RA+ CCR7+ CD8 T cells", "sorted A2-NS4B214 tetramer+ CD8 T cells", 
               "CD8pos_T", "Naive_CD8_T")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "CD8 T cell"
}


# CD4 T cells
for(alias in c("Naive CD4+ T cell", "Naive CD4 T cells", "Total CD4 T cells",  "Primary CD4+ T-cells")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "CD4 T cell"
}
GEO_dat[GEO_dat$source_name == "naïve CD4+ T-cell", "Cell_type_category"] = "CD4 T cell"


# Naive T Regulator
unique(GEO_dat$Cell_type[grep("reg", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("Regulatory_T", "Naive_Tregs", "TREG CD4+ T cell", 
               "Sorted primary Treg cells", "Sorted primary Treg cells TREG CD4+ T cell ")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Tregs"
}

# Naive T Effector
unique(GEO_dat$Cell_type[grep("eff", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("Effector_CD4pos_T", "Naive_Teffs")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Teffs"
}

# Memory cells
unique(GEO_dat$Cell_type[grep("memory", GEO_dat$Cell_type, ignore.case = T)])
for(alias in unique(GEO_dat$Cell_type[grep("memory", GEO_dat$Cell_type, ignore.case = T)])){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = alias
}


# Th1 
unique(GEO_dat$Cell_type[grep("TH1", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("Th1_precursors", "TH1 CD4+ T cell", "TH1-17 CD4+ T cell")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Th1 T cell"
}


# Th2
unique(GEO_dat$Cell_type[grep("TH2", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("Th2_precursors", "TH2 CD4+ T cell")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Th2 T cell"
}

# Th17
unique(GEO_dat$Cell_type[grep("TH17", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("Th17_precursors", "TH17 CD4+ T cell", "Sorted primary Th17 cells")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Th17 T cell"
}

# Follicular
unique(GEO_dat$Cell_type[grep("Follicular", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("Follicular_T_Helper")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Follicular_T_Helper"
}

# TEMRA
for(alias in c("TEMRA CD4 T cells")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "TEMRA CD4 T cells"
}

# Gamma_delta_T
for(alias in c("Gamma_delta_T")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Gamma_delta_T"
}



## Myeloid
unique(GEO_dat$Cell_type[grep("myeloid", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("Myeloid_DCs")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Myeloid_DCs"
}

# Plasmacytoid
unique(GEO_dat$Cell_type[grep("pDC", GEO_dat$Cell_type, ignore.case = T)])
for(alias in c("Plasmablasts", "pDC", "pDCs")){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Plasmacytoid dendritic cells"
}

unique(GEO_dat$Cell_type[grep("cDC", GEO_dat$Cell_type, ignore.case = T)])
for(alias in unique(GEO_dat$Cell_type[grep("cDC", GEO_dat$Cell_type, ignore.case = T)])){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "circulating dendritic cells"
}

# Monocyte
unique(GEO_dat$Cell_type[grep("mono", GEO_dat$Cell_type, ignore.case = T)])
for(alias in unique(GEO_dat$Cell_type[grep("mono", GEO_dat$Cell_type, ignore.case = T)])){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Monocyte"
}

GEO_dat[GEO_dat$source_name == "THP-1 Monocyte", "Cell_type_category"] = "Monocyte THP-1"

# granulocyte/macrophage progenitor
unique(GEO_dat$Cell_type[grep("GMP", GEO_dat$Cell_type, ignore.case = T)])
for(alias in unique(GEO_dat$Cell_type[grep("GMP", GEO_dat$Cell_type, ignore.case = T)])){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "GMP"
}

# normal density neutrophils
GEO_dat[GEO_dat$Cell_type == "NDN", "Cell_type_category"] = "normal density neutrophils"

# low-density granulocytes
GEO_dat[GEO_dat$Cell_type == "LDG", "Cell_type_category"] = "low-density granulocytes"


## Mega
unique(GEO_dat$Cell_type[grep("Mega", GEO_dat$Cell_type, ignore.case = T)])
for(alias in unique(GEO_dat$Cell_type[grep("Mega", GEO_dat$Cell_type, ignore.case = T)])){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Mega"
}

# Erythroid
unique(GEO_dat$Cell_type[grep("Erythroid", GEO_dat$Cell_type, ignore.case = T)])
for(alias in unique(GEO_dat$Cell_type[grep("Erythroid", GEO_dat$Cell_type, ignore.case = T)])){
  GEO_dat[GEO_dat$Cell_type == alias, "Cell_type_category"] = "Erythroid"
}

# Umbilical cord blood
GEO_dat[GEO_dat$source_name == "Umbilical cord blood", "Cell_type_category"] = "Umbilical cord blood"
GEO_dat[GEO_dat$source_name == "Human cord blood CD34+", "Cell_type_category"] = "Umbilical cord blood"

# CAR T cell
GEO_dat[GEO_dat$Cell_type == "CAR T cell", "Cell_type_category"] = "CAR T cell"

unique(GEO_dat[GEO_dat$Cell_type_category == 'NA', "Cell_type"])

as.data.frame(GEO_dat %>% subset(Cell_type_category == "NA") %>% subset(source_name != 'CD4 T cells') )


cell_type_count = GEO_dat %>% group_by(Cell_type_category) %>% count
cell_type_count$leukemia = "False"
cell_type_count[grep("leukemia", cell_type_count$Cell_type_category), "leukemia"] = "True"

g = ggplot(data = cell_type_count) + 
  geom_bar(aes(x = reorder(Cell_type_category,n), y = n, fill = leukemia), stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) + 
  scale_fill_manual(values=c("black", "red")) + 
  xlab("Cell type") + 
  ylab("Number of samples")
ggsave(paste0('results/SraRunTable_Haematopoiesis_Cell_type_category.png'), height = 4, width = 8)



## Large branch

GEO_dat$higher_branch = GEO_dat$Cell_type_category
for(alias in c("Teffs", "Follicular_T_Helper", "Memory_Teffs", "Memory_Tregs", "Tregs", "Th1 T cell", 
               "Th17 T cell", "Th2 T cell", "CD4 T cell", "TEMRA CD4 T cells", "Central memory CD4 T cells", "Effector memory CD4 T cells")){
  GEO_dat[GEO_dat$Cell_type_category == alias, "higher_branch"] = "CD4 T cell"
}

for(alias in c("Central_memory_CD8pos_T", "Effector_memory_CD8pos_T", "CD8 T cell")){
  GEO_dat[GEO_dat$Cell_type_category == alias, "higher_branch"] = "CD8 T cell"
}

for(alias in c("NK cell", "Memory_NK")){
  GEO_dat[GEO_dat$Cell_type_category == alias, "higher_branch"] = "NK cell"
}

for(alias in c("Monocyte THP-1", "Monocyte", "Monocytes")){
  GEO_dat[GEO_dat$Cell_type_category == alias, "higher_branch"] = "Monocyte"
}

for(alias in unique(GEO_dat$Cell_type_category[grep("leukemia", GEO_dat$Cell_type_category)])){
  GEO_dat[GEO_dat$Cell_type_category == alias, "higher_branch"] = "leukemia cell line"
}

for(alias in c("Myeloid_DCs", "Plasmacytoid dendritic cells", "circulating dendritic cells")){
  GEO_dat[GEO_dat$Cell_type_category == alias, "higher_branch"] = "dendritic cells"
}


cell_type_count = GEO_dat %>% group_by(higher_branch) %>% count
cell_type_count$leukemia = "False"
cell_type_count[grep("leukemia", cell_type_count$higher_branch), "leukemia"] = "True"

g = ggplot(data = cell_type_count) + 
  geom_bar(aes(x = reorder(higher_branch,n), y = n, fill = leukemia), stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) + 
  scale_fill_manual(values=c("black", "red")) + 
  xlab("Cell type") + 
  ylab("Number of samples")
ggsave(paste0('results/SraRunTable_Haematopoiesis_Cell_type_category_higher.png'), height = 4, width = 6)


## Even higher

GEO_dat$top_category = GEO_dat$higher_branch
for(alias in c("circulating dendritic cells", "GMP", "low-density granulocytes", "Monocyte","dendritic cells",
               "normal density neutrophils", "low-density granulocytes")){
  GEO_dat[GEO_dat$Cell_type_category == alias, "top_category"] = "myeloid cells"
}

cell_type_count = GEO_dat %>% group_by(top_category) %>% count
cell_type_count$leukemia = "False"
cell_type_count[grep("leukemia", cell_type_count$top_category), "leukemia"] = "True"

g = ggplot(data = cell_type_count) + 
  geom_bar(aes(x = reorder(top_category,n), y = n, fill=leukemia), stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5)) + 
  scale_fill_manual(values=c("black", "red")) + 
  xlab("Cell type") + 
  ylab("Number of samples")
ggsave(paste0('results/SraRunTable_Haematopoiesis_Cell_type_category_top.png'), height = 4, width = 6)

write.table(GEO_dat, 'data/SraRunTable_Haematopoiesis_categorized.csv', sep=',', row.names=F, quote = F)
