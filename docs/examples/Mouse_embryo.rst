
Mouse embryo scNMT-seq data
=================================


Here we use a multi-omics scNMT-seq dataset to demonstrate the simple usage of MultiSpace. The scNMT-seq data was derived from the study of early mouse embryo (`Argelaguet et al., Nature, 2019 <https://www.nature.com/articles/s41586-019-1825-8>`_). From E4.5 to E7.5 are used in this example. The spatial data was collected from E8.5 mouse embryo (`Lohoff et al., Nature, 2022 <https://www.nature.com/articles/s41587-021-01006-2>`_). 


Using MultiSpace to analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Take E8.5 scNMT-seq data for example. Users should better separate data in absolute folder if they have more than one stage.

Step 1 Run MultiSpace Pipelineinit to initialize snakemake
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The first step of running MultiSpace pipeline is to initialize pipeline config file (also initiation sample in same file)and a working directory. All these steps are implemented by :bash:`MultiSpace Pipelineinit` function. 

.. code:: shell

   MultiSpace Pipelineinit --species mm10 \
   --samplesheet metasheet.csv \
   --directory ~/Project/scNMT_WolfReik/ \
   --fasta ~/Reference/mm10.fa \
   --fasta_fai ~/Reference/mm10.fa.fai \
   --lambda_fasta ~/Reference/mm10_lambda.fa \
   --star_annotation ~/Reference/UCSC_mm10/mm10.refGene.gtf \
   --star_index ~/Reference/UCSC_mm10/


The results of :bash:`MultiSpace Pipelineinit` are shown as below.

+---------------------------------------------------+---------------------------------------------------------------------------+
| File                                              | Description                                                               |
+===================================================+===========================================================================+
| config.yaml                                       | Initialized config.yaml generated in user input directory.                |
+---------------------------------------------------+---------------------------------------------------------------------------+
| Snakefile                                         | Snakefile in directory used for running snakemake.                        |
+---------------------------------------------------+---------------------------------------------------------------------------+
| modules/                                          | Snakemake rules stored in subdir modules of directory.                    |
+---------------------------------------------------+---------------------------------------------------------------------------+


Step 2 Run snakemake to preprocess raw data
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

After initialize config file and working directory, :bash:`snakemake -j 5` could be used to run snakemake in working directory. Users can customize cores number ,but be careful not to make the number too large in case limited memory.

.. code:: shell

   cd ~/Project/scNMT_WolfReik/
   snakemake -j 5



The results of :bash:`snakemake -j 5` are shown as below. DNA methylation(WCG, W=A or T). Chromatin accessibility(GCH, H=A, C or T)

+---------------------------------------------------+---------------------------------------------------------------------------+
| File                                              | Description                                                               |
+===================================================+===========================================================================+
| GCH.bin_peak.h5                                   | GCH bin by cell marix stored in H5 format.                                |
+---------------------------------------------------+---------------------------------------------------------------------------+
| WCG.bin_peak.h5                                   | WCG bin by cell marix stored in H5 format.                                |
+---------------------------------------------------+---------------------------------------------------------------------------+
| GCH.bin.merge.peak                                | GCH bin features. Rowname of GCH.bin_peak.h5                              |
+---------------------------------------------------+---------------------------------------------------------------------------+
| WCG.bin.merge.peak                                | WCG bin features. Rowname of WCG.bin_peak.h5                              |
+---------------------------------------------------+---------------------------------------------------------------------------+
| GCH.chr1_chr6.site_peak.h5                        | GCH site by cell matrix separated by chromatin stored in H5 format.       |
| GCH.chr7_chr12.site_peak.h5                       |                                                                           |
| GCH.chr13_chrY.site_peak.h5                       |                                                                           |
+---------------------------------------------------+---------------------------------------------------------------------------+
| WCG.site_peak.h5                                  | WCG site by cell matrix stored in H5 format.                              |
+---------------------------------------------------+---------------------------------------------------------------------------+
| WCG.uniq.peak                                     | WCG site features. Rowname of WCG.site_peak.h5                            |
+---------------------------------------------------+---------------------------------------------------------------------------+
| GCH.uniq.peak                                     | GCH site features. Rowname of GCH.site_peak.h5(merge chromtain together)  |
+---------------------------------------------------+---------------------------------------------------------------------------+
| QCtable.txt                                       | Snakemake rules stored in subdir modules of directory.                    |
+---------------------------------------------------+---------------------------------------------------------------------------+
| usecell.txt                                       | Snakemake rules stored in subdir modules of directory.                    |
+---------------------------------------------------+---------------------------------------------------------------------------+
| QCtable.txt                                       | Snakemake rules stored in subdir modules of directory.                    |
+---------------------------------------------------+---------------------------------------------------------------------------+



Step 3 Run MultiSpace Scorematrix to calculate WCG/GCH score matrix
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

MultiSpace Scorematrix can calculate genebody or gene promoter methylation ratio using snakemake output file.

.. code:: shell

   MultiSpace Scorematrix --species mm10 --cell_barcode 04.WCG.GCH/usecell.txt \
   --file_path 04.WCG.GCH/ --outdir . --matrixtype WCG --region promoter --distance 2000


The results of :bash:`MultiSpace Scorematrix` are gene by cell matrix stored in TXT format.



MultiSpace Scorematrix can calculate gene activity score using RP(regulatory potential) model.

.. code:: shell

   MultiSpace Scorematrix --species mm10 --cell_barcode 04.WCG.GCH/usecell.txt \
   --file_path 04.WCG.GCH/ --outdir . --matrixtype GCH --distance 10000


The results of :bash:`MultiSpace Scorematrix` are gene by cell matrix stored in TXT format.



Step 4 Run MultiSpace Mappingcell to map single cell to spatial
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

MultiSpace Mappingcell can map single cell to spatial location, and get spatially epigenetic signal.
Users can take :bash:`snakemake` output single cell gene expression matrix, bin by cell matrix and bin features as input.
Additionally, users should offer a spatial gene count matrix and cell type file. The count matrix could be tab-separated plain-text file with genes as rows and spots as columns. The celltype file should be a tab-separated plain-text file without header. The first column should be the cell name, and the second column should be the corresponding celltype labels.

.. code:: shell

   MultiSpace Mappingcell --sc_count_file 05.Spatial/RNA_normalized.txt --sc_celltype_file celltype.txt \
   --st_count_file Spatial/seqFISH_scRNA/RNA_st_normalized.txt --spatial_location Spatial/seFISH_scRNA/loc_EM1.txt \
   --epi_binfile WCG.bin_peak.h5 --epi_feature WCG.bin.merge.peak --out_dir . --out_prefix WCG


Users can use :bash:`MultiSpace Mappingcell --help` to see help message.
The results are showed below.


+---------------------------------------------------+---------------------------------------------------------------------------+
| File                                              | Description                                                               |
+===================================================+===========================================================================+
| WCG.signal_mat.npz                                | DNA methylation signal in spatila location.                               |
|                                                   | Bin feature by spot matrix stored in .npz format.                         |
+---------------------------------------------------+---------------------------------------------------------------------------+
| WCG.signal_mat_rowname.txt                        | Rownames of bin feature by spot matrix after filtering.                   |
|                                                   | Colnames of bin feature by spot matrix is colnames of st_count_file.      |
+---------------------------------------------------+---------------------------------------------------------------------------+


Validate mapping accuracy:

.. image:: ../_static/img/thumbnail/validate.png
   :height: 350px


Mapping E7.5 scNMT-seq data to E8.5 spatial location:

.. image:: ../_static/img/thumbnail/expr_spat.png
   :height: 350px
   :align: center



MultiSpace output file downstream analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users can use :bash:`snakemake` output file to do downstream analysis.

Single omic clustering
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


Using Seurat to cluster RNA gene count matrix by stage and celltype.


Mouse embryo gene count matrix cluster by stage(from E4.5 to E7.5)

.. code:: r

   library(Seurat)
   library(ggplot2)
   library(patchwork)
   library(dplyr)
   library(data.table)
   library(stringr)
   samplemeta = read.table("allsamplemeta.txt",sep = " ", header = T)
   RNA_mat <- as.data.frame(read.table("RNA_normalized.txt",header = T,row.names = 1, check.names=FALSE))

   scseurat <- CreateSeuratObject(
    counts = RNA_mat,
    project = "RNA",
    assay = "RNA",
    min.cells = 5
   )
   scseurat@meta.data$type <- "rna"
   scseurat@meta.data$sample <- rownames(scseurat@meta.data)
   scseurat@meta.data = merge(samplemeta,scseurat@meta.data,on = "sample")
   rownames(scseurat@meta.data) = scseurat@meta.data$sample

   scseurat <- NormalizeData(scseurat) %>% ScaleData() 
   scseurat <- SCTransform(scseurat, assay = "RNA",  verbose = FALSE)
   scseurat <- RunPCA(scseurat, dims = 1:30)
   scseurat <- RunUMAP(scseurat, dims = 1:30)
   scseurat <- FindNeighbors(scseurat, dims = 1:30)
   scseurat <- FindClusters(scseurat, resolution = 0.5, verbose = FALSE)

   DimPlot(scseurat,reduction = "umap",group.by = "stage")


.. image:: ../_static/img/thumbnail/clusterbystage.png
   :height: 350px
   :align: center


.. code:: r

   e75samplemeta = samplemeta[which(samplemeta$stage == "E7.5"),]
   e75RNA_mat = RNA_mat[,which(colnames(RNA_mat) %in% e75samplemeta$sample)]

   e75scseurat <- CreateSeuratObject(
    counts = e75RNA_mat,
    project = "RNA",
    assay = "RNA",
    min.cells = 3
   )
   e75scseurat@meta.data$type <- "rna"
   e75scseurat@meta.data$orig.ident <- "E7.5"
   e75scseurat@meta.data$sample <- rownames(e75scseurat@meta.data)
   e75scseurat@meta.data = merge(samplemeta,e75scseurat@meta.data,on = "sample")
   rownames(e75scseurat@meta.data) = e75scseurat@meta.data$sample

   e75scseurat <- SCTransform(e75scseurat, assay = "RNA",  verbose = FALSE)
   e75scseurat <- RunPCA(e75scseurat, dims = 1:30)
   e75scseurat <- RunUMAP(e75scseurat, dims = 1:30)
   e75scseurat <- FindNeighbors(e75scseurat, dims = 1:30)
   e75scseurat <- FindClusters(e75scseurat, resolution = 0.5, verbose = FALSE)

   DimPlot(e75scseurat,reduction = "umap",group.by = "celltype")


.. image:: ../_static/img/thumbnail/clusterbycelltype.png
   :height: 350px
   :align: center


Using Signac to cluster WCG/GCH bin count matrix by stage.
Take WCG bin matrix for example.

.. code:: r
   library(Signac)
   library(Seurat)
   library(ggplot2)
   library(stringr)
   library(reticulate)
   library(Matrix)
   np <- import("numpy")
   scipy <-import("scipy")

   feature = read.csv("WCG.bin.merge.peak",header = F)
   usecell <- read.table("usecells.txt")
   mydata <- h5read("WCG.bin_peak.h5", "Mcsc")
   WCG_mat = sparseMatrix(x = as.numeric(mydata$data),j = as.numeric(mydata$indices),p = as.numeric(mydata$indptr),dims = c(4114260,985),index1 = FALSE)
   colnames(WCG_mat) = usecell$V1
   rownames(WCG_mat) = feature$V1
   meta <- read.table("sample_metadata.txt",sep = "\t", header = T)
   samplemeta = merge(meta,usecell, by.x = "sample",by.y = "V1")[,c('sample','stage','lineage10x')]
   rownames(samplemeta) = samplemeta$sample
   samplemeta$celltype = samplemeta$lineage10x

   WCG <- CreateSeuratObject(
    counts = WCG_mat,
    assay = "peaks",
    min.cells = 5,
    meta = samplemeta
   )
   WCG@meta.data$celltype = samplemeta$celltype

   WCG <- RunTFIDF(WCG)
   WCG <- FindTopFeatures(WCG, min.cutoff = "q90")
   WCG <- RunSVD(WCG)
   WCG <- RunUMAP(
     object = WCG,
     reduction = 'lsi',
     dims = 2:20
   )

   DimPlot(object = WCG, label = TRUE, reduction = "umap", group.by = "stage")


.. image:: ../_static/img/thumbnail/wcgclusterbystage.png
   :width: 50 %
.. image:: ../_static/img/thumbnail/gchclusterbystage.png
   :width: 50 %


Using Signac to cluster WCG/GCH bin count matrix by celltype.

.. code:: r
   e75samplemeta = samplemeta[which(samplemeta$stage == "E7.5"),]
   WCG_mat = WCG_mat[,e75samplemeta$sample]
   e75WCG <- CreateSeuratObject(
    counts = WCG_mat,
    assay = "peaks",
    min.cells = 3,
    meta = e75samplemeta
   )
   e75WCG@meta.data$celltype = e75samplemeta$celltype

   e75WCG <- RunTFIDF(e75WCG)
   e75WCG <- FindTopFeatures(e75WCG, min.cutoff = "q90")
   e75WCG <- RunSVD(e75WCG)
   e75WCG <- RunUMAP(
     object = e75WCG,
     reduction = 'lsi',
     dims = 2:20
   )

   DimPlot(object = e75WCG, label = TRUE, reduction = "umap", group.by = "celltype")

.. image:: ../_static/img/thumbnail/wcgclusterbycelltype.png
   :width: 50%
.. image:: ../_static/img/thumbnail/gchclusterbycelltype.png
   :width: 50%



Multi omics clustering
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code:: r




Spatial multi-omics analysis
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

.. code:: python



