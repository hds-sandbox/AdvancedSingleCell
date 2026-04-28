install.packages('BiocManager', repos='http://cran.us.r-project.org', lib="ENVIRONMENT_FOLDER/lib/R/library/" )

remotes::install_github("mojaveazure/seurat-disk", upgrade="never", lib="ENVIRONMENT_FOLDER/lib/R/library/")

remotes::install_github("SamueleSoraggi/DoubletFinder", lib="ENVIRONMENT_FOLDER/lib/R/library/", upgrade="never")

BiocManager::install(update=FALSE, lib="ENVIRONMENT_FOLDER/lib/R/library/")

BiocManager::install(c("WGCNA", "igraph", "devtools", "GeneOverlap", 'ggrepel'),
                    update=FALSE, lib="ENVIRONMENT_FOLDER/lib/R/library/")

remotes::install_github("NightingaleHealth/ggforestplot", upgrade="never", lib="ENVIRONMENT_FOLDER/lib/R/library/")

remotes::install_github('smorabit/hdWGCNA', ref='dev', upgrade="never", lib="ENVIRONMENT_FOLDER/lib/R/library/")