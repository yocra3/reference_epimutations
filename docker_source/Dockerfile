# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_16

# Install required Bioconductor packages
RUN R -e 'BiocManager::install(c("minfi", "Gviz", "a4Base", "tidyverse", "robustbase", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "UpSetR"))'

# Install github packages (last commit when creating docker)
RUN R -e 'devtools::install_github("perishky/meffil")' ## meffil
