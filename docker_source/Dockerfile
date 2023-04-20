# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_16

# Install required Bioconductor packages
RUN R -e 'BiocManager::install(c("minfi", "Gviz", "a4Base", "tidyverse", "robustbase", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "UpSetR"))'

# Install github packages (last commit when creating docker)
RUN R -e 'devtools::install_github("perishky/meffil")' ## meffil
# RUN R -e 'BiocManager::install("epimutacions")'


## Solve bug with preprocessCore
RUN R -e 'BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE, upgrade = "never")'

## Install extra packages for plotting
RUN R -e 'BiocManager::install(c("cowplot", "VennDetail"))'

## Add devel version  of epimutacions with quantile reference
RUN R -e 'devtools::install_github("isglobal-brge/epimutacions@d82c9e2", upgrade = "never")' ## epimutacions - commit of Dolors before quantile changes
