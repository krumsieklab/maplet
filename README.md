# **maplet**: **M**etabolomics **A**nalysis **P**ipe**L**inE **T**oolbox

<img src="/Images/maplet_hexagon.png" margin-left = "10" align= "right" width = "90"/>
maplet is an R package for statistical data analysis with a special focus on metabolomics datasets. It allows users to create self-contained analytical pipelines. The toolbox builds upon the bioconductor package <a href="https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html">SummarizedExperiment (SE)</a>, which serves as a central repository for each pipeline’s data, analysis steps, and results. maplet provides a suite of functions for interacting with this container including but not limited to data loading, annotation, statistical analysis, visualization, and reporting. maplet is designed to work with the pipe operator (%>%) from the <a href="https://magrittr.tidyverse.org/">magrittr</a> package. This operator allows for smooth connections between pipeline steps, without the need for temporary variables or multiple assignments. The combination of these elements allows for the creation of pipelines which are simple to follow, highly modular, and easily reproducible.

## Installation
The latest stable version of maplet can be easily installed using the following command:
```{r}
devtools::install_github(repo="krumsieklab/maplet@v1.0.0", subdir="maplet")
```

To install from the latest commit:
```{r}
devtools::install_github(repo="krumsieklab/maplet", subdir="maplet")
```

**Note:** maplet is in active development. Any commit without a release tag is not guaranteed to be stable. 

## Getting Started
Users should review the available examples in the [maplet/examples](/examples) folder. A large example pipeline demonstrating the general setup of a maplet pipeline and most of the functions provided by maplet is provided [here](/examples/example_pipeline.R). The mt/example folder also contains several stand-alone examples for more specialized functions not included in the main example.  
  
Users are also encouraged to review the maplet Reference Guide sections 1 and 2.

## Want to Get Involved?
Anyone is welcome to contribute to the maplet R package. For users wishing to contribute a new function, we recommend first reviewing the code and documentation of a few existing functions. We also recommend reviewing section 3 of the maplet Reference Guide.  

**Note**: maplet functions follow strict naming conventions (refer to section 3.2.1 of the maplet Reference Guide for details). A submitted function may initially be rejected if it does not follow the function naming rules.If you have questions about naming and function development, we recommend reaching out to the maintainer, [Kelsey Chetnik](https://krumsieklab.org/).

## Core Contributors
The following people have made significant contributions to the development of maplet:  
Kelsey Chetnik, Jan Krumsiek, Elisa Benedetti, Richa Batra, Daniel Gomari, Mustafa Buyukozkan, Annalise Schweickart, Zeyu Wang,
Karsten Suhre, Matthias Arnold, and Jonas Zierer.
