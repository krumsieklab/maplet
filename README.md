# **maplet**: **M**etabolomics **A**nalysis **P**ipe**L**inE **T**oolbox

<img src="/Images/maplet_hexagon.png" margin-left = "10" align= "right" width = "90"/>
maplet is an R package for statistical data analysis with a special focus on metabolomics datasets. It allows users to create self-contained analytical pipelines. The toolbox builds upon the bioconductor package <a href="https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html">SummarizedExperiment (SE)</a>, which serves as a central repository for each pipeline’s data, analysis steps, and results. maplet provides a suite of functions for interacting with this container including but not limited to data loading, annotation, statistical analysis, visualization, and reporting. maplet is designed to work with a pipe operator - either the popular %>% operator from the <a href="https://magrittr.tidyverse.org/">magrittr</a> package or the recently introduced |> from base R. This operator allows for smooth connections between pipeline steps, without the need for temporary variables or multiple assignments. The combination of these elements allows for the creation of pipelines which are simple to follow, highly modular, and easily reproducible.

## Reference
Chetnik et al. "maplet: An extensible for modular and reproducible omics pipelines". *Bioinformatics*, 2021. [Link to publication.](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab741/6409851)

## Installation
The latest stable version of maplet can be easily installed using the following command:
```{r}
devtools::install_github(repo="krumsieklab/maplet@v1.2.1", subdir="maplet")
```

To install from the latest commit:
```{r}
devtools::install_github(repo="krumsieklab/maplet", subdir="maplet")
```

**Note:** maplet is in active development. Any commit without a release tag is not guaranteed to be stable. 

## Getting Started
Users should review the available examples in the [examples](/examples) folder. A large example pipeline demonstrating the general setup of a maplet pipeline and most of the functions provided by maplet is provided [here](/examples/use_case_examples/example_preprocessing_pipeline.R). The mt/example folder also contains several stand-alone examples for more specialized functions not included in the main example.  
  
Users are also encouraged to review the [maplet Reference Guide](/guide/maplet_Reference_Guide_markdown.md) sections 1 and 2.

## Case Studies
A list of case studies utilizing maplet can be found [here](CaseStudies.md).

## Want to Get Involved?
Anyone is welcome to contribute to the maplet R package. For users wishing to contribute a new function, we recommend first reviewing the code and documentation of a few existing functions. We also recommend reviewing section 3 of the [maplet Reference Guide](/guide/maplet_Reference_Guide_markdown.md).  

**Note**: maplet functions follow strict naming conventions (refer to section 3.2.1 of the [maplet Reference Guide](/guide/maplet_Reference_Guide_markdown.md) for details). A submitted function may initially be rejected if it does not follow the function naming rules. If you have questions about naming and function development, we recommend reaching out to the maintainer, [Kelsey Chetnik](https://krumsieklab.org/).

## Testing Framework
We have developed a testing framework to ensure maplet functions continue work as expected as the package is updated. The testing framework and documentation can be found in the [tests](/tests) folder.

## Core Contributors
The following people have made significant contributions to the development of maplet:  
Kelsey Chetnik, Elisa Benedetti, Daniel P. Gomari, Annalise Schweickart, Richa Batra, Mustafa Buyukozkan, Zeyu Wang, Matthias Arnold, Jonas Zierer, Karsten Suhre, and Jan Krumsiek.
