# spDDB
A Comprehensive Benchmarking of Spatial Deconvolution and Domain Detection Methods across Diverse Tissues and Spatial Transcriptomic Technologies

The Github provides installation instructions, set up and runnable files used in the benchmarking. The overview of the benchmarking is as follows.

<img src='overview.png'>

## Installation

We recommend users to directly clone our stable `main` branch and create conda environment of benchmarking methods using yml files provided in ./ENVIRONMENTS/. Executing these yml files will install all required packages and dependencies. Below is an example showing how to create environments for SynthST and a benchmarking method.

```
git clone https://github.com/Zafar-Lab/spDDB.git
cd spDDB/Environments

conda env create -f SynthST.yml
conda activate SynthST

conda env create -f method_name.yml
conda activate method_name

```

## What Computational tasks can spDDB be used for?

`spDDB` can be used for:
1. Benchmarking study of spatial deconvolution methods
2. Benchmarking study of domain detection methods
3. Providing a suite of evaluation metrics for spatial transcriptomics, including bivariate spatial metrics, cell-type shape characterization metrics, and rare cell-type metrics
4. Simulating synthetic spatial transcriptomics data and synthetic cell-type proportions using `SynthST`
5. Rich spatial dataset repository spanning brain, cancer and organs across tissue, species and technologies.

## spDDB Website
The synthetic datasets are available for download from: [https://zafar-lab.github.io/spDDB_datasets.github.io/](https://zafar-lab.github.io/spDDB_datasets.github.io/)

## Tutorials
Check out the following Colab notebook tutorials to run SynthST and spDDB's evaluation framework. 

1. [SynthST for generation of synthetic cell type proportions - DLPFC 151508](https://colab.research.google.com/drive/17QzIYbda3c1HtHDWnyQcpvpABhAO63UZ?usp=sharing)
2. [SynthST for generation of synthetic spatial gene expression - DLPFC 151508](https://colab.research.google.com/drive/1k-QbrHe6Gq_jmBqKSPawF0JhJlkuU54L?usp=sharing)
3. [Generation of datasets using simulation strategy 2 on MERFISH Lung Cancer dataset](https://colab.research.google.com/drive/1Aqnbtlxy1LXFRNpyXgDXF37RyNA0ammG?usp=sharing)
4. [spDDB's Bi-variate Spatial and Non-spatial evaluation metrics - Autogenes method on DLPFC 151508](https://colab.research.google.com/drive/1wCoktLHypk-NUOMBxgQgZ7-PELYDMddB?usp=sharing)
5. [Identification of Regionally Rare and Rare cell types - DLPFC 151508](https://colab.research.google.com/drive/1M27X8mf14t2JaUFp718v8pEXdZLDFHxq?usp=sharing)
6. [Identification of High Curl, High Elongation, Low Clongation, High Linearity and Low Linearity cell types - DLPFC 151508](https://colab.research.google.com/drive/1iuLaS99PXzK5n8lFs9EHPMs1YcY2MTLN?usp=sharing)
   
## Contributing
In case of any bug reports, enhancement requests, general questions, and other contributions, please create an issue. For more substantial contributions, please fork this repo, push your changes to your fork, and submit a pull request with a good commit message.

## Cite this article
Ajita Shree, Aditya V*, Tanush Kumar* and Hamim Zafar, A Comprehensive Benchmarking of Spatial Deconvolution and Domain Detection Methods across Diverse Tissues and Spatial Transcriptomic Technologies, * equal contribution, https://doi.org/10.64898/2026.05.11.724248
