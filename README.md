# REMIND-Cancer Filtering Pipeline

## Beyond Recurrence: A Novel Workflow to Identify Activating Promoter Mutations in Cancer Genomes
![GitHub Release](https://img.shields.io/github/v/release/nicholas-abad/remind-cancer) ![GitHub last commit](https://img.shields.io/github/last-commit/nicholas-abad/remind-cancer) ![GitHub top language](https://img.shields.io/github/languages/top/nicholas-abad/remind-cancer) ![GitHub repo size](https://img.shields.io/github/repo-size/nicholas-abad/remind-cancer)

Authors: [Nicholas Abad](https://github.com/nicholas-abad)<sup>1,2</sup>, Irina Glas<sup>1,3</sup>, [Chen Hong](https://github.com/chenhong-dkfz)<sup>1,4</sup>, Annika Small<sup>3</sup>, [Yoann Pageaud](https://github.com/YoannPa)<sup>1,5</sup>, Ana Maia<sup>3</sup>, Dieter Weichenhan<sup>6</sup>, Christoph Plass<sup>6</sup>, Barbara Hutter<sup>7</sup>, Benedikt Brors<sup>1,8,9,10</sup>, Cindy Körner<sup>3</sup>, Lars Feuerbach<sup>1</sup>
<br>
![Github Number of Dependencies](https://img.shields.io/badge/repository_maintainer-nicholas_abad_(nicholas.a.abad@gmail.com)-lightgreen)

<sup>1 Division of Applied Bioinformatics, German Cancer Research Center (DKFZ), Heidelberg, Germany</sup><br>
<sup>2 Faculty of Engineering Sciences, Heidelberg University, Heidelberg, Germany </sup><br>
<sup>3 Division of Molecular Genome Analysis, German Cancer Research Center (DKFZ), Heidelberg, Germany</sup><br>
<sup>4 Division of Molecular Genetics, German Cancer Research Center (DKFZ), Heidelberg, Germany.</sup><br>
<sup>5 Faculty of Biosciences, Heidelberg University, Heidelberg, Germany</sup><br>
<sup>6 Division of Cancer Epigenomics, German Cancer Research Center (DKFZ), Heidelberg, Germany</sup><br>
<sup>7 Computational Oncology Group, Molecular Diagnostics Program at the NCT and German Cancer Research Center (DKFZ), Heidelberg, Germany</sup><br>
<sup>8 German Cancer Consortium (DKTK), Core Center Heidelberg, Im Neuenheimer Feld 280, 69120 Heidelberg, Germany</sup><br>
<sup>9 Medical Faculty Heidelberg and Faculty of Biosciences, Heidelberg University, 69120 Heidelberg, Germany</sup><br>
<sup>10 National Center for Tumor Diseases (NCT), Im Neuenheimer Feld 410, 69120 Heidelberg, Germany</sup><br>

## Abstract
Cancer is a heterogeneous disease caused by genetic alterations. The computational analysis of cancer genomes led to the expansion of the catalog of functional mutations. While individual high-impact mutations have been discovered also in gene promoters, frequency-based approaches have only characterized a few candidates so far. To facilitate the identification of rare activating promoter mutations in cancer, we developed a filtering-based computational workflow and applied it to the Pan Cancer Analysis of Whole genomes (PCAWG) dataset. Predicted mutations were investigated using our new visualization framework, pSNV Hunter and prioritized for functional validation by luciferase assay. Here, we positively validated seven candidate pSNVs in vitro, including mutations within the promoters of ANKRD53 and MYB. Our analysis indicates that co-alterations, such as the overexpression or activation of the transcription factors, impact the effectiveness of functional pSNVs. Our analysis more than doubles the number of validated activating promoter mutations in cancer and demonstrates the effectiveness of our filtering pipeline, as well as, pSNV Hunter.

## Additional Repositories
The publication references three additional tools that can be found at the following links:
- [pSNV Hunter](https://github.com/nicholas-abad/pSNV-hunter): Comprehensive visualization tool / dashboard to investigate and select Promoter SNVs (pSNVs) for downstream validation
- [Deep Pileup](https://github.com/nicholas-abad/deep-pileup-wrapper): A quality control approach for evaluating individual genomic loci for potential signal noise
- [Genome Tornado Plots Wrapper](https://github.com/nicholas-abad/genome-tornado-plot-wrapper): Analyzing Copy Number Variation (CNV) Events within the PCAWG dataset via GenomeTornadoPlot

## Contact:
- Please contact Nicholas Abad (nicholas.a.abad@gmail.com) if you have any questions, comments or concerns.
