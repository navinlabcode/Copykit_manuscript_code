# CopyKit Manuscript Code

This repository contains the main data processing workflows, benchmarking scripts, and supplementary code associated with the manuscript: **"Delineating the copy-number substructure of metastatic tumors with CopyKit"**.

While the [CopyKit R package](https://github.com/navinlabcode/copykit) is the primary software tool, this repository houses the specific scripts used for the sample processing, performance benchmarks, simulations and visualizations described in the paper.

## 📂 Repository Structure

### 1. `benchmark/`
This directory contains the scripts used to evaluate the performance of CopyKit and compare it against other single-cell copy number analysis tools.
* **benchmark_accuracy:** Code for generating synthetic single-cell datasets to test accuracy and resolution.
* **benchmark_time:** Scripts to run CopyKit alongside other methods to benchmark runtime, memory usage, and segmentation accuracy.
* **downsampling:** Scripts to run CopyKit alongside other methods to benchmark runtime, memory usage, and segmentation accuracy.

### 2. `sample_process/`
This directory contains the primary pipeline scripts used to process the single-cell DNA sequencing data (scDNA-seq) of the 7 tumor and 3 cell line samples. This directory also includes the scripts for statistical tests and visualizations presented in the study.

### 3. `singularity/`
Contains definition files for building the Singularity container to run the scripts in the `sample_process/` directory and reproduce the results and plots in the paper.

---

## 🔗 Related Repositories

* **Main Software Package:** [navinlabcode/copykit](https://github.com/navinlabcode/copykit)
*The user-friendly R package for single cell copy number analysis.*
* **CopyKit User Guide:** [navinlabcode/CopyKit-UserGuide](https://navinlabcode.github.io/CopyKit-UserGuide/)
*CopyKit complete documentation for general users.*
* **Step-by-step sample process:** [navinlabcode/CopyKit_paper](https://navinlabcode.github.io/CopyKit_paper/)
* [Depreciated/outdated] R markdown tutorial for the tumor samples shown in the manuscript.*

## 📖 Citation

If you use the code or workflows in this repository, please cite the CopyKit manuscript:


## 📞 Contact

For questions regarding the code in this repository, please open an issue or contact the **Navin Lab** at MD Anderson Cancer Center.

```

```
