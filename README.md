<img src="keto.png" width="150" height="50"> 

# A modified Mediterranean ketogenic diet mitigates modifiable risk factors of Alzheimer’s disease: a serum and CSF-based metabolic analysis

This repository contains scripts for analysis and figure generation of the paper<br/> **Schweickart\*, Batra\* et. al.** A modified Mediterranean ketogenic diet mitigates modifiable risk factors of Alzheimer’s disease: a serum and CSF-based metabolic analysis (2023), *medRxiv*. [link]()

# Content
| Script name | Description |
| :--- | :--- |
| custom_functions.R  | Internal customized functions. Will be sourced by other scripts |
| 1_serum_dataprocessing.R| BEAM serum metabolomics preprocessing |
| 2_serum_association_analysis.R  | BEAM serum association analysis |
| 3_csf_dataprocessing.R | BEAM CSF metabolomics preprocessing|
| 4_csf_association_analysis.R | BEAM CSF association analysis |
| 5_metabolome_figures.R  | Figures 2a, 2b, 2d, 3a, 4 (Figures 1, 3b, and 4 (partially) were manually generated) |
| 6a_joint-tensor-nightingale-microbiome.ipynb | Microbiome analysis |
| 6b_joint_tensor_factorization_microbiome_nightingale.R | Microbiome-metabolme joint analysis |
| 6c_microbiome_figures.R  | Figure 2c |
|||
| **Folder** | **Description** |
| input | Contains the input data to be used by the scripts |
| results | Will contain output files |
|||
| **Output files in results/** | **Description** |
| tmp* | Intermediate files used in follow-up scripts |
| supplementary* | Supplementary files from the paper |
| Figure* | Figure panels from the paper |

**Input files for scripts:**  tmp_annotations_sig_serum_diet.csv is provided with content in the input folder. 

# How to

Run scripts in sequence 1 to 5 for metabolomics analysis and figures. 
Run scripts 6a, 6b, 6c for microbiome analysis and figures.