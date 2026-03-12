## Description
This project introduces a data-driven approach to solar radiation quality control. 
By leveraging iForest and a new climatic regime classification, we significantly 
reduce both false alarms and missed detections in BSRN data...

# RadiationClimQC
## Advanced Quality Control for Global Solar Irradiance Data

This repository provides the R implementation of a state-of-the-art framework for the quality control of baseline solar radiation measurements. By combining a novel, data-driven radiation climate classification with the Isolation Forest (iForest) anomaly detection algorithm, this project offers a more robust and precise alternative to traditional QC methods.

### Overview
Quality Control (QC) is critical for ensuring the accuracy of climate datasets and solar energy assessments. Traditional QC tests, such as the standard "Extremely Rare Limits" (ERL) test, often suffer from high rates of false negatives (missed detections) in polluted or high-aerosol environments and false positives (false alarms) in clear-sky regions.

**This project addresses these challenges through two primary innovations:**
1.  **Data-Driven Radiation Climate Classification**: We derived a new classification system for 77 global BSRN stations, allowing for localized and regime-specific QC limit configurations.
2.  **iForest Integration**: The iForest algorithm is utilized for multi-component anomaly detection, effectively identifying outliers that evade standard component-level tests.

The result is a significantly more effective QC process that minimizes both type I and type II errors across diverse climatic regimes.

### Project Structure
```text
.
├── Code/                   # R and Python scripts for processing/analysis
│   ├── 1.1.Read_BSRN.R
│   ├── 1.2.Download_McClear.py
│   ├── 1.3.Feature_extraction.R
│   ├── 1.4.PCA_hclust.R
│   ├── 2.1.Fit_new_limits.R
│   ├── 2.2.Compare_reject.R
│   ├── 2.3.Cases.R
│   ├── 3.1.Metadata.R
│   ├── 3.2.Fitted_parameters.R
│   ├── 3.3.Reject_rate.R
│   ├── 4.1.ERLexample.R
│   ├── 4.2.Harmonic.R
│   ├── 4.3.CompareFit.R
│   ├── 4.4.ClusterVis.R
│   ├── 4.5.WorldMap.R
│   ├── 4.6.Multiplot.R
│   ├── 4.7.ScatterVis.R
│   ├── 4.8.FalsePositive.R
│   └── 4.9.FalseNegative.R
├── Data/                   # Full dataset (Excluded from Git)
├── SampleData/             # Representative datasets for demonstration
│   └── Data_Recreation.md  # Guide to recreating the data environment
├── tex/                    # Manuscript figures (PDF) and tables (TeX)
├── README.md
└── .gitignore
```

> [!IMPORTANT]
> **Git Upload Note**: Not all files are uploaded with Git. Specifically, the `Data/` directory and certain large binary or temporary files are excluded to maintain repository efficiency. Please refer to [SampleData/Data_Recreation.md](SampleData/Data_Recreation.md) for instructions on how to recreate the necessary data structure.


### Key Features
- **iForest (Isolation Forest)**: Utilizes the iForest algorithm for anomaly detection to identify and isolate multi-component solar radiation outliers.
- **Improved QC Limits**: Configuration and optimization of new ERL (Extreme Rate of Change) limits tailored to different climatic regimes.
- **Climatic Regime Segmentation**: Analysis tailored to different radiation regimes (Clear sky vs. Polluted/Aerosol-heavy).
- **Visualization**: High-quality `ggplot2` plots for identifying physical inconsistencies in radiation data.

### Getting Started
1. Clone the repository:
   ```bash
   git clone https://github.com/dazhiyang/clim-limit-fit-solar.git
   ```
2. Set the working directory to the project root in R.
3. To recreate the required `Data/` folder using provided samples, please follow the instructions in [SampleData/Data_Recreation.md](SampleData/Data_Recreation.md).

### Data and Contact
**Authors**: Dazhi Yang  
**University**: Harbin Institute of Technology  
**Contact**: dazhiyang.nus@gmail.com

---
*Generated for GitHub: dazhiyang/clim-limit-fit-solar*
