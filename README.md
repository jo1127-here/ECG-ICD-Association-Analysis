
# ECG–ICD Association Analysis

**Python Programming Portfolio | Yuxuan Zhou**

This project investigates associations between electrocardiogram (ECG) T-wave characteristics and ICD-9 diagnoses using patient-level clinical data. It demonstrates an end-to-end analytical workflow, from data cleaning and feature engineering to statistical modeling, multiple-testing correction, and visualization.

> **Data Availability:** This project uses clinical data accessed through PhysioNet under credentialed-access requirements. Raw data and patient-level processed datasets are not publicly distributed.

## 1. Research Objective

The primary objective is to investigate whether ICD-9 diagnoses are associated with two ECG T-wave characteristics:

- Average T-wave peak amplitude (`avg_t_peak_val`)
- T-wave duration (`t_duration`)

Because individual patients may contribute multiple ECG observations, Linear Mixed-Effects Models (LMMs) are used to account for within-patient correlation.

The analysis focuses on statistical associations rather than causal inference.

## 2. Technical Skills Demonstrated

### Data Processing and Cleaning
- Loaded and validated clinical CSV datasets using pandas.
- Standardized patient identifiers and ICD-9 diagnosis codes.
- Converted invalid measurements to missing values.
- Removed duplicate patient-diagnosis combinations.
- Handled missing values and validated dataset joins.

### Feature Engineering
- Identified the 20 most prevalent ICD-9 diagnoses.
- Retained ICD-9 code 2449 as a prespecified target.
- Converted long-format diagnosis records into patient-level binary indicators using pivot tables.
- Merged diagnosis features with repeated ECG observations.

### Statistical Modeling
- Developed a reusable Python function for fitting LMMs.
- Standardized ECG outcomes before model fitting.
- Included patient-specific random intercepts.
- Implemented convergence checks and model diagnostics.
- Extracted standardized coefficients, standard errors, confidence intervals, and p-values.

### Multiple-Testing Correction
- Applied Benjamini–Hochberg FDR correction separately to each ECG outcome.
- Distinguished nominal significance from FDR-adjusted significance.

### Visualization
- Selected the ten largest absolute standardized coefficients among models without flagged diagnostic warnings.
- Created a forest plot with 95% confidence intervals.
- Distinguished positive and negative associations using color.
- Marked associations passing the FDR significance threshold.

## 3. Statistical Framework

For each ICD-9 diagnosis and ECG outcome, the following random-intercept model was fitted:

```text
Standardized ECG Outcome ~ ICD Indicator + (1 | Patient ID)
```

The ICD indicator is defined at the patient level, while the random intercept accounts for repeated ECG observations within patients.

The resulting coefficient describes a standardized association between the diagnosis indicator and the ECG measurement.

## 4. Example Finding

For ICD-9 code 2449 (Acquired Hypothyroidism), the model estimated a positive association with standardized T-wave peak amplitude.

| Statistic | Estimate |
|---|---:|
| Standardized coefficient | 0.2401 |
| 95% confidence interval | [0.0065, 0.4736] |
| Unadjusted p-value | 0.04397 |
| BH-FDR adjusted q-value | 0.08794 |

Although the association reached nominal statistical significance, it did not remain significant after BH-FDR correction at the 0.05 threshold.

This example illustrates why multiple-testing correction is important when screening numerous clinical associations.

## 5. Repository Structure

```text
ECG-ICD-Association/
├── README.md
├── ECG_ICD_Analysis.ipynb
├── requirements.txt
├── .gitignore
└── figures/
    └── top10_icd_forest.png
```

The `figures/` directory is optional and should contain only aggregate visualizations approved for public sharing.

## 6. Data Requirements

The analysis requires three input files:

| File | Description |
|---|---|
| 01_twave_features.csv | ECG-derived T-wave measurements |
| icd_ecg_filtered.csv | Patient diagnosis records |
| d_icd_diagnoses.csv | ICD-9 diagnosis descriptions |

The first two are project-specific analysis inputs. Their upstream preprocessing is not included in this repository.

The input data must be placed in a local `data/` directory and must not be committed to GitHub.

## 7. Reproducing the Analysis

1. Obtain authorized access to the relevant PhysioNet clinical datasets.
2. Prepare the three required input files.
3. Place them in the local `data/` directory.
4. Install the required Python libraries:

```bash
pip install -r requirements.txt
```

5. Open `ECG_ICD_Analysis.ipynb` in Jupyter Notebook or VS Code.
6. Run the cells sequentially.

Generated results are stored locally in `icd_screening_results/`.

## 8. Data Privacy and Responsible Sharing

The underlying clinical data are subject to PhysioNet credentialed-access requirements and applicable data use agreements.

To comply with these restrictions, this repository does not contain:

- Raw ECG or clinical records
- Patient-level diagnosis matrices
- Patient identifiers
- Restricted intermediate datasets
- Notebook outputs containing patient-level information

Only analysis code, documentation, and appropriately reviewed aggregate outputs may be publicly shared.

Researchers wishing to reproduce the analysis must obtain their own authorized access through [PhysioNet](https://physionet.org/).

## 9. Limitations

- The study is observational and cannot establish causality.
- Diagnosis indicators do not necessarily establish the timing of diagnoses relative to ECG acquisition.
- The analysis screens selected prevalent diagnoses and a prespecified target rather than all ICD codes.
- Model diagnostics cannot eliminate all potential statistical limitations.
- Upstream ECG waveform feature extraction is outside the scope of this repository.

## Author

**Yuxuan Zhou**

Python Programming Portfolio — Graduate Application
