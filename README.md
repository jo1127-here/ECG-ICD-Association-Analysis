# ECG–ICD Association Analysis


This project examines exploratory associations between ECG-derived **average T-wave peak amplitude** (`avg_t_peak_val`) and patient-level ICD-9 diagnosis history. The Python workflow integrates repeated ECG observations with diagnosis indicators, fits linear mixed-effects models, applies multiple-testing correction, and visualizes effect estimates.

> **Data access:** This analysis uses credentialed-access clinical data from PhysioNet. Neither raw records nor patient-level derived datasets are distributed in this repository.

## Research question

Are common ICD-9 diagnoses associated with differences in standardized T-wave peak amplitude? Patients may contribute multiple ECG observations, so the analysis uses a patient-specific random intercept to account for within-patient correlation. The findings are **associations, not causal effects**.

## Workflow

1. **Validate and clean:** Check required columns, normalize `subject_id` and ICD codes, restrict to ICD-9, and convert invalid ECG measurements to missing values.
2. **Define the cohort:** Retain ECG patients with available ICD-9 diagnosis records.
3. **Engineer diagnosis features:** Select the 20 most prevalent ICD-9 codes in the matched cohort; deduplicate patient–diagnosis pairs; construct binary patient-level indicators and merge them onto ECG observations.
4. **Fit mixed-effects models:** For each ICD code, fit a random-intercept linear mixed model with standardized T-wave peak amplitude as the outcome. Check minimum exposed/unexposed patient counts, convergence, model warnings, and coefficient validity.
5. **Correct for multiple testing:** Apply Benjamini–Hochberg false discovery rate (BH-FDR) correction to the valid, non-warning-flagged models.
6. **Visualize:** Display the 10 largest absolute standardized coefficients in a forest plot with 95% confidence intervals. Color indicates association direction; an asterisk indicates BH-FDR q < 0.05.

### Statistical model

```text
standardized(avg_t_peak_val) ~ ICD_indicator + (1 | subject_id)
```

The ICD indicator describes whether a diagnosis appears in a patient's recorded ICD-9 history. It does not establish whether the diagnosis was active at the time of an ECG.

## Illustrative result

For ICD-9 **2449 (acquired hypothyroidism)**, an earlier run of this analysis reported:

| Measure | Estimate |
|---|---:|
| Standardized coefficient | +0.2401 |
| 95% confidence interval | [0.0065, 0.4736] |
| Unadjusted p-value | 0.04397 |
| BH-FDR q-value | 0.08794 |

This association was nominally significant but **did not meet the BH-FDR q < 0.05 threshold**. These values are illustrative results from the prior analysis, not a fresh rerun of the public notebook.

## Repository contents

```text
ECG-ICD-Association-Analysis/
├── README.md
├── ECG_ICD_Analysis_Public.ipynb
├── requirements.txt
├── .gitignore
└── figures/                       # Optional, reviewed aggregate figures only
    └── top10_icd_forest.png
```

The repository contains analysis code and documentation, **not** the clinical input files or patient-level outputs. An optional aggregate forest plot should be reviewed before publication.

## Running locally

Obtain the required clinical data through authorized access and prepare the following project-specific CSV inputs (upstream preparation is not included):

| Local file in `data/` | Purpose |
|---|---|
| `01_twave_features.csv` | ECG-derived T-wave features |
| `icd_ecg_filtered.csv` | Patient-level diagnosis records linked to ECG cohort |
| `d_icd_diagnoses.csv` | ICD-9 diagnosis labels |

Install dependencies with `pip install -r requirements.txt`, then open `ECG_ICD_Analysis_Public.ipynb` and run its cells. The notebook reads from the local `data/` directory and writes model-level summaries to `icd_screening_results/`.

**Do not commit** `data/`, `icd_screening_results/`, raw records, patient identifiers, or executed notebooks containing sensitive outputs. The public notebook has been cleared of execution results. Because the underlying data are restricted, an outside reader cannot run the notebook without their own authorized access and preparation of the required inputs.

## Limitations

- **Temporal alignment:** Diagnosis history is not reliably aligned with each ECG acquisition; the analysis cannot establish whether a diagnosis was active at recording time.
- **Confounding:** The models do not adjust for age, sex, medications, or comorbidities. A patient random intercept accounts for repeated observations but does not remove these confounders.
- **Selection:** The cohort includes only ECG patients with ICD-9 records, potentially limiting generalizability.
- **Exploratory screening:** The 20 most prevalent diagnoses are screened, with BH-FDR applied to valid models. Associations should not be interpreted as causal or clinically validated biomarkers.

## My contributions

I developed the Python workflow for ECG–ICD integration, patient-level feature engineering, mixed-effects modeling, model diagnostics, BH-FDR correction, and forest plot visualization.

## Data availability and responsible sharing

Clinical source data are subject to PhysioNet credentialed-access requirements and applicable data use agreements. Patient-level data and derived matrices are intentionally excluded from this public repository. Access to underlying data must be obtained independently through [PhysioNet](https://physionet.org/).

## Author

**Yuxuan Zhou** — Graduate Application Programming Portfolio
