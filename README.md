# qPCR Standard Curve & Absolute Copy Number Builder

## Overview

This Streamlit application guides users through quantitative PCR (qPCR) data analysis, including:

1. Defining the undiluted stock copy number from either nucleotide sequence or concentration and sequence length
2. Generating a standard curve using a series of 10‑fold dilutions and user‑provided Ct values
3. Performing linear regression (log₁₀(copy number) vs. Ct) and visualizing the standard curve with regression parameters
4. Calculating absolute copy numbers for unknown samples based on Ct values and dilution factors

## Features

* Interactive input forms for sequence, concentration, and Ct values
* Automatic default Ct values spaced by three cycles starting from 9
* Dynamic table display of standards and unknown sample results
* High‑quality matplotlib plots embedded in the Streamlit interface

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/Zubair2021/qPCR_copy_number_calc
   cd qPCR_copy_number_calc
   ```
2. Install dependencies:

   ```bash
   pip install streamlit pandas numpy scipy matplotlib biopython
   ```

## Usage

Run the Streamlit app:

```bash
streamlit run app.py
```

Follow the on‑screen instructions in four sections:

1. **Define Base Copy Number**

   * Choose between computing copy number from concentration & sequence length or from a nucleotide sequence.
   * Enter stock concentration and either sequence length or paste the nucleotide sequence.
2. **Build Standard Curve**

   * Select the number of dilution points.
   * Enter names and Ct values for each 10‑fold dilution (defaults provided but adjustable).
   * Review the table of dilutions and calculated copy numbers.
3. **Fit & Visualize**

   * View the standard curve plot (log₁₀(copy number) vs. Ct).
   * Regression slope, intercept, and R² are displayed on the plot.
4. **Calculate Unknowns**

   * Specify the number of unknown samples.
   * Enter Ct and dilution factor for each sample.
   * View absolute copy number estimates in the results section.

## Methods

* **Copy Number Calculation (ng & length)**
  Uses: copies = (concentration\_ng × 10⁻⁹ g/ng) ÷ (length\_bp × 650 g/mol/bp) × Avogadro’s number.

* **Copy Number Calculation (sequence)**
  Uses Biopython to compute molecular weight of RNA (T→U conversion), then:
  copies = (concentration\_ng × 10⁻⁹ g/ng) ÷ molecular\_weight\_g/mol × Avogadro’s number.

* **Standard Curve Regression**
  Linear regression of Ct (y) against log₁₀(copy number) (x) using SciPy’s `linregress`.
  Outputs slope (m), intercept (b), and coefficient of determination (R²).

* **Unknown Sample Estimation**
  Back‑calculates log₁₀(copy number) = (Ct – b) / m, then applies dilution factor to yield absolute copy number.


