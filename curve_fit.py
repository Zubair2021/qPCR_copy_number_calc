# Streamlit qPCR Standard Curve & Absolute Copy Number App
# Requirements:
#   pip install streamlit pandas numpy scipy matplotlib biopython

import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from Bio.SeqUtils import molecular_weight

# Constants
AVOGADRO = 6.022e23  # molecules per mole
AVERAGE_BP_WEIGHT = 650  # g/mol per base pair (fallback)

# Function Definitions
# --------------------
def calc_base_copies_ng_length(conc_ng: float, length_bp: int) -> float:
    """Calculate copies from concentration (ng) and sequence length (bp)."""
    grams = conc_ng * 1e-9
    moles = grams / (length_bp * AVERAGE_BP_WEIGHT)
    return moles * AVOGADRO


def calc_base_copies_sequence(seq: str, conc_ng: float) -> (float, float):
    """Calculate molecular weight and copies from a sequence and concentration."""
    seq_rna = seq.replace('T', 'U').replace('\n', '').strip().upper()
    mw = molecular_weight(seq_rna, seq_type='RNA')
    grams = conc_ng * 1e-9
    moles = grams / mw
    copies = moles * AVOGADRO
    return copies, mw

# App Configuration
# -----------------
st.set_page_config(
    page_title='qPCR Copy Number Builder',
    layout='wide',
    page_icon='🧬'
)
st.title('🧬 qPCR Standard Curve & Copy Number Builder')

st.markdown("""
This interactive app walks you through:

1. **Defining a base copy number** (from stock concentration & length or sequence).
2. **Building a standard curve** using a series of 10-fold dilutions and Ct values.
3. **Fitting and visualizing** the standard curve (log₁₀ copies vs. Ct).
4. **Calculating unknown sample copy numbers** from Ct and dilution.

Follow each section below and enter the required inputs.
"""
)

# Section 1: Base Copy Number Setup
# -------------------------------
st.markdown('---')
st.header('1️⃣ Define Base Copy Number')
st.markdown(
    'Choose how to calculate the **undiluted** stock copy number. '
    'You only need to enter these details once for your standards.'
)
base_copies = None
method_base = st.radio(
    'Compute base copies by:',
    ['Concentration & Length', 'Sequence & Concentration']
)
if method_base == 'Concentration & Length':
    conc_base = st.number_input(
        'Enter stock concentration (ng)',
        min_value=0.0, value=100.0, step=0.1
    )
    length_base = st.number_input(
        'Enter sequence length (bp)',
        min_value=1, value=1000, step=1
    )
    base_copies = calc_base_copies_ng_length(conc_base, length_base)
    st.markdown(f"### Undiluted stock contains **{base_copies:,.0f} copies**")
else:
    seq_base = st.text_area(
        'Paste nucleotide sequence (DNA/RNA)', height=150
    )
    conc_base = st.number_input(
        'Enter stock concentration (ng)',
        min_value=0.0, value=100.0, step=0.1
    )
    if seq_base:
        base_copies, mw = calc_base_copies_sequence(seq_base, conc_base)
        st.write(f'- Molecular weight: **{mw:,.1f} g/mol**')
        st.markdown(f"### Undiluted stock contains **{base_copies:,.0f} copies**")
    else:
        st.warning('Please paste a valid sequence above.')
        st.stop()

# Section 2: Build Standard Curve
# -------------------------------
st.markdown('---')
st.header('2️⃣ Build Your Standard Curve')
st.markdown(
    'This step auto-generates default Ct values 3 cycles apart starting from 9, '
    'but you can adjust them if needed. Enter names for each standard as well.'
)
num_std = st.slider(
    'Number of dilution points', min_value=2, max_value=12, value=5
)
# Default Ct parameters
start_ct = 9.0
ct_interval = 3.0

def get_default_name(i):
    return f'Standard {i+1}'

records = []
for i in range(num_std):
    cols = st.columns([2, 1, 1])
    name = cols[0].text_input(
        f'Name for point {i+1}', value=get_default_name(i), key=f'name_{i}'
    )
    default_ct = start_ct + ct_interval * i
    ct = cols[1].number_input(
        f'Ct value for {name}', min_value=0.0,
        value=default_ct, step=0.1, key=f'ct_{i}'
    )
    dilution = 10 ** i
    copies = base_copies / dilution
    cols[2].markdown(f'1:{dilution}')
    cols[2].markdown(f'**{int(copies):,} copies**')
    records.append({'Standard': name, 'Dilution': f'1:{dilution}', 'Ct': ct, 'Copies': copies})
df_std = pd.DataFrame(records)
st.table(df_std.style.format({'Copies': '{:,.0f}'}))

# Section 3: Fit & Plot Curve
# ---------------------------
st.markdown('---')
st.header('3️⃣ Fit & Visualize Standard Curve')
st.markdown(
    'This plot shows log₁₀(copy number) vs. Ct. Regression line & R² appear at bottom-right.'
)
# Regression prep
log_copies = np.log10(df_std['Copies'])
ct_vals = df_std['Ct']
slope, intercept, r_val, _, _ = stats.linregress(log_copies, ct_vals)

def plot_curve(x, y, slope, intercept, r_sq):
    fig, ax = plt.subplots(figsize=(6,4))
    ax.scatter(x, y, s=80, edgecolor='k', alpha=0.8)
    x_fit = np.linspace(x.min(), x.max(), 100)
    y_fit = slope*x_fit + intercept
    ax.plot(x_fit, y_fit, linestyle='--', linewidth=2)
    ax.set(xlabel='log₁₀(Copy Number)', ylabel='Ct Value', title='qPCR Standard Curve')
    ax.grid(True, linestyle=':')
    ax.tick_params(labelsize=10)
    ax.title.set_fontsize(14)
    stats_text = f'm={slope:.2f}\nb={intercept:.2f}\nR²={r_sq:.3f}'
    ax.text(0.98, 0.95, stats_text, transform=ax.transAxes,
            ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    st.pyplot(fig)

plot_curve(log_copies, ct_vals, slope, intercept, r_val**2)

# Section 4: Compute Unknowns
# --------------------------
st.markdown('---')
st.header('4️⃣ Calculate Unknown Sample Copies')
st.markdown("Enter each unknown sample's Ct and dilution to get absolute copy numbers.")
num_unk = st.number_input('How many unknown samples?', min_value=1, max_value=12, value=3)
unk_records = []
for j in range(num_unk):
    cols = st.columns(3)
    sample = cols[0].text_input(f'Sample {j+1} name', value=f'Unknown {j+1}', key=f'uname_{j}')
    ct_u = cols[1].number_input(f'Ct for {sample}', min_value=0.0, value=25.0, step=0.1, key=f'ctu_{j}')
    dil_u = cols[2].number_input(f'Dilution factor', min_value=0.0, value=1.0, step=0.1, key=f'dilu_{j}')
    log_n = (ct_u - intercept) / slope
    copies_u = (10**log_n) * dil_u
    unk_records.append({'Sample': sample, 'Ct': ct_u, 'Dilution': dil_u, 'Copies': copies_u})
if unk_records:
    st.markdown('**Results for Unknown Samples:**')
    for rec in unk_records:
        st.markdown(f"- **{rec['Sample']}**: **{rec['Copies']:,.0f} copies** (Ct={rec['Ct']}, Dilution={rec['Dilution']})")

st.markdown('---')
st.info('You have completed all steps. Review the tables and plots above for your results.')
