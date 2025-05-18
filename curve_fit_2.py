# Streamlit qPCR Standard Curve & Absolute Copy Number Builder
# Requirements:
#   pip install streamlit pandas numpy scipy matplotlib biopython plotly openpyxl

import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from Bio.SeqUtils import molecular_weight
import plotly.express as px
import json
from io import BytesIO
from datetime import datetime

# Constants
AVOGADRO = 6.022e23  # molecules per mole
AVERAGE_BP_WEIGHT = 650  # g/mol per base pair

# Initialize audit log
if 'audit' not in st.session_state:
    st.session_state.audit = []

def log_event(msg):
    st.session_state.audit.append(f"{datetime.now().isoformat()} - {msg}")

# Calculation functions
def calc_ng_length(conc_ng, length_bp):
    grams = conc_ng * 1e-9
    moles = grams / (length_bp * AVERAGE_BP_WEIGHT)
    return moles * AVOGADRO

def calc_sequence(seq, conc_ng):
    seq_clean = seq.replace('T','U').replace('\n','').strip().upper()
    mw = molecular_weight(seq_clean, seq_type='RNA')
    grams = conc_ng * 1e-9
    moles = grams / mw
    return moles * AVOGADRO, mw

# App config
st.set_page_config(page_title='qPCR Copy Number Builder', layout='wide')
st.title('qPCR Standard Curve & Absolute Copy Number Builder')

# Section 0: Session management
with st.expander('Session Load/Save'):
    uploaded = st.file_uploader('Load session (JSON)', type=['json'])
    if uploaded:
        data = json.load(uploaded)
        for k,v in data.items():
            st.session_state[k] = v
        st.success('Session loaded')
        log_event('Session loaded')
    if st.button('Save session'):
        tosave = {k:v for k,v in st.session_state.items() if k!='audit'}
        buf = BytesIO()
        buf.write(json.dumps(tosave).encode())
        buf.seek(0)
        st.download_button('Download JSON', buf, 'session.json', 'application/json')
        log_event('Session saved')

# Section 1: Define base copy
st.markdown('---')
st.header('1. Define Base Copy Number')
method = st.selectbox('Method', ['Sequence', 'Concentration & Length'])
if method=='Concentration & Length':
    conc = st.number_input('Stock conc (ng)', min_value=0.0, value=100.0)
    length = st.number_input('Seq length (bp)', min_value=1, value=1000)
    base_copies = calc_ng_length(conc, length)
    st.metric('Undiluted copies', f"{base_copies:,.0f}")
    log_event('Base copies via ng-length')
else:
    seq = st.text_area('Sequence (DNA/RNA)')
    conc2 = st.number_input('Stock conc (ng)', min_value=0.0, value=100.0, key='conc_seq')
    if seq:
        base_copies, mw = calc_sequence(seq, conc2)
        st.write('Molecular weight:', f"{mw:.1f} g/mol")
        st.metric('Undiluted copies', f"{base_copies:,.0f}")
        log_event('Base copies via sequence')
    else:
        st.warning('Enter sequence')
        st.stop()

# Section 2: Build standards
st.markdown('---')
st.header('2. Build Standard Curve')
factor = st.number_input('Dilution factor', min_value=2, value=10)
points = st.slider('Points', 2, 12, 5)
start_ct = st.number_input('Start Ct', min_value=0.0, value=9.0)
interval = st.number_input('Ct interval', min_value=0.1, value=3.0)
records = []
for i in range(points):
    cols = st.columns([2, 1, 1])
    name = cols[0].text_input(f'Name {i+1}', value=f'Std{i+1}', key=f'n{i}')
    default_ct = start_ct + i * interval
    ct = cols[1].number_input('Ct', min_value=0.0, value=default_ct, step=0.1, key=f'ct{i}')
    dil = factor ** i
    copies = base_copies / dil
    cols[2].write(f'1:{dil}')
    cols[2].write(f'**{int(copies):,}**')
    records.append({'Name': name, 'Dilution': dil, 'Ct': ct, 'Copies': copies})
df_std = pd.DataFrame(records)
st.dataframe(df_std)
log_event('Standards defined')

# Section 3: Analytics & plots
st.markdown('---')
st.header('3. Analytics & Visualization')
x = np.log10(df_std['Copies'])
y = df_std['Ct']
slope, intercept, rval, _, _ = stats.linregress(x, y)
eff = (10 ** (-1 / slope) - 1) * 100
lod = df_std['Copies'].min()
col1, col2, col3, col4, col5 = st.columns(5)
col1.metric('Slope', f"{slope:.3f}")
col2.metric('Intercept', f"{intercept:.3f}")
col3.metric('R²', f"{rval**2:.3f}")
col4.metric('Efficiency (%)', f"{eff:.1f}")
col5.metric('LOD', f"{lod:,.0f}")
# Plotly curve
fig = px.scatter(df_std, x=np.log10(df_std['Copies']), y='Ct', text='Name', title='Standard Curve')
fig.add_traces(px.line(x=x, y=slope * x + intercept).data)
st.plotly_chart(fig)
# Residuals
resid = y - (slope * x + intercept)
fig2 = px.scatter(x=x, y=resid, labels={'x': 'log10(C)', 'y': 'Residual Ct'}, title='Residuals')
st.plotly_chart(fig2)
# Bootstrap CI
boot = []
for _ in range(200):
    samp = df_std.copy()
    # sample with replacement
    samp = samp.sample(len(df_std), replace=True)
    # filter valid copies
    valid = samp['Copies'] > 0
    if valid.sum() < 2:
        continue
    x_s = np.log10(samp.loc[valid, 'Copies'])
    y_s = samp.loc[valid, 'Ct']
    try:
        m, b, _, _, _ = stats.linregress(x_s, y_s)
        boot.append((m, b))
    except ValueError:
        # skip this bootstrap sample
        continue
# compute confidence intervals if we have enough samples
if boot:
    barr = np.array(boot)
    ci_m = np.percentile(barr[:, 0], [2.5, 97.5])
    ci_b = np.percentile(barr[:, 1], [2.5, 97.5])
    st.write(f"Slope CI: {ci_m[0]:.3f}-{ci_m[1]:.3f}")
    st.write(f"Intercept CI: {ci_b[0]:.3f}-{ci_b[1]:.3f}")
else:
    st.write("Insufficient data for bootstrap confidence intervals.")

# Section 4: Quality Control
st.markdown('---')
st.header('4. Quality Control')
# User-configurable thresholds
min_slope = st.number_input('Min acceptable slope', value=-3.6, step=0.01)
max_slope = st.number_input('Max acceptable slope', value=-3.1, step=0.01)
r2_thresh = st.number_input('Min acceptable R²', value=0.98, step=0.01)
# Check metrics against thresholds
if (rval**2) < r2_thresh:
    st.warning(f'Low R²: {rval**2:.3f} below threshold {r2_thresh}')
if not (min_slope <= slope <= max_slope):
    st.warning(f'Slope {slope:.3f} outside expected range ({min_slope} to {max_slope})')

# Section 5: Calculate Unknowns
st.markdown('---')
st.header('5. Calculate Unknown Sample Copies')
st.markdown('Enter each unknown sample’s Ct and dilution to get absolute copy numbers.')
num_unk = st.number_input('Unknown count', min_value=1, max_value=12, value=3)
unk = []
for j in range(int(num_unk)):
    a, b, c = st.columns(3)
    nm = a.text_input('Name', value=f'U{j+1}', key=f'unm{j}')
    ctu = b.number_input('Ct', min_value=0.0, value=25.0, step=0.1, key=f'ctu{j}')
    du = c.number_input('Dilution', min_value=0.0, value=1.0, step=0.1, key=f'diu{j}')
    logn = (ctu - intercept) / slope
    cu = 10 ** logn * du
    if ctu < df_std['Ct'].min() or ctu > df_std['Ct'].max():
        b.warning('Ct out of standards range')
    unk.append({'Name': nm, 'Ct': ctu, 'Dilution': du, 'Copies': cu})
df_unk = pd.DataFrame(unk)
st.dataframe(df_unk)


