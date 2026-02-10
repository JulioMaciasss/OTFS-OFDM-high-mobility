# 📡 OTFS vs OFDM in High-Mobility Channels (BSc Thesis, UMA 2025)

This repository contains the MATLAB simulator and results for my BSc thesis at the **University of Málaga (UMA)**:  
**“Simulation and Assessment of Orthogonal Time Frequency Space (OTFS) Modulation in High-Mobility Scenarios”**.  
The work compares **OTFS** and **OFDM** under multiple channel models, with a focus on high mobility and delay–Doppler behavior.

## 🧠 Motivation
OFDM suffers from inter-carrier interference (ICI) when Doppler is high, especially in time-varying multipath channels.  
OTFS operates in the **delay–Doppler domain**, and is often proposed as a more robust waveform in high-mobility scenarios.

## ✅ What this project does
- Implements Monte Carlo simulations to estimate **BER/BLER/SINR vs SNR**.
- Evaluates OTFS vs OFDM under:
  - **AWGN** (sanity check vs theoretical QPSK curve)
  - **Sparse delay–Doppler** multipath (Doppler channel model)
  - **3GPP TDL-A** channel
  - **3GPP CDL-A** channel
- Studies OTFS detector sensitivity to the channel-gain threshold **β**.

Key conclusion (from the thesis): OTFS tends to outperform OFDM in high-mobility environments when the channel is sparse in the delay–Doppler domain, and performance is sensitive to **β**.  
More detailed results depend on the specific Doppler spectrum/channel model (e.g., TDL/CDL).  
(See thesis PDF / presentation for the full discussion.)

## 📂 Repository structure
- `matlab/` — MATLAB scripts and helper functions
- `docs/` — thesis PDF, slides, and exported figures
- `docs/figures/` — key plots (PNG)
- `results/` — optional output files (keep small; avoid large binaries)

## ⚙️ Requirements
-MATLAB R2024b
- Toolboxes:
  - For 3GPP TDL/CDL (`nrTDLChannel`, `nrCDLChannel`): **5G Toolbox** (or equivalent) is typically required (please confirm your setup).

## 🚀 Quickstart
1) Clone the repo and open MATLAB in the repository root.
2) Run the AWGN validation:
   - `matlab/OFDM_validacion.m`
   - `matlab/OTFS_validacion.m`
3) Run the main studies:
   - Sparse delay–Doppler model:
     - `matlab/OFDM_SNRvar.m`
     - `matlab/OTFS_SNRvariable.m`
   - 3GPP TDL-A:
     - `matlab/OFDM_tdl.m`
     - `matlab/OTFS_TDL.m`
   - 3GPP CDL-A:
     - `matlab/OTFS_CDL.m`

> Tip: If you want a single entry point, run `matlab/run_all.m` (provided in this repo) and edit `matlab/config.m`.

## 🧪 Reproducibility notes
Some scripts require changing parameters (SNR range, number of slots, delay spread DS, Doppler DF, guard length G/padLen, β threshold, etc.).  
To avoid manual edits across files, this repo includes:
- `matlab/config.m` — central configuration
- `matlab/run_all.m` — launches the main experiments in a consistent way

## 📈 Results (figures)
Place exported plots here and link them:
- `docs/figures/awgn_validation.png`
- `docs/figures/doppler_high_mobility.png`
- `docs/figures/tdl_results.png`
- `docs/figures/cdl_results.png`
- `docs/figures/beta_sensitivity.png`

## 📚 Citation
If you use this repository, please use the GitHub citation panel (CITATION.cff), or cite the thesis:
Julio Macías Macías, *Simulation and Assessment of OTFS Modulation in High-Mobility Scenarios*, University of Málaga, 2025.

## 🙌 Acknowledgements
Supervisor: Francisco Javier Martín Vega (UMA).

## 📄 License
- Code: TODO (MIT/BSD-3 recommended)
- Thesis document/slides: TODO (see `docs/`)
