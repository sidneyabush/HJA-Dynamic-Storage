
# Folder details

- **[`Create_Master_Hydrometric_Dataset/`](Create_Master_Hydrometric_Dataset/)**  
  A *comprehensive* workflow that harmonizes meteorology and streamflow for ET calculations and downstream analyses.

- **[`Discharge_metrics/`](Discharge_metrics/)**  
  Q-based indicators including:
  - **RBI (flashiness)**
  - **Recession curve slope** (log–log slope of −dQ/dt vs. Q)
  - **Dynamic (active) storage** using the storage–discharge (Kirchner–Staudinger) approach

- **[`ET_calculations/`](ET_calculations/)**  
  Evapotranspiration estimation using multiple methods and combined summaries
  
- **[`Plots/`](Plots/)**  
  Scripts to generate figures for exploration and publication.

- **[`Prior_analyses/`](Prior_analyses/)**  
  Archived scripts from early project stages (e.g., proposal writing, initial concentration–discharge exploration).

- **[`deprecated/`](deprecated/)**  
  Old or replaced scripts retained for reference.

# Workflow details
## Hydromet harmonization & interpolation (highlights)

Builds **daily, site-level meteorology** for each watershed and attaches **discharge**. Focus here is on how multi-station data are combined.

### How data are combined

**Station → Watershed mapping**
- Each watershed lists one or more stations per variable: **temp (T_C)**, **precip (P_mm_d)**, **RH (RH_d_pct)**, **net radiation (NR_Wm2_d)**.

**Single-station case**
- Use the mapped station directly.

**Pairs (two stations for a variable)**
- Fit a simple **OLS** relationship on overlapping days (require **≥ 5 complete cases**).
- Use the fitted line to **predict missing days** at one station from the other.
- If predicting “back” (site1 from site2), use the **inverse** of the fitted line.
- Keep original observations; only fill gaps.

**Triplets (three stations for a variable)**
- Fit **three multiple regressions** (each station as target, other two as predictors) on overlapping days (require **≥ 10 complete cases**).
- Predict missing days **only** for the single missing station in each case.
- If triplet overlap is insufficient, **fall back** to pairwise interpolation.

**Scale-up to watershed**
- After interpolation, if multiple stations are mapped for a watershed/variable, take the **daily average** across mapped stations to produce the **site-ready** series.

**Composite site (GSLOOK_FULL)**
- **Met variables:** daily average of component watersheds (`GSWS01`, `GSWS06`, `LONGER`, `COLD`).
- **Discharge:** computed from raw `GSLOOK` records and drainage area, then joined to the composite.

### Assumptions & thresholds

- **Minimum overlap:** pairs **≥ 5** days; triplets **≥ 10** days.
- **Duplicates:** duplicate `DATE × SITECODE` rows are **averaged** before modeling.
- **Physical constraints:** cap **RH ≤ 100%**; set **precip ≥ 0**.
- **VPD:** computed **after** interpolation from T and RH so it reflects gap-filled inputs.
- **Drainage areas:** read from `drainage_area.csv`; consistent station recodes (e.g., `GSWSMC → GSMACK`) applied prior to joins.

### Outputs 
- Master, site-level table: `.../05_Outputs/MET/data/watersheds_met_data_q.csv`
- QA plots (pair/triplet fits, time series): `.../05_Outputs/MET/plots/`

### Outputs
- **Master table:** watershed-level daily `P`, `Q`, `ET` (mm/day).  
- **QA figures:** pair/triplet fit plots (with R² and 1:1 line) and per-site time series.

## ET_calculations

Produce **daily, watershed-level ET** and a **joined water-balance table** for downstream storage analysis.

### Purpose
- Compute ET from gap-filled meteorology (station → watershed).
- Use a **Hamon** formulation with a **Zhang-style calibrated coefficient** (smoothed/interpolated).
- Assemble `P (mm d⁻¹)`, `Q (mm d⁻¹)`, and `ET (mm d⁻¹)` into a single daily table.

### Inputs
- Harmonized, gap-filled met at watershed level: `T_C`, `P_mm_d`, `RH_d_pct`, `NR_Wm2_d`.
- Discharge converted to depth over drainage area: `Q_mm_d`.
- Site-station mapping from the harmonization step (e.g., single / pair / triplet).

### Workflow 
1. **Load gap-filled met** (already mapped + averaged from stations to watershed).
2. **Compute VPD (kPa)** from filled `T_C` and `RH_d_pct`.
3. **PET (Hamon + Zhang calibration)**  
   - Calculate daily Hamon PET.  
   - Apply a **calibration coefficient** (Zhang-style) and **interpolate** that coefficient across time to avoid step changes.
4. **Select ET series**  
   - Use the calibrated Hamon PET as **ET (mm d⁻¹)** for the master table.
5. **Assemble daily water balance**  
   - Join `P_mm_d`, `Q_mm_d`, and `ET_mm_d` by `DATE × SITECODE`.
6. *(Optional diagnostic)* Export QA plots comparing ET and showing per-site series.

### Key ET calculation:
- Hamon PET scaled by a **Zhang-calibrated coefficient** (time-interpolated), producing **ET (mm d⁻¹)**.

### Assumptions
- Station-watershed combining happens **before** ET:  
  **pairs ≥ 5 overlap (OLS with inverse back-prediction), triplets ≥ 10 overlap (multiple regression);** then daily **averaging** across mapped stations.  
- **Physical caps:** `RH ≤ 100%`, `P ≥ 0`.  
- **VPD** computed **after** interpolation so ET uses the same filled T & RH.  

### Outputs
- **Daily water-balance with ET:** `.../DynamicStorage/daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv`  
- **QA plots (optional):** `.../05_Outputs/MET/plots/` (pair/triplet fits, per-site time series)




