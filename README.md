
### Folder details

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
  Archived notebooks/scripts from early project stages (e.g., proposal writing, initial concentration–discharge exploration).

- **[`deprecated/`](deprecated/)**  
  Old or replaced scripts retained for traceability.

## Workflow details
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

**Roll-up to watershed**
- After interpolation, if multiple stations are mapped for a watershed/variable, take the **daily average** across mapped stations to produce the **site-ready** series.

**Composite site (GSLOOK_FULL)**
- **Met variables:** daily average of component watersheds (`GSWS01`, `GSWS06`, `LONGER`, `COLD`).
- **Discharge:** computed from raw `GSLOOK` records and drainage area, then joined to the composite.

**Discharge units**
- Convert to depth over area:
  \[
  Q_{\text{mm d}^{-1}} = \text{MEAN\_Q} \times 0.0283168 \times 86400 \,/\, \text{DA\_M2} \times 1000
  \]
  *Assumes `MEAN_Q` is in **cfs**. Adjust the factor if in m³ s⁻¹.*

### Assumptions & thresholds

- **Date parsing:** tolerant to ISO/US/EU and **Excel serials**.
- **Minimum overlap:** pairs **≥ 5** days; triplets **≥ 10** days.
- **Duplicates:** duplicate `DATE × SITECODE` rows are **averaged** before modeling.
- **Physical constraints:** cap **RH ≤ 100%**; set **precip ≥ 0**.
- **VPD:** computed **after** interpolation from T and RH so it reflects gap-filled inputs.
- **Drainage areas:** read from `drainage_area.csv`; consistent station recodes (e.g., `GSWSMC → GSMACK`) applied prior to joins.

### Outputs (where to look)
- Master, site-level table: `.../05_Outputs/MET/data/watersheds_met_data_q.csv`
- QA plots (pair/triplet fits, time series): `.../05_Outputs/MET/plots/`


