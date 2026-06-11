# Flow Measurement Techniques Plot Export

This branch starts from the calibrated `main` MATLAB workflow and adds a
lightweight export/plotting layer.

The intention is:

- Keep the HWA calibration and processing logic from `HW_Post.m`.
- Keep the PIV DaVis data convention used by `PIV.m`.
- Export clean CSV tables from MATLAB.
- Use Python only for report-style figures.
- Add an HWA-vs-PIV comparison at the wake profile station.

## Run

From the repository root:

```powershell
matlab -batch "run(fullfile('analysis','export_main_outputs.m'))"
python analysis\plot_main_outputs.py
```

The MATLAB step writes CSV files to:

```text
outputs/tables/
```

The Python step clears any old PNG/PDF files and writes the requested figures to:

```text
outputs/figures/
```

`outputs/` is ignored because it is regenerated from tracked scripts and data.

## Requested Figures

- `01_hwa_calibration_curve.*`
- `02_hwa_autocorrelation_convergence.*`
- `03_piv_mean_fields.*`
- `04_instantaneous_piv_aoa15.*`
- `05_piv_window_size_sensitivity.*`
- `06_piv_overlap_multipass_sensitivity.*`
- `07_piv_ensemble_size_effect.*`
- `08_piv_delta_t_effect.*`
- `09_self_made_piv_vs_davis.*`
- `10_hwa_mean_rms_wake_profiles.*`
- `11_hwa_energy_spectra.*`
- `12_piv_hwa_profile_comparison_xc12.*`

## HWA/PIV Comparison

The comparison uses the HWA wake traverse and the PIV column closest to
`x = -95 mm` in the DaVis coordinate system. This corresponds to about 20 mm
downstream of the trailing edge, because the exported PIV wake runs toward
negative x and the trailing edge is near `x = -75 mm`.

The HWA and PIV vertical coordinates are not the same raw coordinate. HWA uses
the traverse coordinate relative to the zero-AoA trailing-edge reference. DaVis
uses its calibrated image/world coordinate. For profile comparison the PIV
coordinate is therefore mapped as:

```text
y_relative_for_HWA_comparison = -(Y_DaVis - 62 mm)
```

Without this mapping the HWA/PIV plot has an artificial vertical offset. After
the mapping, AoA 0 and AoA 5 align closely; AoA 15 still differs because the
flow is separated and the two methods sample/average different regions.

DaVis `Vx` is negative in the tunnel streamwise direction for this dataset, so
the exported comparison uses:

```text
u_streamwise = -Vx
```

The raw DaVis `Vx` is still exported in `piv_mean_fields.csv`.
