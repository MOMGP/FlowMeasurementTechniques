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

The Python step writes figures to:

```text
outputs/figures/
```

`outputs/` is ignored because it is regenerated from tracked scripts and data.

## Main Outputs

- `main_hwa_calibration.*`
- `main_hwa_profiles.*`
- `main_piv_mean_fields.*`
- `main_hwa_piv_comparison.*`
- `main_piv_processing_summary.*`
- `main_piv_window_size_screenshots.*`

## HWA/PIV Comparison

The comparison uses the HWA wake traverse and the PIV column closest to
`x = -95 mm` in the DaVis coordinate system. This corresponds to about 20 mm
downstream of the trailing edge, because the exported PIV wake runs toward
negative x and the trailing edge is near `x = -75 mm`.

DaVis `Vx` is negative in the tunnel streamwise direction for this dataset, so
the exported comparison uses:

```text
u_streamwise = -Vx
```

The raw DaVis `Vx` is still exported in `piv_mean_fields.csv`.
