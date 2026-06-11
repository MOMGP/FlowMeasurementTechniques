from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "outputs" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "figure.dpi": 130,
        "savefig.dpi": 300,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "lines.linewidth": 1.6,
    }
)

COLORS = {0: "#1f77b4", 5: "#2ca02c", 15: "#d62728"}


def read_table(name: str) -> pd.DataFrame:
    path = TABLES / name
    if not path.exists():
        raise FileNotFoundError(f"Run analysis/export_main_outputs.m first; missing {path}")
    return pd.read_csv(path)


def clear_figures() -> None:
    for pattern in ("*.png", "*.pdf"):
        for path in FIGURES.glob(pattern):
            path.unlink()


def valid_mask(df: pd.DataFrame) -> pd.Series:
    if "valid" not in df.columns:
        return pd.Series(True, index=df.index)
    if df["valid"].dtype == bool:
        return df["valid"]
    return df["valid"].astype(str).str.lower().isin({"true", "1"})


def save(fig: plt.Figure, name: str) -> None:
    if not fig.get_constrained_layout():
        fig.tight_layout()
    fig.savefig(FIGURES / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def value_grid(df: pd.DataFrame, value: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    work = df.copy()
    work.loc[~valid_mask(work), value] = np.nan
    grid = work.pivot(index="y_mm", columns="x_mm", values=value).sort_index().sort_index(axis=1)
    x = grid.columns.to_numpy(float)
    y = grid.index.to_numpy(float)
    X, Y = np.meshgrid(x, y)
    Z = grid.to_numpy(float)
    return X, Y, Z


def vector_grid(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    X, Y, U = value_grid(df, "u_streamwise_ms")
    work = df.copy()
    work.loc[~valid_mask(work), "v_normal_ms"] = np.nan
    v_grid = work.pivot(index="y_mm", columns="x_mm", values="v_normal_ms").sort_index().sort_index(axis=1)
    V = v_grid.to_numpy(float)
    return X, Y, U, V


def field_scale(df: pd.DataFrame, value: str = "u_streamwise_ms") -> tuple[float, float]:
    values = df.loc[valid_mask(df), value].to_numpy(float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return -1.0, 1.0
    vmin, vmax = np.nanpercentile(values, [2, 98])
    if np.isclose(vmin, vmax):
        pad = max(1.0, 0.1 * abs(vmin))
        return float(vmin - pad), float(vmax + pad)
    return float(vmin), float(vmax)


def draw_vector_field(ax: plt.Axes, df: pd.DataFrame, title: str, vmin: float, vmax: float):
    X, Y, U, V = vector_grid(df)
    contour = ax.contourf(X, Y, U, levels=np.linspace(vmin, vmax, 24), cmap="viridis", extend="both")
    skip_y = max(1, U.shape[0] // 16)
    skip_x = max(1, U.shape[1] // 18)
    ax.quiver(
        X[::skip_y, ::skip_x],
        Y[::skip_y, ::skip_x],
        U[::skip_y, ::skip_x],
        V[::skip_y, ::skip_x],
        color="white",
        scale=260,
        width=0.0024,
        alpha=0.85,
    )
    ax.set_title(title)
    ax.set_xlabel("x [mm]")
    ax.set_xticks([-100, -50, 0, 50])
    ax.set_aspect("equal", adjustable="box")
    return contour


def plot_01_hwa_calibration_curve() -> None:
    points = read_table("hwa_calibration_points.csv")
    curve = read_table("hwa_calibration_curve.csv")
    fig, ax = plt.subplots(figsize=(5.6, 3.7))
    ax.scatter(points["signal_v"], points["velocity_ms"], s=22, color="black", label="Calibration points")
    ax.plot(curve["signal_v"], curve["velocity_ms"], color="#1f77b4", label="PCHIP interpolation")
    ax.set_xlabel("CTA output signal [V]")
    ax.set_ylabel("Velocity [m/s]")
    ax.set_title("HWA calibration curve")
    ax.legend()
    save(fig, "01_hwa_calibration_curve")


def plot_02_hwa_autocorrelation_convergence() -> None:
    autocorr = read_table("hwa_autocorrelation.csv")
    convergence = read_table("hwa_convergence.csv").iloc[0]
    first_zero = float(convergence["first_zero_crossing_s"])
    integral_time = float(convergence["integral_time_s"])
    lag_limit = min(float(autocorr["lag_s"].max()), max(0.01, 30.0 * first_zero))
    visible = autocorr[autocorr["lag_s"] <= lag_limit]

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4), constrained_layout=True)

    axes[0].plot(visible["lag_s"] * 1000.0, visible["rho"], color="#1f77b4", lw=1.4)
    axes[0].axhline(0, color="black", lw=0.8)
    axes[0].axvline(1000.0 * first_zero, color="#d62728", ls="--", label="First zero")
    axes[0].axvline(1000.0 * integral_time, color="#2ca02c", ls=":", label="Integral time")
    axes[0].set_xlabel("Lag [ms]")
    axes[0].set_ylabel("Autocorrelation coefficient")
    axes[0].set_title("Autocorrelation")
    axes[0].legend(loc="upper right")

    mean_velocity = abs(float(convergence["mean_velocity_ms"]))
    rms_velocity = float(convergence["rms_velocity_ms"])
    coverage_factor = float(convergence["coverage_factor"])
    target = 100.0 * float(convergence["target_uncertainty_fraction"])
    required_time = float(convergence["sampling_time_required_s"])
    minimum_time = max(2.0 * integral_time, 1.0 / float(convergence["sample_rate_hz"]))
    maximum_time = max(0.1, 30.0 * required_time)
    averaging_time = np.geomspace(minimum_time, maximum_time, 220)
    uncertainty = (
        100.0
        * coverage_factor
        * rms_velocity
        / mean_velocity
        / np.sqrt(averaging_time / (2.0 * integral_time))
    )

    axes[1].plot(averaging_time, uncertainty, color="#1f77b4", lw=1.4)
    axes[1].axhline(target, color="black", ls="--", lw=0.9, label=f"{target:.1f}% target")
    axes[1].axvline(required_time, color="#d62728", ls="--", lw=1.0, label=f"T = {required_time:.4f} s")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Averaging time [s]")
    axes[1].set_ylabel("Estimated uncertainty [%]")
    axes[1].set_title("Mean convergence")
    axes[1].legend(loc="upper right")
    save(fig, "02_hwa_autocorrelation_convergence")


def plot_03_piv_mean_fields() -> None:
    fields = read_table("piv_mean_fields.csv")
    profiles = read_table("piv_profile_xc12.csv")
    profile_x = float(np.nanmedian(profiles["x_profile_mm"]))
    aoas = sorted(fields["aoa_deg"].unique())
    vmin, vmax = field_scale(fields)

    fig, axes = plt.subplots(1, len(aoas), figsize=(12.0, 3.8), sharex=True, sharey=True, constrained_layout=True)
    contour = None
    for ax, aoa in zip(axes, aoas):
        subset = fields[fields["aoa_deg"] == aoa]
        contour = draw_vector_field(ax, subset, f"Mean PIV, AoA {int(aoa)} deg", vmin, vmax)
        ax.axvline(profile_x, color="white", lw=1.0, ls="--", alpha=0.75)
    axes[0].set_ylabel("y [mm]")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, "03_piv_mean_fields")


def plot_04_instantaneous_piv_aoa15() -> None:
    fields = read_table("piv_instantaneous_fields.csv")
    vmin, vmax = field_scale(fields)
    fig, ax = plt.subplots(figsize=(5.6, 4.2), constrained_layout=True)
    contour = draw_vector_field(ax, fields, "Instantaneous PIV field, AoA 15 deg", vmin, vmax)
    ax.set_ylabel("y [mm]")
    cbar = fig.colorbar(contour, ax=ax, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, "04_instantaneous_piv_aoa15")


def plot_05_piv_window_size_sensitivity() -> None:
    fields = read_table("piv_window_size_fields.csv")
    windows = [16, 32, 64]
    vmin, vmax = field_scale(fields)
    fig, axes = plt.subplots(1, len(windows), figsize=(12.0, 3.8), sharex=True, sharey=True, constrained_layout=True)
    contour = None
    for ax, window in zip(axes, windows):
        subset = fields[fields["window_size_px"] == window]
        contour = draw_vector_field(ax, subset, f"{window} x {window} px", vmin, vmax)
    axes[0].set_ylabel("y [mm]")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, "05_piv_window_size_sensitivity")


def plot_06_piv_overlap_multipass_sensitivity() -> None:
    fields = read_table("piv_processing_fields.csv")
    labels = ["0% overlap, single pass", "50% overlap, single pass", "50% overlap, 3-pass"]
    vmin, vmax = field_scale(fields)
    fig, axes = plt.subplots(1, len(labels), figsize=(12.0, 3.8), sharex=True, sharey=True, constrained_layout=True)
    contour = None
    for ax, label in zip(axes, labels):
        subset = fields[fields["processing"] == label]
        contour = draw_vector_field(ax, subset, label, vmin, vmax)
    axes[0].set_ylabel("y [mm]")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, "06_piv_overlap_multipass_sensitivity")


def plot_07_piv_ensemble_size_effect() -> None:
    fields = read_table("piv_ensemble_fields.csv")
    labels = ["10 sequential samples", "10 separated samples", "100 samples"]
    vmin, vmax = field_scale(fields)
    fig, axes = plt.subplots(1, len(labels), figsize=(12.0, 3.8), sharex=True, sharey=True, constrained_layout=True)
    contour = None
    for ax, label in zip(axes, labels):
        subset = fields[fields["processing"] == label]
        contour = draw_vector_field(ax, subset, label, vmin, vmax)
    axes[0].set_ylabel("y [mm]")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, "07_piv_ensemble_size_effect")


def plot_08_piv_delta_t_effect() -> None:
    fields = read_table("piv_delta_t_fields.csv")
    labels = ["Original Delta t", "Short Delta t"]
    vmin, vmax = field_scale(fields)
    fig, axes = plt.subplots(1, len(labels), figsize=(9.4, 3.8), sharex=True, sharey=True, constrained_layout=True)
    contour = None
    for ax, label in zip(axes, labels):
        subset = fields[fields["processing"] == label]
        contour = draw_vector_field(ax, subset, label, vmin, vmax)
    axes[0].set_ylabel("y [mm]")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, "08_piv_delta_t_effect")


def plot_09_self_made_piv_vs_davis() -> None:
    fields = read_table("piv_self_vs_davis_fields.csv")
    valid = valid_mask(fields)
    values = fields.loc[valid, ["davis_u_streamwise_ms", "self_u_streamwise_ms"]].to_numpy(float).ravel()
    vmin, vmax = np.nanpercentile(values[np.isfinite(values)], [2, 98])
    diff_lim = np.nanpercentile(np.abs(fields.loc[valid, "difference_u_ms"].to_numpy(float)), 98)

    panels = [
        ("davis_u_streamwise_ms", "DaVis u", "viridis", vmin, vmax),
        ("self_u_streamwise_ms", "Self-made u", "viridis", vmin, vmax),
        ("difference_u_ms", "Self-made - DaVis", "coolwarm", -diff_lim, diff_lim),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharex=True, sharey=True, constrained_layout=True)
    for ax, (column, title, cmap, lo, hi) in zip(axes, panels):
        X, Y, Z = value_grid(fields, column)
        contour = ax.contourf(X, Y, Z, levels=np.linspace(lo, hi, 24), cmap=cmap, extend="both")
        ax.set_title(title)
        ax.set_xlabel("x [mm]")
        ax.set_aspect("equal", adjustable="box")
        fig.colorbar(contour, ax=ax, shrink=0.82)
    axes[0].set_ylabel("y [mm]")
    save(fig, "09_self_made_piv_vs_davis")


def plot_10_hwa_mean_rms_wake_profiles() -> None:
    profiles = read_table("hwa_profiles.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.1), sharey=True)
    for aoa, group in profiles.groupby("aoa_deg"):
        group = group.sort_values("y_relative_mm")
        color = COLORS.get(int(aoa))
        label = f"AoA {int(aoa)} deg"
        axes[0].plot(group["mean_velocity_ms"], group["y_relative_mm"], "-o", ms=3.5, color=color, label=label)
        axes[1].plot(group["rms_velocity_ms"], group["y_relative_mm"], "-o", ms=3.5, color=color, label=label)
    axes[0].set_title("Mean velocity")
    axes[1].set_title("RMS fluctuation")
    axes[0].set_xlabel("Mean velocity [m/s]")
    axes[1].set_xlabel("RMS velocity [m/s]")
    axes[0].set_ylabel("y relative to trailing edge [mm]")
    for ax in axes:
        ax.axhline(0, color="0.25", lw=0.8)
    axes[0].legend()
    save(fig, "10_hwa_mean_rms_wake_profiles")


def plot_11_hwa_energy_spectra() -> None:
    spectra = read_table("hwa_spectra.csv")
    peaks = read_table("hwa_spectral_peaks.csv")
    spectra = spectra[(spectra["frequency_hz"] >= 1) & (spectra["frequency_hz"] <= 5000)]

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for label, group in spectra.groupby("case_label"):
        aoa = int(group["aoa_deg"].iloc[0])
        ax.loglog(group["frequency_hz"], group["phi_uu"], color=COLORS.get(aoa), label=label)
    for _, peak in peaks.dropna(subset=["frequency_hz"]).iterrows():
        if 1 <= peak["frequency_hz"] <= 5000:
            ax.axvline(peak["frequency_hz"], color=COLORS.get(int(peak["aoa_deg"]), "0.5"), alpha=0.12)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Power spectral density [arb. units]")
    ax.set_title("HWA energy spectra")
    ax.legend()
    save(fig, "11_hwa_energy_spectra")


def plot_12_piv_hwa_profile_comparison_xc12() -> None:
    comparison = read_table("method_comparison_profiles.csv")
    comparison = comparison[valid_mask(comparison)]
    differences = read_table("method_comparison_at_hwa_points.csv")

    aoas = sorted(comparison["aoa_deg"].unique())
    fig, axes = plt.subplots(2, len(aoas), figsize=(12.0, 6.4), sharey=True)
    markers = {"HWA": "o", "PIV": "s"}
    for col, aoa in enumerate(aoas):
        subset = comparison[comparison["aoa_deg"] == aoa]
        for method, group in subset.groupby("method"):
            group = group.sort_values("y_relative_mm")
            marker = markers.get(method, "o")
            axes[0, col].plot(group["mean_velocity_ms"], group["y_relative_mm"], marker=marker, ms=3.5, label=method)
            axes[1, col].plot(group["rms_velocity_ms"], group["y_relative_mm"], marker=marker, ms=3.5, label=method)
        diff = differences[differences["aoa_deg"] == aoa]
        if not diff.empty:
            rmse = np.sqrt(np.nanmean(diff["piv_minus_hwa_mean_ms"] ** 2))
            axes[0, col].text(0.03, 0.04, f"Mean RMSE: {rmse:.2f} m/s", transform=axes[0, col].transAxes)
        axes[0, col].set_title(f"AoA {int(aoa)} deg")
        axes[1, col].set_xlabel("Velocity [m/s]")
        axes[0, col].axhline(0, color="0.25", lw=0.8)
        axes[1, col].axhline(0, color="0.25", lw=0.8)
    axes[0, 0].set_ylabel("y relative to trailing edge [mm]")
    axes[1, 0].set_ylabel("y relative to trailing edge [mm]")
    axes[0, 0].legend()
    axes[0, 0].text(0.02, 1.05, "Mean velocity", transform=axes[0, 0].transAxes, weight="bold")
    axes[1, 0].text(0.02, 1.05, "RMS fluctuation", transform=axes[1, 0].transAxes, weight="bold")
    save(fig, "12_piv_hwa_profile_comparison_xc12")


def main() -> None:
    clear_figures()
    plot_01_hwa_calibration_curve()
    plot_02_hwa_autocorrelation_convergence()
    plot_03_piv_mean_fields()
    plot_04_instantaneous_piv_aoa15()
    plot_05_piv_window_size_sensitivity()
    plot_06_piv_overlap_multipass_sensitivity()
    plot_07_piv_ensemble_size_effect()
    plot_08_piv_delta_t_effect()
    plot_09_self_made_piv_vs_davis()
    plot_10_hwa_mean_rms_wake_profiles()
    plot_11_hwa_energy_spectra()
    plot_12_piv_hwa_profile_comparison_xc12()
    print(f"12 requested figures written to {FIGURES}")


if __name__ == "__main__":
    main()
