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


def plot_hwa_calibration() -> None:
    points = read_table("hwa_calibration_points.csv")
    curve = read_table("hwa_calibration_curve.csv")

    fig, ax = plt.subplots(figsize=(5.6, 3.7))
    ax.scatter(points["signal_v"], points["velocity_ms"], s=22, color="black", label="Calibration points")
    ax.plot(curve["signal_v"], curve["velocity_ms"], color="#1f77b4", label="PCHIP interpolation")
    ax.set_xlabel("CTA output signal [V]")
    ax.set_ylabel("Velocity [m/s]")
    ax.set_title("HWA calibration used by main MATLAB code")
    ax.legend()
    save(fig, "main_hwa_calibration")


def plot_hwa_profiles() -> None:
    profiles = read_table("hwa_profiles.csv")
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.1), sharey=True)

    for aoa, group in profiles.groupby("aoa_deg"):
        group = group.sort_values("y_relative_mm")
        color = COLORS.get(int(aoa))
        label = f"AoA {int(aoa)} deg"
        axes[0].plot(group["mean_velocity_ms"], group["y_relative_mm"], "-o", ms=3.5, color=color, label=label)
        axes[1].plot(group["rms_velocity_ms"], group["y_relative_mm"], "-o", ms=3.5, color=color, label=label)

    axes[0].set_title("HWA mean velocity")
    axes[1].set_title("HWA RMS fluctuations")
    axes[0].set_xlabel("Mean velocity [m/s]")
    axes[1].set_xlabel("RMS velocity [m/s]")
    axes[0].set_ylabel("y relative to trailing edge [mm]")
    for ax in axes:
        ax.axhline(0, color="0.25", lw=0.8)
    axes[0].legend()
    save(fig, "main_hwa_profiles")


def field_grid(df: pd.DataFrame, value: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    work = df.copy()
    invalid = ~valid_mask(work)
    work.loc[invalid, value] = np.nan
    if "v_normal_ms" in work.columns:
        work.loc[invalid, "v_normal_ms"] = np.nan

    z_grid = work.pivot(index="y_mm", columns="x_mm", values=value).sort_index().sort_index(axis=1)
    x = z_grid.columns.to_numpy(float)
    y = z_grid.index.to_numpy(float)
    X, Y = np.meshgrid(x, y)
    Z = z_grid.to_numpy(float)

    if "v_normal_ms" in work.columns:
        v_grid = work.pivot(index="y_mm", columns="x_mm", values="v_normal_ms").sort_index().sort_index(axis=1)
        V = v_grid.to_numpy(float)
    else:
        V = np.zeros_like(Z)
    return X, Y, Z, V


def plot_piv_mean_fields() -> None:
    fields = read_table("piv_mean_fields.csv")
    profiles = read_table("piv_profile_xc12.csv")
    profile_x = float(np.nanmedian(profiles["x_profile_mm"]))
    aoas = sorted(fields["aoa_deg"].unique())
    valid_values = fields.loc[valid_mask(fields), "u_streamwise_ms"].to_numpy(float)
    vmin, vmax = np.nanpercentile(valid_values, [2, 98])

    fig, axes = plt.subplots(1, len(aoas), figsize=(12.0, 3.8), sharex=True, sharey=True, constrained_layout=True)
    if len(aoas) == 1:
        axes = [axes]

    contour = None
    for ax, aoa in zip(axes, aoas):
        subset = fields[fields["aoa_deg"] == aoa]
        X, Y, U, V = field_grid(subset, "u_streamwise_ms")
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
        ax.set_title(f"PIV mean field, AoA {int(aoa)} deg")
        ax.set_xlabel("x [mm]")
        ax.set_xticks([-100, -50, 0, 50])
        ax.axvline(profile_x, color="white", lw=1.0, ls="--", alpha=0.75)
        ax.set_aspect("equal", adjustable="box")
    axes[0].set_ylabel("y [mm]")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.86)
    cbar.set_label("Streamwise velocity [m/s]")
    save(fig, "main_piv_mean_fields")


def plot_method_comparison() -> None:
    comparison = read_table("method_comparison_profiles.csv")
    comparison = comparison[valid_mask(comparison)]
    differences = read_table("method_comparison_at_hwa_points.csv")

    aoas = sorted(comparison["aoa_deg"].unique())
    fig, axes = plt.subplots(2, len(aoas), figsize=(12.0, 6.4), sharey=True)
    if len(aoas) == 1:
        axes = np.array([[axes[0]], [axes[1]]])

    markers = {"HWA": "o", "PIV": "s"}
    for col, aoa in enumerate(aoas):
        subset = comparison[comparison["aoa_deg"] == aoa]
        for method, group in subset.groupby("method"):
            group = group.sort_values("y_relative_mm")
            marker = markers.get(method, "o")
            axes[0, col].plot(
                group["mean_velocity_ms"],
                group["y_relative_mm"],
                marker=marker,
                ms=3.5,
                label=method,
            )
            axes[1, col].plot(
                group["rms_velocity_ms"],
                group["y_relative_mm"],
                marker=marker,
                ms=3.5,
                label=method,
            )

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
    save(fig, "main_hwa_piv_comparison")


def plot_processing_summary() -> None:
    summary = read_table("piv_processing_summary.csv")
    wanted = summary[
        summary["processing"].str.contains("Overlap0SinglePass|Overlap50SinglePass|Overlap50MP3", regex=True)
    ].copy()
    wanted = wanted[wanted["case_id"] == "aoa15"].copy()
    if wanted.empty:
        return

    label_map = {
        "Overlap0SinglePass": "0% overlap\nsingle pass",
        "Overlap50SinglePass": "50% overlap\nsingle pass",
        "Overlap50MP3": "50% overlap\n3-pass",
    }
    wanted["label"] = wanted["processing"].map(label_map).fillna(wanted["processing"])
    wanted = wanted.drop_duplicates("label")

    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.8))
    x = np.arange(len(wanted))
    axes[0].bar(x, wanted["valid_fraction"], color="#4c78a8")
    axes[0].set_ylabel("Valid-vector fraction")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(wanted["label"], rotation=20, ha="right")

    axes[1].bar(x, wanted["spatial_std_u_streamwise_ms"], color="#f58518")
    axes[1].set_ylabel("Spatial std of u [m/s]")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(wanted["label"], rotation=20, ha="right")
    axes[1].set_title("AoA 15 PIV processing comparison")
    save(fig, "main_piv_processing_summary")


def plot_window_size_screenshots() -> None:
    screenshot_dir = ROOT / "data" / "PIV processed" / "Screenshots"
    files = [
        ("16 x 16 px", screenshot_dir / "16x16Test.png"),
        ("32 x 32 px", screenshot_dir / "32x32Test.png"),
        ("64 x 64 px", screenshot_dir / "64x64Test.png"),
    ]
    files = [(label, path) for label, path in files if path.exists()]
    if not files:
        return

    fig, axes = plt.subplots(1, len(files), figsize=(12.2, 4.2), constrained_layout=True)
    if len(files) == 1:
        axes = [axes]
    for ax, (label, path) in zip(axes, files):
        ax.imshow(plt.imread(path))
        ax.set_title(f"DaVis window-size test\n{label}")
        ax.axis("off")
    save(fig, "main_piv_window_size_screenshots")


def main() -> None:
    plot_hwa_calibration()
    plot_hwa_profiles()
    plot_piv_mean_fields()
    plot_method_comparison()
    plot_processing_summary()
    plot_window_size_screenshots()
    print(f"Figures written to {FIGURES}")


if __name__ == "__main__":
    main()
