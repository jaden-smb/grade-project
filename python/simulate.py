#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy import ndimage
from scipy.optimize import curve_fit
import time
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import lbm_shan_chen
except ImportError:
    print("Error: Could not import lbm_shan_chen module.")
    print("Please compile the C++ extension first using the build instructions.")
    sys.exit(1)


def create_density_plot(rho, nx, ny, title="Density Field", vmin=None, vmax=None):
    """2D density heatmap."""
    if isinstance(rho, list) or (isinstance(rho, np.ndarray) and rho.ndim == 1):
        rho_2d = np.array(rho).reshape((ny, nx))
    else:
        rho_2d = rho

    if vmin is None: vmin = rho_2d.min()
    if vmax is None: vmax = rho_2d.max()

    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(rho_2d, cmap='viridis', origin='lower',
                   interpolation='bilinear', vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('X (lattice units)', fontsize=12)
    ax.set_ylabel('Y (lattice units)', fontsize=12)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Density', fontsize=12)
    plt.tight_layout()
    return fig, ax, im


def create_3d_plot(rho, nx, ny, title="Density Field — Surface View", vmin=None, vmax=None):
    """Surface plot of the 2D density field."""
    if isinstance(rho, list) or (isinstance(rho, np.ndarray) and rho.ndim == 1):
        rho_2d = np.array(rho).reshape((ny, nx))
    else:
        rho_2d = rho

    if vmin is None: vmin = rho_2d.min()
    if vmax is None: vmax = rho_2d.max()

    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, rho_2d, cmap='viridis',
                          linewidth=0, antialiased=True,
                          vmin=vmin, vmax=vmax, alpha=0.9)
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('X (lattice units)', fontsize=12)
    ax.set_ylabel('Y (lattice units)', fontsize=12)
    ax.set_zlabel('Density', fontsize=12)
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=20, label='Density')
    ax.view_init(elev=45, azim=45)
    plt.tight_layout()
    return fig, ax, surf


def compute_metrics(rho_2d, rho_threshold):
    total_mass = float(rho_2d.sum())
    max_density = float(rho_2d.max())
    liquid_mask = rho_2d > rho_threshold
    area = int(liquid_mask.sum())
    radius = float(np.sqrt(area / np.pi)) if area > 0 else 0.0

    # Circularity: 4*pi*area / perimeter^2 (1.0 = perfect circle, ~0.785 = square)
    if area > 20:
        # Perimeter = count of liquid cells adjacent to at least one gas cell
        eroded = ndimage.binary_erosion(liquid_mask)
        perimeter = float((liquid_mask & ~eroded).sum())
        if perimeter > 0:
            circularity = 4.0 * np.pi * area / (perimeter * perimeter)
            circularity = min(circularity, 1.0)  # cap at 1.0
        else:
            circularity = 1.0
    else:
        circularity = 0.0

    if area > 20:
        ys, xs = np.where(liquid_mask)
        dx, dy = xs - xs.mean(), ys - ys.mean()
        Ixx = float((dx ** 2).mean())
        Iyy = float((dy ** 2).mean())
        Ixy = float((dx * dy).mean())
        disc = np.sqrt(((Ixx - Iyy) / 2) ** 2 + Ixy ** 2)
        lam1 = (Ixx + Iyy) / 2 + disc
        lam2 = (Ixx + Iyy) / 2 - disc
        aspect_ratio = float(lam2 / lam1) if lam1 > 1e-10 else 1.0
    else:
        aspect_ratio = 0.0

    return total_mass, radius, max_density, aspect_ratio, circularity


def plot_metrics(metrics_history, output_dir):
    steps = [m[0] for m in metrics_history]
    masses = [m[1] for m in metrics_history]
    radii = [m[2] for m in metrics_history]
    max_dens = [m[3] for m in metrics_history]
    aspects = [m[4] for m in metrics_history]
    circs = [m[5] for m in metrics_history]
    clamps = [m[6] for m in metrics_history]

    fig, axes = plt.subplots(3, 2, figsize=(14, 14))

    axes[0, 0].plot(steps, radii, 'b-', linewidth=1.5)
    axes[0, 0].set_title('Effective Radius - size change', fontweight='bold')
    axes[0, 0].set_xlabel('Step'); axes[0, 0].set_ylabel('Radius (lu)')

    axes[0, 1].plot(steps, masses, 'g-', linewidth=1.5)
    axes[0, 1].set_title('Total Mass - conservation check', fontweight='bold')
    axes[0, 1].set_xlabel('Step'); axes[0, 1].set_ylabel('Mass')

    axes[1, 0].plot(steps, max_dens, 'r-', linewidth=1.5)
    axes[1, 0].set_title('Max Density - evaporation indicator', fontweight='bold')
    axes[1, 0].set_xlabel('Step'); axes[1, 0].set_ylabel('Density')

    axes[1, 1].plot(steps, aspects, 'm-', linewidth=1.5)
    axes[1, 1].axhline(1.0, color='k', linestyle='--', linewidth=0.8, label='circular')
    axes[1, 1].set_title('Aspect Ratio - deformation / oscillation', fontweight='bold')
    axes[1, 1].set_xlabel('Step'); axes[1, 1].set_ylabel('Ratio (1= circular)')

    axes[2, 0].plot(steps, circs, 'c-', linewidth=1.5)
    axes[2, 0].axhline(1.0, color='k', linestyle='--', linewidth=0.8, label='perfect circle')
    axes[2, 0].axhline(0.785, color='gray', linestyle=':', linewidth=0.8, label='square')
    axes[2, 0].set_title('Circularity - shape quality', fontweight='bold')
    axes[2, 0].set_xlabel('Step'); axes[2, 0].set_ylabel('Circularity (1=circle)')
    axes[2, 0].legend(fontsize=9)

    axes[2, 1].plot(steps, clamps, 'orange', linewidth=1.5)
    axes[2, 1].set_title('Velocity Clamp Count - stability indicator', fontweight='bold')
    axes[2, 1].set_xlabel('Step'); axes[2, 1].set_ylabel('Clamped cells per interval')

    plt.suptitle('Droplet Assessment Metrics', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/metrics.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def run_simulation(nx=300, ny=300, tau=1.0, G=-5.0, G_final=None,
                   rho_liquid=2.35, rho_gas=0.05,
                   center_x=None, center_y=None, radius=40.0,
                   num_steps=5000, save_plots=True, animate=True,
                   view_surface=True, output_dir="output"):
    """Run LBM simulation and visualize density evolution."""

    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}/")

    if center_x is None: center_x = nx // 2
    if center_y is None: center_y = ny // 2

    print("=" * 60)
    print("Two-Phase LBM Simulation - Shan-Chen Model")
    print("=" * 60)
    print(f"Grid size: {nx} x {ny}")
    print(f"tau={tau}, G={G}, rho_liquid={rho_liquid}, rho_gas={rho_gas}")
    print(f"Droplet: center=({center_x},{center_y}), radius={radius}")
    print(f"Steps: {num_steps}")
    print("=" * 60)

    print("\nInitializing LBM simulator...")
    lbm = lbm_shan_chen.LBMShanChen(nx, ny, tau, G, rho_liquid, rho_gas)
    lbm.initialize_droplet(center_x, center_y, radius)
    rho_threshold = (rho_liquid + rho_gas) / 2.0
    metrics_history = []
    G_current = G

    rho_initial = np.array(lbm.get_density()).reshape((ny, nx))
    density_min = rho_initial.min()
    density_max = rho_initial.max()

    print("\nPlotting initial state...")
    fig_init, ax_init, im_init = create_density_plot(
        rho_initial, nx, ny,
        f"Initial State (t=0)\nDroplet at ({center_x}, {center_y})"
    )
    if save_plots:
        plt.savefig(f'{output_dir}/density_initial.png', dpi=150, bbox_inches='tight')
        print(f"Saved: {output_dir}/density_initial.png")
        if view_surface:
            fig_surf_init, _, _ = create_3d_plot(
                rho_initial, nx, ny, "Initial State — Surface View (t=0)")
            plt.savefig(f'{output_dir}/density_initial_surface.png', dpi=150, bbox_inches='tight')
            plt.close(fig_surf_init)
            print(f"Saved: {output_dir}/density_initial_surface.png")
    plt.close(fig_init)

    if animate:
        frames = [rho_initial.copy()]
        frames_surface = [rho_initial.copy()] if view_surface else None

    print(f"\nRunning simulation for {num_steps} steps...")
    print("Progress: ", end="", flush=True)

    save_interval = max(1, num_steps // 20)
    frame_interval = max(1, num_steps // 200)

    lbm.reset_clamp_count()
    t_start = time.perf_counter()

    for step in range(1, num_steps + 1):
        if G_final is not None:
            G_current = G + (G_final - G) * step / num_steps
            lbm.set_G(G_current)
        lbm.step()

        if step % (num_steps // 20) == 0:
            print(".", end="", flush=True)

        rho = np.array(lbm.get_density()).reshape((ny, nx))
        density_min = min(density_min, rho.min())
        density_max = max(density_max, rho.max())

        if save_plots and step % save_interval == 0:
            fig, ax, im = create_density_plot(rho, nx, ny, f"Time Step {step}",
                                              vmin=density_min, vmax=density_max)
            plt.savefig(f'{output_dir}/density_step_{step:05d}.png', dpi=150, bbox_inches='tight')
            plt.close(fig)

            if view_surface:
                fig_surf, _, _ = create_3d_plot(rho, nx, ny, f"Time Step {step} — Surface View",
                                                vmin=density_min, vmax=density_max)
                plt.savefig(f'{output_dir}/density_step_{step:05d}_surface.png', dpi=150, bbox_inches='tight')
                plt.close(fig_surf)

        if animate and step % frame_interval == 0:
            clamp = lbm.get_clamp_count()
            lbm.reset_clamp_count()
            if clamp > 0:
                print(f"\n  [Warning] Step {step}: velocity clamped in {clamp} cells")
            frames.append(rho.copy())
            m = compute_metrics(rho, rho_threshold)
            metrics_history.append((step, *m, clamp))

            if view_surface and frames_surface is not None:
                frames_surface.append(rho.copy())

    t_end = time.perf_counter()
    elapsed = t_end - t_start

    print(" Done!")

    rho_final = np.array(lbm.get_density()).reshape((ny, nx))

    print("\nPlotting final state...")
    fig_final, _, _ = create_density_plot(
        rho_final, nx, ny, f"Final State (t={num_steps})",
        vmin=density_min, vmax=density_max
    )
    if save_plots:
        plt.savefig(f'{output_dir}/density_final.png', dpi=150, bbox_inches='tight')
        print(f"Saved: {output_dir}/density_final.png")
        if view_surface:
            fig_surf_final, _, _ = create_3d_plot(rho_final, nx, ny,
                                                   f"Final State — Surface View (t={num_steps})",
                                                   vmin=density_min, vmax=density_max)
            plt.savefig(f'{output_dir}/density_final_surface.png', dpi=150, bbox_inches='tight')
            plt.close(fig_surf_final)
            print(f"Saved: {output_dir}/density_final_surface.png")
    plt.close(fig_final)

    if animate and len(frames) > 1:
        print("\nCreating 2D animation...")
        fig_anim, ax_anim = plt.subplots(figsize=(10, 10))
        im_anim = ax_anim.imshow(frames[0], cmap='viridis', origin='lower',
                                 interpolation='bilinear',
                                 vmin=density_min, vmax=density_max)
        ax_anim.set_title('Density Evolution', fontsize=14, fontweight='bold')
        ax_anim.set_xlabel('X (lattice units)', fontsize=12)
        ax_anim.set_ylabel('Y (lattice units)', fontsize=12)
        cbar_anim = plt.colorbar(im_anim, ax=ax_anim)
        cbar_anim.set_label('Density', fontsize=12)

        def animate_frame(frame_num):
            im_anim.set_array(frames[frame_num])
            if frame_num < len(metrics_history):
                step_n, tm, r, md, ar, circ, clamp_n = metrics_history[frame_num]
                ax_anim.set_title(
                    f'Step {step_n} | R={r:.1f} lu | '
                    f'M={tm:.0f} | pmax={md:.3f} | AR={ar:.3f} | C={circ:.3f}',
                    fontsize=11, fontweight='bold')
            return [im_anim]

        anim = animation.FuncAnimation(fig_anim, animate_frame,
                                      frames=len(frames),
                                      interval=50, blit=True, repeat=True)
        if save_plots:
            print("Saving 2D animation...")
            anim.save(f'{output_dir}/density_evolution.gif', writer='pillow', fps=20)
            print(f"Saved: {output_dir}/density_evolution.gif")
        plt.show()

        if view_surface and frames_surface is not None and len(frames_surface) > 1:
            print("\nCreating surface view animation...")
            X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
            fig_surf_anim = plt.figure(figsize=(12, 10))
            ax_surf_anim = fig_surf_anim.add_subplot(111, projection='3d')
            surf_anim = ax_surf_anim.plot_surface(X, Y, frames_surface[0], cmap='viridis',
                                                   linewidth=0, antialiased=True,
                                                   vmin=density_min, vmax=density_max, alpha=0.9)
            ax_surf_anim.set_title('Density Evolution — Surface View', fontsize=14, fontweight='bold', pad=20)
            ax_surf_anim.set_xlabel('X (lattice units)', fontsize=12)
            ax_surf_anim.set_ylabel('Y (lattice units)', fontsize=12)
            ax_surf_anim.set_zlabel('Density', fontsize=12)
            ax_surf_anim.view_init(elev=45, azim=45)
            fig_surf_anim.colorbar(surf_anim, ax=ax_surf_anim, shrink=0.5, aspect=20, label='Density')

            def animate_frame_surface(frame_num):
                ax_surf_anim.clear()
                surf = ax_surf_anim.plot_surface(X, Y, frames_surface[frame_num], cmap='viridis',
                                                  linewidth=0, antialiased=True,
                                                  vmin=density_min, vmax=density_max, alpha=0.9)
                ax_surf_anim.set_title(
                    f'Density Evolution — Surface View - Step {frame_num * frame_interval} (Frame {frame_num}/{len(frames_surface)-1})',
                    fontsize=14, fontweight='bold', pad=20)
                ax_surf_anim.set_xlabel('X (lattice units)', fontsize=12)
                ax_surf_anim.set_ylabel('Y (lattice units)', fontsize=12)
                ax_surf_anim.set_zlabel('Density', fontsize=12)
                ax_surf_anim.view_init(elev=45, azim=45)
                return [surf]

            anim_surf = animation.FuncAnimation(fig_surf_anim, animate_frame_surface,
                                                frames=len(frames_surface),
                                                interval=50, blit=False, repeat=True)
            if save_plots:
                print("Saving surface view animation...")
                anim_surf.save(f'{output_dir}/density_evolution_surface.gif', writer='pillow', fps=20)
                print(f"Saved: {output_dir}/density_evolution_surface.gif")
            plt.show()

        if metrics_history:
            plot_metrics(metrics_history, output_dir)

    print("\n" + "=" * 60)
    print("Simulation Statistics:")
    print("=" * 60)
    print(f"Initial density range: [{rho_initial.min():.3f}, {rho_initial.max():.3f}]")
    print(f"Final density range:   [{rho_final.min():.3f}, {rho_final.max():.3f}]")
    print(f"Density contrast: {rho_final.max() / rho_final.min():.2f}")
    initial_mass = float(rho_initial.sum())
    final_mass = float(rho_final.sum())
    mass_drift_pct = 100.0 * (final_mass - initial_mass) / initial_mass
    print(f"Mass drift:            {mass_drift_pct:+.4f}%")
    print(f"Wall-clock time:       {elapsed:.2f} s")
    print(f"Time per step:         {elapsed/num_steps*1000:.2f} ms")
    print(f"Throughput:            {nx*ny*num_steps/elapsed/1e6:.2f} MLUPS")
    print("=" * 60)

    return lbm, rho_final, elapsed


def run_coexistence_sweep(nx=150, ny=150, tau=1.0, radius=35, num_steps=2000,
                          rho_liquid=2.0, rho_gas=0.1,
                          G_values=None, output_dir="output/scenario_c_coexistence"):
    """Run a parameter sweep over G values and plot the coexistence curve."""
    if G_values is None:
        G_values = [-4.5, -5.0, -5.5, -6.0, -6.5]

    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}/")
    print(f"Grid: {nx}x{ny}, tau={tau}, radius={radius}, steps={num_steps}")
    print(f"G values: {G_values}")

    eq_rho_liquid = []
    eq_rho_gas = []

    for G in G_values:
        print(f"\n  Running G = {G} ...")
        lbm = lbm_shan_chen.LBMShanChen(nx, ny, tau, G, rho_liquid, rho_gas)
        lbm.initialize_droplet(nx // 2, ny // 2, radius)

        t_start = time.perf_counter()
        for _ in range(num_steps):
            lbm.step()
        elapsed = time.perf_counter() - t_start

        rho_final = np.array(lbm.get_density()).reshape((ny, nx))
        eq_rho_liquid.append(float(rho_final.max()))
        eq_rho_gas.append(float(rho_final.min()))
        print(f"    rho_liq={rho_final.max():.4f}, rho_gas={rho_final.min():.4f}, "
              f"t={elapsed:.2f}s ({nx*ny*num_steps/elapsed/1e6:.2f} MLUPS)")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(G_values, eq_rho_liquid, 'b-o', linewidth=1.5, markersize=7, label='Liquid phase (max ρ)')
    ax.plot(G_values, eq_rho_gas, 'r-o', linewidth=1.5, markersize=7, label='Gas phase (min ρ)')
    ax.set_xlabel('G (cohesion parameter)', fontsize=12)
    ax.set_ylabel('Equilibrium Density', fontsize=12)
    ax.set_title('Coexistence Curve: Liquid and Gas Densities vs. G',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/coexistence_curve.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {output_dir}/coexistence_curve.png")

    return G_values, eq_rho_liquid, eq_rho_gas


def plot_density_profile(rho_2d, nx, ny, output_dir):
    """Horizontal cross-section through droplet centre with optional tanh fit."""
    profile = rho_2d[ny // 2, :]
    x = np.arange(nx)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x, profile, 'b-', linewidth=1.5, label='ρ(x)')

    # Attempt tanh fit: rho(x) = rho_avg + rho_diff * tanh((x - x0) / w)
    # Use two-interface model: droplet has left and right edges
    try:
        rho_liq = float(profile.max())
        rho_gas = float(profile.min())
        rho_avg = 0.5 * (rho_liq + rho_gas)
        rho_diff = 0.5 * (rho_liq - rho_gas)
        x0_left  = float(x[profile > rho_avg][0])  if (profile > rho_avg).any() else nx * 0.25
        x0_right = float(x[profile > rho_avg][-1]) if (profile > rho_avg).any() else nx * 0.75

        def two_interface(xx, w):
            return rho_gas + rho_diff * (
                np.tanh((xx - x0_left)  / w) -
                np.tanh((xx - x0_right) / w)
            )

        popt, _ = curve_fit(two_interface, x, profile, p0=[5.0], bounds=(0.5, 50.0))
        interface_width = float(popt[0])
        ax.plot(x, two_interface(x, *popt), 'r--', linewidth=1.2,
                label=f'tanh fit (w={interface_width:.1f} lu)')
        ax.legend(fontsize=10)
    except Exception:
        pass

    ax.set_xlabel('X (lattice units)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Density Cross-Section (y = Ny/2)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/density_profile.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_dir}/density_profile.png")


def plot_velocity_field(lbm, nx, ny, output_dir):
    """Density heatmap overlaid with velocity quiver — visualises spurious currents."""
    rho_2d = np.array(lbm.get_density()).reshape((ny, nx))
    ux_2d  = np.array(lbm.get_velocity_x()).reshape((ny, nx))
    uy_2d  = np.array(lbm.get_velocity_y()).reshape((ny, nx))

    speed = np.sqrt(ux_2d**2 + uy_2d**2)
    max_speed = float(speed.max())
    print(f"  Max spurious velocity magnitude: {max_speed:.6f} lu/step")

    step = 8  # subsample every 8 cells for denser coverage
    xs = np.arange(0, nx, step)
    ys = np.arange(0, ny, step)
    X, Y = np.meshgrid(xs, ys)
    Ux = ux_2d[::step, ::step]
    Uy = uy_2d[::step, ::step]
    Sp = speed[::step, ::step]

    # Normalize arrows to unit length; color encodes magnitude
    mag = np.sqrt(Ux**2 + Uy**2)
    mag_safe = np.where(mag > 1e-15, mag, 1.0)
    Ux_norm = Ux / mag_safe
    Uy_norm = Uy / mag_safe

    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(rho_2d, cmap='viridis', origin='lower', interpolation='bilinear')
    quiv = ax.quiver(X, Y, Ux_norm, Uy_norm, Sp,
                     cmap='hot', alpha=0.8, scale=30, width=0.003,
                     headwidth=3, headlength=4)
    ax.set_title(f'Velocity Field (spurious currents)\nmax|u|={max_speed:.2e} lu/step',
                 fontsize=13, fontweight='bold')
    ax.set_xlabel('X (lattice units)', fontsize=12)
    ax.set_ylabel('Y (lattice units)', fontsize=12)
    plt.colorbar(im, ax=ax, label='Density', location='left', pad=0.08)
    plt.colorbar(quiv, ax=ax, label='|u| (lu/step)')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/velocity_field.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_dir}/velocity_field.png")


def run_laplace_test(tau=0.8, G=-5.0, rho_liquid=2.0, rho_gas=0.1,
                     radii=None, nx=200, ny=200, num_steps=4000,
                     output_dir="output/scenario_d_laplace"):
    """Laplace pressure test: Δp vs 1/R should be linear with slope = surface tension σ."""
    if radii is None:
        radii = [20, 30, 40, 50, 60]

    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}/")
    print(f"Grid: {nx}x{ny}, tau={tau}, G={G}, steps={num_steps}")
    print(f"Radii: {radii}")

    inv_R_list = []
    delta_p_list = []

    for R in radii:
        print(f"\n  Radius R={R} ...")
        lbm = lbm_shan_chen.LBMShanChen(nx, ny, tau, G, rho_liquid, rho_gas)
        lbm.initialize_droplet(nx // 2, ny // 2, float(R))

        t_start = time.perf_counter()
        for _ in range(num_steps):
            lbm.step()
        elapsed = time.perf_counter() - t_start

        rho_arr = np.array(lbm.get_density()).reshape((ny, nx))
        rho_max = float(rho_arr.max())
        rho_min = float(rho_arr.min())

        # Effective radius from liquid area
        rho_threshold = (rho_max + rho_min) / 2.0
        area = float((rho_arr > rho_threshold).sum())
        R_eff = float(np.sqrt(area / np.pi)) if area > 0 else float(R)

        # Shan-Chen EoS: p = rho*cs^2 + G*cs^2/2 * psi(rho)^2
        cs2 = 1.0 / 3.0
        def psi(rho):
            return 1.0 - np.exp(-1.5 * rho)
        p_liq = rho_max * cs2 + G * cs2 / 2.0 * psi(rho_max)**2
        p_gas = rho_min * cs2 + G * cs2 / 2.0 * psi(rho_min)**2
        delta_p = p_liq - p_gas

        print(f"    rho_liq={rho_max:.4f}, rho_gas={rho_min:.4f}, "
              f"R_eff={R_eff:.1f}, Δp={delta_p:.5f}, t={elapsed:.1f}s")

        if R_eff > 0:
            inv_R_list.append(1.0 / R_eff)
            delta_p_list.append(delta_p)

    if len(inv_R_list) < 2:
        print("Not enough data points for Laplace fit.")
        return

    inv_R = np.array(inv_R_list)
    delta_p = np.array(delta_p_list)

    # Linear fit: Δp = σ / R  →  Δp = σ * (1/R)
    coeffs = np.polyfit(inv_R, delta_p, 1)
    sigma = coeffs[0]
    print(f"\n  Fitted surface tension σ = {sigma:.4f} (lattice units)")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(inv_R, delta_p, color='b', s=60, zorder=3, label='Simulated Δp')
    x_fit = np.linspace(0, inv_R.max() * 1.1, 100)
    ax.plot(x_fit, np.polyval(coeffs, x_fit), 'r--', linewidth=1.5,
            label=f'Linear fit  σ = {sigma:.4f}')
    ax.set_xlabel('1 / R_eff  (lu⁻¹)', fontsize=12)
    ax.set_ylabel('Δp = p_liq − p_gas', fontsize=12)
    ax.set_title('Laplace Pressure Test\n(Young–Laplace: Δp = σ/R)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/laplace_test.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {output_dir}/laplace_test.png")

    return inv_R_list, delta_p_list, sigma


def main():
    # ------------------------------------------------------------------ #
    # Scenario A — Steady-state droplet equilibrium (primary validation)  #
    # ------------------------------------------------------------------ #
    print("\n" + "=" * 60)
    print("SCENARIO A — Steady-State Droplet Equilibrium")
    print("=" * 60)
    lbm_a, rho_a, elapsed_a = run_simulation(
        nx=200, ny=200,
        tau=1.0,
        G=-5.0,
        rho_liquid=2.0, rho_gas=0.1,
        radius=40,
        num_steps=5000,
        save_plots=True,
        animate=True,
        view_surface=True,
        output_dir="output/scenario_a_equilibrium",
    )
    rho_a_2d = np.array(lbm_a.get_density()).reshape((200, 200))
    plot_density_profile(rho_a_2d, 200, 200, "output/scenario_a_equilibrium")
    plot_velocity_field(lbm_a, 200, 200, "output/scenario_a_equilibrium")

    # ------------------------------------------------------------------ #
    # Scenario B — G-ramp evaporation analogy                             #
    # ------------------------------------------------------------------ #
    print("\n" + "=" * 60)
    print("SCENARIO B — G-Ramp Evaporation Analogy")
    print("=" * 60)
    lbm_b, rho_b, elapsed_b = run_simulation(
        nx=200, ny=200,
        tau=1.0,
        G=-5.0, G_final=-3.7,   # ramp past critical G_c ≈ -4.0
        rho_liquid=2.0, rho_gas=0.1,
        radius=40,
        num_steps=8000,          # more steps to see full dissolution
        save_plots=True,
        animate=True,
        view_surface=True,
        output_dir="output/scenario_b_evaporation",
    )
    rho_b_2d = np.array(lbm_b.get_density()).reshape((200, 200))
    plot_density_profile(rho_b_2d, 200, 200, "output/scenario_b_evaporation")
    plot_velocity_field(lbm_b, 200, 200, "output/scenario_b_evaporation")

    # ------------------------------------------------------------------ #
    # Scenario C — Parameter sweep for coexistence curve                  #
    # ------------------------------------------------------------------ #
    print("\n" + "=" * 60)
    print("SCENARIO C — Parameter Sweep for Coexistence Curve")
    print("=" * 60)
    run_coexistence_sweep(
        nx=150, ny=150,
        tau=1.0,
        radius=35,
        num_steps=3000,
        rho_liquid=2.0, rho_gas=0.1,
        G_values=[-4.0, -4.5, -5.0, -5.5, -6.0],
        output_dir="output/scenario_c_coexistence",
    )

    # ------------------------------------------------------------------ #
    # Scenario D — Laplace pressure test (surface tension measurement)    #
    # ------------------------------------------------------------------ #
    print("\n" + "=" * 60)
    print("SCENARIO D — Laplace Pressure Test")
    print("=" * 60)
    run_laplace_test(
        tau=1.0,
        G=-5.0,
        rho_liquid=2.0, rho_gas=0.1,
        radii=[20, 30, 40, 50, 60],
        nx=200, ny=200,
        num_steps=5000,
        output_dir="output/scenario_d_laplace",
    )

    print("\nAll scenarios complete!")


if __name__ == "__main__":
    main()
