import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os

from lbm.core import LBMSimulator
from lbm.metrics import compute_metrics
from lbm.visualization import create_density_plot, create_3d_plot, plot_metrics


def run_simulation(nx=300, ny=300, tau=1.0, G=-5.0, G_final=None,
                   rho_liquid=2.35, rho_gas=0.05,
                   center_x=None, center_y=None, radius=40.0,
                   num_steps=5000, save_plots=True, animate=True,
                   view_surface=True, export_obj=False, output_dir="output"):

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
    lbm = LBMSimulator(nx, ny, tau, G, rho_liquid, rho_gas)
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

        if export_obj and len(frames) > 1:
            from lbm.export.obj import export_obj_sequence
            print("\nExporting OBJ sequence...")
            export_obj_sequence(
                frames,
                output_dir=f'{output_dir}/obj_sequence',
                z_scale=50.0,
                subsample=1,
                colormap='viridis',
                rho_min=density_min,
                rho_max=density_max,
            )

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
