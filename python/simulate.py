#!/usr/bin/env python3
"""Two-Phase LBM Shan-Chen simulation with density visualization."""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
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


def create_3d_plot(rho, nx, ny, title="3D Density Field", vmin=None, vmax=None):
    """3D surface plot of the density field."""
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


def run_simulation(nx=200, ny=200, tau=1.0, G=-5.5, 
                   rho_liquid=2.0, rho_gas=0.1,
                   center_x=None, center_y=None, radius=30.0,
                   num_steps=2000, save_plots=True, animate=True,
                   view_3d=True, output_dir="output"):
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
        if view_3d:
            fig_3d_init, _, _ = create_3d_plot(rho_initial, nx, ny, "Initial State 3D (t=0)")
            plt.savefig(f'{output_dir}/density_initial_3d.png', dpi=150, bbox_inches='tight')
            plt.close(fig_3d_init)
            print(f"Saved: {output_dir}/density_initial_3d.png")
    plt.close(fig_init)
    
    if animate:
        frames = [rho_initial.copy()]
        frames_3d = [rho_initial.copy()] if view_3d else None
    
    print(f"\nRunning simulation for {num_steps} steps...")
    print("Progress: ", end="", flush=True)
    
    save_interval = max(1, num_steps // 20)
    frame_interval = max(1, num_steps // 200)
    
    for step in range(1, num_steps + 1):
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
            
            if view_3d:
                fig_3d, _, _ = create_3d_plot(rho, nx, ny, f"Time Step {step} (3D)",
                                              vmin=density_min, vmax=density_max)
                plt.savefig(f'{output_dir}/density_step_{step:05d}_3d.png', dpi=150, bbox_inches='tight')
                plt.close(fig_3d)
        
        if animate and step % frame_interval == 0:
            frames.append(rho.copy())
            if view_3d and frames_3d is not None:
                frames_3d.append(rho.copy())
    
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
        if view_3d:
            fig_3d_final, _, _ = create_3d_plot(rho_final, nx, ny,
                                                f"Final State 3D (t={num_steps})",
                                                vmin=density_min, vmax=density_max)
            plt.savefig(f'{output_dir}/density_final_3d.png', dpi=150, bbox_inches='tight')
            plt.close(fig_3d_final)
            print(f"Saved: {output_dir}/density_final_3d.png")
    plt.show(block=False)
    
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
            ax_anim.set_title(
                f'Density Evolution - Step {frame_num * frame_interval} (Frame {frame_num}/{len(frames)-1})',
                fontsize=14, fontweight='bold')
            return [im_anim]
        
        anim = animation.FuncAnimation(fig_anim, animate_frame, 
                                      frames=len(frames), 
                                      interval=50, blit=True, repeat=True)
        if save_plots:
            print("Saving 2D animation...")
            anim.save(f'{output_dir}/density_evolution.gif', writer='pillow', fps=20)
            print(f"Saved: {output_dir}/density_evolution.gif")
        plt.show()
        
        if view_3d and frames_3d is not None and len(frames_3d) > 1:
            print("\nCreating 3D animation...")
            X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
            fig_3d_anim = plt.figure(figsize=(12, 10))
            ax_3d_anim = fig_3d_anim.add_subplot(111, projection='3d')
            surf_3d_anim = ax_3d_anim.plot_surface(X, Y, frames_3d[0], cmap='viridis',
                                                   linewidth=0, antialiased=True,
                                                   vmin=density_min, vmax=density_max, alpha=0.9)
            ax_3d_anim.set_title('3D Density Evolution', fontsize=14, fontweight='bold', pad=20)
            ax_3d_anim.set_xlabel('X (lattice units)', fontsize=12)
            ax_3d_anim.set_ylabel('Y (lattice units)', fontsize=12)
            ax_3d_anim.set_zlabel('Density', fontsize=12)
            ax_3d_anim.view_init(elev=45, azim=45)
            fig_3d_anim.colorbar(surf_3d_anim, ax=ax_3d_anim, shrink=0.5, aspect=20, label='Density')
            
            def animate_frame_3d(frame_num):
                ax_3d_anim.clear()
                surf = ax_3d_anim.plot_surface(X, Y, frames_3d[frame_num], cmap='viridis',
                                              linewidth=0, antialiased=True,
                                              vmin=density_min, vmax=density_max, alpha=0.9)
                ax_3d_anim.set_title(
                    f'3D Density Evolution - Step {frame_num * frame_interval} (Frame {frame_num}/{len(frames_3d)-1})',
                    fontsize=14, fontweight='bold', pad=20)
                ax_3d_anim.set_xlabel('X (lattice units)', fontsize=12)
                ax_3d_anim.set_ylabel('Y (lattice units)', fontsize=12)
                ax_3d_anim.set_zlabel('Density', fontsize=12)
                ax_3d_anim.view_init(elev=45, azim=45)
                return [surf]
            
            anim_3d = animation.FuncAnimation(fig_3d_anim, animate_frame_3d,
                                             frames=len(frames_3d),
                                             interval=50, blit=False, repeat=True)
            if save_plots:
                print("Saving 3D animation...")
                anim_3d.save(f'{output_dir}/density_evolution_3d.gif', writer='pillow', fps=20)
                print(f"Saved: {output_dir}/density_evolution_3d.gif")
            plt.show()
    
    print("\n" + "=" * 60)
    print("Simulation Statistics:")
    print("=" * 60)
    print(f"Initial density range: [{rho_initial.min():.3f}, {rho_initial.max():.3f}]")
    print(f"Final density range:   [{rho_final.min():.3f}, {rho_final.max():.3f}]")
    print(f"Density contrast: {rho_final.max() / rho_final.min():.2f}")
    print("=" * 60)
    
    return lbm, rho_final


def main():
    # Grid size (larger = more detail but slower)
    nx, ny = 200, 200
    
    # tau: viscosity nu = (tau - 0.5)/3, range 0.5–2.0
    tau = 1.0
    
    # G must be negative; critical value ~-4, typical range -4.5 to -7.0
    G = -5.5
    
    rho_liquid = 2.0
    rho_gas    = 0.1
    
    center_x = None   # None = grid center
    center_y = None
    radius   = 30.0

    num_steps  = 2000
    save_plots = True
    animate    = True
    view_3d    = True
    output_dir = "output"
    
    lbm, rho_final = run_simulation(
        nx=nx, ny=ny, tau=tau, G=G,
        rho_liquid=rho_liquid, rho_gas=rho_gas,
        center_x=center_x, center_y=center_y, radius=radius,
        num_steps=num_steps, save_plots=save_plots, animate=animate,
        view_3d=view_3d, output_dir=output_dir
    )
    
    print("\nSimulation complete!")


if __name__ == "__main__":
    main()
