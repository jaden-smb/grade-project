#!/usr/bin/env python3
"""
Two-Phase Lattice Boltzmann Simulation using Shan-Chen Model

This script runs a liquid-gas phase separation simulation and visualizes
the density evolution over time. The simulation demonstrates:
- Formation of a liquid droplet with high density
- Formation of surrounding gas with low density
- Stable or evolving liquid-gas interface
- Clear phase separation

Usage:
    python simulate.py

Parameters can be modified in the main() function to control:
- Grid size (nx, ny)
- Relaxation time (tau) - controls viscosity
- Cohesion parameter (G) - controls phase separation strength
- Initial densities (rho_liquid, rho_gas)
- Droplet size and position
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import sys
import os

# Add parent directory to path to import the compiled module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import lbm_shan_chen
except ImportError:
    print("Error: Could not import lbm_shan_chen module.")
    print("Please compile the C++ extension first using the build instructions.")
    sys.exit(1)


def create_density_plot(rho, nx, ny, title="Density Field", vmin=None, vmax=None):
    """
    Create a 2D density plot using matplotlib
    
    Args:
        rho: 1D density array (row-major order) or 2D array
        nx: Grid width
        ny: Grid height
        title: Plot title
        vmin, vmax: Color scale limits (auto if None)
    """
    # Reshape to 2D if needed
    if isinstance(rho, list) or (isinstance(rho, np.ndarray) and rho.ndim == 1):
        rho_2d = np.array(rho).reshape((ny, nx))
    else:
        rho_2d = rho
    
    # Auto-determine color scale if not provided
    if vmin is None:
        vmin = rho_2d.min()
    if vmax is None:
        vmax = rho_2d.max()
    
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(rho_2d, cmap='viridis', origin='lower', 
                   interpolation='bilinear', vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('X (lattice units)', fontsize=12)
    ax.set_ylabel('Y (lattice units)', fontsize=12)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Density', fontsize=12)
    
    plt.tight_layout()
    return fig, ax, im


def create_3d_plot(rho, nx, ny, title="3D Density Field", vmin=None, vmax=None):
    """
    Create a 3D surface plot of the density field
    
    Args:
        rho: 1D density array (row-major order) or 2D array
        nx: Grid width
        ny: Grid height
        title: Plot title
        vmin, vmax: Color scale limits (auto if None)
    """
    # Reshape to 2D if needed
    if isinstance(rho, list) or (isinstance(rho, np.ndarray) and rho.ndim == 1):
        rho_2d = np.array(rho).reshape((ny, nx))
    else:
        rho_2d = rho
    
    # Auto-determine color scale if not provided
    if vmin is None:
        vmin = rho_2d.min()
    if vmax is None:
        vmax = rho_2d.max()
    
    # Create coordinate grids
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)
    
    # Create 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create surface plot
    surf = ax.plot_surface(X, Y, rho_2d, cmap='viridis', 
                          linewidth=0, antialiased=True,
                          vmin=vmin, vmax=vmax, alpha=0.9)
    
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('X (lattice units)', fontsize=12)
    ax.set_ylabel('Y (lattice units)', fontsize=12)
    ax.set_zlabel('Density', fontsize=12)
    
    # Add colorbar
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=20, label='Density')
    
    # Set viewing angle
    ax.view_init(elev=45, azim=45)
    
    plt.tight_layout()
    return fig, ax, surf


def run_simulation(nx=200, ny=200, tau=1.0, G=-120.0, 
                   rho_liquid=2.0, rho_gas=0.1,
                   center_x=None, center_y=None, radius=30.0,
                   num_steps=2000, save_plots=True, animate=True,
                   view_3d=True, output_dir="output"):
    """
    Run the LBM simulation and visualize results
    
    Parameters:
    -----------
    nx, ny : int
        Grid dimensions
    tau : float
        Relaxation time (controls viscosity, typically 0.5-2.0)
        Lower tau = higher viscosity
    G : float
        Cohesion parameter (must be negative for phase separation)
        More negative = stronger phase separation
    rho_liquid : float
        Initial liquid density
    rho_gas : float
        Initial gas density
    center_x, center_y : int
        Droplet center position (default: center of grid)
    radius : float
        Initial droplet radius in lattice units
    num_steps : int
        Number of simulation steps
    save_plots : bool
        Whether to save density plots at intervals
    animate : bool
        Whether to create an animated visualization
    view_3d : bool
        Whether to create 3D surface plots
    output_dir : str
        Directory to save output files
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}/")
    
    # Default center to grid center
    if center_x is None:
        center_x = nx // 2
    if center_y is None:
        center_y = ny // 2
    
    print("=" * 60)
    print("Two-Phase LBM Simulation - Shan-Chen Model")
    print("=" * 60)
    print(f"Grid size: {nx} x {ny}")
    print(f"Relaxation time (tau): {tau}")
    print(f"Cohesion parameter (G): {G}")
    print(f"Liquid density: {rho_liquid}")
    print(f"Gas density: {rho_gas}")
    print(f"Droplet center: ({center_x}, {center_y})")
    print(f"Droplet radius: {radius}")
    print(f"Number of steps: {num_steps}")
    print("=" * 60)
    
    # Create LBM simulator
    print("\nInitializing LBM simulator...")
    lbm = lbm_shan_chen.LBMShanChen(nx, ny, tau, G, rho_liquid, rho_gas)
    
    # Initialize droplet
    print("Initializing liquid droplet...")
    lbm.initialize_droplet(center_x, center_y, radius)
    
    # Get initial state
    rho_initial = np.array(lbm.get_density()).reshape((ny, nx))
    
    # Track density range for consistent visualization
    density_min = rho_initial.min()
    density_max = rho_initial.max()
    
    # Plot initial state
    print("\nPlotting initial state...")
    fig_init, ax_init, im_init = create_density_plot(
        rho_initial, nx, ny, 
        f"Initial State (t=0)\nDroplet at ({center_x}, {center_y})"
    )
    if save_plots:
        plt.savefig(f'{output_dir}/density_initial.png', dpi=150, bbox_inches='tight')
        print(f"Saved: {output_dir}/density_initial.png")
        if view_3d:
            fig_3d_init, ax_3d_init, surf_3d_init = create_3d_plot(
                rho_initial, nx, ny,
                f"Initial State 3D (t=0)"
            )
            plt.savefig(f'{output_dir}/density_initial_3d.png', dpi=150, bbox_inches='tight')
            plt.close(fig_3d_init)
            print(f"Saved: {output_dir}/density_initial_3d.png")
    plt.close(fig_init)
    
    # Storage for animation
    if animate:
        frames = []
        frames.append(rho_initial.copy())
        frames_3d = [] if view_3d else None
        if view_3d:
            frames_3d.append(rho_initial.copy())
    
    # Run simulation
    print(f"\nRunning simulation for {num_steps} steps...")
    print("Progress: ", end="", flush=True)
    
    # Save more frames for better visibility - save every N steps
    save_interval = max(1, num_steps // 20)  # Save ~20 plots
    frame_interval = max(1, num_steps // 200)  # Store ~200 frames for animation
    
    for step in range(1, num_steps + 1):
        # Perform one LBM step
        lbm.step()
        
        # Progress indicator
        if step % (num_steps // 20) == 0:
            print(".", end="", flush=True)
        
        # Get current density
        rho = np.array(lbm.get_density()).reshape((ny, nx))
        
        # Update density range tracking
        density_min = min(density_min, rho.min())
        density_max = max(density_max, rho.max())
        
        # Save intermediate plots
        if save_plots and step % save_interval == 0:
            # Use consistent color scale based on all densities seen so far
            fig, ax, im = create_density_plot(
                rho, nx, ny,
                f"Time Step {step}",
                vmin=density_min, vmax=density_max
            )
            plt.savefig(f'{output_dir}/density_step_{step:05d}.png', dpi=150, bbox_inches='tight')
            plt.close(fig)
            
            if view_3d:
                fig_3d, ax_3d, surf_3d = create_3d_plot(
                    rho, nx, ny,
                    f"Time Step {step} (3D)",
                    vmin=density_min, vmax=density_max
                )
                plt.savefig(f'{output_dir}/density_step_{step:05d}_3d.png', dpi=150, bbox_inches='tight')
                plt.close(fig_3d)
        
        # Store frame for animation (more frequent for smoother animation)
        if animate and step % frame_interval == 0:
            frames.append(rho.copy())
            if view_3d and frames_3d is not None:
                frames_3d.append(rho.copy())
    
    print(" Done!")
    
    # Get final state
    rho_final = np.array(lbm.get_density()).reshape((ny, nx))
    
    # Plot final state
    print("\nPlotting final state...")
    fig_final, ax_final, im_final = create_density_plot(
        rho_final, nx, ny,
        f"Final State (t={num_steps})",
        vmin=density_min, vmax=density_max
    )
    if save_plots:
        plt.savefig(f'{output_dir}/density_final.png', dpi=150, bbox_inches='tight')
        print(f"Saved: {output_dir}/density_final.png")
        if view_3d:
            fig_3d_final, ax_3d_final, surf_3d_final = create_3d_plot(
                rho_final, nx, ny,
                f"Final State 3D (t={num_steps})",
                vmin=density_min, vmax=density_max
            )
            plt.savefig(f'{output_dir}/density_final_3d.png', dpi=150, bbox_inches='tight')
            plt.close(fig_3d_final)
            print(f"Saved: {output_dir}/density_final_3d.png")
    plt.show(block=False)
    
    # Create animation if requested
    if animate and len(frames) > 1:
        print("\nCreating 2D animation...")
        fig_anim, ax_anim = plt.subplots(figsize=(10, 10))
        # Use consistent color scale for all frames
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
            step_num = frame_num * frame_interval
            ax_anim.set_title(f'Density Evolution - Step {step_num} (Frame {frame_num}/{len(frames)-1})',
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
        
        # Create 3D animation if requested
        if view_3d and frames_3d is not None and len(frames_3d) > 1:
            print("\nCreating 3D animation...")
            x = np.arange(nx)
            y = np.arange(ny)
            X, Y = np.meshgrid(x, y)
            
            fig_3d_anim = plt.figure(figsize=(12, 10))
            ax_3d_anim = fig_3d_anim.add_subplot(111, projection='3d')
            
            # Initial surface
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
                step_num = frame_num * frame_interval
                ax_3d_anim.set_title(f'3D Density Evolution - Step {step_num} (Frame {frame_num}/{len(frames_3d)-1})',
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
    
    # Print statistics
    print("\n" + "=" * 60)
    print("Simulation Statistics:")
    print("=" * 60)
    print(f"Initial density range: [{rho_initial.min():.3f}, {rho_initial.max():.3f}]")
    print(f"Final density range: [{rho_final.min():.3f}, {rho_final.max():.3f}]")
    print(f"Final max density: {rho_final.max():.3f}")
    print(f"Final min density: {rho_final.min():.3f}")
    print(f"Density contrast: {rho_final.max() / rho_final.min():.2f}")
    print("=" * 60)
    
    return lbm, rho_final


def main():
    """
    Main function - modify parameters here to control simulation behavior
    """
    
    # ========== SIMULATION PARAMETERS ==========
    
    # Grid size (larger = more detail but slower)
    nx = 200
    ny = 200
    
    # Relaxation time (tau)
    # - Controls viscosity: nu = (tau - 0.5) / 3
    # - Typical range: 0.5 < tau < 2.0
    # - Lower tau = higher viscosity = more stable but slower dynamics
    # - Higher tau = lower viscosity = faster dynamics but may be unstable
    tau = 1.0
    
    # Cohesion parameter (G)
    # - MUST be negative for liquid-gas phase separation
    # - More negative = stronger phase separation
    # - Typical range: -200 < G < -50
    # - Too negative = very sharp interface, may be unstable
    # - Too close to zero = weak phase separation
    G = -120.0
    
    # Initial densities
    rho_liquid = 2.0   # Liquid phase density
    rho_gas = 0.1      # Gas phase density
    
    # Droplet initialization
    center_x = None    # None = center of grid
    center_y = None    # None = center of grid
    radius = 30.0      # Initial droplet radius
    
    # Simulation control
    num_steps = 2000   # Number of time steps (increased for longer simulation)
    save_plots = True  # Save density plots
    animate = True     # Create animation
    view_3d = True     # Create 3D visualizations
    output_dir = "output"  # Directory to save output files
    
    # ===========================================
    
    # Run simulation
    lbm, rho_final = run_simulation(
        nx=nx, ny=ny, tau=tau, G=G,
        rho_liquid=rho_liquid, rho_gas=rho_gas,
        center_x=center_x, center_y=center_y, radius=radius,
        num_steps=num_steps, save_plots=save_plots, animate=animate,
        view_3d=view_3d, output_dir=output_dir
    )
    
    print("\nSimulation complete!")
    print("\nTo modify parameters, edit the main() function in simulate.py")


if __name__ == "__main__":
    main()

