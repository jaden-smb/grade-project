import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit


def create_density_plot(rho, nx, ny, title="Density Field", vmin=None, vmax=None):

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


def plot_density_profile(rho_2d, nx, ny, output_dir):

    profile = rho_2d[ny // 2, :]
    x = np.arange(nx)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x, profile, 'b-', linewidth=1.5, label='ρ(x)')

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

    rho_2d = np.array(lbm.get_density()).reshape((ny, nx))
    ux_2d  = np.array(lbm.get_velocity_x()).reshape((ny, nx))
    uy_2d  = np.array(lbm.get_velocity_y()).reshape((ny, nx))

    speed = np.sqrt(ux_2d**2 + uy_2d**2)
    max_speed = float(speed.max())
    print(f"  Max spurious velocity magnitude: {max_speed:.6f} lu/step")

    step = 8
    xs = np.arange(0, nx, step)
    ys = np.arange(0, ny, step)
    X, Y = np.meshgrid(xs, ys)
    Ux = ux_2d[::step, ::step]
    Uy = uy_2d[::step, ::step]
    Sp = speed[::step, ::step]

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
