#!/usr/bin/env python3
"""
Export 2D LBM density fields as Wavefront OBJ heightmap meshes.

Supports two color modes:
  'texture'  — UV-mapped PNG texture (Option A, default, best Blender/Maya compatibility)
  'vertex'   — per-vertex color via unofficial "v x y z r g b" extension (Option B)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mimage


def export_density_to_obj(rho_2d, filepath, z_scale=50.0,
                          subsample=1, colormap='viridis',
                          rho_min=None, rho_max=None,
                          color_mode='texture'):
    """
    Export a 2D density field as a heightmap OBJ mesh.

    Parameters:
        rho_2d      : np.ndarray of shape (ny, nx) — the density field
        filepath    : str — output path WITHOUT extension (e.g. "output/frame_00001")
                      Creates filepath.obj, filepath.mtl, filepath_texture.png
        z_scale     : float — vertical exaggeration factor (z = rho * z_scale)
        subsample   : int — take every Nth point to reduce mesh size (1 = full res)
        colormap    : str — matplotlib colormap name for vertex/texture colors
        rho_min     : float — lower density bound for color normalization (None = auto)
        rho_max     : float — upper density bound for color normalization (None = auto)
        color_mode  : 'texture' (UV-mapped PNG, Option A) or
                      'vertex'  (per-vertex r g b in OBJ, Option B)
    """
    rho_s = rho_2d[::subsample, ::subsample]
    ny_s, nx_s = rho_s.shape

    rmin = float(rho_s.min()) if rho_min is None else float(rho_min)
    rmax = float(rho_s.max()) if rho_max is None else float(rho_max)
    rng = rmax - rmin if rmax != rmin else 1.0
    rho_norm = np.clip((rho_s - rmin) / rng, 0.0, 1.0)

    cmap = plt.get_cmap(colormap)
    colors = cmap(rho_norm)  # (ny_s, nx_s, 4) RGBA

    x_coords = (np.arange(nx_s) - nx_s / 2.0) * subsample
    y_coords = (np.arange(ny_s) - ny_s / 2.0) * subsample
    xx, yy = np.meshgrid(x_coords, y_coords)
    zz = rho_s * z_scale

    mtl_base = os.path.basename(filepath)
    obj_path = f"{filepath}.obj"
    mtl_path = f"{filepath}.mtl"
    tex_path = f"{filepath}_texture.png"

    parts = []
    parts.append(f"# LBM density heightmap  {nx_s}x{ny_s} vertices  z_scale={z_scale}")
    parts.append(f"# color_mode={color_mode}  colormap={colormap}")
    parts.append("")

    if color_mode == 'texture':
        parts.append(f"mtllib {mtl_base}.mtl")
        parts.append("")

    vx = xx.ravel()
    vy = yy.ravel()
    vz = zz.ravel()

    # Output vertices in Blender's Y-up OBJ convention:
    #   OBJ X = grid column (vx)
    #   OBJ Y = height / density (vz)  — Blender maps OBJ Y → Blender Z (up)
    #   OBJ Z = grid row (vy)          — Blender maps OBJ Z → Blender -Y
    if color_mode == 'vertex':
        vr = colors[:, :, 0].ravel()
        vg = colors[:, :, 1].ravel()
        vb = colors[:, :, 2].ravel()
        vertex_lines = [
            f"v {x:.4f} {z:.4f} {y:.4f} {r:.4f} {g:.4f} {b:.4f}"
            for x, y, z, r, g, b in zip(vx, vy, vz, vr, vg, vb)
        ]
    else:
        vertex_lines = [
            f"v {x:.4f} {z:.4f} {y:.4f}"
            for x, y, z in zip(vx, vy, vz)
        ]

    parts.extend(vertex_lines)
    parts.append("")

    if color_mode == 'texture':
        u_coords = np.arange(nx_s) / max(nx_s - 1, 1)
        v_coords = np.arange(ny_s) / max(ny_s - 1, 1)
        uu, vv = np.meshgrid(u_coords, v_coords)
        uv_flat = np.stack([uu.ravel(), vv.ravel()], axis=1)
        uv_lines = [f"vt {u:.6f} {v:.6f}" for u, v in uv_flat]
        parts.extend(uv_lines)
        parts.append("")
        parts.append("usemtl density_material")
        parts.append("")

    JJ, II = np.mgrid[0:ny_s - 1, 0:nx_s - 1]
    v1 = (JJ * nx_s + II + 1).ravel()
    v2 = (JJ * nx_s + II + 2).ravel()
    v3 = ((JJ + 1) * nx_s + II + 2).ravel()
    v4 = ((JJ + 1) * nx_s + II + 1).ravel()

    if color_mode == 'texture':
        face_lines = [
            f"f {a}/{a} {b}/{b} {c}/{c} {d}/{d}"
            for a, b, c, d in zip(v1, v2, v3, v4)
        ]
    else:
        face_lines = [
            f"f {a} {b} {c} {d}"
            for a, b, c, d in zip(v1, v2, v3, v4)
        ]

    parts.extend(face_lines)

    with open(obj_path, 'w') as f:
        f.write('\n'.join(parts))
        f.write('\n')

    if color_mode == 'texture':
        tex_img = colors[::-1, :, :3]
        mimage.imsave(tex_path, tex_img)

        mtl_lines = [
            "newmtl density_material",
            "Ka 1.0 1.0 1.0",
            "Kd 1.0 1.0 1.0",
            "Ks 0.0 0.0 0.0",
            "Ns 0.0",
            f"map_Kd {mtl_base}_texture.png",
        ]
        with open(mtl_path, 'w') as f:
            f.write('\n'.join(mtl_lines))
            f.write('\n')


def export_obj_sequence(frames, output_dir, z_scale=50.0,
                        subsample=1, colormap='viridis',
                        rho_min=None, rho_max=None,
                        color_mode='texture'):
    """
    Export a list of density frames as a numbered OBJ sequence.

    Parameters:
        frames      : list of np.ndarray — density fields from the simulation
        output_dir  : str — directory for the OBJ files
        z_scale, subsample, colormap, color_mode : same as export_density_to_obj
        rho_min, rho_max : density range for consistent color scaling across all
                           frames (None = computed globally from all frames)

    Creates per-frame:  frame_NNNNN.obj  frame_NNNNN.mtl  frame_NNNNN_texture.png
    """
    os.makedirs(output_dir, exist_ok=True)

    n = len(frames)
    if n == 0:
        print("export_obj_sequence: no frames to export.")
        return

    if rho_min is None or rho_max is None:
        all_min = min(float(f.min()) for f in frames)
        all_max = max(float(f.max()) for f in frames)
        if rho_min is None:
            rho_min = all_min
        if rho_max is None:
            rho_max = all_max

    digits = max(5, len(str(n - 1)))

    for idx, frame in enumerate(frames):
        if (idx + 1) % max(1, n // 20) == 0 or idx == 0 or idx == n - 1:
            print(f"  Exporting OBJ {idx + 1}/{n}...", flush=True)
        name = f"frame_{idx:0{digits}d}"
        filepath = os.path.join(output_dir, name)
        export_density_to_obj(
            frame, filepath,
            z_scale=z_scale,
            subsample=subsample,
            colormap=colormap,
            rho_min=rho_min,
            rho_max=rho_max,
            color_mode=color_mode,
        )

    print(f"  OBJ sequence complete: {n} frames written to '{output_dir}/'")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Export LBM density field(s) to Wavefront OBJ heightmap mesh(es)."
    )
    parser.add_argument(
        "input",
        help="Path to a single .npy density file, or a directory containing .npy files.",
    )
    parser.add_argument("-o", "--output", default="obj_output",
                        help="Output directory (default: obj_output)")
    parser.add_argument("--z-scale", type=float, default=50.0,
                        help="Vertical exaggeration factor (default: 50.0)")
    parser.add_argument("--subsample", type=int, default=1,
                        help="Take every Nth grid point (default: 1 = full resolution)")
    parser.add_argument("--colormap", default="viridis",
                        help="Matplotlib colormap name (default: viridis)")
    parser.add_argument("--color-mode", default="texture", choices=["texture", "vertex"],
                        help="Color encoding: 'texture' (UV PNG, default) or 'vertex'")
    args = parser.parse_args()

    input_path = args.input

    if os.path.isfile(input_path) and input_path.endswith(".npy"):
        rho = np.load(input_path)
        if rho.ndim == 1:
            side = int(np.sqrt(rho.size))
            rho = rho.reshape((side, side))
        os.makedirs(args.output, exist_ok=True)
        base = os.path.splitext(os.path.basename(input_path))[0]
        filepath = os.path.join(args.output, base)
        print(f"Exporting {input_path} → {filepath}.obj ...")
        export_density_to_obj(rho, filepath, z_scale=args.z_scale,
                              subsample=args.subsample, colormap=args.colormap,
                              color_mode=args.color_mode)
        print("Done.")

    elif os.path.isdir(input_path):
        npy_files = sorted(
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if f.endswith(".npy")
        )
        if not npy_files:
            print(f"No .npy files found in '{input_path}'.")
            raise SystemExit(1)
        print(f"Found {len(npy_files)} .npy files in '{input_path}'.")
        frames = []
        for p in npy_files:
            rho = np.load(p)
            if rho.ndim == 1:
                side = int(np.sqrt(rho.size))
                rho = rho.reshape((side, side))
            frames.append(rho)
        export_obj_sequence(frames, args.output, z_scale=args.z_scale,
                            subsample=args.subsample, colormap=args.colormap,
                            color_mode=args.color_mode)
    else:
        print(f"Input '{input_path}' is not a .npy file or a directory.")
        raise SystemExit(1)
