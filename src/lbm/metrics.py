import numpy as np
from scipy import ndimage


def compute_metrics(rho_2d, rho_threshold):
    total_mass = float(rho_2d.sum())
    max_density = float(rho_2d.max())
    liquid_mask = rho_2d > rho_threshold
    area = int(liquid_mask.sum())
    radius = float(np.sqrt(area / np.pi)) if area > 0 else 0.0

    if area > 20:
        eroded = ndimage.binary_erosion(liquid_mask)
        perimeter = float((liquid_mask & ~eroded).sum())
        if perimeter > 0:
            circularity = 4.0 * np.pi * area / (perimeter * perimeter)
            circularity = min(circularity, 1.0)
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
