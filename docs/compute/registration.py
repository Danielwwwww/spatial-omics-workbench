"""Orientation search, landmark registration, masked image refinement, and QA."""
from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import least_squares, minimize
from scipy.spatial.distance import cdist


def xy_where(mask):
    return np.column_stack(np.nonzero(mask)[::-1]).astype(float)


def sample(image, xy, order=1, mode="grid-constant"):
    return ndi.map_coordinates(image, xy.T[::-1], order=order, mode=mode)


def apply(matrix, xy):
    return xy @ matrix[:2, :2].T + matrix[:2, 2]


def affine(linear, shift):
    result = np.eye(3)
    result[:2] = np.column_stack([linear, shift])
    return result


def largest(mask):
    labels, n = ndi.label(mask)
    if n == 0:
        return mask.copy()
    counts = np.bincount(labels.ravel())
    return labels == (1 + np.argmax(counts[1:]))


def otsu(image):
    values = image[np.isfinite(image)]
    if np.ptp(values) == 0:
        return float(values[0])
    hist, bins = np.histogram(values, bins=128)
    centers = (bins[1:] + bins[:-1]) / 2
    count = hist.cumsum()[:-1]
    moments = (hist * centers).cumsum()
    score = (moments[:-1] * hist.sum() - count * moments[-1]) ** 2
    score /= np.maximum(count * (hist.sum() - count), 1)
    return float(bins[np.argmax(score) + 1])


def signed(mask):
    return ndi.distance_transform_edt(~mask) - ndi.distance_transform_edt(mask)


def outline_fit(source, target):
    def boundary(mask):
        points = xy_where(mask & ~ndi.binary_erosion(mask))
        return points[::max(1, len(points) // 1600)]

    qs, qt = boundary(source), boundary(target)
    if min(len(qs), len(qt)) < 8:
        raise ValueError("Tissue or capture mask is too small; supply masks or landmarks")
    cs, ct = xy_where(source).mean(0), xy_where(target).mean(0)
    ds, dt = signed(source), signed(target)
    scale = np.sqrt(target.sum() / source.sum())
    params = np.r_[(np.eye(2) * scale).ravel(), 0, 0]
    margin = max(target.shape) * 0.2
    lo = np.r_[scale * np.array([.65, -.38, -.38, .65]), -margin, -margin]
    hi = np.r_[scale * np.array([1.5, .38, .38, 1.5]), margin, margin]
    rows = []
    for swap, a, b in product([np.eye(2), np.eye(2)[::-1]], [1, -1], [1, -1]):
        orientation = np.diag([a, b]) @ swap

        def residual(p):
            linear = p[:4].reshape(2, 2) @ orientation
            offset = ct + p[4:] - linear @ cs
            forward = qs @ linear.T + offset
            backward = (qt - offset) @ np.linalg.inv(linear).T
            return np.r_[sample(dt, forward, mode="nearest"),
                         sample(ds, backward, mode="nearest") * scale * np.sqrt(len(qs) / len(qt))]

        fit = least_squares(residual, params, bounds=(lo, hi), loss="soft_l1", f_scale=4, max_nfev=90)
        linear = fit.x[:4].reshape(2, 2) @ orientation
        matrix = affine(linear, ct + fit.x[4:] - linear @ cs)
        moved = warp(source, matrix, target.shape, order=0)
        dice = 2 * np.count_nonzero(moved & target) / (moved.sum() + target.sum())
        rows.append((float(dice), matrix))
    rows.sort(key=lambda v: v[0], reverse=True)
    return rows[0][1], [v[0] for v in rows]


def landmark_fit(src, dst, weights, model="affine"):
    w = np.sqrt(weights)[:, None]
    if model == "affine":
        coef, _, rank, _ = np.linalg.lstsq(np.c_[src, np.ones(len(src))] * w, dst * w, rcond=None)
        if rank < 3 or abs(np.linalg.det(coef[:2])) < 1e-12:
            raise ValueError("Affine landmarks must be noncollinear in both views")
        return affine(coef[:2].T, coef[2])
    ca, cb = np.average(src, axis=0, weights=weights), np.average(dst, axis=0, weights=weights)
    xa, xb = (src - ca) * w, (dst - cb) * w
    u, sv, vt = np.linalg.svd(xb.T @ xa)
    linear = (u @ vt) * (sv.sum() / np.sum(xa**2))
    return affine(linear, cb - linear @ ca)


@dataclass
class Field:
    centers: np.ndarray
    width: float
    coeff: np.ndarray

    def inverse(self, points):
        return points + self.basis(points) @ self.coeff

    def forward(self, points):
        result = points.copy()
        for start in range(0, len(points), 20_000):
            batch = slice(start, start + 20_000)
            estimate = result[batch]
            target = points[batch]
            for _ in range(32):
                step = target - self.inverse(estimate)
                estimate += step
                if np.max(np.abs(step)) < 1e-4:
                    break
        return result

    def basis(self, points):
        return np.exp(-cdist(points, self.centers, 'sqeuclidean') / (2 * self.width**2))

    def to_json(self):
        return {"centers_xy": self.centers.tolist(), "width_px": self.width, "inverse_coeff_dxdy": self.coeff.tolist()}


def warp(source, matrix, shape, field=None, order=1):
    points = xy_where(np.ones(shape, dtype=bool))
    if field is not None:
        points = field.inverse(points)
    return sample(source, apply(np.linalg.inv(matrix), points), order).reshape(shape)


def descriptor(image, radius):
    # Local patch self-similarity, inspired by MIND. Values are contrast-insensitive.
    image = ndi.gaussian_filter(image.astype(float), max(1.1, radius / 3))
    shifts = [(0, radius), (0, -radius), (radius, 0), (-radius, 0),
              (radius, radius), (radius, -radius), (-radius, radius), (-radius, -radius)]
    diffs = np.stack([ndi.gaussian_filter((image - ndi.shift(image, s, order=0, mode="nearest"))**2,
                                         radius / 2) for s in shifts], axis=0)
    local = diffs[:4].mean(axis=0)
    floor = max(np.quantile(local, .55) * .15, 1e-12)
    return np.exp(-(diffs - diffs.min(0)) / np.maximum(local, floor))


def image_refine(source, fixed, source_mask, fixed_mask, matrix, train_targets=None, source_at_targets=None):
    """Fit regularized inverse displacement on shared interior, withholding spatial image blocks."""
    moving = warp(source, matrix, fixed.shape)
    coverage = warp(source_mask, matrix, fixed.shape, order=0)
    radius = max(3, round(max(fixed.shape) / 110))
    common = ndi.binary_erosion(coverage & fixed_mask, iterations=3 * radius)
    step = max(2, round(max(fixed.shape) / 150))
    clipped = np.zeros_like(common)
    clipped[::step, ::step] = common[::step, ::step]
    xy = xy_where(clipped)
    if len(xy) < 150:
        return None, {"reason": "insufficient_shared_interior"}
    n = max(fixed.shape)
    fg, mg = descriptor(fixed, radius), descriptor(moving, radius)
    signal = ndi.gaussian_gradient_magnitude(fixed.astype(float), sigma=radius / 2)
    weight = np.clip(sample(signal, xy) / max(float(np.quantile(signal[common], .8)), 1e-8), .15, 2.0)
    block = (xy // max(16, int(n / 17))).astype(int)
    held = ((block[:, 0] + 2 * block[:, 1]) % 5) == 0
    if held.sum() < 30 or (~held).sum() < 100:
        return None, {"reason": "insufficient_held_out_image_blocks"}
    span = np.ptp(xy, axis=0)
    grid = np.meshgrid(*[np.linspace(a, b, max(2, min(5, int(s / (n * .19)) + 1)))
                         for a, b, s in zip(xy.min(0), xy.max(0), span)])
    field = Field(np.column_stack([v.ravel() for v in grid]), n * .125, np.empty((grid[0].size, 2)))
    basis = field.basis(xy)
    targets = np.stack([sample(channel, xy) for channel in fg], axis=1)
    deriv = [np.gradient(ch) for ch in mg]
    train = ~held
    anchor_basis = field.basis(train_targets) if train_targets is not None else None
    max_coeff = n * .045
    displacement_scale = n * .022

    def loss(coef, subset, gradient):
        beta = coef.reshape(-1, 2)
        moved = xy[subset] + basis[subset] @ beta
        residual = np.stack([sample(ch, moved) for ch in mg], axis=1) - targets[subset]
        scale = .18
        z = np.sqrt(1 + (residual / scale)**2)
        w = weight[subset, None] / (weight[subset].sum() * mg.shape[0])
        value = np.sum(scale**2 * (z - 1) * w)
        if not gradient:
            return float(value)
        score = residual / z * w
        gd = np.zeros_like(moved)
        for k, (dy, dx) in enumerate(deriv):
            gd[:, 0] += score[:, k] * sample(dx, moved)
            gd[:, 1] += score[:, k] * sample(dy, moved)
        penalty = .0011 / displacement_scale**2
        value += penalty * np.mean(beta**2)
        grad = basis[subset].T @ gd + (2 * penalty / beta.size) * beta
        if anchor_basis is not None:
            anch = train_targets + anchor_basis @ beta - source_at_targets
            anchor_weight = .001 / displacement_scale**2
            value += anchor_weight * np.mean(anch**2)
            grad += 2 * anchor_weight / anch.size * anchor_basis.T @ anch
        return value, grad.ravel()

    x0 = np.zeros(field.coeff.size)
    fit = minimize(loss, x0, args=(train, True), jac=True, method="L-BFGS-B",
                   bounds=[(-max_coeff, max_coeff)] * len(x0), options={"maxiter": 100, "ftol": 1e-10})
    field.coeff = fit.x.reshape(-1, 2)
    # Global Lipschitz upper bound from the summed derivatives of Gaussian basis functions.
    # For multiple centers, use a center-independent bound; shrink before using the inverse.
    bound = np.linalg.norm(field.coeff, axis=1).sum() / (np.sqrt(np.e) * field.width)
    if bound > .66:
        field.coeff *= .66 / bound
        bound = .66
    train_score = loss(field.coeff.ravel(), train, False)
    base_score = loss(x0, held, False)
    validation_score = loss(field.coeff.ravel(), held, False)
    roundtrip = np.linalg.norm(field.forward(field.inverse(xy[::9])) - xy[::9], axis=1)
    info = {"held_image_patch_error_base": base_score, "held_image_patch_error_refined": validation_score,
            "fit_patch_error": train_score, "held_image_samples": int(held.sum()),
            "displacement_lipschitz_bound": float(bound), "roundtrip_p95_px": float(np.quantile(roundtrip, .95)),
            "max_displacement_px": float(np.linalg.norm(field.inverse(xy) - xy, axis=1).max())}
    return field, info


def register(image, mask, he, he_mask, landmarks, refine=True, allow_refined=True,
             evaluation_scale=(1, 1), use_stalign=False, tre_priority='balanced'):
    rows = [r for r in landmarks if r["role"] == "fit"]
    check = [r for r in landmarks if r["role"] == "check"]
    info = {"fit_pairs": len(rows), "check_pairs": len(check), "models": {}}
    if tre_priority not in ('balanced', 'max'):
        raise ValueError('tre_priority must be balanced or max')
    info["tre_priority"] = tre_priority
    src = dst = None
    start = "affine"
    linear_candidates = {}
    if len(rows) >= 3:
        src = np.array([r["source"] for r in rows])
        dst = np.array([r["target"] for r in rows])
        weight = np.array([r.get("weight", 1) for r in rows])
        linear_candidates = {name:landmark_fit(src, dst, weight, name) for name in ("affine", "similarity")}
        groups = {r.get("group", r.get("id", str(i))) for i, r in enumerate(rows)}
        group_errors = {"affine": [], "similarity": []}
        for group in sorted(groups):
            train = np.array([r.get("group",r.get("id",str(i))) != group for i,r in enumerate(rows)])
            if train.sum() < 3:
                continue
            for model in group_errors:
                try:
                    lm = landmark_fit(src[train], dst[train], weight[train], model)
                except ValueError:
                    continue
                err = np.linalg.norm((apply(lm, src[~train]) - dst[~train]) * evaluation_scale, axis=1)
                group_errors[model].append(float(np.mean(err)))
        info["loo_group_mean_errors"] = group_errors
        errors_aff, errors_sim = group_errors["affine"], group_errors["similarity"]
        if len(errors_aff) >= 4 and len(errors_aff) == len(errors_sim) and np.mean(errors_sim) < .95*np.mean(errors_aff):
            start = "similarity"
        matrix = linear_candidates[start]
        info["initialization"] = f"{start}_landmarks"
    else:
        if rows:
            raise ValueError("Supply at least 3 fit pairs, or leave fit pairs empty for automatic contour initialization")
        matrix, scores = outline_fit(mask, he_mask)
        info["orientation_dice_candidates"] = scores
        info["orientation_ambiguity"] = scores[0] - scores[1] < .02
        info["initialization"] = "contour_orientation_search"
        linear_candidates["affine"] = matrix
    field, refqa = image_refine(image, he, mask, he_mask, matrix, dst,
                                apply(matrix, src) if src is not None else None) if refine else (None, {})
    info["image_refinement"] = refqa
    info["refinement_base_model"] = start
    chosen = start
    candidates = [(name, mat, None) for name, mat in linear_candidates.items()]
    if field is not None:
        candidates.append(("image_refined", matrix, field))
    stalign_field = None
    if use_stalign:
        from .stalign_adapter import stalign_candidate

        stalign_field, metrics = stalign_candidate(image, mask, he, he_mask, matrix, src, dst)
        info["stalign_refinement"] = metrics
        if stalign_field is not None:
            candidates.append(("stalign", matrix, stalign_field))
    for label, mat, f in candidates:
        mapped_mask = warp(mask, mat, he_mask.shape, f, order=0)
        intersection = int(np.count_nonzero(mapped_mask & he_mask))
        measurements = []
        if check:
            tc = np.array([r["target"] for r in check])
            estimate = apply(mat, np.array([r["source"] for r in check]))
            if f is not None:
                estimate = f.forward(estimate)
            errors = np.linalg.norm((estimate - tc) * evaluation_scale, axis=1)
            diagonal = np.linalg.norm(np.asarray(he_mask.shape[::-1]) * evaluation_scale)
            measurements = [{"id": r.get("id", str(i)), "tre": float(e), "target": tc[i].tolist(),
                             "relative_tre": float(e / diagonal), "predicted": estimate[i].tolist()}
                            for i, (r, e) in enumerate(zip(check, errors))]
        else:
            errors = np.array([])
        info["models"][label] = {
            "tre": measurements, "mean_tre": float(np.mean(errors)) if len(errors) else None,
            "median_tre": float(np.median(errors)) if len(errors) else None,
            "max_tre": float(np.max(errors)) if len(errors) else None,
            "dice": float(2 * intersection / max(1, mapped_mask.sum() + he_mask.sum())),
            "target_coverage": float(intersection / he_mask.sum()),
            "source_containment": float(intersection / max(1, mapped_mask.sum())),
            "source_raster_to_he_view_affine": mat.tolist()}
    base = info["models"][start]
    refined = info["models"].get("image_refined")
    if (refined is not None and allow_refined and len(check) >= 3
            and refqa["roundtrip_p95_px"] < .1
            and refqa["held_image_patch_error_refined"] < .98 * refqa["held_image_patch_error_base"]
            and refined["mean_tre"] < base["mean_tre"]
            and refined["median_tre"] < base["median_tre"]
            and refined["max_tre"] < base["max_tre"]
            and refined["dice"] >= base["dice"] - .02):
        chosen = "image_refined"
    stalign = info["models"].get("stalign")
    if stalign is not None and allow_refined and len(check) >= 3:
        metrics = info["stalign_refinement"]
        current = info["models"][chosen]
        if tre_priority == 'max':
            better = stalign["max_tre"] < current["max_tre"]
        else:
            better = stalign["mean_tre"] < current["mean_tre"] and stalign["max_tre"] <= current["max_tre"]
        if (metrics["jacobian_min"] > 0 and metrics["roundtrip_p95_px"] < .1
                and metrics["max_displacement_px"] <= .05 * max(he.shape)
                and better
                and stalign["dice"] >= base["dice"] - .02):
            chosen = "stalign"
    info["selected"] = chosen
    info["quality_status"] = "heldout_landmarks_available" if check else "unverified_geometry"
    fields = {"image_refined": field, "stalign": stalign_field}
    return matrix, fields.get(chosen), info, stalign_field if use_stalign else field
