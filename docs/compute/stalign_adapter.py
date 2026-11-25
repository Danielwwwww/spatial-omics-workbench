"""Optional native STalign backend for the registration model comparison."""
from dataclasses import dataclass

import numpy as np

from .registration import descriptor, warp


@dataclass
class STalignField:
    package: object
    state: dict

    def _transform(self, points, transform):
        batches = []
        args = [self.state[key] for key in ('xv', 'v', 'A')]
        for start in range(0, len(points), 32768):
            yx = np.ascontiguousarray(points[start:start + 32768, ::-1])
            batches.append(transform(*args, yx).detach().cpu().numpy()[:, ::-1])
        return np.concatenate(batches) if batches else points.copy()

    def forward(self, points):
        return self._transform(points, self.package.transform_points_source_to_target)

    def inverse(self, points):
        return self._transform(points, self.package.transform_points_target_to_source)

    def to_json(self):
        return {'method': 'STalign LDDMM', 'initialization': 'workbench linear model',
                'affine_in_yx': self.state['A'].detach().cpu().numpy().tolist()}


def stalign_candidate(image, mask, fixed, fixed_mask, matrix, source, target):
    try:
        import STalign
    except ImportError as exc:
        raise ValueError('STalign 需要在本地 Python 环境安装，请参阅 README。') from exc

    step = max(2, int(np.ceil(max(fixed.shape) / 160)))
    radius = max(3, round(max(fixed.shape) / 110))
    moving = warp(image, matrix, fixed.shape)
    moving_mask = warp(mask, matrix, fixed.shape, order=0)
    axes = [np.arange(0, size, step, dtype=float) for size in fixed.shape]
    moving_features = (descriptor(moving, radius)[:4] * moving_mask)[..., ::step, ::step].copy()
    fixed_features = (descriptor(fixed, radius)[:4] * fixed_mask)[..., ::step, ::step].copy()
    kwargs = {}
    if source is not None:
        from .registration import apply
        kwargs['pointsI'] = np.ascontiguousarray(apply(matrix, source)[:, ::-1])
        kwargs['pointsJ'] = np.ascontiguousarray(target[:, ::-1])
    # Freeze the initial affine; assess the additional deformation independently.
    result = STalign.LDDMM(axes, moving_features, axes, fixed_features,
                           A=np.eye(3), prefit=False, epL=0., epT=0.,
                           a=max(3 * step, .07 * max(fixed.shape)), sigma_disp=5.,
                           sigmaP=4., niter=180, draw=False, num_threads=2, **kwargs)
    field = STalignField(STalign, result)
    yy, xx = np.mgrid[:fixed.shape[0]:step, :fixed.shape[1]:step]
    probes = np.column_stack((xx.ravel(), yy.ravel()))
    mapped = field.inverse(probes)
    if not np.isfinite(mapped).all():
        return None, {'reason': 'nonfinite_transform'}
    gy, gx = np.gradient(mapped.reshape((*xx.shape, 2)), step, axis=(0, 1))
    jacobian = gx[..., 0] * gy[..., 1] - gx[..., 1] * gy[..., 0]
    support = fixed_mask[::step, ::step]
    roundtrip = np.linalg.norm(field.forward(mapped) - probes, axis=1)
    metrics = {'jacobian_min': float(jacobian[support].min()),
               'roundtrip_p95_px': float(np.quantile(roundtrip, .95)),
               'max_displacement_px': float(np.linalg.norm(mapped - probes, axis=1).max())}
    if not np.isfinite(list(metrics.values())).all():
        return None, {'reason': 'nonfinite_transform'}
    return field, metrics
