"""Compact display of registered ST counts at H&E preview-pixel resolution."""

import base64

import numpy as np

LOG_CP10K_MAX = np.log1p(10000)


class STDisplay:
    """Collect measured bins on display pixels without filling unsampled positions."""

    def __init__(self, xy, region, umi, width, height):
        coordinates = np.rint(xy).astype(int)
        inside = ((coordinates >= 0) & (coordinates < [width, height])).all(axis=1)
        self.valid = inside & (region > 0) & (umi > 0)
        self.area = width * height
        x, y = coordinates[self.valid].T
        key = region[self.valid].astype(np.int64) * self.area + y * width + x
        labels, groups = np.unique(key, return_inverse=True)
        self.labels = labels
        self.groups = np.full(len(xy), -1, dtype=np.int64)
        self.groups[self.valid] = groups
        self.umi = np.bincount(groups, weights=umi[self.valid])

    def gene(self, rows, counts):
        groups = self.groups[rows]
        use = groups >= 0
        total = np.bincount(groups[use], weights=counts[use], minlength=len(self.labels))
        selected = total > 0
        if not selected.any():
            return {"points": "", "p98": 0, "n": 0, "scale": "fixed-log1p-cp10k"}
        values = np.log1p(10000 * total[selected] / self.umi[selected])
        scale = float(np.quantile(values, .98))
        packed = np.empty(len(values), dtype=[("pixel", "<u4"), ("level", "u1"), ("region", "u1")])
        packed["pixel"] = self.labels[selected] % self.area
        packed["level"] = np.maximum(1, np.rint(255 * np.clip(values / LOG_CP10K_MAX, 0, 1))).astype(np.uint8)
        packed["region"] = self.labels[selected] // self.area
        return {"points": base64.b64encode(packed.tobytes()).decode("ascii"),
                "p98": round(scale, 5), "n": len(values), "scale": "fixed-log1p-cp10k"}
