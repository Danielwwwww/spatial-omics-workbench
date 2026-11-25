"""Validate hand-drawn ROIs in original H&E pixel coordinates."""

import numpy as np


def validate_regions(regions, size, existing=()):
    if not isinstance(regions, list) or len(regions) > 32:
        raise ValueError("手绘脑区应为列表，最多 32 个。")
    used, output = set(existing), []
    width, height = size
    for region in regions:
        name = region["name"].strip()
        xy = np.asarray(region["points"], dtype=float)
        if not name or len(name) > 64 or name in used:
            raise ValueError("请给手绘脑区填写唯一名称，最多 64 字。")
        if (xy.ndim != 2 or xy.shape[1] != 2 or not 3 <= len(xy) <= 512 or not np.isfinite(xy).all()
                or (xy < 0).any() or (xy > [width - 1, height - 1]).any()):
            raise ValueError("请在 H&E 原图范围内标记 3 至 512 个有效顶点。")
        x, y = xy.T
        area = abs(np.sum(x * np.roll(y, -1) - y * np.roll(x, -1))) / 2
        if area < 4:
            raise ValueError("请绘制有面积的区域，按边界顺序选择顶点。")
        output.append({"name": name, "points": xy.tolist()})
        used.add(name)
    return output
