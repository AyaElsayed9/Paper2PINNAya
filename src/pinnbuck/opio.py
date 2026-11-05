"""Reading and writing operator files (`*.pkl`).

From step B1 on, operator files are stored as zlib-compressed pickles (lossless, about 2.3x smaller, 0.05 s to
decompress) because ~1100 designs per side condition do not fit on the disk uncompressed.  `load` reads both
the old plain pickles and the compressed ones, so every older set still works.
"""
from __future__ import annotations

import os
import pickle
import zlib
from pathlib import Path

_MAGIC = b"PBZ1"                     # header of a compressed operator file


def load(path):
    raw = Path(path).read_bytes()
    if raw[:4] == _MAGIC:
        raw = zlib.decompress(raw[4:])
    return pickle.loads(raw)


def save(obj, path, level: int = 6):
    """Atomic write (temporary file, then rename): a stopped job never leaves a partial file."""
    path = Path(path)
    tmp = path.with_suffix(".tmp")
    tmp.write_bytes(_MAGIC + zlib.compress(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL), level))
    os.replace(tmp, path)


def is_compressed(path) -> bool:
    with open(path, "rb") as f:
        return f.read(4) == _MAGIC
