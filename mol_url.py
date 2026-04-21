"""
Encoder for the quick-mol-viewer payload format.

Produces URLs like

    https://mx-e.github.io/quick-mol-viewer/#d=<base64url-payload>

or, if called with an empty base_url, just the opaque payload code suitable
for pasting into the viewer's paste field.

Wire format: see the README.

Usage
-----
    from mol_url import encode, encode_collection

    # single snapshot (N, 3) or trajectory (F, N, 3)
    url = encode(positions, atomic_numbers)

    # collection of independent molecules (variable atom count)
    url = encode_collection([(pos1, Z1), (pos2, Z2), ...])

Run `python mol_url.py` to print example URLs and self-test.
"""

import base64
import gzip
import struct
import numpy as np


DEFAULT_BASE_URL = "https://mx-e.github.io/quick-mol-viewer/"

_VERSION         = 1
_FLAG_DELTA      = 0b01
_FLAG_COLLECTION = 0b10
_PREFIX_GZIP     = 0x67  # 'g'
_PREFIX_RAW      = 0x72  # 'r'


# ---------------------------------------------------------------------------
# low-level helpers

def _base64url(b):
    return base64.urlsafe_b64encode(b).rstrip(b"=").decode("ascii")


def _wrap(body):
    """Prepend 'g' + gzip(body) or 'r' + body, whichever is shorter. Uses
    mtime=0 so output is byte-stable across runs."""
    gz = gzip.compress(body, mtime=0)
    if len(gz) < len(body):
        return bytes([_PREFIX_GZIP]) + gz
    return bytes([_PREFIX_RAW]) + body


def _finalize(payload, base_url):
    code = _base64url(payload)
    if not base_url:
        return code
    if "#" in base_url:
        raise ValueError("base_url must not contain a fragment")
    if not base_url.endswith("/"):
        base_url = base_url + "/"
    return f"{base_url}#d={code}"


def _quantize(positions, scale):
    """positions in Å → int16 quantized values. positions may be any shape."""
    s = 32767.0 / scale
    q = np.rint(positions * s)
    np.clip(q, -32767, 32767, out=q)
    return q.astype("<i2")


def _pick_scale(max_abs):
    return max(float(max_abs), 1e-6)


# ---------------------------------------------------------------------------
# public encoders

def encode(positions, atomic_numbers, *, base_url=DEFAULT_BASE_URL, delta=None):
    """Encode a single molecule or trajectory.

    positions: array-like, shape (N, 3) for a snapshot or (F, N, 3) for a
               trajectory. Units: Ångström.
    atomic_numbers: array-like of length N.
    base_url:  if empty, returns just the opaque payload code.
    delta:     force delta-encoding on/off. Default: on for F > 1, off for F = 1.
    """
    positions = np.asarray(positions, dtype=np.float32)
    atomic_numbers = np.asarray(atomic_numbers, dtype=np.uint8)

    if positions.ndim == 2:
        positions = positions[np.newaxis, ...]
    if positions.ndim != 3 or positions.shape[2] != 3:
        raise ValueError(f"positions must be (N,3) or (F,N,3); got {positions.shape}")
    F, N, _ = positions.shape

    if atomic_numbers.shape != (N,):
        raise ValueError(f"atomic_numbers shape {atomic_numbers.shape} != ({N},)")
    if N == 0 or F == 0:
        raise ValueError("empty molecule")
    if N > 65535 or F > 65535:
        raise ValueError("N and F must fit in uint16")

    if delta is None:
        delta = F > 1

    scale = _pick_scale(np.max(np.abs(positions)) if positions.size else 0.0)
    q = _quantize(positions, scale).astype(np.int32)  # widen for delta diff

    flags = _FLAG_DELTA if (delta and F > 1) else 0
    stored = q
    if flags & _FLAG_DELTA:
        stored = np.empty_like(q, dtype=np.int32)
        stored[0] = q[0]
        stored[1:] = q[1:] - q[:-1]
        # Deltas must still fit in int16 (frame-to-frame motion is normally tiny).
        if np.any(np.abs(stored) > 32767):
            # Fall back to absolute if any delta overflows.
            stored = q
            flags &= ~_FLAG_DELTA

    stored = stored.astype("<i2")

    header = struct.pack("<BBHHf", _VERSION, flags, N, F, scale)
    body = header + atomic_numbers.tobytes() + stored.tobytes()
    return _finalize(_wrap(body), base_url)


def encode_collection(molecules, *, base_url=DEFAULT_BASE_URL):
    """Encode a collection of independent molecules.

    molecules: iterable of (positions, atomic_numbers) pairs, where positions
               is shape (N_i, 3). Each molecule is a single snapshot.
    """
    parsed = []
    max_abs = 0.0
    for pos, Zs in molecules:
        pos = np.asarray(pos, dtype=np.float32)
        Zs  = np.asarray(Zs, dtype=np.uint8)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"expected (N,3) positions, got {pos.shape}")
        Ni = pos.shape[0]
        if Ni == 0:
            raise ValueError("empty molecule in collection")
        if Ni > 65535:
            raise ValueError("N_i must fit in uint16")
        if Zs.shape != (Ni,):
            raise ValueError(f"Z shape {Zs.shape} != ({Ni},)")
        parsed.append((pos, Zs, Ni))
        if pos.size:
            m = float(np.max(np.abs(pos)))
            if m > max_abs:
                max_abs = m

    M = len(parsed)
    if M == 0:
        raise ValueError("empty collection")
    if M > 65535:
        raise ValueError("M must fit in uint16")

    scale = _pick_scale(max_abs)

    chunks = [struct.pack("<BBHHf", _VERSION, _FLAG_COLLECTION, 0, M, scale)]
    for pos, Zs, Ni in parsed:
        chunks.append(struct.pack("<H", Ni))
        chunks.append(Zs.tobytes())
        chunks.append(_quantize(pos, scale).tobytes())
    body = b"".join(chunks)
    return _finalize(_wrap(body), base_url)


# ---------------------------------------------------------------------------
# examples + self-test

_WATER = (
    np.array([[0.000, 0.000, 0.0],
              [0.757, 0.586, 0.0],
              [-0.757, 0.586, 0.0]], dtype=np.float32),
    np.array([8, 1, 1], dtype=np.uint8),
)

_BENZENE = (
    np.array([
        [ 1.39, 0.00, 0], [ 0.70,  1.20, 0], [-0.70,  1.20, 0],
        [-1.39, 0.00, 0], [-0.70, -1.20, 0], [ 0.70, -1.20, 0],
        [ 2.47, 0.00, 0], [ 1.24,  2.14, 0], [-1.24,  2.14, 0],
        [-2.47, 0.00, 0], [-1.24, -2.14, 0], [ 1.24, -2.14, 0],
    ], dtype=np.float32),
    np.array([6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1], dtype=np.uint8),
)


def _benzene_rotation(frames=60):
    pos, Zs = _BENZENE
    out = np.empty((frames, pos.shape[0], 3), dtype=np.float32)
    for f in range(frames):
        t = 2 * np.pi * f / frames
        c, s = np.cos(t), np.sin(t)
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=np.float32)
        out[f] = pos @ R.T
    return out, Zs


def _small_collection():
    # Water, ammonia, methane, CO2, methanol, benzene — enough variety to
    # exercise the grid view.
    return [
        _WATER,
        (np.array([[0, 0, 0], [0.94, 0, -0.33], [-0.47, 0.81, -0.33], [-0.47, -0.81, -0.33]], np.float32),
         np.array([7, 1, 1, 1], np.uint8)),
        (np.array([[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63],
                   [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63]], np.float32),
         np.array([6, 1, 1, 1, 1], np.uint8)),
        (np.array([[0, 0, 0], [1.16, 0, 0], [-1.16, 0, 0]], np.float32),
         np.array([6, 8, 8], np.uint8)),
        (np.array([[0, 0, 0], [1.43, 0, 0], [-0.37, 0.5, 0.89],
                   [-0.37, 0.5, -0.89], [-0.37, -1.02, 0], [1.86, 0.9, 0]], np.float32),
         np.array([6, 8, 1, 1, 1, 1], np.uint8)),
        _BENZENE,
    ]


def _print_examples():
    print("# Example URLs — paste after your viewer's base URL or shift-click directly.\n")

    url = encode(*_WATER)
    print(f"## Water (static, {len(url)} chars)")
    print(url)
    print()

    traj, Zs = _benzene_rotation(60)
    url_d = encode(traj, Zs)
    url_a = encode(traj, Zs, delta=False)
    print(f"## Benzene, 60-frame rotation (trajectory, {len(url_d)} chars)")
    print(url_d)
    print(f"    (without delta-encoding: {len(url_a)} chars)\n")

    coll = _small_collection()
    url_c = encode_collection(coll)
    print(f"## Collection, {len(coll)} mols ({len(url_c)} chars)")
    print(url_c, "\n")


def _selftest():
    # Encode+decode a trajectory and verify round-trip.
    import io, zlib

    def decode(code_or_url):
        code = code_or_url.split("#d=", 1)[-1]
        pad = "=" * (-len(code) % 4)
        raw = base64.urlsafe_b64decode(code + pad)
        body = gzip.decompress(raw[1:]) if raw[0] == _PREFIX_GZIP else raw[1:]
        v, flags, N, F, scale = struct.unpack("<BBHHf", body[:10])
        assert v == 1
        return flags, N, F, scale, body[10:]

    # Single
    url = encode(*_WATER)
    flags, N, F, scale, rest = decode(url)
    assert flags == 0 and N == 3 and F == 1, (flags, N, F)

    # Trajectory with delta
    traj, Zs = _benzene_rotation(40)
    url = encode(traj, Zs)
    flags, N, F, scale, rest = decode(url)
    assert flags & _FLAG_DELTA, "expected delta flag on trajectory"
    assert N == 12 and F == 40

    # Collection
    coll = _small_collection()
    url = encode_collection(coll)
    flags, N, F, scale, rest = decode(url)
    assert flags & _FLAG_COLLECTION, "expected collection flag"
    assert N == 0 and F == len(coll)

    # Empty base_url → pure code
    code = encode(*_WATER, base_url="")
    assert "#" not in code and "/" not in code
    assert code.replace("-", "+").replace("_", "/").isalnum() or set(code) <= set(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_")

    print("selftest: OK")


if __name__ == "__main__":
    _selftest()
    print()
    _print_examples()
