"""Design space of Paper 2: sections, parameter ranges, lay-ups, sampling.

Version-1 ranges are the ranges validated against Abaqus (see docs/DESIGN_SPACE.md).
Lay-ups are sampled ply by ply from ANGLE_SET, which contains every validated lay-up
(cross-ply, [45/-45], [60/-30], [0/30/60/90]); the network sees them through the
12 lamination parameters of each laminate (plus its thickness), so the inputs have a
fixed size whatever the number of plies.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, asdict

import numpy as np
from scipy.stats import qmc

ANGLE_SET = (0.0, 30.0, -30.0, 45.0, -45.0, 60.0, -60.0, 90.0)

CFRP_OMEGA = dict(type="orthotropic", E1=150000.0, E2=9000.0, G12=5300.0, G13=5300.0, G23=5300.0, nu12=0.33)
CFRP_MEMBER = dict(type="orthotropic", E1=171400.0, E2=9080.0, G12=5290.0, G13=5290.0, G23=2645.0, nu12=0.32)

# ------------------------------------------------------------------ families
FAMILIES = {
    # omega-stiffened panel (Abaqus model of Paper 1 / thesis), validated ranges
    "omega": dict(
        continuous={"bs": (150.0, 285.5), "beta": (0.75, 3.0)},
        fixed=dict(b3=27.75, h=30.0, theta=25.0, material=CFRP_OMEGA, t_ply=0.125,
                   n_skin=8, n_stringer=7),
        categorical={"bc": ("SS", "CC"), "side": ("symmetry", "periodic")},
    ),
    # thin-walled members (Abaqus parametric scripts): web H = 50, flange B = 20, 4 plies
    **{sec: dict(
        continuous={"L": (200.0, 1000.0)},
        fixed=dict(H=50.0, B=20.0, material=CFRP_MEMBER, t_ply=0.25, n_plies=4),
        categorical={"bc": ("SS", "CC")},
    ) for sec in ("c", "z", "l", "t", "i", "box")},
}


# ------------------------------------------------------------------ lay-ups
def sample_layup(n: int, rng: np.random.Generator, symmetric: bool) -> list[float]:
    """n ply angles from ANGLE_SET; symmetric=True mirrors the first half (B = 0)."""
    if symmetric:
        half = list(rng.choice(ANGLE_SET, size=n // 2))
        mid = [float(rng.choice(ANGLE_SET))] if n % 2 else []
        return [float(a) for a in half + mid + half[::-1]]
    return [float(a) for a in rng.choice(ANGLE_SET, size=n)]


INDUSTRIAL_ANGLES = (0.0, 45.0, -45.0, 90.0)


def industrial_layups(n: int) -> list[tuple[float, ...]]:
    """All stacking sequences of n plies that follow the usual industrial design rules:
    symmetric (B = 0), balanced (as many +45 as -45 plies), angles 0 / +-45 / 90 only,
    and at least 10 % of the plies in each of the four directions.
    n = 8 (skin): 24 sequences; n = 7 (stringer): 12 sequences."""
    from itertools import product
    half, odd = n // 2, n % 2
    out = set()
    for h in product(INDUSTRIAL_ANGLES, repeat=half):
        for mid in (INDUSTRIAL_ANGLES if odd else ((),)):
            seq = tuple(h) + ((mid,) if odd else ()) + tuple(h[::-1])
            cnt = {a: seq.count(a) for a in INDUSTRIAL_ANGLES}
            if cnt[45.0] == cnt[-45.0] and all(c >= 0.1 * n for c in cnt.values()):
                out.add(tuple(float(a) for a in seq))
    return sorted(out)


def lamination_parameters(angles_deg, t_ply: float) -> np.ndarray:
    """12 lamination parameters [xi_A(1..4), xi_B(1..4), xi_D(1..4)] (Tsai-Pagano).

    xi_A = int f dz~,  xi_B = 4 int f z~ dz~,  xi_D = 12 int f z~^2 dz~ over z~ in [-1/2, 1/2],
    f = [cos 2t, sin 2t, cos 4t, sin 4t]; ply 1 at the bottom (z~ = -1/2).
    """
    a = np.radians(np.asarray(angles_deg, float)); n = len(a)
    z = np.linspace(-0.5, 0.5, n + 1)
    f = np.stack([np.cos(2 * a), np.sin(2 * a), np.cos(4 * a), np.sin(4 * a)], axis=1)
    i0 = np.diff(z); i1 = np.diff(z ** 2) / 2; i2 = np.diff(z ** 3) / 3
    return np.concatenate([f.T @ i0, 4 * (f.T @ i1), 12 * (f.T @ i2)])


# ------------------------------------------------------------------ sample
@dataclass
class Sample:
    family: str
    params: dict                       # continuous + categorical values
    layups: dict                       # name -> list of ply angles [deg]
    fixed: dict = field(default_factory=dict)

    @property
    def sid(self) -> str:
        """Stable identifier (hash of the full definition)."""
        s = json.dumps(dict(f=self.family, p=self.params, l=self.layups, x=self.fixed), sort_keys=True, default=str)
        return hashlib.sha1(s.encode()).hexdigest()[:12]

    def features(self) -> np.ndarray:
        """Network input vector: scaled continuous parameters, categorical flags,
        lamination parameters and thickness of every laminate."""
        fam = FAMILIES[self.family]; out = []
        for k, (lo, hi) in fam["continuous"].items():
            out.append(2 * (self.params[k] - lo) / (hi - lo) - 1)
        for k, opts in fam["categorical"].items():
            out.append(float(opts.index(self.params[k])))
        tp = self.fixed["t_ply"]
        for name in sorted(self.layups):
            out += list(lamination_parameters(self.layups[name], tp))
        return np.asarray(out, float)

    def to_dict(self):
        return dict(sid=self.sid, **asdict(self))


def sample_design(family: str, n: int, seed: int = 0, symmetric_share: float = 0.5,
                  categorical: dict | None = None, layups: str = "general") -> list[Sample]:
    """n samples: scrambled Sobol for the continuous parameters, independent
    seeded draws for categorical values and lay-ups.  `categorical` fixes some
    categorical parameters (e.g. {"bc": "SS"}).  layups = "general" (ply by ply from
    ANGLE_SET) or "industrial" (uniform over `industrial_layups` of every laminate)."""
    fam = FAMILIES[family]; keys = list(fam["continuous"])
    sob = qmc.Sobol(d=len(keys), scramble=True, seed=seed)
    m = int(np.ceil(np.log2(max(n, 2))))
    u = sob.random_base2(m)[:n]
    lo = np.array([fam["continuous"][k][0] for k in keys]); hi = np.array([fam["continuous"][k][1] for k in keys])
    x = qmc.scale(u, lo, hi)
    rng = np.random.default_rng(seed + 1)
    out = []
    for i in range(n):
        p = {k: float(round(x[i, j], 4)) for j, k in enumerate(keys)}
        for k, opts in fam["categorical"].items():
            p[k] = (categorical or {}).get(k) or str(rng.choice(opts))
        sym = bool(rng.random() < symmetric_share)
        fx = dict(fam["fixed"])
        if layups == "industrial":
            pick = lambda k: [float(a) for a in (lambda L_: L_[rng.integers(len(L_))])(industrial_layups(k))]
            lay = ({"skin": pick(fx["n_skin"]), "stringer": pick(fx["n_stringer"])} if family == "omega"
                   else {"wall": pick(fx["n_plies"])})
            out.append(Sample(family, p, lay, fx))
            continue
        if family == "omega":
            lay = {"skin": sample_layup(fx["n_skin"], rng, sym), "stringer": sample_layup(fx["n_stringer"], rng, sym)}
        else:
            lay = {"wall": sample_layup(fx["n_plies"], rng, sym)}
        out.append(Sample(family, p, lay, fx))
    return out


# ------------------------------------------------------------------ validated cases
VALIDATED_LAYUPS = {
    "omega": {"Cross-ply": ([0, 90, 0, 90, 90, 0, 90, 0], [0, 90, 0, 90, 0, 90, 0]),
              "Sym1": ([45, -45, -45, 45, 45, -45, -45, 45], [45, -45, 45, -45, 45, -45, 45]),
              "Sym2": ([60, -30, -30, 60, 60, -30, -30, 60], [60, -30, 60, -30, 60, -30, 60]),
              "Unsym": ([0, 30, 60, 90, 0, 30, 60, 90], [0, 30, 60, 90, 0, 30, 60])},
    "member": {"Sym1": [45, -45, -45, 45], "Sym2": [60, -30, -30, 60], "Unsym": [0, 30, 60, 90]},
}


def validated_omega_sample(layup: str, bs: float, beta: float, bc: str, side: str) -> Sample:
    sk, st = VALIDATED_LAYUPS["omega"][layup]
    return Sample("omega", dict(bs=bs, beta=beta, bc=bc, side=side),
                  {"skin": [float(a) for a in sk], "stringer": [float(a) for a in st]},
                  dict(FAMILIES["omega"]["fixed"]))


def validated_member_sample(section: str, layup: str, L: float, bc: str) -> Sample:
    return Sample(section, dict(L=L, bc=bc), {"wall": [float(a) for a in VALIDATED_LAYUPS["member"][layup]]},
                  dict(FAMILIES[section]["fixed"]))
