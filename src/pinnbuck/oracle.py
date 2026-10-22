"""Reference solutions: the spline finite strip (oracle) and full GBT.

`solve(sample, solver="strip", mesh="fast")` returns one record; `run_batch`
evaluates many samples in parallel, appends every record to a JSON-lines file as
soon as it is done, and skips samples already in that file (resumable).

Units: omega panel -> critical line load N_cr [N/mm]; members -> P_cr [N].
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import time
import traceback
from pathlib import Path

from .designspace import Sample
from .solver import configure, solver, solver_hash


# ---------------------------------------------------------------- inputs
def omega_panel(sample: Sample):
    m = solver(); fx = sample.fixed; t = fx["t_ply"]
    spec = dict(material=fx["material"], b3=fx["b3"], h=fx["h"], theta=fx["theta"],
                pskin=[(a, t) for a in sample.layups["skin"]],
                pstiff=[(a, t) for a in sample.layups["stringer"]])
    return m.make_omega_panel(sample.params["bs"], None, [], [], spec=spec)


def omega_length(sample: Sample) -> float:
    return sample.params["beta"] * sample.params["bs"]


def member_cfg(sample: Sample):
    m = solver(); fx = sample.fixed
    return m._study_cfg(sample.family, fx["material"], list(sample.layups["wall"]),
                        H=fx["H"], B=fx["B"], TPLY=fx["t_ply"])


# ---------------------------------------------------------------- solves
def _omega(sample: Sample, kind: str, mesh: str, nmodes: int):
    p = sample.params
    with configure(OMEGA_SIDE_BC=p["side"], OMEGA_MESH=mesh) as m:
        pan = omega_panel(sample); L = omega_length(sample)
        if kind == "strip":
            modes, info = m.omega_spline_eigs(pan, L, p["bc"], nmodes=nmodes)
        else:
            modes, info = m.omega_gbt_buckling(pan, L, p["bc"], nmodes=nmodes)
    return dict(value=modes[0][0], m=modes[0][1], eigs=[float(v) for v, _ in modes],
                mode=info.get("mode_type"), ndof=info.get("ndof"), unit="N/mm")


def _member(sample: Sample, kind: str, mesh: str, nmodes: int):
    p = sample.params; L = p["L"]
    mesh_sw = {} if mesh == "standard" else dict(MESH_STRIP=2.5, MESH_SPAN=10.0)
    with configure(**mesh_sw) as m:
        cfg = member_cfg(sample)
        if kind == "strip":
            r = m.unified_member_buckling(cfg, L, p["bc"], nev=nmodes)
            return dict(value=r["P"], m=r.get("m"), eigs=r.get("eigs"), mode=r.get("kind"),
                        ndof=r.get("ndof"), unit="N")
        solver_, fam, meth, seg, hx = m._gbt_opts({})
        sec = m.GBTSection(cfg, seg, fam, m.GBT_NMIN)
        r = m.gbt_fsm_solve(sec, float(L), p["bc"], method="spline", hx=min(hx, L / 20.0), nev=nmodes)
        return dict(value=r["P"], m=None, eigs=r.get("eigs"), mode=None, ndof=r.get("ndof"), unit="N")


def solve(sample: Sample, solver_kind: str = "strip", mesh: str = "fast", nmodes: int = 3) -> dict:
    """One reference solution; never raises (errors are recorded)."""
    t0 = time.perf_counter()
    rec = dict(sid=sample.sid, family=sample.family, solver=solver_kind, mesh=mesh,
               solver_hash=solver_hash(), sample=sample.to_dict())
    try:
        fn = _omega if sample.family == "omega" else _member
        rec.update(fn(sample, solver_kind, mesh, nmodes)); rec["ok"] = True
    except Exception as e:                                   # record, do not stop a batch
        rec.update(ok=False, error=f"{type(e).__name__}: {e}", trace=traceback.format_exc(limit=3))
    rec["time_s"] = round(time.perf_counter() - t0, 3)
    return rec


# ---------------------------------------------------------------- batches
def _worker(args):
    d, solver_kind, mesh = args
    s = Sample(d["family"], d["params"], d["layups"], d["fixed"])
    return solve(s, solver_kind, mesh)


def done_ids(path: Path, solver_kind: str, mesh: str) -> set:
    if not Path(path).exists():
        return set()
    out = set()
    for line in open(path):
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        if r.get("ok") and r.get("solver") == solver_kind and r.get("mesh") == mesh:
            out.add(r["sid"])
    return out


def run_batch(samples, out_path, solver_kind="strip", mesh="fast", workers=1, threads_per_worker=2,
              progress=True):
    """Evaluate samples not yet in out_path; append records as they finish."""
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    skip = done_ids(out_path, solver_kind, mesh)
    todo = [s for s in samples if s.sid not in skip]
    if progress:
        print(f"{len(samples)} samples, {len(skip)} already done, {len(todo)} to run "
              f"({solver_kind}, mesh={mesh}, workers={workers})", flush=True)
    jobs = [(dict(family=s.family, params=s.params, layups=s.layups, fixed=s.fixed), solver_kind, mesh) for s in todo]
    t0 = time.perf_counter(); n = 0
    with open(out_path, "a") as f:
        if workers <= 1:
            it = map(_worker, jobs)
            pool = None
        else:
            os.environ.setdefault("OMP_NUM_THREADS", str(threads_per_worker))
            os.environ.setdefault("MKL_NUM_THREADS", str(threads_per_worker))
            pool = mp.get_context("spawn").Pool(workers)       # spawn: same behaviour on Windows
            it = pool.imap_unordered(_worker, jobs)
        try:
            for rec in it:
                f.write(json.dumps(rec, default=float) + "\n"); f.flush(); n += 1
                if progress:
                    v = f"{rec['value']:.4f} {rec['unit']}" if rec.get("ok") else "FAILED " + rec.get("error", "")
                    print(f"  [{n}/{len(jobs)}] {rec['sid']} {v}  ({rec['time_s']:.1f} s, "
                          f"total {time.perf_counter() - t0:.0f} s)", flush=True)
        finally:
            if pool is not None:
                pool.close(); pool.join()
    return out_path


def load_records(path, solver_kind=None, mesh=None, ok_only=True) -> list[dict]:
    recs = []
    for line in open(path):
        r = json.loads(line)
        if ok_only and not r.get("ok"):
            continue
        if solver_kind and r["solver"] != solver_kind:
            continue
        if mesh and r["mesh"] != mesh:
            continue
        recs.append(r)
    return recs
