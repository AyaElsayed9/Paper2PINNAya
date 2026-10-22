"""Generate reference data for one family from a YAML config.

    python scripts/generate_oracle_data.py configs/omega_pilot.yaml

Resumable: records already in the output file are skipped.  Every record carries
the sample definition, the solver hash, the mesh and the run time.
"""
import argparse
import json
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pinnbuck.designspace import sample_design            # noqa: E402
from pinnbuck.oracle import run_batch                      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--workers", type=int, default=None, help="override the config")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    root = Path(a.config).resolve().parents[1]
    samples = sample_design(cfg["family"], cfg["n_samples"], seed=cfg.get("seed", 0),
                            symmetric_share=cfg.get("symmetric_share", 0.5),
                            categorical=cfg.get("fix", None))
    out_dir = root / cfg["out_dir"]; out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "samples.jsonl", "w") as f:
        for s in samples:
            f.write(json.dumps(s.to_dict(), default=float) + "\n")
    for job in cfg["jobs"]:
        run_batch(samples, out_dir / f"{job['solver']}_{job['mesh']}.jsonl", job["solver"], job["mesh"],
                  workers=a.workers or cfg.get("workers", 1), threads_per_worker=cfg.get("threads_per_worker", 2))


if __name__ == "__main__":          # required for multiprocessing on Windows
    main()
