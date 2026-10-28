# -*- coding: utf-8 -*-
# =============================================================================
#  export_abaqus_mode1.py  --  buckling mode 1 of finished Abaqus jobs -> CSV
#
#  Reads existing .odb files (nothing is re-run) and writes, for each job,
#  fe_modes/MODE1_<jobname>.csv with every node of all instances:
#      instance, node, x, y, z, U1, U2, U3      (header also gives the eigenvalue)
#  and one summary file fe_modes/job_summary.csv with, per job, the critical
#  eigenvalue and the times / model size read from the .dat and .msg files.
#
#  RUN (PowerShell), with the folder that contains the job folders, e.g.
#     abaqus python export_abaqus_mode1.py "$HOME\Downloads\omega_modes_runs" val_bs150_beta1p5_Unsym_PBC
#  arguments after the folder = job-name filters (substring); none = all jobs.
# =============================================================================
import os
import re
import sys

from odbAccess import openOdb


def eigval(desc):
    low = desc.lower()
    if "eigenvalue" in low:
        m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", desc[low.index("eigenvalue") + 10:])
        if m:
            return float(m.group(0))
    return None


def dat_info(base):
    """CPU / wallclock time and number of variables from <job>.dat / <job>.msg."""
    info = {"total_cpu_s": "", "wallclock_s": "", "variables": ""}
    for ext in (".dat", ".msg"):
        f = base + ext
        if not os.path.exists(f):
            continue
        txt = open(f, "r").read()
        for key, pat in (("total_cpu_s", r"TOTAL CPU TIME \(SEC\)\s*=\s*([\d.Ee+-]+)"),
                         ("wallclock_s", r"WALLCLOCK TIME \(SEC\)\s*=\s*([\d.Ee+-]+)"),
                         ("variables", r"TOTAL NUMBER OF VARIABLES IN THE MODEL\s+([\d]+)")):
            if not info[key]:
                m = re.search(pat, txt)
                if m:
                    info[key] = m.group(1)
    return info


root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.expanduser("~"), "Downloads", "omega_modes_runs")
filters = sys.argv[2:]
out = os.path.join(os.path.dirname(os.path.abspath(root)), "fe_modes")
if not os.path.isdir(out):
    os.makedirs(out)
summary = open(os.path.join(out, "job_summary.csv"), "w")
summary.write("job,Ncr_mode1,total_cpu_s,wallclock_s,variables\n")
for dirpath, _dirs, files in os.walk(root):
    for fn in sorted(files):
        if not fn.endswith(".odb"):
            continue
        job = fn[:-4]
        if filters and not any(f in job for f in filters):
            continue
        odb = openOdb(os.path.join(dirpath, fn), readOnly=True)
        step = odb.steps[list(odb.steps.keys())[-1]]
        best = None
        for fr in step.frames:
            ev = eigval(fr.description)
            if ev is not None and ev > 1.0e-3 and (best is None or ev < best[0]):
                best = (ev, fr)
        if best is None:
            print("%s: no positive eigenvalue" % job); odb.close(); continue
        with open(os.path.join(out, "MODE1_%s.csv" % job), "w") as fm:
            fm.write("instance,node,x,y,z,U1,U2,U3,eigenvalue=%.6f\n" % best[0])
            for ikey in odb.rootAssembly.instances.keys():
                ins = odb.rootAssembly.instances[ikey]
                cxyz = dict((n.label, n.coordinates) for n in ins.nodes)
                if not cxyz:
                    continue
                for v in best[1].fieldOutputs["U"].getSubset(region=ins).values:
                    c = cxyz.get(v.nodeLabel)
                    if c is not None:
                        fm.write("%s,%d,%.5f,%.5f,%.5f,%.8e,%.8e,%.8e\n" % (ikey, v.nodeLabel, c[0], c[1], c[2],
                                                                       v.data[0], v.data[1], v.data[2]))
        odb.close()
        inf = dat_info(os.path.join(dirpath, job))
        summary.write("%s,%.5f,%s,%s,%s\n" % (job, best[0], inf["total_cpu_s"], inf["wallclock_s"], inf["variables"]))
        print("%s: N_cr = %.4f  cpu = %s s  wall = %s s  variables = %s"
              % (job, best[0], inf["total_cpu_s"], inf["wallclock_s"], inf["variables"]))
summary.close()
print("written to", out)
