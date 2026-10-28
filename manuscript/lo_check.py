"""Render the docx with LibreOffice (math module installed) for a visual check of the equations.
The section-level <w:bidi/> of the MDPI template makes LibreOffice draw the equation tables
right-to-left (number on the left); it is removed in the check copy only."""
import os, re, shutil, subprocess, sys, zipfile, glob
from pathlib import Path
H = Path(__file__).parent
src = H / "out" / "GBT_periodic_omega_manuscript.docx"
d = H / "out" / "lo3"; shutil.rmtree(d, ignore_errors=True); d.mkdir()
z = zipfile.ZipFile(src); files = {n: z.read(n) for n in z.namelist()}
files["word/document.xml"] = re.sub(rb"<w:bidi ?/>", b"", files["word/document.xml"])
with zipfile.ZipFile(d / "check.docx", "w", zipfile.ZIP_DEFLATED) as zo:
    for n, b in files.items():
        zo.writestr(n, b)
soff = "/root/.claude/skills/synced/89684747-8288-4727-bcdd-0f93eb064259_c99e516e-a235-41ea-88e1-b79caa917733/docx/scripts/office/soffice.py"
subprocess.run([sys.executable, soff, "--headless", "--convert-to", "pdf", "check.docx"], cwd=d, capture_output=True)
first, last = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("1", "99")
subprocess.run(["pdftoppm", "-r", "75", "-png", "-f", first, "-l", last, "check.pdf", "p"], cwd=d)
print(subprocess.run(["pdfinfo", "check.pdf"], cwd=d, capture_output=True, text=True).stdout.split("Pages:")[1].split()[0], "pages")
print(sorted(os.path.basename(f) for f in glob.glob(str(d / "p-*.png"))))
