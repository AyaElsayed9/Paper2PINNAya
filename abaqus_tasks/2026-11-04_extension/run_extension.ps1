# Abaqus runs for Paper 2, extended design space (4 Nov 2026)
# Put this file next to Omegastiffener_PBC_v2.py, open PowerShell in that folder and run:
#     powershell -ExecutionPolicy Bypass -File .\run_extension.ps1
# Each set writes into its own folder in Downloads. Send back the two CSV files of every folder.

$base = Join-Path $HOME "Downloads\abaqus_paper2_extension"

# Set A: periodic side edges, SS ends (9 cells, the first is a control)  -- needed
$env:PANEL_SIDES = "PBC"; $env:PANEL_BC = "SS"
$env:PANEL_CASES = "150:1.5:Sym1,143.12:1.53:Cross-ply,131.65:3.12:Sym1,129.3:1:Sym2,143.58:2.16:Unsym,147.12:1.34:Cross-ply,126:2.5:Sym1,267.14:3.12:Unsym,210.84:3.35:Sym2"
$env:PANEL_OUT = Join-Path $base "A_periodic_SS"
abaqus cae noGUI=Omegastiffener_PBC_v2.py

# Set B: symmetry side edges (as Paper 1), SS ends (7 cells, the first is a control)  -- needed
$env:PANEL_SIDES = "SYM"; $env:PANEL_BC = "SS"
$env:PANEL_CASES = "150:1.5:Sym1,149.3:2.19:Sym1,127.91:1.58:Sym2,134.4:3.06:Unsym,142.91:1.05:Sym1,139.72:3.46:Sym2,131.15:1.27:Unsym"
$env:PANEL_OUT = Join-Path $base "B_symmetry_SS"
abaqus cae noGUI=Omegastiffener_PBC_v2.py

# Set C: periodic side edges, CC ends (4 cells)  -- optional, for the clamped network later
$env:PANEL_SIDES = "PBC"; $env:PANEL_BC = "CC"
$env:PANEL_CASES = "126.57:1.31:Cross-ply,137.84:3.15:Sym1,145.03:1.58:Sym2,133.76:2.49:Unsym"
$env:PANEL_OUT = Join-Path $base "C_periodic_CC"
abaqus cae noGUI=Omegastiffener_PBC_v2.py

Remove-Item Env:PANEL_CASES, Env:PANEL_OUT, Env:PANEL_SIDES, Env:PANEL_BC
