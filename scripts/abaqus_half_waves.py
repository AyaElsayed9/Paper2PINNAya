"""18 Nov 2026: half-wave number of Abaqus mode files (MODE1_*.csv from export_abaqus_mode1.py), per instance and component.
    cd data/abaqus_sym1_modes && gunzip -k *.gz && python ../../scripts/abaqus_half_waves.py"""
import glob, numpy as np, pandas as pd
def dom(z, Y, L, J=40):
    """Y (nz, ncols) sampled at z (uniform grid assumed); dominant harmonic share"""
    tot=(Y**2).sum(); best=(0,0)
    for j in range(1,J+1):
        B=np.stack([np.sin(j*np.pi*z/L),np.cos(j*np.pi*z/L)],1)
        c,*_=np.linalg.lstsq(B,Y,rcond=None); sh=((B@c)**2).sum()/tot
        if sh>best[1]: best=(j,sh)
    return best
for f in sorted(glob.glob("*MODE1*.csv")):
    d=pd.read_csv(f); ev=float(d.columns[-1].split("=")[1]); d.columns=list(d.columns[:-1])+["U3"] if False else d.columns
    d=d.rename(columns={d.columns[-1]:"U3"})
    L=d.z.max(); out=[f.split("MODE1_")[1][:-4], "N %.3f"%ev, "L %.0f"%L]
    for inst in ["INST_SKIN","INST_STIFF"]:
        s=d[d.instance==inst]
        for comp in ["U1","U2","U3"]:
            P=s.pivot_table(index="z",columns=["x","y"],values=comp).dropna(axis=1)
            z=P.index.values; Y=P.values
            out.append("%s.%s m=%d(%.2f) rms=%.2g"%(inst[5:9],comp,*dom(z,Y,L),np.sqrt((Y**2).mean())))
    print(" | ".join(out))
