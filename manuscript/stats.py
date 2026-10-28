import pandas as pd, numpy as np
d=pd.read_csv('src/omega_results_V32.csv')
e=lambda a,b:100*(a-b)/b
d['eg']=e(d.gbt_sym,d.abq_sym); d['es']=e(d.strip_sym,d.abq_sym); d['egs']=e(d.gbt_sym,d.strip_sym)
d['egp']=e(d.gbt_per,d.abq_per); d['esp']=e(d.strip_per,d.abq_per); d['egsp']=e(d.gbt_per,d.strip_per)
d['sh_a']=e(d.abq_per,d.abq_sym); d['sh_g']=e(d.gbt_per,d.gbt_sym); d['sh_s']=e(d.strip_per,d.strip_sym)
L=['Isotropic','Cross-ply','Sym1','Sym2','Unsym']
f=lambda s:f"{s.abs().mean():.2f} | {s.abs().max():.2f}"
print("TABLE sym")
for l in L:
  for bc in ['SS','CC']:
    x=d[(d.layup==l)&(d.bc==bc)]
    print(f"| {l} | {bc} | {f(x.eg)} | {f(x.es)} | {f(x.egs)} |")
print(f"| all 60 | | {f(d.eg)} | {f(d.es)} | {f(d.egs)} |")
for bs in [150,285.5]:
  x=d[d.bs==bs]; print(bs,'sym gbt',f(x.eg),x.eg.min(),x.eg.max(),'strip',x.es.min(),x.es.max())
s=d[d.bc=='SS']
print("TABLE per by layup")
for l in L:
  x=s[s.layup==l]
  print(f"| {l} | {f(x.egp)} | {f(x.esp)} | {x.sh_a.min():+.2f} to {x.sh_a.max():+.2f} | {x.sh_g.min():+.2f} to {x.sh_g.max():+.2f} | {x.sh_s.min():+.2f} to {x.sh_s.max():+.2f} |")
print('all30',f(s.egp),f(s.esp))
for bs in [150,285.5]:
  x=s[s.bs==bs]; print(bs,'per',f(x.egp), x.egp.min(), x.egp.max())
dv=(s.sh_g-s.sh_a); print('shift dev gbt',dv.abs().mean(),dv.abs().max(),np.corrcoef(s.sh_g,s.sh_a)[0,1])
dv=(s.sh_s-s.sh_a); print('shift dev strip',dv.abs().mean(),dv.abs().max(),np.corrcoef(s.sh_s,s.sh_a)[0,1])
print('gbt vs strip per all60',f(d.egsp))
c=d[d.bc=='CC']; print('CC per gbt-strip',f(c.egsp))
print("TABLE CC shift")
for l in L:
  x=c[c.layup==l]; print(f"| {l} | {x.sh_g.min():+.2f} to {x.sh_g.max():+.2f} | {x.sh_s.min():+.2f} to {x.sh_s.max():+.2f} | {f(x.egsp)} |")
print("TABLE full SS periodic")
for _,r in s.iterrows():
  print(f"| {r.layup} | {r.bs:g} | {r.beta:g} | {r.abq_sym:.3f} | {r.abq_per:.3f} | {r.gbt_per:.3f} | {r.strip_per:.3f} | {r.egp:+.2f} | {r.esp:+.2f} | {r.sh_a:+.2f} | {r.sh_g:+.2f} |")
print("TABLE full CC")
for _,r in c.iterrows():
  print(f"| {r.layup} | {r.bs:g} | {r.beta:g} | {r.abq_sym:.3f} | {r.gbt_sym:.3f} | {r.strip_sym:.3f} | {r.eg:+.2f} | {r.gbt_per:.3f} | {r.strip_per:.3f} | {r.sh_g:+.2f} | {r.sh_s:+.2f} |")
print("TABLE full SS sym")
for _,r in s.iterrows():
  print(f"| {r.layup} | {r.bs:g} | {r.beta:g} | {r.abq_sym:.3f} | {r.gbt_sym:.3f} | {r.strip_sym:.3f} | {r.eg:+.2f} | {r.es:+.2f} |")
