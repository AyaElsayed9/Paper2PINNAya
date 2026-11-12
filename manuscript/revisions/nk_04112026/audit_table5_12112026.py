import docx, json, csv, numpy as np
M='/home/claude/paper2pinnaya/manuscript/'
d=docx.Document('../aya12/clean.docx'); T=d.tables[41]
num=lambda s: float(s.replace('−','-').replace('+',''))
pc=[json.loads(l) for l in open(M+'results/periodic_curves.jsonl')]
def gbt(lay,bs,beta,bc):
    r=[x for x in pc if x['layup']==lay and abs(x['bs']-bs)<1e-6 and abs(x['beta']-beta)<1e-9 and x['bc']==bc]
    return r[-1]['N'], len(r)
jobs={r['job']:float(r['Ncr_mode1']) for r in csv.DictReader(open(M+'src/abaqus_job_summary.csv'))}
cc={(r['layup'],float(r['bs']),float(r['beta'])):float(r['Ncr']) for r in csv.DictReader(open(M+'src/abaqus_cc_periodic.csv'))}
v32={(r['layup'],r['bc'],float(r['bs']),float(r['beta'])):r for r in csv.DictReader(open(M+'src/omega_results_V32.csv'))}
lay=None; es={'SS':[],'CC':[]}; bad=0; rows=[]
for r in T.rows[1:]:
    c=[x.text.strip() for x in r.cells]
    lay=c[0] or lay; bs=float(c[1]); beta=float(c[2])
    for bc,o in (('SS',3),('CC',6)):
        A,G,E=num(c[o]),num(c[o+1]),num(c[o+2])
        g,n=gbt(lay,bs,beta,bc)
        a_src = float(v32[(lay,'SS',bs,beta)]['abq_per']) if bc=='SS' else cc[(lay,bs,beta)]
        tag=lay.replace('-','m'); b=f"{beta:g}".replace('.','p')
        job=f"val_bs{bs:g}_beta{b}_{tag}{'_CC' if bc=='CC' else ''}_PBC".replace('285.5','285p5')
        aj=jobs.get(job)
        e=100*(g-a_src)/a_src
        ok = abs(A-a_src)<=6e-4*max(1,a_src/10) and abs(G-g)<0.0006*max(1,g/10) and abs(E-round(e,2))<0.0051
        if aj is not None and abs(aj-a_src)>1e-3: ok=False
        if not ok or n!=1:
            bad+=1; print('MISMATCH',lay,bs,beta,bc,'doc',A,G,E,'src',a_src,aj,round(g,4),round(e,3),'nrec',n)
        es[bc].append(e); rows.append((lay,bs,beta,bc,e))
print('bad',bad, 'jobs found', sum(1 for x in rows))
for bc in es:
    x=np.abs(es[bc]); print(bc,'mean %.3f max %.3f'%(x.mean(),x.max()))
import collections
for bc in ('SS',):
  for l in ['Isotropic','Cross-ply','Sym1','Sym2','Unsym']:
    x=[abs(e) for (a,b,c_,d_,e) in rows if a==l and d_==bc]; print(l,'%.3f'%np.mean(x))
  for bsv in (150,285.5):
    x=[e for (a,b,c_,d_,e) in rows if b==bsv and d_=='SS']; print(bsv,'mean %.3f max %.3f nneg %d'%(np.mean(np.abs(x)),np.max(np.abs(x)),sum(v<0 for v in x)))
# CC < SS count
ss={(a,b,c_):None for a,b,c_,d_,e in rows}
n=0;nn=0
for r in T.rows[1:]:
    pass
print('job matches', sum(1 for (a,b,c_,d_,e) in rows if jobs.get(f"val_bs{b:g}_beta{f'{c_:g}'.replace('.','p')}_{a.replace('-','m')}{'_CC' if d_=='CC' else ''}_PBC".replace('285.5','285p5')) is not None))
print([k for k in jobs][:8])
n=0
for (a,b,c_,d_,e) in rows:
    if d_=='SS':
        g1,_=gbt(a,b,c_,'SS'); g2,_=gbt(a,b,c_,'CC')
        a1=float(v32[(a,'SS',b,c_)]['abq_per']); a2=cc[(a,b,c_)]
        if g2<g1 or a2<a1: n+=1; print('CC<SS',a,b,c_, g2<g1, a2<a1)
print('n CC<SS',n)
g1,_=gbt('Isotropic',150,0.5,'SS'); g2,_=gbt('Isotropic',150,0.5,'CC'); print('raise b0.5', g2/g1-1)
