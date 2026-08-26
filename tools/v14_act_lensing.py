#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path
import numpy as np
from pixell import enmap
RA0=193.309380; DEC0=2.514639; Z0=.822175
RNG=np.random.default_rng(2608261402)

def angdist(dec,ra,dec0,ra0):
 return np.arccos(np.clip(np.sin(dec)*math.sin(dec0)+np.cos(dec)*math.cos(dec0)*np.cos(ra-ra0),-1,1))

def aper(m,ra0,dec0,rdeg=1.55,outer=2.8):
 dra=outer/max(math.cos(math.radians(dec0)),.2)
 box=np.deg2rad([[dec0-outer,ra0-dra],[dec0+outer,ra0+dra]])
 try:s=enmap.submap(m,box)
 except Exception:return None
 if s.size==0:return None
 if s.ndim>2:s=s.reshape((-1,)+s.shape[-2:])[0]
 pos=s.posmap(); d=angdist(pos[0],pos[1],math.radians(dec0),math.radians(ra0)); v=np.asarray(s,float)
 good=np.isfinite(v)
 inn=good&(d<=math.radians(rdeg)); ring=good&(d>math.radians(rdeg))&(d<=math.radians(outer))
 if inn.sum()<20:return None
 return {'inner_mean':float(np.mean(v[inn])),'inner_median':float(np.median(v[inn])),'ring_mean':float(np.mean(v[ring])) if ring.any() else None,'npix_inner':int(inn.sum()),'finite_fraction':float(good.mean())}

def score(p):
 n=p.name.lower(); s=0
 for k,v in [('kappa',10),('convergence',10),('baseline',5),('mv',4),('map',2),('data',2)]:
  if k in n:s+=v
 for k,v in [('mask',-20),('sim',-20),('noise',-10),('curl',-10),('meanfield',-10),('mf_',-10),('ivar',-10)]:
  if k in n:s+=v
 return s

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--out',default='actout');a=ap.parse_args();o=Path(a.out);o.mkdir(parents=True,exist_ok=True)
 fs=sorted(Path(a.root).rglob('*.fits'),key=score,reverse=True); print('FITS',len(fs));
 for f in fs[:50]:print(score(f),f)
 chosen=None; m=None; cand=None
 for f in fs:
  if score(f)<0:continue
  try:
   mm=enmap.read_map(str(f)); cc=aper(mm,RA0,DEC0)
   if cc is not None and cc['finite_fraction']>.1:
    chosen=f;m=mm;cand=cc;break
  except Exception as e: print('SKIP',f,type(e).__name__,str(e)[:160])
 if chosen is None:
  out={'status':'NO_USABLE_ACT_MAP_AT_V14','files':[str(x) for x in fs[:30]]};(o/'act_v14.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));return
 vals=[]
 ras=RNG.uniform(0,360,500)
 for ra in ras:
  q=aper(m,float(ra),DEC0)
  if q is not None and q['finite_fraction']>.5:vals.append(q['inner_mean'])
 vals=np.asarray(vals,float)
 if len(vals):
  z=(cand['inner_mean']-vals.mean())/max(vals.std(ddof=1),1e-30);p=(np.sum(vals<=cand['inner_mean'])+1)/(len(vals)+1)
 else:z=np.nan;p=np.nan
 out={'status':'OK','map':str(chosen),'candidate':cand,'control_n':int(len(vals)),'control_mean':float(vals.mean()) if len(vals) else None,'control_sd':float(vals.std(ddof=1)) if len(vals)>1 else None,'control_q05':float(np.quantile(vals,.05)) if len(vals) else None,'control_q50':float(np.quantile(vals,.5)) if len(vals) else None,'control_q95':float(np.quantile(vals,.95)) if len(vals) else None,'candidate_zscore':float(z),'lower_tail_empirical_p':float(p),'aperture_deg':1.55,'interpretation':'negative z supports projected matter deficit; null/positive disfavors lensing counterpart'}
 (o/'act_v14.json').write_text(json.dumps(out,indent=2));np.savetxt(o/'act_control_values.txt',vals);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
