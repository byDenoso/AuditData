#!/usr/bin/env python3
import argparse,json,math
from pathlib import Path
import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord
import astropy.units as u
RA0=193.309380;DEC0=2.514639
RNG=np.random.default_rng(2608261403)

def aper(m,nside,vec,rdeg=1.55):
 pix=hp.query_disc(nside,vec,np.deg2rad(rdeg),inclusive=False)
 v=m[pix];v=v[np.isfinite(v)&(v!=hp.UNSEEN)]
 return float(np.mean(v)) if len(v) else np.nan

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--out',default='planckout');a=ap.parse_args();o=Path(a.out);o.mkdir(parents=True,exist_ok=True)
 root=Path(a.root); klms=list(root.rglob('dat_klm.fits')); masks=list(root.rglob('mask.fits.gz'))
 if not klms:
  out={'status':'NO_KLM'};(o/'planck_v14.json').write_text(json.dumps(out,indent=2));print(out);return
 f=klms[0]
 try:
  kalm=hp.read_alm(str(f),hdu=1)
 except Exception:
  kalm=hp.read_alm(str(f))
 lmax_in=hp.Alm.getlmax(len(kalm)); lmax=min(512,lmax_in)
 # dat_klm is the delivered Planck minimum-variance convergence K_LM estimator.
 # Truncate harmonics for degree-scale targeted diagnostic, do not apply an extra l(l+1)/2.
 if lmax<lmax_in:
  cut=np.zeros(hp.Alm.getsize(lmax),dtype=complex)
  for ell in range(lmax+1):
   for m in range(ell+1): cut[hp.Alm.getidx(lmax,ell,m)]=kalm[hp.Alm.getidx(lmax_in,ell,m)]
  kalm=cut
 knside=512; kappa=hp.alm2map(kalm,knside,lmax=lmax,verbose=False)
 mask=None
 if masks:
  try: mask=hp.read_map(str(masks[0]),verbose=False)
  except Exception: mask=None
 c=SkyCoord(ra=RA0*u.deg,dec=DEC0*u.deg,frame='icrs').galactic;l=float(c.l.deg);b=float(c.b.deg);vec=hp.ang2vec(l,b,lonlat=True)
 obs=aper(kappa,knside,vec)
 vals=[];tries=0
 while len(vals)<1000 and tries<30000:
  tries+=1;ll=float(RNG.uniform(0,360));bb=float(RNG.uniform(max(-85,b-10),min(85,b+10)));v=hp.ang2vec(ll,bb,lonlat=True)
  if mask is not None:
   p0=hp.vec2pix(hp.get_nside(mask),*v)
   if not np.isfinite(mask[p0]) or mask[p0]<=0:continue
  q=aper(kappa,knside,v)
  if np.isfinite(q):vals.append(q)
 vals=np.asarray(vals,float); z=float((obs-vals.mean())/vals.std(ddof=1));p=float((np.sum(vals<=obs)+1)/(len(vals)+1))
 out={'status':'OK','file':str(f),'lmax_input':int(lmax_in),'lmax_used':int(lmax),'kappa_nside':knside,'gal_l_deg':l,'gal_b_deg':b,'aperture_deg':1.55,'candidate_kappa':obs,'control_n':int(len(vals)),'control_mean':float(vals.mean()),'control_sd':float(vals.std(ddof=1)),'control_q05':float(np.quantile(vals,.05)),'control_q50':float(np.quantile(vals,.5)),'control_q95':float(np.quantile(vals,.95)),'candidate_zscore':z,'lower_tail_empirical_p':p,'interpretation':'negative kappa supports projected matter deficit; this is a targeted relative map diagnostic, not a Planck lensing likelihood'}
 (o/'planck_v14.json').write_text(json.dumps(out,indent=2));np.savetxt(o/'control.txt',vals);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
