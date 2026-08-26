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
 fs=list(Path(a.root).rglob('*.fits'));print('FITS',fs)
 if not fs:
  (o/'planck_v14.json').write_text(json.dumps({'status':'NO_FITS'},indent=2));return
 f=fs[0]
 try:
  phi,hdr=hp.read_map(str(f),field=0,h=True,verbose=False); hd=dict(hdr)
  try:mask=hp.read_map(str(f),field=1,verbose=False)
  except Exception:mask=np.ones_like(phi)
 except Exception as e:
  out={'status':'READ_FAIL','error':repr(e),'file':str(f)};(o/'planck_v14.json').write_text(json.dumps(out,indent=2));print(out);return
 nside=hp.get_nside(phi);nest=str(hd.get('ORDERING','RING')).upper().startswith('NEST')
 if nest:
  phi=hp.reorder(phi,n2r=True);mask=hp.reorder(mask,n2r=True)
 good=np.isfinite(phi)&(phi!=hp.UNSEEN)&(mask>0);phi=np.where(good,phi,0.0)
 lmax=512; alm=hp.map2alm(phi,lmax=lmax,iter=0); ell=np.arange(lmax+1); fl=.5*ell*(ell+1); kalm=hp.almxfl(alm,fl); knside=512;kappa=hp.alm2map(kalm,knside,lmax=lmax,verbose=False)
 # equatorial -> galactic because Planck lensing product is Galactic
 c=SkyCoord(ra=RA0*u.deg,dec=DEC0*u.deg,frame='icrs').galactic;l=float(c.l.deg);b=float(c.b.deg);vec=hp.ang2vec(l,b,lonlat=True)
 obs=aper(kappa,knside,vec)
 # controls conditioned on similar Galactic latitude, avoiding masked regions in original map
 vals=[];tries=0
 while len(vals)<1000 and tries<20000:
  tries+=1;ll=float(RNG.uniform(0,360));bb=float(RNG.uniform(max(-85,b-10),min(85,b+10)));v=hp.ang2vec(ll,bb,lonlat=True);p0=hp.vec2pix(nside,*v)
  if not good[p0]:continue
  q=aper(kappa,knside,v)
  if np.isfinite(q):vals.append(q)
 vals=np.asarray(vals);z=float((obs-vals.mean())/vals.std(ddof=1));p=float((np.sum(vals<=obs)+1)/(len(vals)+1))
 out={'status':'OK','file':str(f),'input_nside':int(nside),'kappa_nside':knside,'lmax':lmax,'gal_l_deg':l,'gal_b_deg':b,'aperture_deg':1.55,'candidate_kappa':obs,'control_n':int(len(vals)),'control_mean':float(vals.mean()),'control_sd':float(vals.std(ddof=1)),'control_q05':float(np.quantile(vals,.05)),'control_q50':float(np.quantile(vals,.5)),'control_q95':float(np.quantile(vals,.95)),'candidate_zscore':z,'lower_tail_empirical_p':p,'interpretation':'negative kappa supports projected matter deficit'}
 (o/'planck_v14.json').write_text(json.dumps(out,indent=2));np.savetxt(o/'control.txt',vals);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
