#!/usr/bin/env python3
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.special import gammaincc,gamma
C=299792.458;H0=70.391;OM=.33
CMB_RA=167.94;CMB_DEC=-6.94;CMB_V=369.82
BINS=np.arange(0.002,0.1221,0.002)
OUTER_WINDOWS=[(.060,.100),(.070,.100),(.075,.100),(.075,.110)]
INNER_WINDOWS=[(.005,.050),(.005,.075),(.010,.050),(.010,.075)]
MODEL_GRID=[
 ('LH11',-23.24,-.86,-2.9),
 ('LH11_noKE',-23.24,-.86,0.0),
 ('Mbright',-23.44,-.86,-2.9),
 ('Mfaint',-23.04,-.86,-2.9),
 ('alpha_shallow',-23.24,-.76,-2.9),
 ('alpha_steep',-23.24,-.96,-2.9),
]

def E(z):return np.sqrt(OM*(1+z)**3+1-OM)
def grids():
 z=np.linspace(1e-5,.15,50000);Ez=E(z);dc=(C/H0)*cumulative_trapezoid(1/Ez,z,initial=0);dl=(1+z)*dc;dvdz=(C/H0)*dc**2/Ez;return z,dc,dl,dvdz
ZG,DCG,DLG,DVDZG=grids()
def interp(arr,z):return np.interp(z,ZG,arr)
def unitvec(ra,dec):
 a=np.deg2rad(np.asarray(ra,float));d=np.deg2rad(np.asarray(dec,float));q=np.cos(d);return np.c_[q*np.cos(a),q*np.sin(a),np.sin(d)]
def apex_vec():return unitvec([CMB_RA],[CMB_DEC])[0]
def read2mrs(p):
 specs=[(0,16),(17,26),(27,36),(37,46),(47,56),(57,63),(78,84),(173,178)]
 names=['ID','RAdeg','DEdeg','GLON','GLAT','Kcmag','Ktmag','cz']
 d=pd.read_fwf(p,colspecs=specs,names=names)
 for c in names[1:]:d[c]=pd.to_numeric(d[c],errors='coerce')
 return d

def prepare(d,frame,magmode,mlim=11.5):
 x=d[np.isfinite(d.cz)&np.isfinite(d.RAdeg)&np.isfinite(d.DEdeg)].copy()
 nh=unitvec(x.RAdeg,x.DEdeg)
 cz=x.cz.to_numpy(float)
 if frame=='cmb':cz=cz+CMB_V*(nh@apex_vec())
 x['z']=np.maximum(cz/C,0)
 m=x.Ktmag if magmode=='total' else x.Kcmag
 x=x[np.isfinite(m)&(m<=mlim)&(x.z<.13)].copy()
 return x

def model_nz(z,Mh,alpha,kcoef,mlim):
 # M* = Mh + 5 log10 h. Additive count-model k+e = kcoef*z.
 h=H0/100;Mstar=Mh+5*np.log10(h)
 dl=interp(DLG,z);dvdz=interp(DVDZG,z)
 DM=5*np.log10(np.maximum(dl,1e-9))+25
 Delta=kcoef*z
 Mlim=mlim-DM-Delta
 xmin=10**(-.4*(Mlim-Mstar))
 s=alpha+1
 if s<=0: raise ValueError('alpha must be > -1 for gammaincc implementation')
 sel=gammaincc(s,xmin)*gamma(s)
 return dvdz*sel

def integrate_model(lo,hi,Mh,a,kc,mlim):
 z=np.linspace(lo,hi,3000);return np.trapz(model_nz(z,Mh,a,kc,mlim),z)
def observed_count(x,lo,hi):return int(((x.z>=lo)&(x.z<hi)).sum())
def test_setup(x,frame,magmode,model,Mh,a,kc,mlim):
 rows=[]
 for olo,ohi in OUTER_WINDOWS:
  no=observed_count(x,olo,ohi); mo=integrate_model(olo,ohi,Mh,a,kc,mlim)
  amp=no/mo if mo>0 else np.nan
  for ilo,ihi in INNER_WINDOWS:
   ni=observed_count(x,ilo,ihi); mi=integrate_model(ilo,ihi,Mh,a,kc,mlim)*amp
   ratio=ni/mi if mi>0 else np.nan
   # Poisson only for a transparent floor; outer normalization uncertainty included approximately.
   fracerr=math.sqrt(1/max(ni,1)+1/max(no,1));err=ratio*fracerr
   rows.append(dict(frame=frame,magmode=magmode,mlim=mlim,model=model,Mh=Mh,alpha=a,kcoef=kc,outer_lo=olo,outer_hi=ohi,N_outer=no,inner_lo=ilo,inner_hi=ihi,N_inner=ni,ratio_selfnorm=ratio,delta_selfnorm=ratio-1,ratio_err_poisson=err,significance=(ratio-1)/err if err>0 else np.nan))
 return rows

def bin_profile(x,model,Mh,a,kc,mlim,olo=.075,ohi=.100):
 mo=integrate_model(olo,ohi,Mh,a,kc,mlim);no=observed_count(x,olo,ohi);amp=no/mo
 out=[]
 for lo,hi in zip(BINS[:-1],BINS[1:]):
  n=observed_count(x,lo,hi);m=integrate_model(lo,hi,Mh,a,kc,mlim)*amp
  out.append(dict(frame='cmb',magmode='total',model=model,zlo=lo,zhi=hi,zmid=(lo+hi)/2,N=n,N_model=m,ratio=n/m if m>0 else np.nan))
 return out

def sector_relative(x,model,Mh,a,kc,mlim,ilo=.005,ihi=.075,olo=.075,ohi=.100):
 # Per longitude sector: compare inner/outer ratio to model inner/outer. Area cancels.
 mi=integrate_model(ilo,ihi,Mh,a,kc,mlim);mo=integrate_model(olo,ohi,Mh,a,kc,mlim);mr=mi/mo
 out=[]
 for sec in range(12):
  l0=sec*30;l1=(sec+1)*30;q=x[(x.GLON>=l0)&(x.GLON<l1)]
  ni=observed_count(q,ilo,ihi);no=observed_count(q,olo,ohi)
  ratio=(ni/max(no,1))/mr
  err=ratio*math.sqrt(1/max(ni,1)+1/max(no,1))
  out.append(dict(sector=sec,l0=l0,l1=l1,N_inner=ni,N_outer=no,ratio_selfnorm=ratio,delta_selfnorm=ratio-1,err_poisson=err))
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--data',required=True);ap.add_argument('--out',default='selfnorm');args=ap.parse_args();o=Path(args.out);o.mkdir(parents=True,exist_ok=True);d=read2mrs(args.data)
 rows=[]
 for frame in ['helio','cmb']:
  for mm in ['iso','total']:
   x=prepare(d,frame,mm,11.5)
   for model,Mh,a,kc in MODEL_GRID:rows+=test_setup(x,frame,mm,model,Mh,a,kc,11.5)
 R=pd.DataFrame(rows);R.to_csv(o/'selfnorm_grid.csv',index=False)
 # Frozen primary: CMB, total K<11.5, LH11 shape incl k+e, outer .075-.10.
 P=R[(R.frame=='cmb')&(R.magmode=='total')&(R.model=='LH11')&(R.outer_lo==.075)&(R.outer_hi==.100)].copy();P.to_csv(o/'primary_selfnorm.csv',index=False)
 x=prepare(d,'cmb','total',11.5);B=pd.DataFrame(bin_profile(x,'LH11',-23.24,-.86,-2.9,11.5));B.to_csv(o/'primary_nz_bins.csv',index=False)
 S=pd.DataFrame(sector_relative(x,'LH11',-23.24,-.86,-2.9,11.5));S.to_csv(o/'primary_sector_selfnorm.csv',index=False)
 # Summaries for inner .005-.075 and .005-.05 across all systematics.
 g75=R[(R.inner_lo==.005)&(R.inner_hi==.075)].copy();g50=R[(R.inner_lo==.005)&(R.inner_hi==.050)].copy()
 p75=P[(P.inner_lo==.005)&(P.inner_hi==.075)].iloc[0];p50=P[(P.inner_lo==.005)&(P.inner_hi==.050)].iloc[0]
 summary={
  'catalog_rows':int(len(d)),
  'primary':{'setup':'CMB,total K<11.5,LH11 selection shape,k+e=-2.9z,normalized at 0.075<z<0.10','z_lt_0p075':{'N_inner':int(p75.N_inner),'N_outer':int(p75.N_outer),'ratio':float(p75.ratio_selfnorm),'delta':float(p75.delta_selfnorm),'poisson_err':float(p75.ratio_err_poisson),'sigma':float(p75.significance)},'z_lt_0p05':{'N_inner':int(p50.N_inner),'N_outer':int(p50.N_outer),'ratio':float(p50.ratio_selfnorm),'delta':float(p50.delta_selfnorm),'poisson_err':float(p50.ratio_err_poisson),'sigma':float(p50.significance)}},
  'systematics_z_lt_0p075':{'n':int(len(g75)),'delta_median':float(g75.delta_selfnorm.median()),'delta_p10':float(g75.delta_selfnorm.quantile(.1)),'delta_p90':float(g75.delta_selfnorm.quantile(.9)),'fraction_delta_le_minus_0p15':float((g75.delta_selfnorm<=-.15).mean()),'fraction_delta_le_minus_0p20':float((g75.delta_selfnorm<=-.20).mean()),'fraction_delta_le_minus_0p25':float((g75.delta_selfnorm<=-.25).mean())},
  'systematics_z_lt_0p05':{'n':int(len(g50)),'delta_median':float(g50.delta_selfnorm.median()),'delta_p10':float(g50.delta_selfnorm.quantile(.1)),'delta_p90':float(g50.delta_selfnorm.quantile(.9))},
  'sector_primary':{'delta_median':float(S.delta_selfnorm.median()),'delta_min':float(S.delta_selfnorm.min()),'delta_max':float(S.delta_selfnorm.max()),'sectors_under':int((S.delta_selfnorm<0).sum()),'n':12}
 }
 (o/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print('\nPRIMARY\n',P.to_string(index=False));print('\nSYSTEMATICS z<.075\n',g75.sort_values('delta_selfnorm').to_string(index=False));print('\nSECTORS\n',S.to_string(index=False))
if __name__=='__main__':main()
