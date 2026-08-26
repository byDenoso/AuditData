#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.special import gammaincc, gamma

C=299792.458
H0=70.391
OM=0.33
FSKY=0.91
MLIM_ISO=11.75
# CMB dipole apex in ICRS, approximate solar speed wrt CMB.
CMB_RA=167.94; CMB_DEC=-6.94; CMB_V=369.82
SHELLS=[(40,70),(70,100),(100,150),(150,200),(200,250),(250,300),(300,350),(350,400)]
CUM_R=[100,150,200,250,300,350,400]
LF_GRID=[
    ('cole_base',-23.44,-0.96),
    ('brightM',-23.65,-0.96),
    ('faintM',-23.20,-0.96),
    ('shallow_alpha',-23.44,-0.85),
    ('steep_alpha',-23.44,-1.05),
]

def E(z):
    return np.sqrt(OM*(1+z)**3 + (1-OM))

def distances(z):
    z=np.asarray(z,float)
    zg=np.linspace(0,max(.16,float(np.nanmax(z))*1.02),30000)
    chi=(C/H0)*cumulative_trapezoid(1/E(zg),zg,initial=0)
    dc=np.interp(z,zg,chi)
    dl=(1+z)*dc
    return dc,dl

def unitvec(ra,dec):
    a=np.deg2rad(np.asarray(ra,float)); d=np.deg2rad(np.asarray(dec,float)); q=np.cos(d)
    return np.c_[q*np.cos(a),q*np.sin(a),np.sin(d)]

def apex_vec():
    return unitvec([CMB_RA],[CMB_DEC])[0]

def read_2mrs(path):
    # Fixed-width columns from CDS J/ApJS/199/26 table3.dat.
    colspecs=[(0,16),(17,26),(27,36),(37,46),(47,56),(57,63),(78,84),(173,178)]
    names=['ID','RAdeg','DEdeg','GLON','GLAT','Kcmag','Ktmag','cz']
    d=pd.read_fwf(path,colspecs=colspecs,names=names,na_values=['','     '])
    for c in names[1:]: d[c]=pd.to_numeric(d[c],errors='coerce')
    return d

def sector_access_frac(l0,l1):
    # Exact for our 30-deg sectors under the published 2MRS latitude mask.
    width=(l1-l0)/360.0
    center=(l0+l1)/2.0
    bmin=8.0 if (center<30 or center>=330) else 5.0
    return width*(1-math.sin(math.radians(bmin)))

def frame_catalog(d,frame,kcorr_on,mag_mode,lf_name,Mstar_h,alpha):
    x=d.copy()
    ok=np.isfinite(x.cz)&np.isfinite(x.RAdeg)&np.isfinite(x.DEdeg)&np.isfinite(x.Kcmag)&(x.Kcmag<=MLIM_ISO)
    x=x.loc[ok].copy()
    nhat=unitvec(x.RAdeg,x.DEdeg)
    if frame=='cmb':
        cz=x.cz.to_numpy(float) + CMB_V*(nhat@apex_vec())
    else:
        cz=x.cz.to_numpy(float)
    z=np.maximum(cz/C,1e-5)
    dc,dl=distances(z)
    x['z']=z; x['dc']=dc; x['dl']=dl
    kcorr=-6*np.log10(1+z) if kcorr_on else np.zeros_like(z)
    h=H0/100.0
    Mstar=Mstar_h+5*np.log10(h)
    # Total magnitude is a closer match to the Cole et al. Kron LF. 2MRS selection is isophotal,
    # so infer an effective total-mag limit from the observed median total-isophotal offset.
    delta_med=float(np.nanmedian(x.Ktmag-x.Kcmag))
    if mag_mode=='total':
        m=x.Ktmag.to_numpy(float); mlim=MLIM_ISO+delta_med
    else:
        m=x.Kcmag.to_numpy(float); mlim=MLIM_ISO
    DM=5*np.log10(np.maximum(dl,1e-8))+25
    M=m-DM-kcorr
    Mlim=mlim-DM-kcorr
    L=10**(-0.4*(M-Mstar))
    xmin=10**(-0.4*(Mlim-Mstar))
    fracL=np.clip(gammaincc(alpha+2,xmin),1e-5,1.0)
    x['Lstar']=L; x['fracL']=fracL; x['Lcorr']=L/fracL
    x['setup']=f'{frame}|kc{int(kcorr_on)}|{mag_mode}|{lf_name}'
    return x,dict(frame=frame,kcorr=kcorr_on,mag_mode=mag_mode,lf=lf_name,Mstar_h=Mstar_h,alpha=alpha,Mstar=Mstar,delta_mag_total_iso=delta_med)

def j_expected(alpha):
    h=H0/100.0
    phistar=0.0108*h**3
    return phistar*gamma(alpha+2)

def summarize_setup(x,meta):
    rows=[]; alpha=meta['alpha']; jexp=j_expected(alpha)
    for a,b in SHELLS:
        q=x[(x.dc>=a)&(x.dc<b)]
        vol=FSKY*4*math.pi/3*(b**3-a**3)
        j=float(q.Lcorr.sum()/vol)
        rows.append({**meta,'kind':'shell','rlo':a,'rhi':b,'N':len(q),'j_Lstar_Mpc3':j,'ratio_external':j/jexp,'delta_external':j/jexp-1})
    for R in CUM_R:
        q=x[(x.dc>=40)&(x.dc<R)]
        vol=FSKY*4*math.pi/3*(R**3-40**3)
        j=float(q.Lcorr.sum()/vol)
        rows.append({**meta,'kind':'cumulative','rlo':40,'rhi':R,'N':len(q),'j_Lstar_Mpc3':j,'ratio_external':j/jexp,'delta_external':j/jexp-1})
    # Internal ratios to an outer reference. Useful because LF normalization cancels.
    for reflo,refhi in [(300,400),(300,450)]:
        qr=x[(x.dc>=reflo)&(x.dc<refhi)]
        vr=FSKY*4*math.pi/3*(refhi**3-reflo**3)
        jr=float(qr.Lcorr.sum()/vr)
        for R in [200,250,300]:
            qi=x[(x.dc>=40)&(x.dc<R)]
            vi=FSKY*4*math.pi/3*(R**3-40**3)
            ji=float(qi.Lcorr.sum()/vi)
            rows.append({**meta,'kind':f'internal_ref_{reflo}_{refhi}','rlo':40,'rhi':R,'N':len(qi),'N_ref':len(qr),'j_Lstar_Mpc3':ji,'j_ref':jr,'ratio_external':np.nan,'delta_external':np.nan,'ratio_internal':ji/jr if jr>0 else np.nan,'delta_internal':ji/jr-1 if jr>0 else np.nan})
    return rows

def sector_test(x,meta,R=300):
    out=[]; jexp=j_expected(meta['alpha'])
    for i in range(12):
        l0=30*i;l1=30*(i+1); f=sector_access_frac(l0,l1)
        q=x[(x.dc>=40)&(x.dc<R)&(x.GLON>=l0)&(x.GLON<l1)]
        vol=f*4*math.pi/3*(R**3-40**3)
        j=float(q.Lcorr.sum()/vol)
        out.append({**meta,'R':R,'sector':i,'l0':l0,'l1':l1,'N':len(q),'access_sky_fraction':f,'ratio_external':j/jexp,'delta_external':j/jexp-1})
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',required=True);ap.add_argument('--out',default='density2mrs');a=ap.parse_args()
    o=Path(a.out);o.mkdir(parents=True,exist_ok=True)
    d=read_2mrs(a.data)
    allrows=[]; sectors=[]
    for frame in ['helio','cmb']:
      for kc in [False,True]:
       for mag_mode in ['total','iso']:
        for lf_name,Mh,alpha in LF_GRID:
            x,meta=frame_catalog(d,frame,kc,mag_mode,lf_name,Mh,alpha)
            allrows += summarize_setup(x,meta)
            if lf_name=='cole_base': sectors += sector_test(x,meta,300)
    R=pd.DataFrame(allrows); S=pd.DataFrame(sectors)
    R.to_csv(o/'density_profiles.csv',index=False);S.to_csv(o/'sector_300mpc.csv',index=False)
    # Frozen primary: CMB-frame, K-correction on, total magnitudes, Cole LF.
    primary=R[(R.frame=='cmb')&(R.kcorr==True)&(R.mag_mode=='total')&(R.lf=='cole_base')].copy()
    primary.to_csv(o/'primary_profile.csv',index=False)
    # Systematics distribution for the exact PEER gate (40-300 Mpc).
    gate=R[(R.kind=='cumulative')&(R.rlo==40)&(R.rhi==300)].copy()
    gate.to_csv(o/'gate_40_300_systematics.csv',index=False)
    # Sector primary for anisotropy.
    sp=S[(S.frame=='cmb')&(S.kcorr==True)&(S.mag_mode=='total')&(S.lf=='cole_base')].copy()
    sp.to_csv(o/'primary_sectors_300.csv',index=False)
    prow=primary[(primary.kind=='cumulative')&(primary.rhi==300)].iloc[0]
    internal=primary[(primary.kind=='internal_ref_300_400')&(primary.rhi==300)].iloc[0]
    # Translate luminosity contrast to a rough matter interval for b_K in [1,1.2].
    dl=float(prow.delta_external)
    dm_b1=dl; dm_b12=dl/1.2
    target=-0.25
    summary={
      'catalog_rows':int(len(d)),
      'primary_setup':'cmb|Kcorr|total|Cole2001/2dF K LF',
      'primary_40_300':{'N':int(prow.N),'ratio_to_external_K_luminosity_density':float(prow.ratio_external),'delta_K_luminosity':dl,'rough_delta_matter_b1':dm_b1,'rough_delta_matter_b1p2':dm_b12},
      'primary_internal_40_300_vs_300_400':{'N_inner':int(internal.N),'N_outer':int(internal.N_ref),'ratio':float(internal.ratio_internal),'delta':float(internal.delta_internal)},
      'systematics_40_300':{'ratio_median':float(gate.ratio_external.median()),'ratio_p10':float(gate.ratio_external.quantile(.1)),'ratio_p90':float(gate.ratio_external.quantile(.9)),'delta_median':float(gate.delta_external.median()),'n_setups':int(len(gate))},
      'sector_primary_300':{'ratio_median':float(sp.ratio_external.median()),'ratio_min':float(sp.ratio_external.min()),'ratio_max':float(sp.ratio_external.max()),'sectors_below_0p75':int((sp.ratio_external<.75).sum()),'n_sectors':int(len(sp))},
      'PEER_required_mass_delta_approx':target,
      'passes_PEER_delta_gate_if_b1':bool(dm_b1<=target),
      'passes_PEER_delta_gate_if_b1p2':bool(dm_b12<=target)
    }
    (o/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    print('\nPRIMARY PROFILE\n',primary.to_string(index=False))
    print('\n40-300 SYSTEMATICS\n',gate[['frame','kcorr','mag_mode','lf','N','ratio_external','delta_external']].sort_values('ratio_external').to_string(index=False))
    print('\nPRIMARY SECTORS\n',sp[['sector','l0','l1','N','ratio_external','delta_external']].to_string(index=False))

if __name__=='__main__': main()
