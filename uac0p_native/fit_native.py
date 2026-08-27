import json, pathlib, yaml
import camb
from cobaya.run import run
from fs_bao_likelihoods.desi_fs_bao_all import desi_fs_bao_all as DESILikelihood
from fs_bao_likelihoods.reptvelocileptors import reptvelocileptors as REPTTheory

ROOT=pathlib.Path('.').resolve()
DATA=str(ROOT/'data/likelihood')
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
TRACERS=['bgs_z0','lrg_z0','lrg_z1','lrg_z2','elg_z1','qso_z0']
official=yaml.safe_load((ROOT/'desi_official/dr1/cobaya/desi_fs_bao_all.yaml').read_text())
nuisance={k:v for k,v in official['params'].items() if k.startswith('pre_')}

def camb_sigma8(model, As):
    common=dict(H0=model['H0'],ombh2=model['ombh2'],omch2=model['omch2'],mnu=0.06,
                tau=0.0544,As=As,ns=model['ns'],num_massive_neutrinos=1,nnu=3.044)
    if model['kind']=='uac':
        p=camb.set_params(**common,dark_energy_model='EarlyQuintessence',use_zc=True,
                          n=3.0,zc=model['zc'],fde_zc=model['fde_zc'],theta_i=model['theta_i'])
    else:
        p=camb.set_params(**common)
    p.set_matter_power(redshifts=[0.0],kmax=2.0)
    return float(camb.get_results(p).get_sigma8_0())

def tune_As(model):
    As=2.1e-9
    for _ in range(3):
        As *= (model['sigma8_target']/camb_sigma8(model,As))**2
    return As,camb_sigma8(model,As)

MODELS={
 'UAC0P':dict(kind='uac',H0=70.84499,ombh2=0.022920,omch2=0.124814,ns=0.98958,
              sigma8_target=0.81853,zc=10**3.81,fde_zc=0.080503,theta_i=2.89155),
 'LCDM':dict(kind='lcdm',H0=68.751,ombh2=0.0226180693,omch2=0.11698460023,ns=0.9669,
             sigma8_target=0.799182)
}

def make_info(name,m,As,tracers,suffix):
    params={'H0':m['H0'],'ombh2':m['ombh2'],'omch2':m['omch2'],'As':As,'ns':m['ns'],'tau':0.0544,'mnu':0.06}
    if m['kind']=='uac': params|={'n':3.0,'zc':m['zc'],'fde_zc':m['fde_zc'],'theta_i':m['theta_i']}
    params|=nuisance
    extra={'num_massive_neutrinos':1,'nnu':3.044,'lens_potential_accuracy':0}
    if m['kind']=='uac': extra|={'dark_energy_model':'EarlyQuintessence','use_zc':True}
    return {
      'theory':{
        'camb':{'extra_args':extra,'speed':2},
        'reptvelocileptors':{'external':REPTTheory}},
      'likelihood':{
        'desi_fs_bao_all':{'external':DESILikelihood,'data_dir':DATA,
          'observable_name':'spectrum-poles-rotated+bao-recon','tracers':tracers,'solve':'marg'}},
      'params':params,
      'sampler':{'minimize':{'method':'scipy','ignore_prior':False,'best_of':1,'seed':20260827,'max_evals':5000}},
      'output':str(OUT/f'{name}_{suffix}'),'force':True,'stop_at_error':True}

def point_dict(pt):
    out={'text':str(pt)}
    for a in ['logpost','logpriors','loglikes','derived']:
        try:
            v=getattr(pt,a); v=v() if callable(v) else v
            if hasattr(v,'tolist'): v=v.tolist()
            out[a]=v
        except Exception as e: out[a+'_error']=str(e)
    return out

def parse(prefix):
    candidates=[OUT/f'{prefix}.minimum.txt',OUT/f'{prefix}.bestfit.txt']
    path=next(p for p in candidates if p.exists())
    ls=[x for x in path.read_text().splitlines() if x.strip()]
    return {k:float(v) for k,v in zip(ls[0].lstrip('#').split(),ls[-1].split())}

summary={'contract':'official DESI DR1 FS+BAO six-tracer likelihood; fixed cosmology; native REPT/velocileptors; official nuisance priors+marginalization','tracers':TRACERS,'models':{},'environment':{'camb':camb.__version__}}
for name,m in MODELS.items():
    As,s=tune_As(m); m['As']=As; m['sigma8_realized']=s
    print(name,'As',As,'sigma8',s,flush=True)
    updated,sampler=run(make_info(name,m,As,TRACERS,'all'),no_mpi=True)
    prod=sampler.products()
    summary['models'][name]={'input':m,'bestfit':point_dict(prod['minimum']),'full_set_of_mins':str(prod.get('full_set_of_mins'))}
    (OUT/f'{name}_updated.json').write_text(json.dumps(updated,default=str,indent=2))
try:
    summary['tables']={'UAC0P':parse('UAC0P_all'),'LCDM':parse('LCDM_all')}
except Exception as e:
    summary['table_parse_error']=repr(e)
for name in ['UAC0P','LCDM']:
    tab=summary.get('tables',{}).get(name,{})
    chis={k:v for k,v in tab.items() if k.startswith('chi2__')}
    summary['models'][name]['chi2_components']=chis
    summary['models'][name]['chi2_likelihood_sum']=sum(chis.values()) if chis else None
u=summary['models']['UAC0P']['chi2_likelihood_sum']; l=summary['models']['LCDM']['chi2_likelihood_sum']
if u is not None and l is not None:
    summary['delta_chi2_UAC0P_minus_LCDM']=u-l
(OUT/'NEXO_UAC0P_DESI_DR1_NATIVE_FULLSHAPE_GATE.json').write_text(json.dumps(summary,indent=2,default=str))
print(json.dumps(summary,indent=2,default=str),flush=True)
