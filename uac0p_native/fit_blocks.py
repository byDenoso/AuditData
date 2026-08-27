import json, pathlib
import camb
from cobaya.run import run

ROOT=pathlib.Path('.').resolve(); OUT=ROOT/'block_results'; OUT.mkdir(exist_ok=True)
THEORY='fs_bao_likelihoods.reptvelocileptors'
COMPONENTS={
  'all_fs':'fs_bao_likelihoods.desi_fs_all_nolya',
  'all_fs_bao':'fs_bao_likelihoods.desi_fs_bao_all_nolya',
  'bgs':'fs_bao_likelihoods.desi_fs_bao_bgs',
  'lrg':'fs_bao_likelihoods.desi_fs_bao_lrg',
  'elg':'fs_bao_likelihoods.desi_fs_bao_elg',
  'qso':'fs_bao_likelihoods.desi_fs_bao_qso',
}
MODELS={
 'UAC0P':dict(kind='uac',H0=70.84499,ombh2=0.022920,omch2=0.124814,ns=0.98958,As=2.1035499892002698e-09,zc=10**3.81,fde_zc=0.080503,theta_i=2.89155),
 'LCDM':dict(kind='lcdm',H0=68.751,ombh2=0.0226180693,omch2=0.11698460023,ns=0.9669,As=2.0925566586443135e-09)
}

def info(model_name,m,block,component):
    params={'H0':m['H0'],'ombh2':m['ombh2'],'omch2':m['omch2'],'As':m['As'],'ns':m['ns'],'tau':0.0544,'mnu':0.06}
    if m['kind']=='uac': params|={'n':3.0,'zc':m['zc'],'fde_zc':m['fde_zc'],'theta_i':m['theta_i']}
    extra={'num_massive_neutrinos':1,'nnu':3.044,'lens_potential_accuracy':0}
    if m['kind']=='uac': extra|={'dark_energy_model':'EarlyQuintessence','use_zc':True}
    return {'theory':{'camb':{'extra_args':extra,'speed':2},THEORY:None},
            'likelihood':{component:None},'params':params,
            'sampler':{'minimize':{'method':'scipy','ignore_prior':False,'best_of':3,'seed':20260827,'max_evals':6000}},
            'output':str(OUT/f'{model_name}_{block}'),'force':True,'stop_at_error':True}

def parse(prefix):
    p=next(p for p in [OUT/f'{prefix}.minimum.txt',OUT/f'{prefix}.bestfit.txt'] if p.exists())
    lines=[x for x in p.read_text().splitlines() if x.strip()]
    return {k:float(v) for k,v in zip(lines[0].lstrip('#').split(),lines[-1].split())}

summary={'contract':'DESI DR1 official generated block likelihoods; fixed cosmology; native REPT; official nuisance priors/marginalization; scipy best_of=3','blocks':{},'models':MODELS}
for block,component in COMPONENTS.items():
    summary['blocks'][block]={'component':component,'models':{}}
    for mn,m in MODELS.items():
        print('RUN',block,mn,flush=True)
        updated,sampler=run(info(mn,m,block,component),no_mpi=True)
        tab=parse(f'{mn}_{block}')
        chis={k:v for k,v in tab.items() if k.startswith('chi2__')}
        summary['blocks'][block]['models'][mn]={'chi2':sum(chis.values()),'chi2_components':chis,'minuslogpost':tab.get('minuslogpost'),'table':tab}
    u=summary['blocks'][block]['models']['UAC0P']['chi2']; l=summary['blocks'][block]['models']['LCDM']['chi2']
    summary['blocks'][block]['delta_chi2_UAC_minus_LCDM']=u-l
    print('DELTA',block,u-l,flush=True)

# Exact no-QSO profile under separable generated blocks = BGS + LRG + ELG.
summary['derived']={
 'no_qso_delta_chi2':sum(summary['blocks'][b]['delta_chi2_UAC_minus_LCDM'] for b in ['bgs','lrg','elg']),
 'block_sum_delta_chi2':sum(summary['blocks'][b]['delta_chi2_UAC_minus_LCDM'] for b in ['bgs','lrg','elg','qso']),
 'fs_only_delta_chi2':summary['blocks']['all_fs']['delta_chi2_UAC_minus_LCDM'],
 'fs_bao_delta_chi2':summary['blocks']['all_fs_bao']['delta_chi2_UAC_minus_LCDM'],
 'incremental_BAO_effect_on_model_ordering':summary['blocks']['all_fs_bao']['delta_chi2_UAC_minus_LCDM']-summary['blocks']['all_fs']['delta_chi2_UAC_minus_LCDM']
}
(OUT/'NEXO_UAC0P_DESI_NATIVE_BLOCK_BATTERY.json').write_text(json.dumps(summary,indent=2,default=str))
print(json.dumps(summary['derived'],indent=2),flush=True)
