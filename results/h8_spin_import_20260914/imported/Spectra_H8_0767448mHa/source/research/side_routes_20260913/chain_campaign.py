"""Bounded independent amplitude-screening campaign; no main-task files changed."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from pathlib import Path
import argparse,json,time,hashlib
from research.side_routes_20260913.positive_chain import bounds


def case(name,m):
    t=[F(1) for _ in range(m-1)]
    if name=='bond_disorder':t=[F(2+(i%3),3) for i in range(m-1)]
    fields=[t[0]]+[t[i-1]+t[i] for i in range(1,m-1)]+[t[-1]]
    if name=='parent':v=[-2*a for a in t]
    elif name=='near_parent':v=[-2*a+F(1,10000) for a in t]
    elif name=='repulsive':v=[F(1)]*(m-1);fields=[F(0)]*m
    elif name=='bond_disorder':v=[-2*a+F(1,4) for a in t]
    elif name=='site_disorder':v=[F(1,2)]*(m-1);fields=[F([-2,1,3,-1,0][i%5],10) for i in range(m)]
    else:raise ValueError('Unknown case')
    return {'modes':m,'particles':m//2,'hopping':list(map(str,t)),'interaction':list(map(str,v)),
            'fields':list(map(str,fields)),'energy_unit':'t_reference=1; abstract lattice model'}


def experiment(task):
    name,m,x,r=task;model=case(name,m)
    amplitude={'sites':[str(r if i%2 else F(1)) for i in range(m)],'bonds':[str(x)]*(m-1)}
    receipt=bounds(model,amplitude)
    return {'case':name,'modes':m,'bond_weight':str(x),'alternating_site_weight':str(r),
            'model':model,'amplitude':amplitude,'claim':receipt}


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    names=['parent','near_parent','repulsive','bond_disorder','site_disorder']
    xs=list(map(F,['1/2','3/5','2/3','3/4','4/5','9/10','1','11/10','5/4','3/2']))
    rs=list(map(F,['4/5','1','5/4','3/2']))
    tasks=[(name,12,x,r) for name in names for x in xs for r in rs]
    with ProcessPoolExecutor(max_workers=4) as pool:screen=list(pool.map(experiment,tasks))
    screened=time.monotonic();(out/'screen.json').write_text(json.dumps(screen,indent=2)+'\n')
    winners={name:min((a for a in screen if a['case']==name),key=lambda a:F(a['claim']['width'])) for name in names}
    ladder=[]
    for name,w in winners.items():
        for m in (4,8,16,32,64):
            result=experiment((name,m,F(w['bond_weight']),F(w['alternating_site_weight'])))
            p=out/f'{name}_m{m}.json';p.write_text(json.dumps(result,indent=2)+'\n');ladder.append(result)
            print(name,m,result['claim']['width_float'],flush=True)
    summary={'screen_experiments':len(screen),'held_out_size_experiments':len(ladder),'screen_size':12,'workers':4,
        'screen_wall_seconds':screened-start,'total_wall_seconds':time.monotonic()-start,
        'sum_experiment_wall_seconds':sum(a['claim']['wall_seconds'] for a in screen+ladder),
        'selection':'Minimum exact interval width on M12, then unchanged amplitude parameters on held-out sizes.',
        'winners':{name:{k:w[k] for k in ('bond_weight','alternating_site_weight')} for name,w in winners.items()},
        'ladder':[{'case':a['case'],'modes':a['modes'],**a['claim']} for a in ladder],
        'scope':'200 parameter screens across five scenarios, not 200 independent mathematical hypotheses. Every width is exact for the declared abstract chain.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('DONE',summary['total_wall_seconds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();run(a.out)
