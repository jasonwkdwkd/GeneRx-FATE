"""Reproduce the fixed-split feature permutation test."""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

parser=argparse.ArgumentParser()
parser.add_argument('--permutations',type=int,default=200)
parser.add_argument('--seed',type=int,default=20261110)
parser.add_argument('--out', type=Path, default=Path('results/permutation'))
parser.add_argument('pairs', type=Path, help='GEO-reconstructed pair arrays')
args=parser.parse_args()
assert args.permutations>=19
out = args.out
out.mkdir(parents=True, exist_ok=True)
import generx_fate as gf
with np.load(args.pairs,allow_pickle=False) as d:
    a,b,y=[np.asarray(d[k],dtype=np.float64) for k in ('a','b','y')]
    pairs=d['pairs'].astype(str)
assert a.shape==b.shape==y.shape and pairs.shape==(len(a),2)
add=a+b; resid=y-add; x=gf.feature_matrix(a,b); mask=gf.panels(a,b,200)
genes=np.unique(pairs.ravel()); holdout_n=max(1,int(round(.25*len(genes))))
assert len(genes)==73 and len(pairs)==131 and holdout_n==18

# Exact same rowwise masked SSE function as the v0.3 ablation.
def sse_rows(pred, te):
    # gf._sse returns one masked SSE per pair.
    vals=np.asarray(gf._sse(y,pred,mask,np.asarray(te,dtype=int)),dtype=np.float64)
    assert vals.shape==(len(te),), (vals.shape,len(te))
    return vals

nperm=args.permutations
sum_sse=np.zeros((nperm+1,len(pairs)),dtype=np.float64)
counts=np.zeros(len(pairs),dtype=np.int64)
mean_sse=np.zeros(len(pairs),dtype=np.float64)
base_sse=np.zeros(len(pairs),dtype=np.float64)
add_sse=np.zeros(len(pairs),dtype=np.float64)
observed_seed_rows=[]
for seed in range(100,120):
    held=set(np.random.default_rng(seed).choice(genes,size=holdout_n,replace=False))
    tr=np.array([i for i,(u,v) in enumerate(pairs) if u not in held and v not in held],dtype=int)
    te=np.array([i for i,(u,v) in enumerate(pairs) if u in held and v in held],dtype=int)
    assert len(tr)>=10 and len(te)>=1
    assert set(pairs[tr].ravel()).isdisjoint(set(pairs[te].ravel()))
    base=gf.Baseline('B3').fit(add[tr],resid[tr],mask[tr])
    pb=add[te]+base.predict(add[te])
    target=np.zeros_like(resid[tr]); k=min(5,len(tr))
    for it,iv in KFold(n_splits=k,shuffle=True,random_state=int(seed)).split(tr):
        bl=gf.Baseline('B3').fit(add[tr[it]],resid[tr[it]],mask[tr[it]])
        target[iv]=resid[tr[iv]]-bl.predict(add[tr[iv]])
    pca=PCA(n_components=1,svd_solver='full').fit(target)
    scores=pca.transform(target)
    scaler=StandardScaler().fit(x[tr]); xt=scaler.transform(x[tr]); xv=scaler.transform(x[te])
    pm=pb+target.mean(axis=0)[None,:]
    sr_mean=sse_rows(pm,te)
    sr_b3=sse_rows(pb,te)
    sr_add=sse_rows(add[te],te)
    mean_sse[te]+=sr_mean; base_sse[te]+=sr_b3; add_sse[te]+=sr_add; counts[te]+=1
    # Each null replicate has its own independent seed-specific permutation.
    for p in range(nperm+1):
        if p==0: training_x=xt
        else:
            rng=np.random.default_rng(np.random.SeedSequence([args.seed,p,seed]))
            training_x=xt[rng.permutation(len(tr))]
        ridge=Ridge(alpha=100.0).fit(training_x,scores)
        pred_scores=np.asarray(ridge.predict(xv)).reshape(len(te),1)
        pf=pb+pca.inverse_transform(pred_scores)
        sum_sse[p,te]+=sse_rows(pf,te)
    observed_seed_rows.append({'seed':seed,'n_train':len(tr),'n_test':len(te)})
    print(f'seed={seed}: train={len(tr)} test={len(te)} finished {nperm} null permutations',flush=True)

valid=counts>0
assert int(counts.sum())==172,(int(counts.sum()),'expected 172')
assert int(valid.sum())==92,(int(valid.sum()),'expected 92')
# Equal weight per unique pair: average each pair's SSE across appearances.
mean_u=mean_sse[valid]/counts[valid]; base_u=base_sse[valid]/counts[valid]; add_u=add_sse[valid]/counts[valid]
full_u=sum_sse[:,valid]/counts[valid][None,:]
null_stat=1-full_u.sum(axis=1)/mean_u.sum()
observed=float(null_stat[0]); null=null_stat[1:]
p_emp=(1+int(np.count_nonzero(null>=observed)))/(1+nperm)
summary={'n_permutations':nperm,'permutation_seed':args.seed,'n_splits':20,
 'n_test_instances':int(counts.sum()),'n_unique_pairs':int(valid.sum()),
 'r2_b3_vs_add':float(1-base_u.sum()/add_u.sum()),
 'r2_mean_vs_add':float(1-mean_u.sum()/add_u.sum()),
 'r2_full_vs_add':float(1-full_u[0].sum()/add_u.sum()),
 'observed_full_vs_mean':observed,
 'null_mean':float(null.mean()),'null_std':float(null.std(ddof=1)),
 'null_q025':float(np.quantile(null,.025)),'null_q975':float(np.quantile(null,.975)),
 'n_null_ge_observed':int(np.count_nonzero(null>=observed)),
 'empirical_one_sided_p':float(p_emp),
 'null_definition':'training feature-row shuffles; heldout features and targets unchanged',
 'weighting':'average SSE per unique pair, then pooled SSE ratio'}
pd.DataFrame({'permutation':np.arange(nperm+1),'r2_full_vs_mean':null_stat,
 'is_observed':np.arange(nperm+1)==0}).to_csv(out/'permutation_distribution.csv',index=False)
pd.DataFrame({'pair_idx':np.flatnonzero(valid),'pair_a':pairs[valid,0],
 'pair_b':pairs[valid,1],'appearances':counts[valid],
 'sse_add':add_u,'sse_b3':base_u,'sse_mean':mean_u,
 'sse_full_observed':full_u[0]}).to_csv(out/'observed_unique_pairs.csv',index=False)
pd.DataFrame(observed_seed_rows).to_csv(out/'split_counts.csv',index=False)
(out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print('\n=== FEATURE PERMUTATION TEST (unique-pair weighted) ===')
print(json.dumps(summary,indent=2,ensure_ascii=False))
print('Saved outputs to',out)
