"""Bounded, inspectable CART regression and nearest-record baselines."""
from __future__ import annotations
import json
import numpy as np

def _data(X, y=None):
    X=np.asarray(X,dtype=float)
    if X.ndim!=2 or X.shape[0]<1 or X.shape[1]<1 or X.shape[1]>20 or X.shape[0]>10000 or not np.isfinite(X).all(): raise ValueError("invalid feature matrix")
    if y is None:return X
    y=np.asarray(y,dtype=float)
    if y.ndim!=1 or len(y)!=len(X) or not np.isfinite(y).all(): raise ValueError("invalid targets")
    return X,y

class CARTRegressor:
    def __init__(self,max_depth=4,min_leaf=8):
        if not isinstance(max_depth,int) or isinstance(max_depth,bool) or not 0<=max_depth<=8: raise ValueError("invalid max_depth")
        if not isinstance(min_leaf,int) or isinstance(min_leaf,bool) or min_leaf<1: raise ValueError("invalid min_leaf")
        self.max_depth=max_depth; self.min_leaf=min_leaf; self.fit_operations=0; self.split_candidates=0; self.sse_evaluations=0; self.featurevalues_scanned=0; self.predict_operations=0; self.tree=None
    def fit(self,X,y):
        X,y=_data(X,y); self.n_features_=X.shape[1]; self.fit_operations=0; self.split_candidates=0; self.sse_evaluations=0; self.featurevalues_scanned=0
        def build(ix,depth):
            vals=y[ix]; node={"value":float(vals.mean()),"count":len(ix)}
            if depth>=self.max_depth or len(ix)<2*self.min_leaf or np.all(vals==vals[0]): return node
            base=float(np.sum((vals-vals.mean())**2)); best=None
            for f in range(X.shape[1]):
                u=np.unique(X[ix,f]); self.featurevalues_scanned += len(u); self.split_candidates += len(u)-1
                for a,b in zip(u[:-1],u[1:]):
                    t=(a+b)/2; left=ix[X[ix,f]<=t]; right=ix[X[ix,f]>t]
                    if len(left)<self.min_leaf or len(right)<self.min_leaf: continue
                    self.sse_evaluations += 1
                    s=float(np.sum((y[left]-y[left].mean())**2)+np.sum((y[right]-y[right].mean())**2))
                    key=(s,f,float(t))
                    if best is None or key<best[0]: best=(key,left,right)
            if best is None or best[0][0]>=base: return node
            node.update(feature=best[0][1],threshold=best[0][2],left=build(best[1],depth+1),right=build(best[2],depth+1)); return node
        self.tree=build(np.arange(len(X)),0); self.fit_operations=self.split_candidates+self.sse_evaluations; return self
    def predict(self,X):
        X=_data(X); 
        if self.tree is None: raise ValueError("fit required")
        if X.shape[1] != self.n_features_: raise ValueError("feature dimension mismatch")
        out=[]; ops=0
        for row in X:
            n=self.tree
            while "feature" in n: ops+=1; n=n["left"] if row[n["feature"]]<=n["threshold"] else n["right"]
            out.append(n["value"])
        self.predict_operations=ops; return np.asarray(out)
    def to_dict(self):
        if self.tree is None: raise ValueError("fit required")
        return {"max_depth":self.max_depth,"min_leaf":self.min_leaf,"n_features":self.n_features_,"tree":self.tree}
    @classmethod
    def from_dict(cls,d):
        o=cls(d["max_depth"],d["min_leaf"]); o.n_features_=d["n_features"]
        if not isinstance(o.n_features_,int) or isinstance(o.n_features_,bool) or not 1<=o.n_features_<=20: raise ValueError("invalid serialized feature count")
        def check(n,depth):
            if not isinstance(n,dict) or not isinstance(n.get("count"),int) or n["count"]<1 or not np.isfinite(n.get("value",np.nan)): raise ValueError("invalid serialized node")
            if "feature" not in n:
                if "threshold" in n or "left" in n or "right" in n: raise ValueError("malformed leaf")
                return
            if depth>=o.max_depth or not isinstance(n["feature"],int) or isinstance(n["feature"],bool) or not 0<=n["feature"]<o.n_features_ or not np.isfinite(n["threshold"]): raise ValueError("invalid serialized split")
            if set(n)!={"value","count","feature","threshold","left","right"}: raise ValueError("malformed split")
            check(n["left"],depth+1); check(n["right"],depth+1)
        check(d["tree"],0); o.tree=d["tree"]; return o

class NearestRecordRegressor:
    def fit(self,X,y): X,y=_data(X,y); self.X,self.y=X.copy(),y.copy(); self.fit_operations=0; return self
    def predict(self,X):
        X=_data(X)
        if not hasattr(self,"X"): raise ValueError("fit required")
        if X.shape[1]!=self.X.shape[1]: raise ValueError("feature dimension mismatch")
        out=[]; ops=0
        for row in X:
            # Per-query vectorization bounds temporary memory at O(N*d), while
            # argmin preserves first-match tie behavior of the scalar loop.
            dist=np.sum((self.X-row)**2,axis=1); ops+=len(self.X)*X.shape[1]
            out.append(self.y[int(np.argmin(dist))])
        self.predict_operations=ops; return np.asarray(out)
