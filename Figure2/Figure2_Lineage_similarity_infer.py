# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:57:09 2024

@author: junyao
"""

import scanpy as sc
import numpy as np
import pandas as pd
import os
import sys
import anndata as an
import cospar as cs
args1 = sys.argv[1]
args2 = sys.argv[2]
clone = pd.read_csv(args1,index_col=0)
adata = an.AnnData(clone.T)
adata
adata.obs['state_info'] = clone.columns
adata.obsm["X_clone"] = clone.T
adata.obs['time_info'] = 'E18.5'
adata.uns["data_des"] = ["hi"]
cs.tl.fate_coupling(
            adata, source="X_clone"
        )
df1=pd.DataFrame(adata.uns['fate_coupling_X_clone']['X_coupling'])
df1.index=adata.uns['fate_coupling_X_clone']['fate_names']
df1.columns=adata.uns['fate_coupling_X_clone']['fate_names']
df1.to_csv(args2)