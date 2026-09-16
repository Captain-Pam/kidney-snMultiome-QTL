#!/usr/bin/env python
"""
plot_atac_interpret.py

Per prioritized eQTL variant, plot ref-allele contribution logos from the two
ATAC models, on the ATAC quantity each score was computed from:
  left : Cerberus ATAC ISM (logSUM, local; 8-fold ensemble), matched ATAC track
  right: ChromBPNet DeepLIFT-SHAP (counts; 5-fold ensemble), matched CT model
Variant at window center (dashed line). Window = +/- WIN bp.
"""
import os, numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logomaker

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figlib import (  # noqa: E402
    ALT_COLOR, ATAC_IDX, CT_COLORS, CT_LABELS, NT_COLORS, REF_COLOR,
    add_ols, add_ref_lines, cfg, model_ct, qtl_ct, select_lead,
    setup_style, style_ax,
)

setup_style()


B   = os.path.join(cfg("WORK_DIR"), "eqtl", "atac_interpret")
PRI = f"{B}/prioritized_variants.tsv"
BZ  = f"{B}/borzoi_atac_ism/ensemble_scores.h5"
OUT = os.path.join(cfg("FIGURE_DIR"), "eqtl_atac_interpretation")
WIN = 20
ATAC_IDX = {"CNT_CD_PC":0,"DCT":2,"DTL_ATL":4,"EC":6,"IC":8,"Immune":10,
            "PEC":12,"PTS":14,"Podocyte":16,"Stromal":18,"TAL":20,"injPT":22}

def ism_to_logo(ism, oh):                    # ism:(L,4) oh:(L,4)
    L=ism.shape[0]; c=np.zeros((L,4))
    for i in range(L):
        r=int(np.argmax(oh[i])); others=[j for j in range(4) if j!=r]
        c[i,r]=-float(np.mean(ism[i,others]))
    return pd.DataFrame(c,columns=list("ACGT"))

pri=pd.read_csv(PRI,sep="\t")
pri["scorer_vid"]=pri.chrom+":"+pri.pos.astype(str)+":"+pri.ref+":"+pri.alt
pri=pri.sort_values(["celltype","chrom","pos"]).reset_index(drop=True)

# Cerberus ISM ensemble
bz=h5py.File(BZ)
bz_labels=[x.decode() for x in bz["label"][:]]
bz_lab2i={v:i for i,v in enumerate(bz_labels)}

# ChromBPNet per-CT shap
def cbp_row(ct,vid_us):
    with h5py.File(f"{B}/chrombpnet_shap/{ct}_ism.h5") as f:
        vids=[x.decode() for x in f["variant_ids"][:]]
        if vid_us not in vids: return None,None
        i=vids.index(vid_us)
        return f["ref_contrib"][i], f.attrs["half_win"]

n=len(pri)
fig,axes=plt.subplots(n,2,figsize=(11,1.5*n))
if n==1: axes=axes[None,:]
for k,row in pri.iterrows():
    ct=row.celltype; vid=row.scorer_vid
    # ---- Cerberus ATAC ISM ----
    axL=axes[k,0]
    bi=bz_lab2i[vid]; tj=ATAC_IDX[ct]
    ism=bz["ref/logSUM"][bi,:,:,tj].astype(np.float32)   # (200,4)
    oh =bz["ref/seqs"][bi].astype(np.float32).T          # (200,4)
    c=100                                                # variant index in 200bp
    lg=ism_to_logo(ism[c-WIN:c+WIN+1], oh[c-WIN:c+WIN+1]); lg.index=range(-WIN,WIN+1)
    logomaker.Logo(lg,color_scheme=NT_COLORS,ax=axL)
    axL.axvline(0,ls="--",lw=0.6,color="grey")
    axL.set_title(f"{row.gene} | {ct} | β={row.eqtl_beta_aligned:+.2f}  "
                  f"Cerberus-ATAC logSUM={row.borzoi_atac_logSUM:+.1f}",fontsize=6)
    # ---- ChromBPNet SHAP ----
    axR=axes[k,1]
    vid_us=vid.replace(":","_")
    rc,hw=cbp_row(ct,vid_us)
    if rc is not None:
        cc=int(hw)                                       # variant index in 2*hw
        sub=rc[cc-WIN:cc+WIN+1]
        lg2=pd.DataFrame(sub,columns=list("ACGT")); lg2.index=range(-WIN,WIN+1)
        logomaker.Logo(lg2,color_scheme=NT_COLORS,ax=axR)
        axR.axvline(0,ls="--",lw=0.6,color="grey")
    axR.set_title(f"{row.gene} | {ct}  ChromBPNet logfc={row.chrombpnet_logfc:+.2f}",fontsize=6)
    for ax in (axL,axR):
        ax.tick_params(labelsize=5); [ax.spines[s].set_visible(False) for s in ["top","right"]]
axes[0,0].annotate("Cerberus ATAC ISM (logSUM, local)",xy=(0.5,1.0),xytext=(0.5,1.45),
                   xycoords="axes fraction",ha="center",fontsize=9,fontweight="bold")
axes[0,1].annotate("ChromBPNet DeepLIFT SHAP (counts)",xy=(0.5,1.0),xytext=(0.5,1.45),
                   xycoords="axes fraction",ha="center",fontsize=9,fontweight="bold")
fig.tight_layout(h_pad=1.4)
fig.savefig(OUT+".pdf",dpi=300,bbox_inches="tight")
fig.savefig(OUT+".png",dpi=150,bbox_inches="tight")
print("Saved",OUT+".pdf")
bz.close()
