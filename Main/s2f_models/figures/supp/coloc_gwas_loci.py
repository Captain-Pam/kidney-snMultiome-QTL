#!/usr/bin/env python
"""
gwas_locus_figure.py  <rs_id>   (or no arg = all in config.json)

Parameterized version of borzoi_gwas_figure.py: for an eGFR credible-set GWAS
variant with eQTL+caQTL coloc, builds the locus figure:
  Manhattan / PIP / model-prioritization lollipops /
  ChromBPNet(left)+Cerberus(right) ATAC coverage+contrib+motif /
  Cerberus RNA pred + gene + RNA ISM + motif.
Motif is auto-selected as the best CWM (TomTom-style) match to the alt ATAC ISM.
Handles + and - strand eGenes.
"""
import os, re, sys, json, glob
import numpy as np, pandas as pd, h5py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
import logomaker

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figlib import (  # noqa: E402
    ALT_COLOR, ATAC_IDX, CT_COLORS, CT_LABELS, NT_COLORS, REF_COLOR,
    add_ols, add_ref_lines, cfg, model_ct, pwm_to_ic_df, qtl_ct, rc_matrix,
    select_lead, setup_style, style_ax,
)

setup_style()

# ── Paths / constants ─────────────────────────────────────────────────────────
SUS=cfg("GWAS_SUSIE_DIR")
ROOT=os.path.join(cfg("WORK_DIR"), "gwas", "loci")
GTF_FILE=cfg("EQTL_GTF")
GENOME=cfg("CERBERUS_FASTA")
MOTIF_DB=cfg("MOTIF_DB")
MACS2=os.path.join(cfg("CHROMBPNET_MACS2_DIR"), "{ct}_peaks.narrowPeak")
OUTDIR=os.path.join(cfg("FIGURE_DIR"), "coloc_gwas_loci"); os.makedirs(OUTDIR,exist_ok=True)
WIDE_FLANK=200_000; ATAC_PRED_FLANK=1_000; CBP_COV_FLANK=500
ATAC_ISM_TIDX_DEFAULT=10
CS_COLOR="#E07B39"; GRAY_DASH=dict(color="#888888",lw=0.6,ls="--"); PEAK_SHADE="#d9d9d9"

# ── motif PWM helpers ─────────────────────────────────────────────────────────
def icw(p):
    ic=np.log2(np.clip(p,1e-9,1))-np.log2(0.25); ic[ic<0]=0; return p*ic
def place_pwm(prob,is_rc,bp,n):
    """Build an IC logo with the motif placed at index bp (clamped), given orientation."""
    use=rc_matrix(prob) if is_rc else prob; L=len(prob)
    bp=max(0,min(int(bp),n-L)); al=np.zeros((n,4)); al[bp:bp+L]=np.clip(use,1e-6,1)
    df=pwm_to_ic_df(np.clip(al,1e-6,1)); df.iloc[:bp]=0.0; df.iloc[bp+L:]=0.0
    return df

# ── Motif helpers, local to this script ───────────────────────────────────────
# This panel searches the whole JASPAR database rather than using a named motif,
# so it needs a scored variant of the alignment that figlib does not provide.

def load_all_motifs():
    """Every motif in the MEME database as {id: (name, probability matrix)}."""
    text = open(MOTIF_DB).read()
    out = {}
    for m in re.finditer(r"MOTIF (\S+)\s+(\S+)?\n(.*?)(?=\nMOTIF|\Z)", text, re.DOTALL):
        mid, name, body = m.group(1), m.group(2), m.group(3)
        rows = [l.split() for l in body.splitlines()
                if re.match(r"^\s*[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+", l)]
        if not rows:
            continue
        p = np.array([[float(x) for x in r] for r in rows])
        if p.shape[1] != 4:
            continue
        out[mid] = (name or mid, p / p.sum(1, keepdims=True))
    return out


def icw(p):
    ic = np.log2(np.clip(p, 1e-9, 1)) - np.log2(0.25)
    ic[ic < 0] = 0
    return p * ic


def align_pwm_to_cwm(cwm, prob, n, require_cover=None):
    """Place a PWM against an attribution profile. Returns (logo, is_rc, offset)."""
    L = len(prob)
    Pf, Pr = icw(prob), icw(rc_matrix(prob))
    best, bp, brc = -np.inf, 0, False
    for is_rc, P in [(False, Pf), (True, Pr)]:
        for pos in range(n - L + 1):
            if require_cover is not None and not (pos <= require_cover <= pos + L - 1):
                continue
            sc = float(np.sum(P * cwm[pos:pos + L]))
            if sc > best:
                best, bp, brc = sc, pos, is_rc
    use = rc_matrix(prob) if brc else prob
    al = np.zeros((n, 4))
    al[bp:bp + L] = np.clip(use, 1e-6, 1)
    df = pwm_to_ic_df(np.clip(al, 1e-6, 1))
    df.iloc[:bp] = 0.0
    df.iloc[bp + L:] = 0.0
    return df, brc, bp


def best_motif(cwm,motifs,var_idx,n):
    # NOTE: this ranks candidates by the third return value of
    # align_pwm_to_cwm, which is the placement OFFSET, not a match score. That
    # looks like a bug -- it selects the motif that happens to sit furthest
    # right rather than the best match. Preserved as-is so this panel keeps
    # reproducing the published figure; see the README before changing it.
    best=(-np.inf,None)
    for mid,(name,prob) in motifs.items():
        L=len(prob)
        if L<6 or L>n: continue
        _,_,sc=align_pwm_to_cwm(cwm,prob,n,require_cover=var_idx)
        if sc>best[0]: best=(sc,(mid,name,prob))
    return best[1]
def draw_logo(ax,df,ylim,vpos,n):
    logomaker.Logo(df,color_scheme=NT_COLORS,ax=ax,flip_below=True,baseline_width=0,font_name="DejaVu Sans")
    for p in ax.patches: p.set_zorder(3)
    ax.axvline(vpos,color="#888",lw=0.5,ls="--",zorder=1); ax.axhline(0,color="#bbb",lw=0.4,zorder=1)
    ax.set_xlim(-0.5,n-0.5); ax.set_ylim(*ylim); ax.set_xticks([])
    ax.tick_params(axis="y",labelsize=5,length=2); ax.spines[["top","right","bottom"]].set_visible(False)
def draw_pwm(ax,df,vpos,n,rc,label):
    logomaker.Logo(df,color_scheme=NT_COLORS,ax=ax,baseline_width=0,font_name="DejaVu Sans")
    ax.axvline(vpos,color="#888",lw=0.5,ls="--",zorder=1)
    ax.set_xlim(-0.5,n-0.5); ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top","right","left","bottom"]].set_visible(False)
    ax.text(0.5,-0.12,f"{label}{' (−)' if rc else ''}",transform=ax.transAxes,ha="center",va="top",
            fontsize=5,style="italic",color="#333")
def gene_exons(name):
    ex=[]; gs=ge=None; strand="+"
    for line in open(GTF_FILE):
        if line.startswith("#"): continue
        f=line.rstrip().split("\t")
        if len(f)<9: continue
        m=re.search(r'gene_name "([^"]+)"',f[8])
        if not m or m.group(1)!=name: continue
        s,e=int(f[3]),int(f[4]); strand=f[6]
        if f[2]=="gene": gs,ge=s,e
        elif f[2]=="exon": ex.append((s,e))
    return gs,ge,strand,sorted(set(ex))
def atac_ism_logo(seqs_oh,ls,tidx):
    L=seqs_oh.shape[1]; ohT=seqs_oh.T; v=ls[:,:,tidx].astype(float); a=np.zeros((L,4))
    for p in range(L):
        i=int(np.argmax(ohT[p])); o=[j for j in range(4) if j!=i]; a[p,i]=v[p,i]-np.mean(v[p,o])
    return pd.DataFrame(a,columns=list("ACGT"))
def rna_ism_logo(seqs_oh,ls,tidx):
    L=seqs_oh.shape[1]; ohT=seqs_oh.T; v=ls[:,:,tidx].astype(float); a=np.zeros((L,4))
    for p in range(L):
        i=int(np.argmax(ohT[p])); o=[j for j in range(4) if j!=i]; a[p,i]=-float(np.mean(v[p,o]))
    return pd.DataFrame(a,columns=list("ACGT"))
def find_peak(ct,chrom,pos):
    path=MACS2.format(ct=ct)
    if not os.path.exists(path): return None
    for line in open(path):
        f=line.split("\t")
        if f[0]==chrom and int(f[1])<=pos<=int(f[2]): return (int(f[1]),int(f[2]))
    return None

ALL_MOTIFS=load_all_motifs()

# ════════════════════════════════════════════════════════════════════════════
def build(cfg):
    rs=cfg["rs"]; ct=cfg["ct"]; chrom=cfg["chrom"]; pos=int(cfg["pos"])
    ref=cfg["ref"]; alt=cfg["alt"]; egene=cfg["egene"]; gene_id=cfg["gene_id"]
    strand=cfg["strand"]; vk=f"{chrom}_{pos}_{ref}_{alt}"
    coloc_cs=int(cfg["coloc_cs"]); role=cfg.get("role","lead"); lead_snp=cfg.get("lead_snp",vk)
    atac_only=cfg.get("atac_only",False)
    ld=glob.glob(f"{SUS}/{rs}_*")[0]
    print(f"\n=== {rs} {ct} {egene} {chrom}:{pos} {ref}>{alt} ({strand}) cs={coloc_cs} role={role} ===")

    # GWAS
    ss=pd.read_csv(f"{ld}/sumstats.tsv",sep="\t"); pipdf=pd.read_csv(f"{ld}/pip.tsv",sep="\t")
    cnum=int(chrom.replace("chr","")); lo,hi=pos-WIDE_FLANK,pos+WIDE_FLANK
    ssw=ss[(ss.hg38_chr==cnum)&ss.hg38_bp.between(lo,hi)].copy()
    ssw["nlp"]=-np.log10(ssw.P.clip(lower=1e-300))
    pw=pipdf[(pipdf.hg38_chr==cnum)&pipdf.hg38_bp.between(lo,hi)].copy()
    # credible set(s) to show (default: the single colocalizing set)
    show_cs=cfg.get("show_cs",[coloc_cs])
    cspos=set(pw.loc[pw.cs.isin(show_cs),"hg38_bp"])

    # CS model scores (restricted to the colocalizing credible set)
    csb=pd.read_csv(f"{ROOT}/{rs}/cs_borzoi_scores.tsv",sep="\t")
    csc=pd.read_csv(f"{ROOT}/{rs}/cs_chrombpnet_scores.tsv",sep="\t")
    css=csb.merge(csc[["SNP","chrombpnet_log2fc"]],on="SNP")
    css=css[css.pos.isin(cspos)].sort_values("pos")
    _pad=max(2500,(css.pos.max()-css.pos.min())*0.10)   # side white-space so dots aren't clipped
    CS_LO,CS_HI=css.pos.min()-_pad,css.pos.max()+_pad

    # Cerberus pred
    bp=np.load(f"{ROOT}/borzoi_pred/{vk}.npz",allow_pickle=True)
    bs=int(bp["bins_start"]); bz=int(bp["bin_size"]); nb=len(bp["atac_ref"])
    bc=bs+(np.arange(nb)+0.5)*bz
    atac_ref,atac_alt=bp["atac_ref"],bp["atac_alt"]; rna_ref,rna_alt=bp["rna_ref"],bp["rna_alt"]
    alo,ahi=pos-ATAC_PRED_FLANK,pos+ATAC_PRED_FLANK
    if not atac_only:
        gs,ge,_,ex=gene_exons(egene)
        rlo=max(min(pos,gs)-5000, bc[0]); rhi=min(max(pos,ge)+5000, bc[-1])
        rmask=(bc>=rlo)&(bc<=rhi)

    # ChromBPNet pred
    cp=np.load(f"{ROOT}/chrombpnet/pred/{vk}.npz",allow_pickle=True)
    co=int(cp["out_start"][0]); cl=int(cp["out_len"][0]); cx=co+np.arange(cl)+0.5
    cbp_ref,cbp_alt=cp["ref"],cp["alt"]; clo,chi=pos-CBP_COV_FLANK,pos+CBP_COV_FLANK
    peak=find_peak(ct,chrom,pos)

    # ATAC ISM (panel_h track)
    with h5py.File(f"{ROOT}/atac_ism/scores.h5") as h:
        labels=list(h["label"][:].astype(str)); si=labels.index(vk)
        a_oh=h["ref/seqs"][si]; aa_oh=h["alt/seqs"][si]
        a_ls=h["ref/cov/logSUM"][si]; aa_ls=h["alt/cov/logSUM"][si]
    aidx=int(cfg["atac_ism_ph"])
    atac_lr=atac_ism_logo(a_oh,a_ls,aidx); atac_la=atac_ism_logo(aa_oh,aa_ls,aidx)
    avar=int(np.where((a_oh!=aa_oh).any(0))[0][0])
    am=max(atac_lr.values.max(),atac_la.values.max()); an=min(atac_lr.values.min(),atac_la.values.min())
    aylim=(an-abs(am-an)*.05,am+abs(am-an)*.05)

    # RNA ISM (lookup CT RNA+/- track in targets_cov order)
    if not atac_only:
        rna_dir=cfg.get("rna_ism_dir","rna_ism")
        tcov=pd.read_csv(f"{ROOT}/{rna_dir}/targets_cov.txt",sep="\t",index_col=0)
        # cov output is strand-collapsed to RNA+ slots; covgene writes sense-strand
        # logSED there for both + and - genes, so always use RNA+:{ct}.
        want=f"RNA+:{ct}:"
        ridx=[i for i,d in enumerate(tcov.description) if d.startswith(want)][0]
        with h5py.File(f"{ROOT}/{rna_dir}/scores.h5") as h:
            rlabels=list(h["label"][:].astype(str))
            gids=h["gene_ids"][:].astype(str)        # sorted unique gene list
            snp_idx=h["snp_idx"][:]; gene_idx=h["gene_idx"][:]; rsi=rlabels.index(vk)
            gbase=np.array([g.split(".")[0] for g in gids]); gtarget=gene_id.split(".")[0]
            # pair p: snp = snp_idx[p], gene = gene_ids[gene_idx[p]]
            pidx=[p for p in range(len(snp_idx))
                  if int(snp_idx[p])==rsi and gbase[int(gene_idx[p])]==gtarget][0]
            r_oh=h["ref/seqs"][rsi]; ra_oh=h["alt/seqs"][rsi]
            r_ls=h["ref/covgene/logSED"][pidx]; ra_ls=h["alt/covgene/logSED"][pidx]
        rna_lr=rna_ism_logo(r_oh,r_ls,ridx); rna_la=rna_ism_logo(ra_oh,ra_ls,ridx)
        rvar=int(np.where((r_oh!=ra_oh).any(0))[0][0])
        rm=max(rna_lr.values.max(),rna_la.values.max()); rn=min(rna_lr.values.min(),rna_la.values.min())
        rylim=(rn-abs(rm-rn)*.05,rm+abs(rm-rn)*.05)

    # ChromBPNet SHAP
    with h5py.File(f"{ROOT}/chrombpnet/ism_contribs_v4/{ct}_ism.h5") as h:
        ids=list(h["variant_ids"][:].astype(str)); ix=ids.index(vk)
        s_oh=h["ref_onehot"][ix]; sa_oh=h["alt_onehot"][ix]
        s_ct=h["ref_contrib"][ix]; sa_ct=h["alt_contrib"][ix]
    shap_lr=pd.DataFrame(s_oh*s_ct,columns=list("ACGT")); shap_la=pd.DataFrame(sa_oh*sa_ct,columns=list("ACGT"))
    svar=int(np.where((s_oh!=sa_oh).any(1))[0][0])
    sm=max(shap_lr.values.max(),shap_la.values.max()); sn=min(shap_lr.values.min(),shap_la.values.min())
    sylim=(sn-abs(sm-sn)*.05,sm+abs(sm-sn)*.05)

    # motif: fixed override (cfg["motif_id"], optionally from cfg["motif_meme"])
    # or auto-selected from alt ATAC ISM CWM
    if cfg.get("motif_id"):
        mid=cfg["motif_id"]; mname=cfg.get("motif_name",mid)
        if cfg.get("motif_meme"):
            _t=open(cfg["motif_meme"]).read()
            _blk=re.search(rf"MOTIF {re.escape(mid)}\b.*?\n(.*?)(?=\nMOTIF|\Z)",_t,re.DOTALL).group(0)
            _rows=[l.split() for l in _blk.splitlines() if re.match(r"^\s*[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+",l)]
            _m=np.array([[float(x) for x in r] for r in _rows]); mprob=_m/_m.sum(1,keepdims=True)
        else:
            mprob=ALL_MOTIFS[mid][1]
        print(f"  fixed motif: {mid} ({mname})")
    else:
        mid,mname,mprob=best_motif(atac_la.values,ALL_MOTIFS,avar,50)
        print(f"  auto motif: {mid} {mname}")
    # align each panel's PWM to its alt CWM (covering the variant) + optional shift
    _,atac_rc,atac_bp=align_pwm_to_cwm(atac_la.values,mprob,50,require_cover=avar)
    _,shap_rc,shap_bp=align_pwm_to_cwm(shap_la.values,mprob,50,require_cover=svar)
    atac_bp+=cfg.get("motif_shift_atac",0); shap_bp+=cfg.get("motif_shift_shap",0)
    atac_pwm=place_pwm(mprob,atac_rc,atac_bp,50)
    shap_pwm=place_pwm(mprob,shap_rc,shap_bp,50)
    def exon_sum(ref,alt,xf,exons):
        r=s=0.0
        for es,ee in exons:
            mk=(xf>=es)&(xf<=ee); r+=ref[mk].sum(); s+=alt[mk].sum()
        return r,s
    if not atac_only:
        # RNA motif: same genomic position/orientation as the Cerberus-ATAC motif
        # (RNA window is 100 bp with variant at rvar; ATAC is 50 bp at avar)
        rna_rc=atac_rc
        rna_bp=atac_bp+(rvar-avar)+cfg.get("motif_shift_rna",0)
        rna_pwm=place_pwm(mprob,rna_rc,rna_bp,100)
        rsum,asum=exon_sum(rna_ref,rna_alt,bc,ex)

    # ── layout ────────────────────────────────────────────────────────────────
    FIG_W=4.2; MTOP=0.20; MBOT=0.34
    FULL_L,FULL_W=0.15,0.82; COL_L_L,COL_W=0.155,0.355; COL_L_R=0.585
    if atac_only:
        CS3_L=[0.155,0.585]; CS3_W=0.355   # 2 model columns
        ROWS=[("manh",0.41,"full",0.05),("pip",0.45,"full",0.34),("csprior",0.31,"cs3",0.30),
              ("atac_cov",0.57,"atac2",0.05),("atac_ref",0.42,"atac2",0.03),("atac_alt",0.42,"atac2",0.03),
              ("atac_pwm",0.18,"atac2",0.0)]
    else:
        CS3_L=[0.155,0.445,0.735]; CS3_W=0.235
        ROWS=[("manh",0.41,"full",0.05),("pip",0.45,"full",0.34),("csprior",0.31,"cs3",0.30),
              ("atac_cov",0.57,"atac2",0.05),("atac_ref",0.42,"atac2",0.03),("atac_alt",0.42,"atac2",0.03),
              ("atac_pwm",0.18,"atac2",0.34),("rna_pred",0.80,"full",0.02),("gene",0.28,"full",0.30),
              ("rna_ism_ref",0.42,"full",0.03),("rna_ism_alt",0.42,"full",0.03),("rna_pwm",0.18,"full",0.0)]
    FIG_H=MTOP+MBOT+sum(h for _,h,_,_ in ROWS)+sum(g for *_,g in ROWS)
    fig=plt.figure(figsize=(FIG_W,FIG_H)); A={}; yc=FIG_H-MTOP
    for k,h,mode,gap in ROWS:
        bot=yc-h
        if mode=="full": A[k]=fig.add_axes([FULL_L,bot/FIG_H,FULL_W,h/FIG_H])
        elif mode=="cs3":
            for j in range(len(CS3_L)): A[f"{k}_{j}"]=fig.add_axes([CS3_L[j],bot/FIG_H,CS3_W,h/FIG_H])
        else:
            A[k+"_l"]=fig.add_axes([COL_L_L,bot/FIG_H,COL_W,h/FIG_H]); A[k+"_r"]=fig.add_axes([COL_L_R,bot/FIG_H,COL_W,h/FIG_H])
        yc=bot-gap
    ct0=A["atac_cov_l"].get_position().y1
    fig.text(COL_L_L+COL_W/2,ct0+0.004,"ChromBPNet",ha="center",va="bottom",fontsize=6.5,color="#333")
    fig.text(COL_L_R+COL_W/2,ct0+0.004,"Cerberus",ha="center",va="bottom",fontsize=6.5,color="#333")

    lead=pos
    # Manhattan
    ax=A["manh"]; ncs=ssw[~ssw.hg38_bp.isin(cspos)]; ocs=ssw[ssw.hg38_bp.isin(cspos)&(ssw.hg38_bp!=lead)]; foc=ssw[ssw.hg38_bp==lead]
    ax.scatter(ncs.hg38_bp,ncs.nlp,s=3,color="#9a9a9a",alpha=.6,rasterized=True,linewidths=0)
    ax.scatter(ocs.hg38_bp,ocs.nlp,s=9,color=CS_COLOR,alpha=.95,zorder=4,edgecolors="white",linewidths=.3)
    ax.scatter(foc.hg38_bp,foc.nlp,s=9,color=CS_COLOR,edgecolors="black",linewidths=.7,zorder=6)
    ax.axhline(-np.log10(5e-8),color="#d62728",lw=.7,ls="--"); ax.axvline(lead,**GRAY_DASH)
    ax.set_ylabel("−log$_{10}$(P)",fontsize=6.5); ax.set_ylim(0,ssw.nlp.max()*1.08); ax.set_xlim(lo,hi)
    ax.tick_params(labelsize=6,length=2); ax.set_xticks([]); ax.spines[["top","right","bottom"]].set_visible(False)
    ax.text(0.98,0.92,"eGFR GWAS",ha="right",va="top",fontsize=6,transform=ax.transAxes,color="#333")
    role_lbl={"pipLead":"highest-PIP variant","modelLead":"highest model-score variant",
              "pip+model":"highest-PIP = highest model-score"}.get(role,role)
    ax.text(0.02,0.92,f"cs{coloc_cs}: {role_lbl}",ha="left",va="top",fontsize=5.5,transform=ax.transAxes,color="#555")
    # PIP
    ax=A["pip"]; ncp=pw[pw.cs==-1]; ocp=pw[(pw.cs!=-1)&(pw.hg38_bp!=lead)]; fp=pw[pw.hg38_bp==lead]
    ax.scatter(ncp.hg38_bp,ncp.pip,s=3,color="#9a9a9a",alpha=.5,rasterized=True,linewidths=0)
    ax.scatter(ocp.hg38_bp,ocp.pip,s=9,color=CS_COLOR,alpha=.95,zorder=4,edgecolors="white",linewidths=.3)
    ax.scatter(fp.hg38_bp,fp.pip,s=9,color=CS_COLOR,edgecolors="black",linewidths=.7,zorder=6)
    ax.axvline(lead,**GRAY_DASH); ax.set_ylabel("PIP",fontsize=6.5); ax.set_ylim(0,1); ax.set_yticks([0,.5,1]); ax.set_xlim(lo,hi)
    ax.tick_params(labelsize=6,length=2)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x/1e6:.2f}")); ax.tick_params(axis="x",labelsize=6,length=2)
    ax.spines[["top","right"]].set_visible(False)
    # CS prioritization
    cs_models=[("csprior_0","chrombpnet_log2fc","ChromBPNet","log$_2$FC"),
               ("csprior_1","borzoi_atac_logSUM_eff","Cerberus ATAC","ΔlogSUM")]
    if not atac_only:
        cs_models.append(("csprior_2","borzoi_rna_logSED_eff","Cerberus RNA","ΔlogSED"))
    for k,col,mn,met in cs_models:
        ax=A[k]; vals=css[col].values; ax.axhline(0,color="#ccc",lw=.4,zorder=1); ax.axvline(lead,**GRAY_DASH)
        for _,r in css.iterrows():
            ld_=(r["pos"]==lead)
            ax.scatter([r["pos"]],[r[col]],s=16 if ld_ else 9,color=CS_COLOR,
                       edgecolors="black" if ld_ else "white",linewidths=.7 if ld_ else .3,zorder=3)
        ax.set_xlim(CS_LO,CS_HI); l0,h0=min(0,vals.min()),max(0,vals.max()); sp=h0-l0 or 1
        ax.set_ylim(l0-sp*.08,h0+sp*.16); ax.set_title(mn,fontsize=6,pad=2); ax.set_ylabel(met,fontsize=5,labelpad=1)
        ax.tick_params(axis="y",labelsize=4.5,width=.4,length=2,pad=1); ax.yaxis.set_major_locator(MaxNLocator(3,prune="both"))
        ax.set_xticks([CS_LO+(CS_HI-CS_LO)*.25,CS_LO+(CS_HI-CS_LO)*.75])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x/1e6:.3f}")); ax.tick_params(axis="x",labelsize=4.3,length=2)
        ax.spines[["top","right"]].set_visible(False); ax.spines["left"].set_linewidth(.4); ax.spines["bottom"].set_linewidth(.4)
    cs_lbl="cs"+"+".join(str(c) for c in show_cs)
    fig.text(CS3_L[0],A["csprior_0"].get_position().y1+0.013,
             f"Predicted effect across the {len(css)} {cs_lbl} variants  ({lead_snp} PIP={cfg.get('lead_pip',0):.2f} outlined)",
             ha="left",va="bottom",fontsize=5.2,color="#555")
    # ATAC coverage panel-G style
    leg=[Line2D([0],[0],color=REF_COLOR,lw=.8,label=f"ref ({ref})"),Line2D([0],[0],color=ALT_COLOR,lw=.8,label=f"alt ({alt})")]
    def cov(ax,x,rf,al,xl,ylabel=None,pk=False):
        if pk and peak: ax.axvspan(peak[0],peak[1],color=PEAK_SHADE,alpha=.8,zorder=0,lw=0)
        ax.axvline(lead,color="#888",lw=.5,ls="--",zorder=2)
        mk=(x>=xl[0])&(x<=xl[1]); xp,rp,ap=x[mk],rf[mk],al[mk]
        if rp.sum()<=ap.sum(): ax.plot(xp,ap,color=ALT_COLOR,lw=.4,zorder=3); ax.plot(xp,rp,color=REF_COLOR,lw=.4,zorder=4)
        else: ax.plot(xp,rp,color=REF_COLOR,lw=.4,zorder=3); ax.plot(xp,ap,color=ALT_COLOR,lw=.4,zorder=4)
        ax.set_xlim(*xl); ax.set_ylim(bottom=0); ax.set_xticks([])
        ax.spines[["top","right","bottom"]].set_visible(False); ax.spines["left"].set_linewidth(.4)
        ax.tick_params(axis="y",labelsize=4.5,width=.4,length=2,pad=1); ax.yaxis.set_major_locator(MaxNLocator(2,prune="both"))
        if ylabel: ax.set_ylabel(ylabel,fontsize=5.5,rotation=0,ha="right",va="center",labelpad=2)
        ax.legend(handles=leg,fontsize=4.5,loc="upper right",frameon=False,handlelength=.8,handleheight=.55,borderpad=0,labelspacing=.2)
    cov(A["atac_cov_l"],cx,cbp_ref,cbp_alt,(clo,chi),ylabel=f"pred\n({ct} ATAC)",pk=True)
    cov(A["atac_cov_r"],bc,atac_ref,atac_alt,(alo,ahi),pk=True)
    if peak: A["atac_cov_l"].text((peak[0]+peak[1])/2,1.0,"peak",transform=A["atac_cov_l"].get_xaxis_transform(),ha="center",va="bottom",fontsize=4.5,color="#777",clip_on=False)
    draw_logo(A["atac_ref_l"],shap_lr,sylim,svar,50); A["atac_ref_l"].set_ylabel(f"contrib\nref ({ref})",fontsize=5.5,rotation=0,ha="right",va="center",labelpad=2)
    draw_logo(A["atac_alt_l"],shap_la,sylim,svar,50); A["atac_alt_l"].set_ylabel(f"contrib\nalt ({alt})",fontsize=5.5,rotation=0,ha="right",va="center",labelpad=2)
    draw_logo(A["atac_ref_r"],atac_lr,aylim,avar,50); draw_logo(A["atac_alt_r"],atac_la,aylim,avar,50)
    draw_pwm(A["atac_pwm_l"],shap_pwm,svar,50,shap_rc,mname); A["atac_pwm_l"].set_ylabel("motif",fontsize=5.5,rotation=0,ha="right",va="center",labelpad=2)
    draw_pwm(A["atac_pwm_r"],atac_pwm,avar,50,atac_rc,mname)
    if not atac_only:
        # RNA pred panel-I style
        ax=A["rna_pred"]; xr,rr,ar2=bc[rmask],rna_ref[rmask],rna_alt[rmask]; ytop=max(rr.max(),ar2.max())*1.15 or 1
        if rr.sum()>=ar2.sum(): ax.plot(xr,rr,color=REF_COLOR,lw=.5,zorder=3); ax.plot(xr,ar2,color=ALT_COLOR,lw=.5,zorder=4)
        else: ax.plot(xr,ar2,color=ALT_COLOR,lw=.5,zorder=3); ax.plot(xr,rr,color=REF_COLOR,lw=.5,zorder=4)
        ax.legend(handles=[Line2D([0],[0],color=REF_COLOR,lw=.8,label=f"ref ({ref})  Σ={rsum:.1f}"),
                           Line2D([0],[0],color=ALT_COLOR,lw=.8,label=f"alt ({alt})  Σ={asum:.1f}")],
                  fontsize=3.8,loc="upper right",frameon=False,handlelength=.8,handleheight=.55,borderpad=0,labelspacing=.2)
        ax.axvline(lead,color="#888",lw=.7,ls="--",zorder=7)
        ax.text(0.01,0.96,f"{chrom}:{pos:,} {ref}→{alt}",transform=ax.transAxes,ha="left",va="top",fontsize=3.8,color="#888")
        ax.set_xlim(rlo,rhi); ax.set_ylim(0,ytop); ax.set_xticks([]); ym=max(rr.max(),ar2.max())
        ax.set_yticks([0,ym]); ax.set_yticklabels(["0",f"{ym:.1f}"]); ax.tick_params(axis="y",labelsize=4.5,width=.4,length=2,pad=1)
        ax.text(0.01,0.03,f"{chrom}:{int(rlo):,}–{int(rhi):,}",transform=ax.transAxes,ha="left",va="bottom",fontsize=3.8,color="#888")
        ax.spines[["top","right","bottom"]].set_visible(False); ax.spines["left"].set_linewidth(.4)
        ax.set_ylabel(f"Cerberus pred\n({ct} RNA{'+' if strand=='+' else '−'})",fontsize=5.5,rotation=0,ha="right",va="center",labelpad=2)
        # gene annotation
        ax=A["gene"]; ax.set_xlim(rlo,rhi); ax.set_ylim(0,1)
        s,e=max(gs,rlo),min(ge,rhi); ax.plot([s,e],[.5,.5],color="#333",lw=1,solid_capstyle="butt")
        for es,ee in ex:
            if ee<rlo or es>rhi: continue
            ax.add_patch(Rectangle((max(es,rlo),.5-.19),min(ee,rhi)-max(es,rlo),.38,color="#333",zorder=3))
        # strand arrows
        span=rhi-rlo; step=span/8
        for ap in np.arange(s+step*0.5,e,step):
            dx=step*0.12*(1 if strand=="+" else -1)
            ax.annotate("",xy=(ap+dx,.5),xytext=(ap-dx,.5),arrowprops=dict(arrowstyle="-|>",color="#333",lw=.4,mutation_scale=3),zorder=4)
        ax.text(s+span*.01,.80,f"{egene} ({strand})",fontsize=5.5,va="bottom",color="#333"); ax.axvline(lead,**GRAY_DASH)
        ax.spines[["top","right","left"]].set_visible(False); ax.set_yticks([])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f"{x/1e6:.3f}")); ax.tick_params(axis="x",labelsize=5.5,length=2)
        # RNA ISM + motif
        draw_logo(A["rna_ism_ref"],rna_lr,rylim,rvar,100); A["rna_ism_ref"].set_ylabel(f"Ref ({ref})",fontsize=5.5,rotation=0,ha="right",va="center",labelpad=2)
        draw_logo(A["rna_ism_alt"],rna_la,rylim,rvar,100); A["rna_ism_alt"].set_ylabel(f"Alt ({alt})",fontsize=5.5,rotation=0,ha="right",va="center",labelpad=2)
        A["rna_ism_ref"].set_title(f"RNA ISM / {egene} exons ({ct})",fontsize=6,pad=2,loc="left")
        draw_pwm(A["rna_pwm"],rna_pwm,rvar,100,rna_rc,mname); A["rna_pwm"].set_ylabel("motif",fontsize=5.5,rotation=0,ha="right",va="center",labelpad=2)

    out=f"{OUTDIR}/{rs}_{egene}_{ct}_{role}_{lead_snp}.pdf"
    fig.savefig(out,dpi=200); fig.savefig(out.replace(".pdf",".png"),dpi=200)
    plt.close(fig); print(f"  saved {out}")

if __name__=="__main__":
    SPECS=json.load(open(f"{ROOT}/figure_specs.json"))
    targets=[c for c in SPECS if (len(sys.argv)<2 or c["rs"]==sys.argv[1])]
    for c in targets: build(c)
