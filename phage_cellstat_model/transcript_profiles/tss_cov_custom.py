import argparse, sys, re, json
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np, pandas as pd, pysam

def parse_attrs(s: str) -> Dict[str,str]:
    d={}
    for part in s.strip().strip(";").split(";"):
        part=part.strip()
        if not part: continue
        if "=" in part:
            k,v = part.split("=",1); d[k]=v
        else:
            m=re.match(r'(\S+)\s+"([^"]+)"', part)
            if m: d[m.group(1)]=m.group(2)
    return d

def load_gene_like_records(gtf: str, feature: str, keys: List[str]) -> Dict[str,dict]:
    """
    Returns a dict id->record with fields: chrom, start, end, strand.
    If feature == 'CDS', we collapse all CDS for a gene into one span and use
    the extreme end as TSS/TES proxies (useful when no 'gene' rows exist).
    """
    rows=[]
    with open(gtf) as f:
        for ln in f:
            if not ln or ln.startswith("#"): continue
            p=ln.rstrip("\n").split("\t")
            if len(p)!=9: continue
            seqid,src,ft,start,end,score,strand,phase,attrs=p
            if feature!="CDS":
                if ft!=feature: continue
                a=parse_attrs(attrs)
                rows.append((seqid,int(start)-1,int(end),strand,a,ft))
            else:
                if ft!="CDS": continue
                a=parse_attrs(attrs)
                rows.append((seqid,int(start)-1,int(end),strand,a,ft))

    idx={}
    if feature!="CDS":
        for chrom,start,end,strand,a,ft in rows:
            gid=None
            for k in keys:
                if k in a: gid=a[k]; break
            if gid is None:
                gid=a.get("ID")
            if gid is None:
                continue
            idx[gid]={"chrom":chrom,"start":start,"end":end,"strand":strand,"attrs":a}
    else:
        by_gene={}
        for chrom,start,end,strand,a,ft in rows:
            gid=None
            for k in ["gene_id","gene_name","Parent"]+keys:
                if k in a: gid=a[k]; break
            if gid is None:
                gid=a.get("transcript_id")
            if gid is None:
                continue
            by_gene.setdefault((gid,chrom,strand), []).append((start,end))
        for (gid,chrom,strand), spans in by_gene.items():
            starts=[s for s,e in spans]; ends=[e for s,e in spans]
            start=min(starts); end=max(ends)
            idx[gid]={"chrom":chrom,"start":start,"end":end,"strand":strand,"attrs":{"gene_like":gid}}
    return idx

def rpm_scale(bam: pysam.AlignmentFile) -> float:
    return 1.0 / max(bam.mapped,1) * 1e6

def read_is_sense(read: pysam.AlignedSegment, gene_strand: str, lib: str) -> bool:
    if lib=="unstranded":
        return True
    frag_on_plus = not read.is_reverse
    sense_on_plus = (lib=="fr-secondstrand")  # fr-secondstrand keeps sense on plus
    sense_gene_plus = (gene_strand=="+")
    return frag_on_plus == sense_gene_plus if sense_on_plus else frag_on_plus != sense_gene_plus

def count_cov(bam: pysam.AlignmentFile, chrom: str, start: int, end: int,
              sense_only: bool, gene_strand: str, lib: str) -> np.ndarray:
    if not sense_only:
        A=bam.count_coverage(chrom,start,end,quality_threshold=0)
        return np.sum(np.vstack(A),axis=0).astype(float)
    def rc(read): return read_is_sense(read,gene_strand,lib)
    A=bam.count_coverage(chrom,start,end,quality_threshold=0,read_callback=rc)
    return np.sum(np.vstack(A),axis=0).astype(float)

def bin_pair(x: np.ndarray, y: np.ndarray, binsize: int) -> Tuple[np.ndarray,np.ndarray]:
    if binsize<=1: return x,y
    n=len(y); pad=(-n)%binsize
    if pad:
        y=np.pad(y,(0,pad),mode="edge")
        x=np.pad(x,(0,pad),mode="edge")
    yb=y.reshape(-1,binsize).mean(axis=1)
    xb=x.reshape(-1,binsize).mean(axis=1).astype(int)
    return xb,yb

def apply_log(y: np.ndarray, mode: str, pseudocount: float) -> np.ndarray:
    if mode=="none":
        return y
    pc = float(pseudocount)
    if mode=="log2":
        return np.log2(y + pc)
    if mode=="log10":
        return np.log10(y + pc)
    if mode=="ln":
        return np.log(y + pc)
    return y  # fallback

def main():
    ap=argparse.ArgumentParser(description="Genome coverage around custom GTF genes/transcripts/CDS.")
    ap.add_argument("--bam",required=True,help="BAM on Drive")
    ap.add_argument("--gtf",required=True,help="Custom GTF with recombinant features")
    ap.add_argument("--genes",nargs="+",help="Gene IDs/names to profile")
    ap.add_argument("--genes-file",help="File with one gene ID/name per line")
    ap.add_argument("--feature",choices=["gene","transcript","CDS"],default="gene",
                    help="Which feature to anchor on (CDS collapses all CDS per gene).")
    ap.add_argument("--anchor",choices=["tss","tes","cds_start","cds_end"],default="tss",
                    help="Anchor point in transcript sense; cds_* require --feature=CDS or transcripts with CDS extents.")
    ap.add_argument("--upstream",type=int,default=1000)
    ap.add_argument("--downstream",type=int,default=1000)
    ap.add_argument("--bin",type=int,default=10)
    ap.add_argument("--norm",choices=["none","rpm"],default="rpm")
    ap.add_argument("--strand",choices=["unstranded","fr-firststrand","fr-secondstrand"],default="unstranded")
    ap.add_argument("--attrs",nargs="+",default=["gene_name","gene_id","Name","locus_tag","ID","RecombinantID"],
                    help="Attribute keys to match your custom IDs (order matters).")
    ap.add_argument("--outdir",default="profiles")
    ap.add_argument("--plot",action="store_true")

    ap.add_argument("--log", choices=["none","log2","log10","ln"], default="none",
                    help="Apply log transform to coverage; adds a second log column to TSV and uses it for plotting.")
    ap.add_argument("--pseudocount", type=float, default=1.0,
                    help="Pseudocount added before log to avoid log(0). Default 1.0.")
    args=ap.parse_args()

    if not args.genes and not args.genes_file:
        sys.exit("Provide --genes or --genes-file")
    gene_list = args.genes or [ln.strip() for ln in open(args.genes_file) if ln.strip()]

    bam=pysam.AlignmentFile(args.bam,"rb")
    ref_len={r:l for r,l in zip(bam.references,bam.lengths)}
    scale = rpm_scale(bam) if args.norm=="rpm" else 1.0

    feats = load_gene_like_records(args.gtf, feature=args.feature, keys=args.attrs)
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True)

    manifest=[]; missing=[]; bad_chr=[]
    for g in gene_list:
        if g not in feats:
            missing.append(g); continue
        rec=feats[g]; chrom=rec["chrom"]; strand=rec["strand"]
        if chrom not in ref_len:
            bad_chr.append((g,chrom)); continue

        # anchors
        if args.anchor in ["tss","tes"]:
            tss = rec["start"] if strand=="+" else rec["end"]-1
            tes = rec["end"]-1 if strand=="+" else rec["start"]
            anchor_pos = tss if args.anchor=="tss" else tes
        else:
            cds_start = rec["start"] if strand=="+" else rec["end"]-1
            cds_end   = rec["end"]-1 if strand=="+" else rec["start"]
            anchor_pos = cds_start if args.anchor=="cds_start" else cds_end

        # window (clamped to contig)
        s = max(0, anchor_pos - args.upstream)
        e_req = anchor_pos + args.downstream + 1
        L=ref_len[chrom]; e=min(e_req, L)
        if e<e_req:
            print(f"[INFO] {g}: window clamped at {chrom} length {L}", file=sys.stderr)

        y = count_cov(bam, chrom, s, e, sense_only=(args.strand!="unstranded"),
                      gene_strand=strand, lib=args.strand) * scale
        n=len(y); x=(np.arange(n)+s)-anchor_pos

        x,y = bin_pair(x,y,args.bin)

        # flip so “downstream” is + for both strands
        if strand=="-":
            x = -x[::-1]; y = y[::-1]

        # Prepare DataFrame with raw and (optional) log coverage
        raw_col = f"coverage_{args.norm}"
        log_suffix = "" if args.log=="none" else f"_{args.log}"
        df_dict = {"rel_bp": x, raw_col: np.round(y, 6)}
        if args.log != "none":
            y_log = apply_log(y, args.log, args.pseudocount)
            df_dict[raw_col + log_suffix] = np.round(y_log, 6)

        df = pd.DataFrame(df_dict)

        # filename reflects log choice
        tsv=outdir/f"{g}_{args.feature}_{args.anchor}_up{args.upstream}_down{args.downstream}_bin{args.bin}_{args.strand}_{args.norm}{('_'+args.log) if args.log!='none' else ''}.tsv"
        df.to_csv(tsv,sep="\t",index=False)

        if args.plot:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(8,3))
            if args.log == "none":
                plt.plot(df["rel_bp"].values, df[raw_col].values, linewidth=1.25, label="raw")
                ylabel = "Coverage" + (" (RPM)" if args.norm=="rpm" else "")
            else:
                plt.plot(df["rel_bp"].values, df[raw_col + log_suffix].values, linewidth=1.25, label=args.log)
                ylabel = f"log coverage ({args.log}, +{args.pseudocount} pc)" + (" of RPM" if args.norm=="rpm" else " of counts")
            plt.axvline(0, linestyle="--")
            plt.title(f"{g} {chrom}:{rec['start']+1}-{rec['end']} ({strand}), anchor={args.anchor.upper()}")
            plt.xlabel(f"Position vs {args.anchor.upper()} (bp)")
            plt.ylabel(ylabel)
            plt.tight_layout(); plt.savefig(outdir/f"{tsv.stem}.png", dpi=200); plt.close()

        manifest.append({"gene":g,"chrom":chrom,"strand":strand,"feature":args.feature,
                         "anchor":args.anchor,"upstream":args.upstream,"downstream":args.downstream,
                         "binsize":args.bin,"norm":args.norm,"libtype":args.strand,
                         "log":args.log,"pseudocount":args.pseudocount,"tsv":str(tsv)})

    with open(outdir/"manifest.json","w") as f: json.dump(manifest,f,indent=2)

    if missing:
        print(f"[WARN] {len(missing)} IDs not found in GTF: {missing[:8]}{'...' if len(missing)>8 else ''}", file=sys.stderr)
        print(f"       Try --attrs to include your custom key (e.g., --attrs RecombinantID gene_id gene_name).", file=sys.stderr)
    if bad_chr:
        print(f"[WARN] {len(bad_chr)} IDs on contigs not in BAM: {bad_chr[:5]}{'...' if len(bad_chr)>5 else ''}", file=sys.stderr)
        print("       Check that BAM was aligned to the same reference (including plasmids/inserts).", file=sys.stderr)

    print(f"Done. Wrote {len(manifest)} profiles to {outdir}")

if __name__=="__main__":
    main()
