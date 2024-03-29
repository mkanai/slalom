# coding: utf-8
import argparse
import os.path
import numpy as np
import scipy as sp
import pandas as pd
import hail as hl
from hail.linalg import BlockMatrix
from hail.utils import new_temp_file

gnomad_latest_versions = {"GRCh37": "2.1.1", "GRCh38": "3.1.2"}
gnomad_pops = {"GRCh37": ["afr", "amr", "eas", "fin", "nfe"], "GRCh38": ["afr", "amr", "eas", "fin", "nfe", "sas"]}
gnomad_ld_variant_indices = {
    "GRCh37": "gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.{pop}.common.adj.ld.variant_indices.ht",
    "GRCh38": "gs://finucane-requester-pays/slalom/gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.{pop}.common.adj.ld.variant_indices.b38.ht",
}


class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            if value.isnumeric():
                value = float(value)
            getattr(namespace, self.dest)[key] = value


# cf. https://github.com/armartin/prs_disparities/blob/master/run_prs_holdout.py
def flip_text(base):
    """
    :param StringExpression base: Expression of a single base
    :return: StringExpression of flipped base
    :rtype: StringExpression
    """
    return hl.switch(base).when("A", "T").when("T", "A").when("C", "G").when("G", "C").default(base)


def align_alleles(ht, ht_gnomad, flip_rows=None):
    ht = ht.annotate(
        **(
            hl.case()
            .when(
                hl.is_defined(ht_gnomad[ht.locus, hl.array([ht.alleles[0], ht.alleles[1]])]),
                hl.struct(alleles=[ht.alleles[0], ht.alleles[1]], flip_row=False),
            )
            .when(
                hl.is_defined(ht_gnomad[ht.locus, hl.array([ht.alleles[1], ht.alleles[0]])]),
                hl.struct(alleles=[ht.alleles[1], ht.alleles[0]], flip_row=True),
            )
            .when(
                hl.is_defined(ht_gnomad[ht.locus, hl.array([flip_text(ht.alleles[0]), flip_text(ht.alleles[1])])]),
                hl.struct(alleles=[flip_text(ht.alleles[0]), flip_text(ht.alleles[1])], flip_row=False),
            )
            .when(
                hl.is_defined(ht_gnomad[ht.locus, hl.array([flip_text(ht.alleles[1]), flip_text(ht.alleles[0])])]),
                hl.struct(alleles=[flip_text(ht.alleles[1]), flip_text(ht.alleles[0])], flip_row=True),
            )
            .default(hl.struct(alleles=[ht.alleles[0], ht.alleles[1]], flip_row=False))
        )
    )

    if flip_rows is not None:
        ht = ht.annotate(**{row: hl.if_else(ht.flip_row, -ht[row], ht[row]) for row in flip_rows})
        ht = ht.drop("flip_row")

    return ht


def get_diag_mat(diag_vec: BlockMatrix):
    x = diag_vec.T.to_numpy()
    diag_mat = np.identity(len(x)) * np.outer(np.ones(len(x)), x)
    return BlockMatrix.from_numpy(diag_mat)


def abf(beta, se, W=0.04):
    z = beta / se
    V = se**2
    r = W / (W + V)
    lbf = 0.5 * (np.log(1 - r) + (r * z**2))
    denom = sp.special.logsumexp(lbf)
    prob = np.exp(lbf - denom)
    return lbf, prob


def get_cs(variant, prob, coverage=0.95):
    ordering = np.argsort(prob)[::-1]
    idx = np.where(np.cumsum(prob[ordering]) > coverage)[0][0]
    cs = variant[ordering][: (idx + 1)]
    return cs


def main(args):
    hl._set_flags(no_whole_stage_codegen="1")
    reference_genome = args.reference_genome
    gnomad_version = gnomad_latest_versions[reference_genome]
    gnomad_ht_path = f"gs://finucane-requester-pays/slalom/gnomad/release/{gnomad_version}/ht/genomes/gnomad.genomes.r{gnomad_version}.sites.most_severe.ht"

    ht_snp = hl.import_table(args.snp, impute=True, types={"chromosome": hl.tstr}, delimiter="\s+")
    ht_snp = ht_snp.annotate(
        locus=hl.parse_locus(
            hl.delimit([ht_snp.chromosome, hl.str(ht_snp.position)], delimiter=":"), reference_genome=reference_genome
        ),
        alleles=[ht_snp.allele1, ht_snp.allele2],
    )
    if args.align_alleles:
        ht_gnomad = hl.read_table(gnomad_ht_path)
        ht_snp = align_alleles(ht_snp, ht_gnomad, flip_rows=["beta"])
        ht_snp = ht_snp.annotate(variant_aligned=hl.variant_str(ht_snp.locus, ht_snp.alleles))

    ht_snp = ht_snp.annotate(variant=hl.variant_str(ht_snp.locus, ht_snp.alleles))
    ht_snp = ht_snp.key_by("locus", "alleles")
    ht_snp = ht_snp.add_index("idx_snp")

    # annotate in novel CUPs and reject
    if args.annotate_cups:
        cup = hl.read_table(
            f"gs://finucane-requester-pays/slalom/cup_files/FASTA_BED.ALL_{reference_genome}.novel_CUPs.ht"
        )
        reject = hl.read_table(
            f"gs://finucane-requester-pays/slalom/cup_files/FASTA_BED.ALL_{reference_genome}.reject_2.ht"
        )
        ht_snp = ht_snp.annotate(in_cups=hl.is_defined(cup[ht_snp.locus]) | hl.is_defined(reject[ht_snp.locus]))

    # annotate vep and freq
    if args.annotate_consequence or args.annotate_gnomad_freq:
        ht_gnomad = hl.read_table(gnomad_ht_path)
        consequences = ["most_severe", "gene_most_severe", "consequence"] if args.annotate_consequence else []
        freq_expr = (
            {f"gnomad_v{gnomad_version[0]}_af_{pop}": ht_gnomad.freq[pop].AF for pop in gnomad_pops[reference_genome]}
            if args.annotate_gnomad_freq
            else {}
        )
        ht_gnomad = ht_gnomad.select(*consequences, **freq_expr)
        ht_snp = ht_snp.join(ht_gnomad, how="left")
    ht_snp = ht_snp.checkpoint(new_temp_file())

    t = new_temp_file()
    df = ht_snp.key_by().drop("locus", "alleles", "idx_snp").export(t)
    df = pd.read_csv(t, sep="\t", dtype={"chromosome": "string"})

    if args.abf:
        lbf, prob = abf(df.beta, df.se, W=args.abf_prior_variance)
        cs = get_cs(df.variant, prob, coverage=0.95)
        cs_99 = get_cs(df.variant, prob, coverage=0.99)
        df["lbf"] = lbf
        df["prob"] = prob
        df["cs"] = df.variant.isin(cs)
        df["cs_99"] = df.variant.isin(cs_99)

    if args.lead_variant is None:
        if args.lead_variant_choice == "p":
            lead_idx_snp = df.p.idxmin()
        elif args.lead_variant_choice == "prob":
            lead_idx_snp = df.prob.idxmax()
        elif args.lead_variant_choice in ["gamma", "gamma-p"]:
            lead_idx_snp = df.index[df.gamma]
            if len(lead_idx_snp) == 0:
                if args.lead_variant_choice == "gamma-p":
                    lead_idx_snp = df.p.idxmin()
                else:
                    raise ValueError("No lead variants found with gamma.")
            elif len(lead_idx_snp) > 1:
                raise ValueError("Multiple lead variants found with gamma.")
            else:
                lead_idx_snp = lead_idx_snp[0]
        args.lead_variant = df.variant[lead_idx_snp]
    else:
        lead_idx_snp = df.index[df.variant == args.lead_variant]

    df["lead_variant"] = False
    df.iloc[lead_idx_snp, df.columns.get_loc("lead_variant")] = True

    # annotate LD
    r2_label = "r2" if not args.export_r else "r"
    if args.ld_reference == "gnomad":
        ld_matrices = [
            f"gs://gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.{pop}.common.ld.bm"
            for pop in gnomad_pops["GRCh37"]
        ]
        ld_variant_indices = [
            gnomad_ld_variant_indices[reference_genome].format(pop=pop) for pop in gnomad_pops["GRCh37"]
        ]
        ld_labels = [f"gnomad_lead_{r2_label}_{pop}" for pop in gnomad_pops["GRCh37"]]
    else:
        ld_matrices = [args.custom_ld_path]
        ld_variant_indices = [args.custom_ld_variant_index_path]
        ld_labels = [f"{args.custom_ld_label}_lead_{r2_label}"]

    for ld_bm_path, ld_ht_path, col in zip(ld_matrices, ld_variant_indices, ld_labels):
        ht = hl.read_table(ld_ht_path)
        ht = ht_snp.join(ht, "inner")
        ht = ht.checkpoint(new_temp_file())

        lead_idx = ht.filter(hl.variant_str(ht.locus, ht.alleles) == args.lead_variant).head(1).idx.collect()

        if len(lead_idx) == 0:
            df[col] = np.nan
            continue

        idx = ht.idx.collect()
        # sorted index - BlockMatrix.filter requires sorted index
        idx2 = sorted(list(set(idx)))
        # relative lead pos in sorted index
        idx3 = np.where(np.array(idx2) == lead_idx[0])[0].tolist()
        # relative lead pos in original index
        idx4 = np.where(np.array(idx) == lead_idx[0])[0].tolist()

        bm = BlockMatrix.read(ld_bm_path)
        bm = bm.filter(idx2, idx2)
        v_row = bm.filter_rows(idx3).to_numpy()[0]
        v_col = bm.filter_cols(idx3).T.to_numpy()[0]
        if not np.all(np.diff(idx) > 0):
            order = np.argsort(idx)
            rank = np.empty_like(order)
            _, inv_idx = np.unique(np.sort(idx), return_inverse=True)
            rank[order] = inv_idx
            v_row = v_row[rank]
            v_col = v_col[rank]

        # re-densify
        r2 = v_row + v_col
        r2[idx4] = r2[idx4] / 2
        if not args.export_r:
            r2 = r2**2

        idx_snp = ht.idx_snp.collect()
        df[col] = np.nan
        df.iloc[idx_snp, df.columns.get_loc(col)] = r2

    if args.weighted_average_r is not None:
        n_samples = []
        ld = []
        for k, v in args.weighted_average_r.items():
            if isinstance(v, str):
                if v not in df.columns:
                    print(f"Column {v} not found.")
                    continue
                n_samples.append(df[v].values)
            else:
                n_samples.append(np.tile(v, len(df.index)))
            ld.append(df[f"gnomad_lead_r_{k}"].values)
        if len(n_samples) == 1:
            df["r"] = ld[0]
        else:
            n_samples = np.array(n_samples).T
            ld = np.array(ld).T
            df["r"] = np.nansum(n_samples * ld, axis=1) / np.nansum(n_samples * ~np.isnan(ld), axis=1)
    elif args.ld_reference == "custom":
        df["r"] = df[ld_labels[0]]
    else:
        df["r"] = df["gnomad_lead_r_nfe"]

    if args.dentist_s:
        lead_z = (df.beta / df.se).iloc[lead_idx_snp]
        df["t_dentist_s"] = ((df.beta / df.se) - df.r * lead_z) ** 2 / (1 - df.r**2)
        df["t_dentist_s"] = np.where(df["t_dentist_s"] < 0, np.inf, df["t_dentist_s"])
        df.iloc[lead_idx_snp, df.columns.get_loc("t_dentist_s")] = np.nan
        df["nlog10p_dentist_s"] = sp.stats.chi2.logsf(df["t_dentist_s"], df=1) / -np.log(10)

    if args.out.startswith("gs://"):
        fopen = hl.hadoop_open
    else:
        fopen = open

    with fopen(args.out, "w") as f:
        df.drop(columns=["variant"]).to_csv(f, sep="\t", na_rep="NA", index=False)

    if args.summary:
        df["r2"] = df.r**2
        n_r2 = np.sum(df.r2 > args.r2_threshold)
        n_na = np.sum(np.isnan(df.r2))
        outlier_idx = (df.r2 > args.r2_threshold) & (df.nlog10p_dentist_s > args.nlog10p_dentist_s_threshold)
        n_dentist_s_outlier = np.sum(outlier_idx)
        max_pip_idx = df.prob.idxmax()
        variant = df.chromosome.str.cat([df.position.astype(str), df.allele1, df.allele2], sep=":")
        expr = {
            "lead_pip_variant": [variant.iloc[max_pip_idx]],
            "n_total": [len(df.index)],
            "n_r2": [n_r2],
            "n_na": [n_na],
            "n_dentist_s_outlier": [n_dentist_s_outlier],
            "fraction": [n_dentist_s_outlier / n_r2 if n_r2 > 0 else 0],
            "max_pip": [np.max(df.prob)],
        }

        if args.annotate_consequence:
            nonsyn_idx = (df.r2 > args.r2_threshold) & df.consequence.isin(["pLoF", "Missense"])
            expr = {
                **expr,
                "n_nonsyn": [np.sum(nonsyn_idx)],
                "n_nonsyn_outlier": [np.sum(nonsyn_idx & outlier_idx)],
                "max_pip_nonsyn": [np.max(df.prob.loc[nonsyn_idx])],
                "cs_nonsyn": [np.any(df.cs.loc[nonsyn_idx])],
                "cs_99_nonsyn": [np.any(df.cs_99.loc[nonsyn_idx])],
                "nonsyn_variants": [",".join(variant.loc[nonsyn_idx].values)],
            }

        if "n_samples" in df.columns:
            if args.case_control:
                if "n_cases" in df.columns:
                    df["n_eff_samples"] = df.n_samples * (df.n_cases / df.n_samples) * (1 - df.n_cases / df.n_samples)
                else:
                    df["n_eff_samples"] = np.nan
            else:
                df["n_eff_samples"] = df.n_samples

            n_eff_r2 = df.n_eff_samples.loc[df.r2 > args.r2_threshold]
            expr = {
                **expr,
                "min_neff_r2": [np.nanmin(n_eff_r2) if n_r2 > 0 else np.nan],
                "max_neff_r2": [np.nanmax(n_eff_r2)] if n_r2 > 0 else np.nan,
            }

        df_summary = pd.DataFrame(expr)
        with fopen(args.out_summary, "w") as f:
            df_summary.to_csv(f, sep="\t", na_rep="NA", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--snp", type=str, required=True, help="Input snp file from fine-mapping")
    parser.add_argument("--out", type=str, required=True, help="Output path")
    parser.add_argument("--out-summary", type=str, help="Output summary path")
    parser.add_argument("--delimiter", type=str, default=" ", help="Delimiter for output ld matrix")
    parser.add_argument("--lead-variant", type=str, help="Lead variant to annotate gnomAD LD")
    parser.add_argument(
        "--lead-variant-choice",
        type=str,
        default="p",
        choices=["p", "prob", "gamma", "gamma-p"],
        help="Strategy for choosing a lead variant",
    )
    parser.add_argument("--align-alleles", action="store_true", help="Whether to align alleles with gnomAD")
    parser.add_argument(
        "--annotate-cups", action="store_true", help="Whether to annotate novel conversion-unstable positions (CUPs)"
    )
    parser.add_argument("--annotate-consequence", action="store_true", help="Whether to annotate VEP consequences")
    parser.add_argument("--annotate-gnomad-freq", action="store_true", help="Whether to annotate gnomAD frequencies")
    parser.add_argument(
        "--ld-reference", type=str, default="gnomad", choices=["gnomad", "custom"], help="Choice of LD reference"
    )
    parser.add_argument("--custom-ld-path", type=str, help="Path of user-provided LD BlockMatrix")
    parser.add_argument("--custom-ld-variant-index-path", type=str, help="Path of user-provided LD variant index table")
    parser.add_argument("--custom-ld-label", type=str, help="Label of user-provided LD")
    parser.add_argument("--export-r", action="store_true", help="Export signed r values instead of r2")
    parser.add_argument("--weighted-average-r", type=str, nargs="+", action=ParseKwargs, help="")
    parser.add_argument("--dentist-s", action="store_true", help="Annotate DENTIST-S statistics")
    parser.add_argument("--abf", action="store_true", help="Run ABF")
    parser.add_argument("--abf-prior-variance", type=float, default=0.04, help="Prior effect size variance for ABF")
    parser.add_argument(
        "--reference-genome",
        type=str,
        default="GRCh37",
        choices=["GRCh37", "GRCh38"],
        help="Reference genome of sumstats",
    )
    parser.add_argument("--summary", action="store_true", help="Whether to output a summary file")
    parser.add_argument("--case-control", action="store_true", help="Whether the input is from a case-control study")
    parser.add_argument(
        "--r2-threshold", type=float, default=0.6, help="r2 threshold of DENTIST-S outlier variants for prediction"
    )
    parser.add_argument(
        "--nlog10p-dentist-s-threshold",
        type=float,
        default=4,
        help="-log10 DENTIST-S P value threshold of DENTIST-S outlier variants for prediction",
    )

    args = parser.parse_args()

    if args.out_summary is None:
        args.out_summary = f"{os.path.splitext(args.out)[0]}.summary.txt"

    if args.ld_reference == "custom" and (
        (args.custom_ld_path is None) or (args.custom_ld_variant_index_path is None) or (args.custom_ld_label is None)
    ):
        raise argparse.ArgumentError(
            "All of --custom-ld-path, --custom-ld-variant-index-path, and --custom-ld-label should be provided"
        )

    main(args)
