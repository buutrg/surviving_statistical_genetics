
library("optparse")
# ==========================================================================

# ==========================================================================
option_list = list(
  make_option(c("--file"), type = "character"),
  make_option(c("--ncores"), type = "numeric"),
  make_option(c("--chr_col"), type = "character"),
  make_option(c("--pos"), type = "character"),
  make_option(c("--rsid"), type = "character"),
  make_option(c("--method"), type = "character", default="auto"),
  make_option(c("--effect_allele"), type = "character"),
  make_option(c("--non_effect_allele"), type = "character"),
  make_option(c("--n_eff"), type = "character"),
  make_option(c("--isbinary"), type = "logical"),
  make_option(c("--se"), type = "character"),
  make_option(c("--beta"), type = "character"),  
  make_option(c("--sdd_file"), type = "character"),
  make_option(c("--out"), type = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



library(data.table)
library(dplyr)
options(datatable.fread.datatable=FALSE)
library(bigsnpr)

# opt = data.frame(
#   file = "jointGwasMc_TG.txt_new",
#   chr_col = "chr",
#   pos = "poshg37",
#   rsid = "rawRSID",
#   effect_allele = "A1",
#   non_effect_allele = "A2",
#   n_eff = "N",
#   se = "se",
#   ncores=4,
#   beta = "b",
#   isbinary = F,
#   out = "TG"
# )

print(opt)

isbinary = opt$isbinary
ncores = opt$ncores
sdd_file = opt$sdd_file
method = opt$method


# if (file.exists(opt$out)) {
# 	q()
# }


sumstat_all = fread(opt$file)

if (method == "auto") sdd = fread(sdd_file)

head(sumstat_all)


sumstats_inp = sumstat_all %>% select(
  opt$chr_col, 
  opt$pos, 
  opt$rsid, 
  opt$effect_allele,
  opt$non_effect_allele,
  opt$n_eff,
  opt$se,
  opt$beta)

names(sumstats_inp) =
    c("chr",
    "pos",
    "rsid",
    "a1",
    "a0",
    "n_eff",
    "beta_se",
    "beta")

sumstats_inp$chr = as.numeric(sumstats_inp$chr)
sumstats_inp$pos = as.numeric(sumstats_inp$pos)
sumstats_inp$beta = as.numeric(sumstats_inp$beta)
sumstats_inp$beta_se = as.numeric(sumstats_inp$beta_se)
# sumstats_inp$MAF = as.numeric(sumstats_inp$MAF)

sumstats_inp$chisq = (sumstats_inp$beta / sumstats_inp$beta_se) ^ 2
sumstats_inp = sumstats_inp %>% filter(chisq > 0 & beta_se > 0)

map_ldref = readRDS("/broad/hptmp/btruong/ldpred2_prsmixGWAS/data/ldref/map.rds")

df_beta = snp_match(sumstats_inp, map_ldref)

in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], map_ldref[, c("chr", "pos")])
df_beta <- df_beta[in_test, ]


if (method == "auto") {
    # df_beta$beta_se <- with(df_beta, 1/sqrt(2*af_UKBB*(1 - af_UKBB)*(n_eff + beta^2)))
    # df_beta$beta <- with(df_beta, beta * beta_se)
    df_beta$oldid = sumstat_all$SNP[match(df_beta$rsid.ss, sumstat_all[,opt$rsid])]
    df_beta$sd_val = sdd$sdd[match(df_beta$oldid, sdd$ID)]

    sd_val = df_beta$sd_val

    if (!isbinary) {
    slope <- median(df_beta$sd_val * df_beta$beta_se * sqrt(df_beta$n_eff))
    print("Slope from lin reg")
    } else {
    slope <- 2
    print("Slope from log reg")
    }

    print(slope)


    sd_ss <- with(df_beta, slope / (beta_se * sqrt(n_eff)))

    # Remove SNPs with large diff between sd_ss-sd_val
    is_bad <- sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05

    print(paste("Filtered", sum(is_bad), "SNPs"))

    print(dim(df_beta))

    df_beta = df_beta[!is_bad,]

}

df_beta = df_beta %>% filter(beta_se > 0)


head(df_beta)

tmp <- tempfile(tmpdir = "tmp-data")

for (chr in 1:22) {

    cat(chr, ".. ", sep = "")

    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ## indices in 'corr_chr'
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

    corr_chr <- readRDS(paste0("/broad/hptmp/btruong/ldpred2_prsmixGWAS/data/ldref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]

    if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
    } else {
    corr$add_columns(corr_chr, nrow(corr))
    }
}

writeLines("")
n_eff = max(df_beta$n_eff, na.rm=T)

# Heritability estimation of LD score regression
# to be used as a starting value in LDpred2-auto
ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                chi2 = (beta / beta_se)^2,
                                sample_size = n_eff,
                                ncores = ncores))
print(ldsc)
h2_est <- ldsc[["h2"]]

ldsc_out = t(as.data.frame(ldsc))
fwrite(ldsc_out, paste0(opt$out, "_ldsc.txt"), row.names=F, quote=F, sep="\t")


if (method == "auto") {
    # LDpred2-auto
    set.seed(1)
    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
        vec_p_init = seq_log(1e-4, 0.9, length.out = 30),
        # vec_p_init = c(0.02),
        burn_in = 800, num_iter = 400,
        allow_jump_sign = FALSE, shrink_corr = 0.95,
        ncores = ncores, verbose = T)

    multi_auto_eff = do.call(cbind, 
    lapply(multi_auto, function(x) return(x$beta_est)))

    rr = function(x,d=5) round(as.numeric(x),d)

    colnames(multi_auto_eff) = unlist(lapply(multi_auto, function(x) {
    return(paste0(
        "h2init_", rr(x$h2_init), "__h2est_", rr(x$h2_est),
        "__pinit_", rr(x$p_init), "__pest_", rr(x$p_est)
        ))
    }))

    df_beta$newSNP = sumstat_all$SNP[match(df_beta$rsid.ss, sumstat_all$rawRSID)]

    multi_auto_eff_out = data.frame(
    df_beta %>% select(newSNP, a1, a0, rsid.ss),
    multi_auto_eff
    )
    head(multi_auto_eff_out)
    colnames(multi_auto_eff_out)[c(1:4)] = c("SNP", "A1", "A2", "rsid")

    fwrite(multi_auto_eff_out, paste0(opt$out, "_snpeff_ldpred2auto.txt"), row.names=F, quote=F, sep="\t")

}

if (method == "inf") {
    beta_inf = snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
    beta_inf_out = data.frame(
        df_beta %>% select(newSNP, a1, a0, rsid.ss),
        beta_inf
    )
    head(beta_inf_out)
    colnames(beta_inf_out)[c(1:4)] = c("SNP", "A1", "A2", "rsid")

    fwrite(beta_inf_out, paste0(opt$out, "_snpeff_ldpred2inf.txt"), row.names=F, quote=F, sep="\t")

}


