#!/usr/bin/env python

import sys
import re
from collections import defaultdict, namedtuple

import numpy as np
import pandas as pd
import limix # requires v2
from limix.stats import pca
from numpy_sugar.linalg import economic_qs
from glimix_core.glmm import GLMMExpFam
from glimix_core.lmm import LMM

def G_to_K(G):
    if G.shape[0] > 10000:
        sys.stderr.write("Very long G - has it been transposed?\n")
    return(np.dot(G, G.T))

def calc_var_comp(y, K_part, K_all, covar=None):
    K_resid = limix.qc.normalise_covariance(K_all - K_part)
    K_part = limix.qc.normalise_covariance(K_part)

    glmm = limix.glmm.GLMMComposer(len(y))
    #glmm._likname = 'bernoulli'
    glmm.y = y
    glmm.fixed_effects.append_offset()
    if covar is not None:
        glmm.fixed_effects.append(covar)

    glmm.covariance_matrices.append(K_part)
    glmm.covariance_matrices.append(K_resid)
    glmm.covariance_matrices.append_iid_noise()
    glmm.fit(verbose=False)

    part_scale = glmm.covariance_matrices[0].scale
    resid_scale = glmm.covariance_matrices[1].scale
    noise_scale = glmm.covariance_matrices[2].scale
    lml = glmm.lml

    return(part_scale, resid_scale, noise_scale, lml)


def main():

    print("Reading input")
    # read in genotype matrix (GT only)
    #G_all = np.loadtxt("nl_fixed.mat.csv", delimiter=",")
    #np.save("G_all.npy", G_all)
    G_all = np.load("G_all.npy")

    # read in metadata for genotype matrix
    # LoF information
    Allele = namedtuple('ref', 'alt')
    lofs = defaultdict(Allele)
    with open('nl_fixed.lof_only.pos.txt', 'r') as lof_file:
        for line in lof_file:
            (pos, ref, alt) = line.rstrip().split("\t")
            lofs[pos] = ref, alt

    # all positions
    positions = []
    lof_loc = []
    with open('nl_fixed.pos.txt', 'r') as pos_file:
        for idx, line in enumerate(pos_file):
            (pos, ref, alt) = line.rstrip().split("\t")
            positions.append(int(pos))

            if pos in lofs.keys():
                if ref == lofs[pos][0] and alt == lofs[pos][1]:
                    lof_loc.append(idx)
    positions = np.array(positions, dtype=int)

    # genes
    Gene = namedtuple('start', 'end')
    genes = defaultdict(Gene)
    with open('23FSpn_CDS.reg', 'r') as gene_file:
        header = gene_file.readline()
        for line in gene_file:
            (chrom, start, end, gene_id) = line.rstrip().split("\t")
            genes[gene_id] = int(start), int(end)


    # Construct phenotype (special for NL data, as pheno can be inferred from sample name)
    pheno  = []
    samples = []
    with open('vcf_sample_order.txt', 'r') as sample_file:
        for line in sample_file:
            samples.append(line.rstrip())
            if re.search("^1", line):
                pheno.append(1)
            else:
                pheno.append(0)
    y = np.array(pheno)

    # set up limix
    print("Set up limix")
    var_counts = np.sum(G_all, axis=1)
    #G_common = G_all[np.where((var_counts > 18) & (var_counts < 1819))[0], :].T
    #np.save("G_common.npy", G_common)
    G_common = np.load("G_common.npy")

    #r = limix.stats.pca(G_common.T, ncomp = 10)
    #np.save("r_components.npy", r['components'])
    #print("10 PCs variance explained: " + str(r['explained_variance_ratio']) + "\n")
    r = np.load("r_components.npy").T

    #K_all = np.dot(G_common, G_common.T)
    #np.save("K_all.npy", K_all)
    K_all_vcf = np.load("K_all.npy")

    K = pd.read_table("nl_C_vcf_fasttree.txt", index_col = 0)
    K = K.reindex(index=samples, columns=samples)
    K_all_tree = K.values

    print('h^2 %.3f' % limix.her.estimate(y, 'normal', K_all_tree, verbose=True))
    print('h^2 %.3f' % limix.her.estimate(y, 'bernoulli', K_all_tree, verbose=True))

    print("PC variance explained")
    # top PCs

    QS = economic_qs(K_all_tree)

    glmm = LMM(y, r, QS)
    glmm.fit(verbose=True)
    print(glmm.fixed_effects_variance/(glmm.v0 + glmm.v1 + glmm.fixed_effects_variance))

    print("Serotype variance explained")
    # serotype: cps locus 302490-323780
    #   w PCs
    #   w/o PCs

    # from alignment
    sero_idx = np.where((positions >= 302490) & (positions <= 323780))[0]
    G_sero = G_all[sero_idx,:].T
    print(calc_var_comp(y, G_to_K(G_sero), K_all_vcf, covar=None))
    print(calc_var_comp(y, G_to_K(G_sero), K_all_vcf, covar=r))

    # from tree
    K_sero_tree = pd.read_table("cps_C.tsv", index_col = 0)
    K_sero_tree = K_sero_tree.reindex(index=samples, columns=samples)
    K_sero_tree = K_sero_tree.values

    # both tree branches at same scale
    print(calc_var_comp(y, K_sero_tree, K_all_tree, covar=None))
    print(calc_var_comp(y, K_sero_tree, K_all_tree, covar=r))

    # categorical
    K_sero_cat = pd.read_table("sero_K_categorical.tsv", index_col = 0)
    found = np.isin(samples, K_sero_cat.index)
    K_sero_cat = K_sero_cat.reindex(index=samples, columns=samples)
    K_sero_cat = K_sero_cat.values
    print('h^2 %.3f' % limix.her.estimate(y[found], 'bernoulli', K_sero_cat[found,:][:,found], verbose=True)) #h^2 = 0.495
    #print('h^2 %.3f' % limix.her.estimate(y[found], 'bernoulli', K_sero_cat[found,:][:,found], M = r[found,:], verbose=True))# fails at iter 71
    print('h^2 %.3f' % limix.her.estimate(y[found], 'bernoulli', K_sero_cat[found,:][:,found], M = r[found,1:3], verbose=True))#h^2 = 0.508

    # fixed effects
    sero_covars_cat = pd.read_table("serotype_covars.txt", index_col = 0)
    assert(np.all(found == np.isin(samples, sero_covars_cat.index)))
    sero_covars_cat = sero_covars_cat.reindex(index=samples)
    sero_covars_cat = sero_covars_cat.values

    # variance decomposition w/ serotype fixed effects
    glmm = limix.glmm.GLMMComposer(len(y[found]))
    glmm.y = y[found]
    glmm.fixed_effects.append_offset()
    glmm.fixed_effects.append(sero_covars_cat[found,:])
    glmm.covariance_matrices.append(limix.qc.normalise_covariance(K_all_tree[found,:][:,found]))
    glmm.covariance_matrices.append_iid_noise()
    glmm.fit(verbose=False)
    print(glmm) #LML: -781.7942874960014

    glmm = limix.glmm.GLMMComposer(len(y[found]))
    glmm.y = y[found]
    glmm.fixed_effects.append_offset()
    glmm.covariance_matrices.append(limix.qc.normalise_covariance(K_all_tree[found,:][:,found]))
    glmm.covariance_matrices.append_iid_noise()
    glmm.fit(verbose=False)
    print(glmm) #LML: -870.4633045267495

    print("Gene burden variance explained")
    # LoF
    #   aggregate LoF burden
    #       w PCs
    #       w/o PCs
    #       <1% MAF only
    #       w/o singletons + doubletons
    #   LoF burden by gene
    lof_mat = G_all[lof_loc, :]
    lof_mat_pos = positions[lof_loc]

    gene_burden = []
    gene_burden_af001 = []
    gene_burden_ac3 = []
    for gene_id in genes.keys():
        start, end = genes[gene_id]
        gene_G = lof_mat[(lof_mat_pos >= start) & (lof_mat_pos <= end), :]
        gene_burden.append(np.where(np.sum(gene_G, axis=0) > 0, 1, 0))

        var_counts = np.sum(gene_G, axis=1)
        gene_G_af001 = gene_G[var_counts < 18, :]
        gene_burden_af001.append(np.where(np.sum(gene_G_af001, axis=0) > 0, 1, 0))
        gene_G_ac3 = gene_G[var_counts > 2, :]
        gene_burden_ac3.append(np.where(np.sum(gene_G_ac3, axis=0) > 0, 1, 0))

    gene_burden = np.array(gene_burden).T
    gene_burden_af001 = np.array(gene_burden_af001).T
    gene_burden_ac3 = np.array(gene_burden_ac3).T

    print(calc_var_comp(y, G_to_K(gene_burden), K_all_vcf, covar=None))
    print(calc_var_comp(y, G_to_K(gene_burden), K_all_vcf, covar=r))
    print(calc_var_comp(y, G_to_K(gene_burden_af001), K_all_vcf, covar=None))
    print(calc_var_comp(y, G_to_K(gene_burden_af001), K_all_vcf, covar=r))
    print(calc_var_comp(y, G_to_K(gene_burden_ac3), K_all_vcf, covar=None))
    print(calc_var_comp(y, G_to_K(gene_burden_ac3), K_all_vcf, covar=r))

    for gene_id, gene_G in zip(genes.keys(), gene_burden.T):
        try:
            print(gene_id)
            print(calc_var_comp(y, G_to_K(gene_G.reshape(-1,1)), K_all_vcf, covar=None))
        except:
            sys.stderr.write("Could not fit model\n")
            sys.exit(0)

if __name__ == '__main__':
    main()
