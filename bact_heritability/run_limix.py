#!/usr/bin/env python

import re
from collections import defaultdict, namedtuple

import numpy as np
import limix
from limix.stats import pca
from numpy_sugar.linalg import economic_qs
from glimix_core.glmm import GLMMExpFam

def calc_var_comp(y, G_part, K_all, covar=None):
    K_part = np.dot(G_part, G_part.T)
    K_resid = limix.qc.normalise_covariance(K_all - K_part)
    K_part = limix.qc.normalise_covariance(K_part)

    glmm = limix.glmm.GLMMComposer(len(y))
    glmm.likname('bernoulli')
    glmm.y = y
    glmm.fixed_effects.append_offset()
    if covar:
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
            positions.append(pos)

            if pos in lofs.keys():
                if ref == lofs[pos][0] and alt == lofs[pos][1]:
                    lof_loc.append(idx)
    positions = np.array(positions)

    # genes
    Gene = namedtuple('start', 'end')
    genes = defaultdict(Gene)
    with open('23FSpn_CDS.reg', 'r') as gene_file:
        header = gene_file.readline()
        for line in gene_file:
            (chrom, start, end, gene_id) = line.rstrip().split("\t")
            genes[gene_id] = start, end


    # Construct phenotype (special for NL data, as pheno can be inferred from sample name)
    pheno = []
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
    var_counts = np.sum(G_all, axis=1)
    G_common = G_all[np.where((var_counts > 18) & (var_counts < 1819))[0], :]

    r = limix.stats.pca(G_common, ncomp = 10)
    np.save("r_components.npy", r['components'])
    print("10 PCs variance explained: " + str(r['explained_variance_ratio']) + "\n")

    K_all = np.dot(G_common, G_common.T)
    np.save("K_all.npy", K_all)
    print('h^2 %.3f' % limix.her.estimate(y, 'bernoulli', K, verbose=False))

    # top PCs
    # for PCs (https://media.nature.com/original/nature-assets/srep/2016/160525/srep26471/extref/srep26471-s1.pdf)
    #   ev(Xa) = ev(Y) - v0 - v1 - psi(X(X'(v0 + v1)X)^-1*X')
    #   ev(Y) = psi(YY')
    #   psi(A) = (1/(n-1))*(tr(A) - 1/n*A)
    QS = economic_qs(K_all)

    glmm = GLMMExpFam(y, 'bernoulli', r['components'], QS)
    glmm.fit(verbose=False)
    rand_var = glmm.v0 + glmm.v1
    psi = lambda A, n: (1/(n - 1)) * (np.trace(A) - (1/n)*A)
    evY = psi(np.dot(y, y.T), len(y))
    adjust = psi(np.dot(r['components'], np.dot(np.linalg.inv(np.dot(r['components'], r['components'])/rand_var), r['components'])), len(y))
    print(evY - rand_var)
    print(evY - rand_var - adjust)

    # serotype: cps locus 302490-323780
    #   w PCs
    #   w/o PCs
    sero_idx = np.where((positions >= 302490) & (positions <= 323780))[0]
    G_sero = G_all[sero_idx,]
    print(calc_var_comp(y, G_sero, K_all, covar=None))
    print(calc_var_comp(y, G_sero, K_all, covar=r['components']))

    # LoF
    #   aggregate LoF burden
    #       w PCs
    #       w/o PCs
    #       w/o singletons + doubletons
    #   LoF burden by gene
    lof_mat = G_all[lof_loc, :]
    lof_mat_pos = positions[lof_loc]

    gene_burden = []
    gene_burden_ac3 = []
    for gene_id in genes.keys():
        start, end = genes[gene_id]
        gene_G = lof_mat[np.which(lof_mat_pos >= start & lof_mat_pos <= end), :]
        gene_burden.append(np.where(np.sum(gene_G, axis=0) > 0, 1, 0))

        var_counts = np.sum(gene_G, axis=1)
        gene_G_ac3 = geneG[:, var_counts > 2]
        gene_burden_ac3.append(np.where(np.sum(gene_G_ac3, axis=0) > 0, 1, 0))

    print(calc_var_comp(y, gene_burden, K_all, covar=None))
    print(calc_var_comp(y, gene_burden, K_all, covar=r['components']))
    print(calc_var_comp(y, gene_burden_ac3, K_all, covar=None))

    for gene_id, gene_G in zip(genes.keys(), gene_burden):
        print(gene_id)
        print(calc_var_comp(y, gene_G, K_all, covar=None))

if __name__ == '__main__':
    main()
