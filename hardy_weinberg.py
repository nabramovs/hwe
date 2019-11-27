import sys
import time
import math as m
from decimal import Decimal

def Hardy_Weinberg_Equilibrium_exact_test_user_Kantale(obs_hets, obs_hom1, obs_hom2, mid_p=True):
    """
	 This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
	 Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
	 Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000

	 Written by Jan Wigginton
    """

    if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
            raise Exception("FATAL ERROR - SNP-HWE: Current genotype configuration (%s  %s %s) includes negative count" % (obs_hets, obs_hom1, obs_hom2))

    obs_homc = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
    obs_homr = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2

    rare_copies = 2 * obs_homr + obs_hets # A
    genotypes   = obs_hets + obs_homc + obs_homr # N

    het_probs = [0.0] * (rare_copies + 1) # list of probabilities for all possible allele combinations

    #start at midpoint
    mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes) # don't understand why this is complicated, rare_copies / 2 works the same
    #mid = rare_copies / 2  

    #check to ensure that midpoint and rare alleles have same parity
    if (rare_copies & 1) ^ (mid & 1):
            mid += 1

    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = genotypes - curr_hets - curr_homr

    het_probs[mid] = 1.0
    sum = float(het_probs[mid]) # "N, na" from the paper 

    for curr_hets in xrange(mid, 1, -2):

            het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
            sum += het_probs[curr_hets - 2];

            # 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
            curr_homr += 1
            curr_homc += 1

    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = genotypes - curr_hets - curr_homr

    for curr_hets in xrange(mid, rare_copies - 1, 2):

            het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
            sum += het_probs[curr_hets + 2]

            #add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
            curr_homr -= 1
            curr_homc -= 1

    for i in xrange(0, rare_copies + 1):
            het_probs[i] /= sum

    '''
    #alternate p-value calculation for p_hi/p_lo
    p_hi = float(het_probs[obs_hets])
    for i in xrange(obs_hets, rare_copies+1):
            p_hi += het_probs[i]

    p_lo = float(het_probs[obs_hets])
    for i in xrange(obs_hets-1, -1, -1):
            p_lo += het_probs[i]

    p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo
    '''
    p_hwe = 0.0
    

    #  p-value calculation for p_hwe
    for i in xrange(0, rare_copies + 1):
            if het_probs[i] > het_probs[obs_hets]:
                    continue;
            p_hwe += het_probs[i]

    p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

    if mid_p:
        p_hwe -= het_probs[obs_hets] / 2

    return p_hwe


def main():
    # EXAMPLE
    '''
    het = 6
    hom1 = 1
    hom2 = 2
      
    print Hardy_Weinberg_Equilibrium_exact_test_user_Kantale(het, hom1, hom2, mid_p=True)
    '''
    
if __name__ == "__main__":
    sys.exit(main())