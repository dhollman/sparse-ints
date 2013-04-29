# -*- coding: utf-8 -*-
from __future__ import division
from itertools import product
import math
from warnings import warn
import numpy as np
import sys
import os
from os.path import exists
import grendel
import urllib
from grendel.chemistry import Molecule
from grendel.gmath import Tensor, IndexRange, IndexingContext, DeclareIndexRange, chopped, Matrix, issubrange, Vector
from grendel.gmath.einsum import EinsumTensor, EinsumSum, EinsumContraction
import psi4
from binfile_loader import ProgressBar, tensor_from_binfile, tensor_and_basis_sets_from_binfile

termwise_printing = EinsumSum.termwise_printing

#--------------------------------------------------------------------------------#

geoms = {
    "h2" : """
        H   0.0  0.0  -0.35
        H   0.0  0.0   0.35
    """,
    "h2o" : """
        H   1.5  0.0  -0.3
        H  -1.5  0.0  -0.3
        O   0.0  0.0   1.0
    """,
    "NeAtom" : """
        Ne  0.0  0.0  0.0
    """
}

basis_map = {
    "cc-pVDZ-F12-CABS" : "cc-pVDZ-F12-OPTRI"
}

def psi_basis_name(mpqc_basis_name):
    if mpqc_basis_name in basis_map:
        return basis_map[mpqc_basis_name]
    else:
        return mpqc_basis_name

#--------------------------------------------------------------------------------#

if len(sys.argv) == 2:
    filename = sys.argv[1]
elif __name__ != "__main__":
    filename = "matrices/h2_cc-pVDZ_cc-pVDZ-F12-CABS_O.bin"
else:
    raise ValueError("Expecting exactly one filename as argument.")

if not exists(filename):
    raise IOError("File not found: {}".format(filename))
directory = os.path.join(filename.split("/")[:-1])
baseext = filename.split("/")[-1]
ext = baseext.split(".")[-1]
base = ".".join(baseext.split(".")[:-1])
nameparts = base.split("_")
matname = nameparts[-1]
molname = nameparts[0]
basis_name = psi_basis_name(nameparts[1])
aux_basis_name = psi_basis_name(nameparts[2])
modifiers = nameparts[3:-1]
if "bsrrrr" in modifiers:
    raise NotImplementedError()
if matname not in ["P", "Q", "O"] and "allints" not in modifiers:
    raise NotImplementedError()

#--------------------------------------------------------------------------------#

if molname not in geoms:
    raise ValueError("Don't know geometry for molecule named '{}'".format(molname))

#--------------------------------------------------------------------------------#

# constants used in the program, put up here for quick changing
correlation_factor_exponent = 1.0
linear_dependency_cutoff = 1e-8

#--------------------------------------------------------------------------------#

# Debugging switches and such
sanity_check = False
print_intermediates = False
load_integrals = True
load_transform = load_integrals
save_integrals = True
save_transform = save_integrals
keep_ao_versions = False
enforce_brillouin = False
debug = False
EinsumContraction.factorize_contraction = True
diagonal_ansatz = True

#--------------------------------------------------------------------------------#
# Create an indexing context unique to this file
val_checking_context = IndexingContext()
IndexRange.begin_indexing_context(val_checking_context)

# Quick print functions
cstr = lambda x: chopped(x, 1e-10).formatted_string()
fstr = lambda x: x.formatted_string(float_digits=10, float_width=14, label_width=14, float_type_string='f')

#--------------------------------------------------------------------------------#

print "{:#^80s}".format(" Begin PSI4 Babbeling ")

# Set up the molecule and options
mol = psi4.geometry("symmetry c1\n" + geoms[molname].strip())
psi4.options.basis = basis_name
psi4.options.df_basis_scf = aux_basis_name
psi4.options.e_convergence = 10
psi4.options.scf_type = 'direct'

# Create the BasisSet objects needed
parser = psi4.Gaussian94BasisSetParser()
ao_basis = psi4.BasisSet.construct(parser, mol, 'BASIS')
aux_basis = psi4.BasisSet.construct(parser, mol, 'DF_BASIS_SCF')

# Run SCF
E_scf = psi4.energy("scf")
print('SCF Energy: {}'.format(E_scf))

# Get a reference to the SCF wavefunction object
wfn = psi4.reference_wavefunction()

# Get some OrbitalSpace objects
docc_space = wfn.alpha_orbital_space('i', 'SO', 'OCC')
virt_space = wfn.alpha_orbital_space('a', 'SO', 'VIR')

# Build the RI space (the CABS+ space)
ri_space = psi4.OrbitalSpace.build_ri_space(aux_basis, ao_basis, linear_dependency_cutoff)
ri_basis = ri_space.basisset()

# Create a couple of MintsHelper objects to let us get stuff
mints_ao = psi4.MintsHelper(ao_basis)
mints_ri = psi4.MintsHelper(ri_basis)

# Get sizes of stuff
ncore = 0 # for now
nao = ao_basis.nbf()
nobs = nmo = wfn.nmo()
ndocc = wfn.doccpi()[0]
nao_ri = ri_basis.nbf()

# Build the CABS space using Psi (just to get the dimensions of the CABS space, so we
#   can construct the index ranges before building the CABS space ourselves)
cabs = psi4.OrbitalSpace.build_cabs_space(docc_space, ri_space, linear_dependency_cutoff)
cabs = psi4.OrbitalSpace.build_cabs_space(virt_space, cabs, linear_dependency_cutoff)

# Get sizes of more stuff
ncabs = cabs.dim()[0]
nri = nmo + ncabs
naux = aux_basis.nbf()
nao_aux = aux_basis.nbf()

#--------------------------------------------------------------------------------#

# Declare the index ranges
# Yeah, these are spinfree orbitals.  But I get tired of using upper case all the time...so deal with it
DeclareIndexRange("p',q',r',s',t',u'", 0, nri, name='Full RI MO space').with_subranges(
    IndexRange('p,q,r,s,t,u', 0, nmo, name='MO Orbital Basis').with_subranges(
        IndexRange('w,x,y,z', 0, ndocc, name='Geminal Generating Space').with_subranges(
            IndexRange('fi,fj,fk,fl', 0, ncore, name='MO Core Orbitals'),
            IndexRange('i,j,k,l,m,n', ncore, ndocc, name='MO Occupied Space')
        ),
        IndexRange('a,b,c,d,e,f', ndocc, nmo, name='MO Virtual Space')
    ),
    IndexRange("p'',q'',r'',s'',t'',u''", nmo, nri, name='Auxiliary MO space')
)
DeclareIndexRange("µ',ν',ρ',δ',σ',τ',κ',λ',γ',ε',ζ'", 0, nao_ri, name='Full RI AO Basis').with_subranges(
    IndexRange("µ,ν,ρ,σ,δ,τ,κ,λ,γ,ε,ζ", 0, nao, name='AO Orbital Space'),
    IndexRange("µ'',ν'',ρ'',σ''", nao, nao_ri, name='Auxiliary AO Space')
)

#--------------------------------------------------------------------------------#

# Get the overlap matrix
S = Tensor("µ',ν'")
Smat = mints_ri.so_overlap()
for mup, pp in product(xrange(nri), xrange(nri)):
    S[mup, pp] = Smat[0, mup, pp]

#--------------------------------------------------------------------------------#

# Get the orthogonalization coefficients
cmatri = ri_space.C()
Cortho = Matrix([[cmatri[0,i,j] for j in xrange(nri)] for i in xrange(nao + nao_aux)], indices="µ',p'")

# Get the orbital basis set coefficients
Cobs = Tensor("µ,p")
Cmat_docc = docc_space.C()
for mu, i in product(xrange(nao), xrange(ndocc)):
    Cobs[mu, i] = Cmat_docc[0, mu, i]
Cmat_virt = virt_space.C()
for mu, a in product(xrange(nao), xrange(ndocc,nmo)):
    Cobs[mu, a] = Cmat_virt[0, mu, a-ndocc]

# Build the CABS space
S12 = Matrix(indices="p,p'")
S12["p,p'"] = Cortho["µ',p'"] * S["µ',ν"] * Cobs["ν,p"]
_, __, V = np.linalg.svd(S12)
V = V.view(Matrix) * Cortho.T
Ccabs = Matrix(V[nmo:,:]).T

# and the orbital energies...
emat = wfn.epsilon_a()
e = Tensor("p")
for p in xrange(nmo):
    e[p] = emat[p]
ef = e[ndocc-1]

#--------------------------------------------------------------------------------#

# Now get the coefficient matrices
C = Tensor("µ',p'")
C["µ,p"] = Cobs["µ,p"]
C["µ',p''"] = Ccabs["µ',p''"]

# and the orbital energies...
emat = wfn.epsilon_a()
e = Tensor("p")
for p in xrange(nmo):
    e[p] = emat[p]

# create a CorrelationFactor object
cf = psi4.FittedSlaterCorrelationFactor(correlation_factor_exponent)

#--------------------------------------------------------------------------------#
# Get the density matrices

P = Tensor("µ,ν", name="P")
P["µ,ν"] = C["µ,i"] * C["ν,i"]
Q = Tensor("µ,ν", name="Q")
Q["µ,ν"] = C["µ,a"] * C["ν,a"]
O = Tensor("µ',ν'", name="O")
O["µ',ν'"] = C["µ',r''"] * C["ν',r''"]

#Sinv_obs = Tensor(S["µ,ν"].t.I, indices="µ,ν", name="Sinv_obs")
#Sinv = Tensor(S.I, indices="µ',ν'", name="Sinv")

#--------------------------------------------------------------------------------#

def psi_memmap_filename(tensname, indices):
    return os.path.join(
        directory,
        ".psi4ints_" + "_".join([molname, basis_name, aux_basis_name, tensname, urllib.quote(indices)]) + ".npmm"
    )

def load_psi_memmap(tensname, indices, force_reload=False, **kwargs):
    fname = psi_memmap_filename(tensname, indices)
    idxs = EinsumTensor.split_indices(indices)
    shape = tuple(len(val_checking_context[idx]) for idx in idxs)
    newkw = dict(
        mode="XX",
        dtype='float64',
        shape=shape
    )
    newkw.update(kwargs)
    if not force_reload and exists(fname):
        if newkw['mode'] == "XX": newkw['mode'] = "r"
        rv = np.memmap(fname, **newkw)
        return rv, True
    else:
        if newkw['mode'] == "XX": newkw['mode'] = "w+"
        rv = np.memmap(fname, **newkw)
        return rv, False

def get_2e_ao_integrals(name, indices, force_reload=False):
    mints_functions = {
        "g" : "ao_eri",
        "F" : "ao_f12",
        "Fg" : "ao_f12g12",
        "Fsq" : "ao_f12_squared",
        "DC" : "ao_f12_double_commutator"
    }
    if name not in mints_functions:
        raise ValueError("Don't know how to compute integrals named '{}'".format(name))
    intmap, needs_compute = load_psi_memmap(name, indices, force_reload=force_reload)
    if needs_compute:
        idxs = EinsumTensor.split_indices(indices)
        idxrngs = [val_checking_context[idx] for idx in idxs]
        if all(rng.name == 'Full RI AO Basis' for rng in idxrngs):
            function = getattr(mints_ri, mints_functions[name])
        elif all(rng.name == 'AO Orbital Space' for rng in idxrngs):
            function = getattr(mints_ao, mints_functions[name])
        else:
            raise NotImplementedError("Can't compute anything but all RI or all OBS integrals currently.")

        # Tell Psi to do the actual computation
        print "Computing tensor '{}'...".format(name),
        sys.stdout.flush()
        if name is "g":
            mat = function()
        else:
            mat = function(cf)
        print "done!"

        # Move the computed integrals into the memmap so we can use them
        print "Moving tensor '{}' data to python tensor structure...".format(name)
        pbar = ProgressBar(np.prod([len(rng) for rng in idxrngs]), indent='  ')
        size_q, size_s = len(idxrngs[1]), len(idxrngs[3])
        for p, q, r, s in product(*[xrange(rng.begin_index, rng.end_index) for rng in idxrngs]):
            pp,qq,rr,ss = tuple(x - rng.begin_index for x, rng in zip((p,q,r,s), idxrngs))
            intmap[pp,qq,rr,ss] = mat[0, p*size_q + q, r*size_s + s]
            pbar+= 1
        print "\ndone!"
    return Tensor(intmap, indices=indices)

def check_same(tens1, tens2):
    diff = tens1 - tens2
    maxabs = diff.max_abs()
    print("Should be zero: {}".format(maxabs))
    if maxabs > 1e-8:
        print('='*80)
        print(diff.zero_structure())
        print('-'*80)
        print(fstr(tens1))
        print('-'*80)
        print(fstr(tens2))
        print('='*80)

#--------------------------------------------------------------------------------#

print "{:#^80s}".format(" End PSI4 Babbeling ")
print "\n"

print "Loading tensor file '{}'".format(filename)
loaded_memmap_tens, bs = tensor_and_basis_sets_from_binfile(filename, verbose=True)
# mapping of the indices from PSI ordering to MPQC ordering
#   e.g. p-orbitals (y,z,x relative to PSI) of p-orbitals
index_map = []
if not all(bs[0] == b for b in bs):
    raise NotImplementedError()
for shell in bs[0]:
    if shell.nfunction == 1:
        rel_order = [0]
    elif shell.nfunction == 3:
        rel_order = [1,2,0]
    elif shell.nfunction == 5:
        rel_order = [4,2,0,1,3]
    else:
        # We probably need to do some reordering, I just haven't figured it out yet
        rel_order = range(shell.nfunction)
    # Now append the indices
    for disp in rel_order:
        index_map.append(shell.first_bf+disp)

print "done!\n"

# Get the psi counterpart of the tensor we want to check
if matname in ["P", "Q", "O"]:
    print "Checking density matrix {}".format(matname)
    psi_mat = globals()[matname]
elif matname in ["g", "F", "Fg", "Fsq", "DC"]:
    print "Checking untransformed integrals {}".format(matname)
    psi_mat = get_2e_ao_integrals(matname, "µ,ν,ρ,σ")
else:
    raise NotImplementedError()

print "Reindexing PSI tensor to match MPQC ordering..."
psi_mat = psi_mat.reindexed(index_map, reverse=True)

nonzero_thresh = 1e-8
nonzero_diff_amount = 1e-5
def maxabs_with_index(tens):
    idxs = tuple(np.unravel_index(np.argmax(abs(tens)), tens.shape))
    return float(tens[idxs]), idxs


def nonzero_diff(tens, reference):
    tmp = abs(tens)
    scaled = reference*abs(nonzero_diff_amount)
    return np.logical_and(tmp > scaled, tmp > nonzero_thresh)

def num_nonzero_diff(tens, reference):
    return np.sum(nonzero_diff(tens, reference))


if not psi_mat.shape == loaded_memmap_tens.shape:
    print "  Failed:  Tensors aren't even the same shape."
    print "    Loaded tensor shape: {}".format(loaded_memmap_tens.shape)
    print "    PSI tensor shape: {}".format(psi_mat.shape)
else:
    tens_to_check = Tensor(loaded_memmap_tens, indices=psi_mat.indices, name="Loaded "+matname)
    diff = psi_mat - tens_to_check
    maxdiff, position = maxabs_with_index(diff)
    avg_mag = (abs(tens_to_check) + abs(psi_mat)) / 2
    diffnz = nonzero_diff(diff, avg_mag)
    nonzero_count = np.sum(diffnz)
    nonzero_count = num_nonzero_diff(diff, avg_mag)
    tot_elements = np.prod(diff.shape)
    nonzero_percent = float(nonzero_count/tot_elements)*100

    if nonzero_percent > 10:
        print "\n" + "="*80 + "\n"
        print "PSI version:\n"
        print fstr(psi_mat)
        print "-"*80 + "\n"
        print "Loaded version:"
        print fstr(tens_to_check)
        print "-"*80 + "\n"
        print "Difference zero structure:"
        print (diffnz*np.sign(diff)).zero_structure(max_width=200)
        print "\n" + "="*80 + "\n"

    print "Summary\n" + "-"*80
    print "  Maximum difference magnitude: {:12.7e}".format(maxdiff)
    print "  Index of maximum difference element: {}".format(position)
    print "    PSI4 tensor value at {}: {:12.7f}".format(position, psi_mat[position])
    print "    MPQC tensor value at {}: {:12.7f}".format(position, tens_to_check[position])

    print "  Number of elements significantly different: {} (out of {})".format(
        nonzero_count, tot_elements
    )
    print "  Percent of elements significantly different: {:5.2f}%".format(
        nonzero_thresh, nonzero_percent
    )
    print "  Average difference magnitude: {:12.7e}".format(np.mean(abs(diff)))





# Get the various two electron integrals we need
#g = get_2e_ao_integrals('g', "µ,ν,ρ,σ")
#F = get_2e_ao_integrals('F', "µ,ν,ρ,σ")


