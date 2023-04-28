Test : Parabands - Zn2O2
Enabled : Yes
TestGroups : hdf5 serial parallel

Unpack: data.tgz
Copy: WFN_ref_no_pseudobands.h5
Copy: parabands.inp

%for l, nb in enumerate(('-1', '200')):
<%
npools = 2
cur_str = '{}'.format(l)
pb_out = 'parabands_{}.out'.format(cur_str)
compare_out = 'compare_wfns_{}.out'.format(cur_str)
wfn_out = 'WFN_out_{}.h5'.format(cur_str)
wfn_ref = 'WFN_ref_no_pseudobands.h5'
compare_nb = ''
%>\
${"###############################################################################"}

Banner: nb=${nb} (${cur_str})

Command: sed -i.bak 's/.*number_pools.*/number_pools ${npools}/;\
s/.*number_bands.*/number_bands ${nb}/;\
s/.*output_wfn_file.*/output_wfn_file ${wfn_out}/' parabands.inp

Executable: parabands.cplx.x
Processors: 4
Output: ${pb_out}
Input: NONE

Executable: compare_wfns.cplx.x
Arguments: ${wfn_ref} ${wfn_out} --terse --tol_deg=1e-10 ${compare_nb} > ${compare_out}

Precision: 1e-10
match; Maximum error in eigenvalues             ; LINE(${compare_out}, 1, 1); 0.0
match; Maximum error in WFN orthonormality      ; LINE(${compare_out}, 3, 1); 0.0
match; Maximum error in cross WFN orthonormality; LINE(${compare_out}, 4, 1); 0.0
Precision: 1e-9
match; Maximum error in reassembled Hamiltonian ; LINE(${compare_out}, 5, 1); 0.0

%endfor
