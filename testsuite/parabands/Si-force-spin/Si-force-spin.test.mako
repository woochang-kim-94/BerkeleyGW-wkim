Test : Parabands - Si, constrained spin pol.
Enabled : Yes
TestGroups : hdf5 serial parallel

Unpack: data.tgz
Copy: WFN_ref_no_pseudobands.h5
Copy: parabands.inp

%for k, npools in enumerate((1, 3)):
<%
cur_str = '{}'.format(k)
pb_out = 'parabands_{}.out'.format(cur_str)
compare_out = 'compare_wfns_{}.out'.format(cur_str)
wfn_out = 'WFN_out_{}.h5'.format(cur_str)
wfn_ref = 'WFN_ref_no_pseudobands.h5'
%>\
${"###############################################################################"}

Banner: npools=${npools} (${cur_str})

Command: sed -i.bak 's/.*number_pools.*/number_pools ${npools}/;\
s/.*output_wfn_file.*/output_wfn_file ${wfn_out}/' parabands.inp

Executable: parabands.cplx.x
Processors: 4
Output: ${pb_out}
Input: NONE

Executable: compare_wfns.cplx.x
Arguments: ${wfn_ref} ${wfn_out} --terse --tol_deg=1e-10 > ${compare_out}

Precision: 1e-10
match; Maximum error in eigenvalues             ; LINE(${compare_out}, 1, 1); 0.0
match; Maximum error in WFN orthonormality      ; LINE(${compare_out}, 3, 1); 0.0
match; Maximum error in cross WFN orthonormality; LINE(${compare_out}, 4, 1); 0.0
Precision: 1e-9
match; Maximum error in reassembled Hamiltonian ; LINE(${compare_out}, 5, 1); 0.0

%endfor
