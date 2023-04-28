Test : Parabands - Benzothiazole
Enabled : Yes
TestGroups : hdf5 serial parallel

Unpack: data.tgz
Copy: WFN_ref.h5
Copy: parabands.inp

%for l, solver in enumerate((-1,1)):
<%
nb = -1
cur_str = '{}'.format(l)
pb_out = 'parabands_{}.out'.format(cur_str)
compare_out = 'compare_wfns_{}.out'.format(cur_str)
wfn_out = 'WFN_out_{}.h5'.format(cur_str)
# No need to find non-deg. subspace if we include all bands and there is a single kpt
no_deg = '--no_deg' if nb=='-1' else ''
%>\
${"###############################################################################"}

Banner: solver=${solver} (${cur_str})

Command: sed -i.bak 's/.*solver_algorithm.*/solver_algorithm ${solver}/;\
s/.*number_bands.*/number_bands ${nb}/;\
s/.*output_wfn_file.*/output_wfn_file ${wfn_out}/' parabands.inp

Executable: parabands.cplx.x
Processors: 4
Output: ${pb_out}
Input: NONE

Executable: compare_wfns.cplx.x
Arguments: WFN_ref.h5 ${wfn_out} --terse --tol_deg=1e-10 ${no_deg} > ${compare_out}

Precision: 1e-10
match; Maximum error in eigenvalues             ; LINE(${compare_out}, 1, 1); 0.0
match; Maximum error in WFN orthonormality      ; LINE(${compare_out}, 3, 1); 0.0
match; Maximum error in cross WFN orthonormality; LINE(${compare_out}, 4, 1); 0.0
Precision: 1e-9
match; Maximum error in reassembled Hamiltonian ; LINE(${compare_out}, 5, 1); 0.0

%endfor
