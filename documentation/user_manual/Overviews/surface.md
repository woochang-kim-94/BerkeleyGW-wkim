# `surface` code


## Input file structure (`surface.inp`)
```
header text                   
inputfilename file            
inputfileformat cube|xsf      
outputfilename file           
outputfileformat povray       
isovalue (0.0,1.0)            
sign positive|negative|both   
power 0|1|2                   
algorithm cube|tetrahedron    
smooth T|F                    
box T|F                       
basis T|F                     
uc T|F                        
uco                           
ucox ucoy ucoz                
ucv                           
ucv1x ucv1y ucv1z             
ucv2x ucv2y ucv2z             
ucv3x ucv3y ucv3z             
ucu bohr|angstrom|latvec      
sc T|F                        
sco                           
scox scoy scoz                
scv                           
scv1x scv1y scv1z             
scv2x scv2y scv2z             
scv3x scv3y scv3z             
scu bohr|angstrom|latvec      
sct T|F                       
```


## Example (HOMO of benzene)
```
header C6H6_band_15
inputfilename C6H6.b_15.cube
inputfileformat cube
outputfilename C6H6.b_15.pov
outputfileformat povray
isovalue 0.9
sign both
power 1
algorithm cube
smooth T
box F
basis F
uc F
uco
0.0 0.0 0.0
ucv
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
ucu latvec
sc T
sco
-0.5 -0.5 -0.5
scv
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
scu latvec
sct T
```


## Documentation

The isosurface is defined by `isovalue * max |psi|` for `power = 0` and by
`integral_v |psi|^power = isovalue * integral |psi|^power` for `power > 0`
where `psi` is the scalar field and `v` is the volume inside the isosurface.

Note that espresso already outputs `|psi|^2`, setting `power = 2` yields
`|psi|^4` in this case, to produce `|psi|^2` isosurface set `power = 1`.

`algorithm` stands for marching cubes and marching tetrahedra
smooth defines whether to average normals for smooth rendering

`box` specifies whether to generate supercell box in povray script

`basis` specifies whether to generate coordinate axes in povray script

Acronyms:

- `sfo` = scalar field origin (read from volumetric file)
- `sfv` = scalar field vectors (read from volumetric file)
- `uco` = unit cell origin
- `ucv` = unit cell vectors
- `ucu` = unit cell units
- `sco` = supercell origin
- `scv` = supercell vectors
- `scu` = supercell units
- `sct` = supercell translational symmetry

Rules:

- if `uc = T` then `uco`/`ucv` are read from parameter file
  else `uco`/`ucv` are set to `sfo`/`sfv`
- if `sc = T` then `sco`/`scv` are read from parameter file
  else `sco`/`scv` are set to` uco`/`ucv`
- if `uc = T` and `ucu = latvec` then `uco`/`ucv` are scaled by `sfv`
- if `sc = T` and `scu = latvec` then `sco`/`scv` are scaled by `ucv`
- if `sct = T` then `scv` is used for translational symmetry
  else `ucv` is used for translational symmetry

Hints: 

- In most cases, set `uc = F`, `uc = T` is only needed if volumetric data do
  not span the whole unit cell and the supercell is being used, so the correct
  unit cell must be defined.

- Use `sc = T` to construct the supercell spanning several unit cells or to
  assemble the unit cell around the origin of the coordinate system as in
  benzene example above (`sco = (-0.5 -0.5 -0.5)` in latvec units).

- The supercell may span a fractional number of unit cells, in this case set
  `sct = F` to produce the correct bonds at the faces of the supercell.
