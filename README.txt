Compiles/runs as in manual in doc/ folder, but new additions wrt LbyL *not* described there. The new inputs for this are provided at the bottom of the input files - I have set up these for the run2 and run3 kinematics described in the overleaf (additional cuts on e.g. acoplanarity not included). 

The actual implementation of the various LbyL amplitudes is in 'src/subamps/lightlightpol.f'

As discussed in the manual, you need to run

./init < input_run2.DAT

(or input_run3.DAT, it doesn't matter - only the rts has to be consistent for this) once. After that you should run

./superchic < input_***.DAT

for the corresponding case.

LbyL specific flags set in lines 149 to 154. The 'loop' flag can be set to:

lep : lepton loop
quark: pQCD quark loop
w: w loop
pion: pion loop
meson: meson exchanges + pion loop
sum: rum rule result
tot_quark: lepton + pQCD quark
tot_meson: meson exchanges + pion loop + lepton
tot_sum: sum rule + lepton

'interpolate' is Fermi distribution interpolation, and can be used for meson, sum rule cases (as well as the corresponding 'tot').

'bottomonium' flag includes (.true.) or excludes (.false.) bottomonium resonances at the amplitude level in the considered lbl process. For particular implementation see lines 431-504 of lightlightpol.f
