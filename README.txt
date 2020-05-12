Compiles/runs as in manual in doc/ folder, but new additions wrt LbyL *not* described there. The new inputs for this are provided at the bottom of the input files - I have set up these for the run2 and run3 kinematics described in the overleaf (additional cuts on e.g. acoplanarity not included). 

The actual implementation of the various LbyL amplitudes is in 'src/subamps/lightlightpol.f'

As discussed in the manual, you need to run

./init < input_run2.DAT

(or input_run3.DAT, it doesn't matter - only the rts has to be consistent for this) once. After that you should run

./superchic < input_***.DAT

for the corresponding case.

-------------------------------------------------------------------------------------------------------------------------
POSSIBLE ISSUES WITH LHAPDF INSTALLATION

Installation instructions and files for LHAPDF library can be found at https://lhapdf.hepforge.org/install.html .
The necessary path modification is shown in /doc/superchic3.pdf guide.

Probably LHAPDF library will need Doxygen. It can be found at http://www.doxygen.nl/download.html .

The library makes use of the PDF set: MMHT2014lo68cl. If LHAPDF cannot find this set, the code terminates with the following error:

    libc++abi.dylib: terminating with uncaught exception of type LHAPDF::ReadError: Info file not found for PDF set 'MMHT2014lo68cl'

    Program received signal SIGABRT: Process abort signal.

This set can be found in the list of the all sets at https://lhapdf.hepforge.org/pdfsets.html .
One should untar it in the /share/LHAPDF directory

-------------------------------------------------------------------------------------------------------------------------
ELASTIC LbL SCATTERING PART OF THE CODE

'src/subamps/lightlightpol.f' - main subroutine path.
'src/subamps/RHint.f'         - source of 'rhint' interpolator of the input forward LbL low-energy amplitudes.
'src/subamps/Msrule.dat'      - forward LbL low-energy amplitudes data file, M^{ansatz}(t=0) in Overleaf note, 
                                depending on rts in GeV from 0.3 to 20 with step 0.05.
