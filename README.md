# SuperChic v3.05 LbyL (by [Lucian Harland-Lang](lucian.harland-lang@physics.ox.ac.uk))

SuperChic 3 is a Fortran based Monte Carlo event generator for central
exclusive production at parton level, for a range of Standard Model final
states. User–defined histograms may be output, as well as unweighted events
in HEPEVT and Les Houches formats. By default the program makes use
of the LHAPDF library, but otherwise the code is completely stand–alone.


Compiles/runs as in manual `superchic3.05LbyL/doc/superchic3.pdf`, but new additions wrt LbyL *not* described there. The new inputs for this are provided at the bottom of the input files - they had been set up for the run2 and run3 kinematics described in the Overleaf (additional cuts on e.g. acoplanarity not included). 

## Installing

General remarks on the installation process are given in the manual `superchic3.05LbyL/doc/superchic3.pdf` in sect.2 and 3.

1. Install [LHAPDF](https://lhapdf.hepforge.org/install.html) library and perform the neccessary path modification according to sect.3 of manual.

Possibly LHAPDF library will need [Doxygen](http://www.doxygen.nl/download.html) to be installed as well.

The final step is to download the pdf set  [MMHT2014lo68cl](http://lhapdfsets.web.cern.ch/lhapdfsets/current/MMHT2014lo68cl.tar.gz) and untar it in the `.../share/LHAPDF` directory

## Running

To run
```
./init < input_run2.DAT
```
(or input_run3.DAT, it doesn't matter - only the rts has to be consistent for this) once. After that you should run
```
./superchic < input_***.DAT
```
for the corresponding case.

### Possible issues

The error
```
    libc++abi.dylib: terminating with uncaught exception of type LHAPDF::ReadError: Info file not found for PDF set 'MMHT2014lo68cl'

    Program received signal SIGABRT: Process abort signal.
```
signalizes that the last step of installing process was skipped, and LHAPDF package cannot find the required PDF set.

## Elastic LbyL scattering part of the code

The actual implementation of the various LbyL amplitudes can be found in the following sources:
- `src/subamps/lightlightpol.f` - main subroutine path.
- `src/subamps/RHint.f`         - source of `rhint` interpolator of the input forward LbL low-energy amplitudes.
- `src/subamps/Msrule.dat`      - forward LbL low-energy amplitudes data file, M^{ansatz}(t=0) in Overleaf note, 
                                depending on `rts` in GeV from 0.3 to 20 with step 0.05.
