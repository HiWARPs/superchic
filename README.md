## SuperChic v3.05 LbyL (by [Lucian Harland-Lang](lucian.harland-lang@physics.ox.ac.uk))

SuperChic 3 is a Fortran based Monte Carlo event generator for central
exclusive production at parton level, for a range of Standard Model final
states. User–defined histograms may be output, as well as unweighted events
in HEPEVT and Les Houches formats. By default the program makes use
of the LHAPDF library, but otherwise the code is completely stand–alone.


Compiles/runs as in manual `doc/superchic3.pdf`, but new additions wrt LbyL *not* described there. The new inputs for this are provided at the bottom of the input files - they had been set up for the run2 and run3 kinematics described in the Overleaf (additional cuts on e.g. acoplanarity not included). 

## Installation
1. Clone the repository

```
git clone https://github.com/HiWARPs/superchic.git
```

2. Install [LHAPDF](https://lhapdf.hepforge.org/install.html) library using Homebrew package manager

```
brew tap davidchall/hep
brew install lhapdf
```

  2a. Alternatively, follow Sections 2 and 3 of the manual `doc/superchic3.pdf' to install the library. In that case, [Doxygen](http://www.doxygen.nl/download.html) might need to be installed as well.

3. Download the PDF set  [MMHT2014lo68cl](http://lhapdfsets.web.cern.ch/lhapdfsets/current/MMHT2014lo68cl.tar.gz) and untar into the `$LHAPDF/share/LHAPDF/` directory, where $LHAPDF is the intallation directory of the library (usually, export $LHAPDF=/usr/local/Cellar/lhapdf/6.2.1 )
```
mv MMHT2014lo68cl.tar.gz $LHAPDF/share/LHAPDF/
cd $LHAPDF/share/LHAPDF
tar xvfz MMHT2014lo68cl.tar.gz
```

4. Compile by running make in the project directory containing the makefile
```
make
```

## Running

In the bin directory run
```
./init < input_run3.DAT
```
(or input_run2.DAT). Next run
```
./superchic < input_run3.DAT
```
(or input_run2.DAT).

### Possible issues


```
    libc++abi.dylib: terminating with uncaught exception of type LHAPDF::ReadError: Info file not found for PDF set 'MMHT2014lo68cl'

    Program received signal SIGABRT: Process abort signal.
```
means the required PDF set is not found, see Step 3 of the Installation.

## Elastic LbyL scattering part of the code

The actual implementation of the various LbyL amplitudes can be found in the following sources:
- `src/subamps/lightlightpol.f` - main subroutine path.
- `src/subamps/RHint.f`         - source of `rhint` interpolator of the input forward LbL low-energy amplitudes.
- `src/subamps/Msrule.dat`      - forward LbL low-energy amplitudes data file, M^{ansatz}(t=0) in Overleaf note, 
                                depending on `rts` in GeV from 0.3 to 20 with step 0.05.
