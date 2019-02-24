
# SigmaTransform-js
Implements the Sigma Transform in JavaScript.

## Contents
- [General Information](#general-information)
- [Usage and Examples](#usage-and-examples)

# General information
This repository shows an exemplary implementation of the "SigmaTransform", as defined in the thesis _"Quantum Frames and Uncertainty Principles arising from Symplectomorphisms"_, written in JavaScript.

Note that this library is not intended to show maximal performance, but show the usability of the universal interface of the "Sigma Transform" to perform well-known signal processing transforms / algorithms like the Short-Time Fourier Transform and Wavelet Transform.

Currently, only the one- and two-dimensional SigmaTransform is implemented and utilizes the jsfft library from
https://github.com/dntj/jsfft.

# Usage and Examples
Clone this repository, or copy the files to a local folder.
## Usage
Perform a Wavelet Transform on a signal "f" and get coefficients and reconstruction
```javascript
// get SigmaTransform class
const st = require('./SigmaTransform.js');
// define a spectral diffeomorphism, e.g. logarithm for WaveletTransform
const sig   = x => Math.log2(x*(x>0));
// define a window, e.g. a standard Gaussian (in the warped Domain)
const gauss = x => Math.exp(-Math.PI * x**2 );
// create WaveletTransform object
const WT = new st.SigmaTransform1D( 
    sig   , // the diffeomorphism
    gauss , // the window
    Fs    , // the sampling frequency
    N     , // the signal-size
    chans   // the channels in warped fourier domain
);
// compute WT of signal "f"...
WT.analyze( f );
// ...get coefficients
var coeff = WT.coeff;
// try to reconstruct (no dual windows involved) and get data
var rec = WT.synthesize().rec;
```
## Examples
The javascript files
```
    Example1D_STFT.js
    Example1D_Wavelet.js
    Example2D_STFT.js
    Example1D_SIM2.js
```
provide more elaborated usage examples and may be run with node.js via
```
    node Example1D_STFT.js
    node Example1D_Wavelet.js
    node Example2D_STFT.js
    node Example1D_SIM2.js
```
