// get SigmaTransform class
const st = require('./SigmaTransform.js');
// get filesystem
const fs = require('fs');
// get lineReader
const lineReader = require('readline').createInterface({
	  input: fs.createReadStream('bat.asc')
});

// read bat signal
var bat = [];
lineReader.on('line', function (line) {
	bat.push( parseFloat( line ) );
});

// done with reading?
lineReader.on('close', function() {
    // get length
	const N  = bat.length;
    // set sampling frequency
	const Fs = 143000;
    // set spectral diffeomorphism
    const sig   = x => Math.log2(x*(x>0));
    // set lower and upper bounds of the axis and stepsize
    const lower = sig(Fs*0.005), upper = sig(Fs/2*1.1), dx = (upper-lower)/N;
    // set width of the window in warped domain and number of channels
    const winwidth = 4, numsteps = 400;
    // create frequency channels
    const chans = new Float32Array( N ).map( (x,i) => lower+i*dx );
    // define a (Gaussian) window 
	const win = x => Math.exp(-Math.PI * ( x / sig(Fs) * numsteps / winwidth )**2 );
    // create WaveletTransform object
	const WaveletT = new st.SigmaTransform1D( sig , win , Fs , N , chans );
    // compute Wavelet Transform of bat signal
	var coeff = WaveletT.analyze( bat ).coeff;
    // try to reconstruct (using the same windows, so no dual windows involved)
	var rec = WaveletT.synthesize().rec;

    // save coefficients to file
    var file = fs.createWriteStream('coeff.asc');
    file.on('error', function(err) { console.log( "error", err ); });
    coeff.map( v => { v.map(x=>{file.write( '' + (x.real**2+x.imag**2) +',' );});file.write('\n'); } );
    file.end();

    // save reconstruction to file
    file = fs.createWriteStream('rec.asc');
    file.on('error', function(err) { console.log( "error", err ); });
    rec.map( v => { file.write( ''+v.real+'\n' ); });
    file.end();

    // save spectra of windows to file
    //file = fs.createWriteStream('wins.asc');
    //file.on('error', function(err) { console.log( "error", err ); });
    //WaveletT.windows.map( v => { v.map(x=>{file.write( '' + (x.real**2+x.imag**2) +',' );});file.write('\n'); } );
    //file.end();
});



