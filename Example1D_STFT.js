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
    // create frequency channels
    const chans = new Float32Array( N ).map( (x,i) => -Fs/2+i*Fs/N );
    // define a window
	const gauss = x => Math.exp(-Math.PI * ( x / Fs * N / 16.0 )**2 );
    // create STFT object
	const STFT = new st.SigmaTransform1D( x=>x , gauss , Fs , N , chans );
    // compute STFT of bat signal
	var coeff = STFT.analyze( bat ).coeff;
    // try to reconstruct (using the same windows, so no dual windows involved)
	var rec = STFT.synthesize().rec;

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
    //STFT.windows.map( v => { v.map(x=>{file.write( '' + (x.real**2+x.imag**2) +',' );});file.write('\n'); } );
    //file.end();
});



