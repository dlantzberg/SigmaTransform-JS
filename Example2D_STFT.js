// get SigmaTransform class
const st = require('./SigmaTransform.js');
// get filesystem
const fs = require('fs');
// get lineReader
const lineReader = require('readline').createInterface({
	  input: fs.createReadStream('lena.asc')
});

// read lena signal
var lena = [];
lineReader.on('line', function (line) {
    line.split(' ').map( x => {
	    lena.push( parseFloat( x ) );
    } );
});

// done with reading?
lineReader.on('close', function() {
    // get length
	const N  = [128,128];
    // set sampling frequency
	const Fs = [128,128];
    // set 2D spectral diffeomorphism, (each x is an array: [x,y] => [x,y] )
    const sigma = xy => xy;
    // set stepwidth for channels
    const stepwidth = 5;
    // create frequency channels
//    console.log( new Float32Array( 5 ).map( (x,i) => -Fs[0]/2+i*stepwidth ) );

    const chans = st.meshgrid(  new Float32Array( 25 ).map( (x,i) => -Fs[0]/2+i*stepwidth ) ,
                                new Float32Array( 25 ).map( (x,i) => -Fs[1]/2+i*stepwidth ) ) ;
    // define a rectangular window
	const win   = xy =>  (xy[0]>=-2.5)*(xy[0]<2.5)
                       * (xy[1]>=-2.5)*(xy[1]<2.5);
    // create STFT object
	const STFT = new st.SigmaTransform2D( sigma , win , Fs , N , chans );
    // compute STFT of bat signal
	var coeff = STFT.analyze( lena ).coeff;
    // try to reconstruct (using the same windows, so no dual windows involved)
	var rec = STFT.synthesize().rec;

    // save coefficients to file
    //var file = fs.createWriteStream('coeff.asc');
    //file.cork();
    //file.on('error', function(err) { console.log( "error", err ); });
    //coeff.map( v => { v.map(x=>{file.write( '' + (x.real**2+x.imag**2).toFixed(4) +',' );});file.write('\n'); } );
    //file.end();

    // save reconstruction to file
    var file = fs.createWriteStream('rec.asc');
    file.on('error', function(err) { console.log( "error", err ); });
    rec.map( v => { file.write( ''+v.real.toFixed(4)+'\n' ); });
    file.end();

    // save spectra of windows to file
    //file = fs.createWriteStream('wins.asc');
    //file.on('error', function(err) { console.log( "error", err ); });
    //STFT.windows.map( v => { v.map(x=>{file.write( '' + (x.real**2+x.imag**2).toFixed(4) +',' );});file.write('\n'); } );
    //file.end();
});



