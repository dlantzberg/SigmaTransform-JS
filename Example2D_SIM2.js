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
    // set 2D spectral diffeomorphism as a (exponential) polar map
    const sigma = xy => [ Math.log2(xy[0]**2+xy[1]**2 +.0000001)/2 , 
                          Math.atan(xy[1]/xy[0]) ];
    // create frequency channels
    const chans = st.meshgrid(  new Float32Array( 7 ).map( (x,i) => i ) ,
                                new Float32Array( 17 ).map( (x,i) => -3.14/2+i*3.14/16 ) );
    // define a rectangular window
	const win   = xy =>  (xy[0]>0)*(xy[0]<=1)
                       * (xy[1]>=0)*(xy[1]<=3.14/16);
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
    file = fs.createWriteStream('wins.asc');
    file.on('error', function(err) { console.log( "error", err ); });
    STFT.windows.map( v => { v.map(x=>{file.write( '' + (x.real**2+x.imag**2).toFixed(4) +',' );});file.write('\n'); } );
    file.end();
});



