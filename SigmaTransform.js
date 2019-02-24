'use strict';
class baseComplexArray {
    constructor(other, arrayType = Float32Array) {
        if (other instanceof ComplexArray) {
            // Copy constuctor.
            this.ArrayType = other.ArrayType;
            this.real = new this.ArrayType(other.real);
            this.imag = new this.ArrayType(other.imag);
        } else {
            this.ArrayType = arrayType;
            // other can be either an array or a number.
            this.real = new this.ArrayType(other);
            this.imag = new this.ArrayType(this.real.length);
        }
        this.length = this.real.length;
    }

    toString() {
        const components = [];
        this.forEach((value, i) => {
            components.push(
                `(${value.real.toFixed(2)}, ${value.imag.toFixed(2)})`
            );
        });
        return `[${components.join(', ')}]`;
    }

    forEach(iterator) {
        const n = this.length;
        // For gc efficiency, re-use a single object in the iterator.
        const value = Object.seal(Object.defineProperties({}, {
            real: {writable: true}, imag: {writable: true},
        }));

        for (let i = 0; i < n; i++) {
            value.real = this.real[i];
            value.imag = this.imag[i];
            iterator(value, i, n);
        }
    }

    // In-place mapper.
    map(mapper) {
        this.forEach((value, i, n) => {
            mapper(value, i, n);
            this.real[i] = value.real;
            this.imag[i] = value.imag;
        });

        return this;
    }

    conjugate() {
        return new baseComplexArray(this).map((value) => {
            value.imag *= -1;
        });
    }

    magnitude() {
        const mags = new this.ArrayType(this.length);

        this.forEach((value, i) => {
            mags[i] = Math.sqrt(value.real*value.real + value.imag*value.imag);
        })

        return mags;
    }
    magsq() {
        const mags = new this.ArrayType(this.length);

        this.forEach((value, i) => {
            mags[i] = value.real*value.real + value.imag*value.imag;
        })

        return mags;
    }
}

// Math constants and functions we need.
const PI = Math.PI;
const SQRT1_2 = Math.SQRT1_2;

function FFT(input) {
    return ensureComplexArray(input).FFT();
};

function InvFFT(input) {
    return ensureComplexArray(input).InvFFT();
};

function frequencyMap(input, filterer) {
    return ensureComplexArray(input).frequencyMap(filterer);
};

class ComplexArray extends baseComplexArray {
    FFT() {
        return fft(this, false);
    }

    InvFFT() {
        return fft(this, true);
    }

    // Applies a frequency-space filter to input, and returns the real-space
    // filtered input.
    // filterer accepts freq, i, n and modifies freq.real and freq.imag.
    frequencyMap(filterer) {
        return this.FFT().map(filterer).InvFFT();
    }
}

function ensureComplexArray(input) {
    return input instanceof ComplexArray && input || new ComplexArray(input);
}

function fft(input, inverse) {
    const n = input.length;

    if (n & (n - 1)) {
        return FFT_Recursive(input, inverse);
    } else {
        return FFT_2_Iterative(input, inverse);
    }
}

function FFT_Recursive(input, inverse) {
    const n = input.length;

    if (n === 1) {
        return input;
    }

    const output = new ComplexArray(n, input.ArrayType);

    // Use the lowest odd factor, so we are able to use FFT_2_Iterative in the
    // recursive transforms optimally.
    const p = LowestOddFactor(n);
    const m = n / p;
    const normalisation = 1 / Math.sqrt(p);
    let recursive_result = new ComplexArray(m, input.ArrayType);

    // Loops go like O(n Σ p_i), where p_i are the prime factors of n.
    // for a power of a prime, p, this reduces to O(n p log_p n)
    for(let j = 0; j < p; j++) {
        for(let i = 0; i < m; i++) {
            recursive_result.real[i] = input.real[i * p + j];
            recursive_result.imag[i] = input.imag[i * p + j];
        }
        // Don't go deeper unless necessary to save allocs.
        if (m > 1) {
            recursive_result = fft(recursive_result, inverse);
        }

        const del_f_r = Math.cos(2*PI*j/n);
        const del_f_i = (inverse ? -1 : 1) * Math.sin(2*PI*j/n);
        let f_r = 1;
        let f_i = 0;

        for(let i = 0; i < n; i++) {
            const _real = recursive_result.real[i % m];
            const _imag = recursive_result.imag[i % m];

            output.real[i] += f_r * _real - f_i * _imag;
            output.imag[i] += f_r * _imag + f_i * _real;

            [f_r, f_i] = [
                f_r * del_f_r - f_i * del_f_i,
                f_i = f_r * del_f_i + f_i * del_f_r,
            ];
        }
    }

    // Copy back to input to match FFT_2_Iterative in-placeness
    // TODO: faster way of making this in-place?
    for(let i = 0; i < n; i++) {
        input.real[i] = normalisation * output.real[i];
        input.imag[i] = normalisation * output.imag[i];
    }

    return input;
}

function FFT_2_Iterative(input, inverse) {
    const n = input.length;

    const output = BitReverseComplexArray(input);
    const output_r = output.real;
    const output_i = output.imag;
    // Loops go like O(n log n):
    //   width ~ log n; i,j ~ n
    let width = 1;
    while (width < n) {
        const del_f_r = Math.cos(PI/width);
        const del_f_i = (inverse ? -1 : 1) * Math.sin(PI/width);
        for (let i = 0; i < n/(2*width); i++) {
            let f_r = 1;
            let f_i = 0;
            for (let j = 0; j < width; j++) {
                const l_index = 2*i*width + j;
                const r_index = l_index + width;

                const left_r = output_r[l_index];
                const left_i = output_i[l_index];
                const right_r = f_r * output_r[r_index] - f_i * output_i[r_index];
                const right_i = f_i * output_r[r_index] + f_r * output_i[r_index];

                output_r[l_index] = SQRT1_2 * (left_r + right_r);
                output_i[l_index] = SQRT1_2 * (left_i + right_i);
                output_r[r_index] = SQRT1_2 * (left_r - right_r);
                output_i[r_index] = SQRT1_2 * (left_i - right_i);

                [f_r, f_i] = [
                    f_r * del_f_r - f_i * del_f_i,
                    f_r * del_f_i + f_i * del_f_r,
                ];
            }
        }
        width <<= 1;
    }

    return output;
}

function BitReverseIndex(index, n) {
    let bitreversed_index = 0;

    while (n > 1) {
        bitreversed_index <<= 1;
        bitreversed_index += index & 1;
        index >>= 1;
        n >>= 1;
    }
    return bitreversed_index;
}

function BitReverseComplexArray(array) {
    const n = array.length;
    const flips = new Set();

    for(let i = 0; i < n; i++) {
        const r_i = BitReverseIndex(i, n);

        if (flips.has(i)) continue;

        [array.real[i], array.real[r_i]] = [array.real[r_i], array.real[i]];
        [array.imag[i], array.imag[r_i]] = [array.imag[r_i], array.imag[i]];

        flips.add(r_i);
    }

    return array;
}

function LowestOddFactor(n) {
    const sqrt_n = Math.sqrt(n);
    let factor = 3;

    while(factor <= sqrt_n) {
        if (n % factor === 0) return factor;
        factor += 2;
    }
    return n;
}


function FFTImageDataRGBA(data, nx, ny) {
    const rgb = splitRGB(data);

    return mergeRGB(
        FFT2D(new ComplexArray(rgb[0], Float32Array), nx, ny),
        FFT2D(new ComplexArray(rgb[1], Float32Array), nx, ny),
        FFT2D(new ComplexArray(rgb[2], Float32Array), nx, ny)
    );
};

function splitRGB(data) {
    const n = data.length / 4;
    const r = new Uint8ClampedArray(n);
    const g = new Uint8ClampedArray(n);
    const b = new Uint8ClampedArray(n);

    for(let i = 0; i < n; i++) {
        r[i] = data[4 * i    ];
        g[i] = data[4 * i + 1];
        b[i] = data[4 * i + 2];
    }

    return [r, g, b];
}

function gray2RGB(data) {
    const n = data.length;
    if( data instanceof ComplexArray ) {
        const output = new ComplexArray(n * 4);

        for(let i = 0; i < n; i++) {
            output.real[4 * i    ] = data.real[i];
            output.imag[4 * i    ] = data.imag[i];
            output.real[4 * i + 1] = data.real[i];
            output.imag[4 * i + 1] = data.imag[i];
            output.real[4 * i + 2] = data.real[i];
            output.imag[4 * i + 2] = data.imag[i];
            output.real[4 * i + 3] = 255;
            output.imag[4 * i + 3] = 255;
        }
        return output;
    } else {
        const output = new Uint8ClampedArray(n * 4);

        for(let i = 0; i < n; i++) {
            output[4 * i    ] = data[i];
            output[4 * i + 1] = data[i];
            output[4 * i + 2] = data[i];
            output[4 * i + 3] = 255;
        }
        return output;
    }
}


function RGB2gray(data) {
    const n = data.length / 4,
        G = new Uint8ClampedArray(n);
    for(let i = 0; i < n; i++)
        G[i]  = (data[4 * i] + data[4 * i + 1] + data[4 * i + 2])/3;
    return G;
}


function mergeRGB(r, g, b) {
    const n = r.length;
    const output = new ComplexArray(n * 4);

    for(let i = 0; i < n; i++) {
        output.real[4 * i    ] = r.real[i];
        output.imag[4 * i    ] = r.imag[i];
        output.real[4 * i + 1] = g.real[i];
        output.imag[4 * i + 1] = g.imag[i];
        output.real[4 * i + 2] = b.real[i];
        output.imag[4 * i + 2] = b.imag[i];
    }

    return output;
}
function mergeRB(r, b) {
    const n = r.length;
    const output = new ComplexArray(n * 4);

    for(let i = 0; i < n; i++) {
        output.real[4 * i    ] = r.real[i];
        output.imag[4 * i    ] = r.imag[i];
        output.real[4 * i + 1] = 0;
        output.imag[4 * i + 1] = 0;
        output.real[4 * i + 2] = b.real[i];
        output.imag[4 * i + 2] = b.imag[i];
    }

    return output;
}

function FFT2D(input, nx, ny, inverse) {
    const transform = inverse ? 'InvFFT' : 'FFT';
    const output = new ComplexArray(input.length, input.ArrayType);
    const row = new ComplexArray(nx, input.ArrayType);
    const col = new ComplexArray(ny, input.ArrayType);

    for(let j = 0; j < ny; j++) {
        row.map((v, i) => {
            v.real = input.real[i + j * nx];
            v.imag = input.imag[i + j * nx];
        });
        row[transform]().forEach((v, i) => {
            output.real[i + j * nx] = v.real;
            output.imag[i + j * nx] = v.imag;
        });
    }

    for(let i = 0; i < nx; i++) {
        col.map((v, j) => {
            v.real = output.real[i + j * nx];
            v.imag = output.imag[i + j * nx];
        });
        col[transform]().forEach((v, j) => {
            output.real[i + j * nx] = v.real;
            output.imag[i + j * nx] = v.imag;
        });
    }

    return output;
}

class SigmaTransform1D {
	constructor( sigma,win,Fs,size,chans ) {
        // store data
		this.win = win;
		this.sigma = sigma;
		this.Fs  = Fs;
		this.N = size;
		this.chans = chans;
        this.numchan = chans.length;
        // make the spectra of windows - could be post-poned
        this.makeWindows();
	}
	analyze(signal) {
        // get spectrum of signal
		const F = new ComplexArray( signal ).FFT();
        this.coeff = [];
		for(var c=0;c<this.numchan;c++) {
            // get new space for one coefficient line
            this.coeff.push( new ComplexArray( this.N ) );
            // calc coefficient
            this.coeff[c].map( (x,i) => {
                // windows in this implementation are always real, so no need for conjugation..
                x.real = F.real[i] * this.windows[c].real[i];
                x.imag = F.imag[i] * this.windows[c].real[i];
            }).InvFFT();
		}
        // return object for chainable api
		return this;
	}
	synthesize() {
        this.rec = new ComplexArray( this.N );
        for(var c=0;c<this.numchan;c++) {
            const C = new ComplexArray( this.coeff[c] ).FFT();
            this.rec.map( (x,i) => {
                x.real += C.real[i]*this.windows[c].real[i];
                x.imag += C.imag[i]*this.windows[c].real[i];
            });
        }
		// ifft
        this.rec.InvFFT();
        // return object for chainable api
		return this;
	}
	makeWindows() {
        var dx = this.Fs/this.N;
        var curr = this.Fs/2 - dx;
        // make frequencies
        var warpDom = new Float32Array( this.N ).map(x => { curr += dx; return this.sigma(( curr % this.Fs ) - this.Fs/2); } );
        this.windows = [];
		for(var c=0;c<this.numchan;c++) {
            this.windows.push( new ComplexArray( this.N ) );
            // currently only real windows ...
            this.windows[c].map( (x,i) => {
                x.real = this.win( warpDom[i] - this.chans[c] );
//              x.imag = 0;
            } );
		}
	}
}


class SigmaTransform2D {
    constructor( sigma,win,Fs,size,chans ) {
        // store data
        this.win = win;
        this.sigma = sigma;
        this.Fs  = Fs;
        this.N = size;
        this.chans = chans;
        this.numchan = chans.length;
        // could be post-poned
        this.makeWindows();
    }
    analyze(signal) {
        const F = FFT2D( new ComplexArray( signal.map(x=>isNaN(x)?0:x)),this.N[0],this.N[1],false );
        this.coeff = [];
        for(var c=0;c<this.numchan;c++) {
            this.coeff.push( new ComplexArray( this.N[0]*this.N[1] ) );
            this.coeff[c] = FFT2D( this.coeff[c].map( (x,i) => {
                // windows in this implementation are always real, so no need for conjugation..
                x.real = F.real[i] * this.windows[c].real[i];
                x.imag = F.imag[i] * this.windows[c].real[i];
            } ) , this.N[0] , this.N[1] , true );
        }
        // return object for chainable api
        return this;
    }
    synthesize() {
        this.rec = new ComplexArray( this.N[0]*this.N[1] );
        for(var c=0;c<this.numchan;c++) {
            const C = FFT2D( new ComplexArray( this.coeff[c] ),this.N[0],this.N[1],false );
            this.rec.map( (x,i) => {
                x.real += C.real[i]*this.windows[c].real[i];
                x.imag += C.imag[i]*this.windows[c].real[i];
            } );
        }
        // ifft
        this.rec = FFT2D( this.rec , this.N[0] , this.N[1] , true);
        return this;
    }
    makeWindows() {
        var dx   = [ this.Fs[0]/this.N[0] , this.Fs[1]/this.N[1] ];
        var curr = [ this.Fs[0]/2 - dx[0] , this.Fs[1]/2 - dx[1] ];
        // make frequencies
        var ax0 = new Float32Array( this.N[0] ).map(x => { 
            curr[0] += dx[0]; return ( curr[0] % this.Fs[0] ) - this.Fs[0]/2; 
        });
        var ax1 = new Float32Array( this.N[1] ).map(x => { 
            curr[1] += dx[1]; return ( curr[1] % this.Fs[1] ) - this.Fs[1]/2; 
        });
        // meshgrid and warp the frequencies
        var warpDom = meshgrid( ax0 , ax1 ).map(x=>this.sigma(x));
        // create spectra of windows
        this.windows = [];
        for(var c=0;c<this.numchan;c++) {
            this.windows.push( new ComplexArray( this.N[0]*this.N[1] ) );
            // currently only real windows ...
            this.windows[c].map( (x,i) => {
                x.real = this.win( [ warpDom[i][0] - this.chans[c][0],
                                     warpDom[i][1] - this.chans[c][1] ] );
                //              x.imag = 0;
            } );
        }
    }
}

function meshgrid( X , Y ) {
//    console.log( "Los: ", Array.from(X).map( x => [x,2] ) );
    return Array.from(X).map( x => Array.from(Y).map( y => [x,y] ) ).flat();
}

module.exports.SigmaTransform1D = SigmaTransform1D;
module.exports.SigmaTransform2D = SigmaTransform2D;
module.exports.meshgrid         = meshgrid;
module.exports.ComplexArray     = ComplexArray;
