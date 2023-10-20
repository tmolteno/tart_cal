'''
    GNSS Satellite acquisition stuff
    
    Author Tim Molteno. 2011-2019 tim@elec.ac.nz
    
    # aptitude install libfftw3-dev
    # pip install pyfftw
'''

import time
import numpy as np
import pyfftw
import scipy
from scipy import optimize

def generateCAcode(PRN):
    # A lookup table from lecture notes from the danish GPS center
    code_delay_table = [5,6,7,8,17,18,139,140,141,251,252,254,255,256,257,258,469,470,471,472,473,474,509,512,513,514,515,516,859,860,861,862,863,950,947,948,950];
    #    g2shift = circular shift of G2 maximal length code relative to the
    #    G1 maximal length code (must be an integer in the range 0:1023)
    g2shift = code_delay_table[PRN - 1];

    lfsr = -1*np.ones(10);

    g1 = np.empty(1023)
    #--- Generate all G1 signal chips based on the G1 feedback polynomial -----
    for i in range(1023):
        g1[i] = lfsr[9]
        saveBit = lfsr[2]*lfsr[9]
        lfsr = np.concatenate(([saveBit], lfsr[0:9]))

    #--- Generate G2 code -----------------------------------------------------

    #--- Initialize g2 output to speed up the function ---
    g2 = np.empty(1023)
    #--- Load shift register ---
    lfsr = -1*np.ones(10);

    #--- Generate all G2 signal chips based on the G2 feedback polynomial -----
    for i in range(1023):
        g2[i] = lfsr[9]
        saveBit = lfsr[1]*lfsr[2]*lfsr[5]*lfsr[7]*lfsr[8]*lfsr[9];
        lfsr = np.concatenate(([saveBit], lfsr[0:9]))

    #--- Shift G2 code --------------------------------------------------------
    g2 = np.roll(g2, g2shift);

    #--- Form single sample C/A code by multiplying G1 and G2 -----------------
    CAcode = -(g1 * g2);
    return CAcode


def gold(samples_per_code, PRN, epochs):
    samples_per_chip = samples_per_code/1023.0;
    CAcode = generateCAcode(PRN);
    code_samples = np.arange(np.floor(samples_per_code*epochs));	# An array of sample indexes
    code_indices = np.floor(code_samples / samples_per_chip).astype('int') % 1023;
    code = CAcode[code_indices];
    return code

def peak_func(x, a, b, c):
    return a*np.exp(-b*(x-c)**2)

def residuals(p, y, x):
    a,b,c = p
    err = y-peak_func(x,a,b,c)
    return err

def peak_fit(xdata,ydata,p0):
    return scipy.optimize.leastsq(residuals, p0, args=(ydata, xdata))

# FFTW variables
fft_in = 0
fft_out = 0
fft_machine = 0
ifft_in = 0
ifft_out = 0
ifft_machine = 0

fft2d_in = 0
fft2d_out = 0
ifft2d_in = 0
ifft2d_out = 0
fft2d_machine = 0
ifft2d_machine = 0


def acquire(x, sampling_freq, center_freq, searchBand, PRN, debug=False):

    pyfftw.interfaces.cache.enable()
    pyfftw.interfaces.cache.set_keepalive_time(50.0)

    sampling_period = 1.0/sampling_freq
    samples_per_ms = sampling_freq/1000.0    # Samples per millisecond
    samples_per_chip = samples_per_ms / 1023
    samples_per_chunk = int(samples_per_ms)

    epochs_available = np.floor(np.size(x)/samples_per_ms)

    freqBinSize= 0.5e3
    numberOfFrqBins = int(round(2*searchBand/freqBinSize)) + 1
    fc = np.linspace(center_freq - searchBand, center_freq + searchBand, numberOfFrqBins)

    # Construct variables for use of FFTW
    align = pyfftw.simd_alignment
    dtype = 'complex64'
    #dtype = 'complex128'
    global fft_in, fft_out, fft_machine, ifft_in, ifft_out, ifft_machine
    global fft2d_in, fft2d_out, fft2d_machine, ifft2d_in, ifft2d_out, ifft2d_machine


    write_wisdom=False
    try:
        import pickle
        wisdom = pickle.load( open( "wisdom.wis", "rb" ) )
        pyfftw.import_wisdom(wisdom)
    except:
        write_wisdom = True
        print('no wisdom file')

    fft_in = pyfftw.empty_aligned(samples_per_chunk, dtype=dtype, n=align)
    fft_out = pyfftw.empty_aligned(samples_per_chunk, dtype=dtype, n=align)
    ifft_in = pyfftw.empty_aligned(samples_per_chunk, dtype=dtype, n=align)
    ifft_out = pyfftw.empty_aligned(samples_per_chunk, dtype=dtype, n=align)

    start = time.time()

    fft_machine = pyfftw.FFTW(fft_in, fft_out, flags=('FFTW_MEASURE',))
    ifft_machine = pyfftw.FFTW(ifft_in, ifft_out, direction='FFTW_BACKWARD', flags=('FFTW_MEASURE',))

    fft2d_in = pyfftw.empty_aligned((numberOfFrqBins,samples_per_chunk), dtype=dtype, n=align)
    fft2d_out = pyfftw.empty_aligned((numberOfFrqBins,samples_per_chunk), dtype=dtype, n=align)
    ifft2d_in = pyfftw.empty_aligned((numberOfFrqBins,samples_per_chunk), dtype=dtype, n=align)
    ifft2d_out = pyfftw.empty_aligned((numberOfFrqBins,samples_per_chunk), dtype=dtype, n=align)
    fft2d_machine = pyfftw.FFTW(fft2d_in, fft2d_out, flags=('FFTW_MEASURE',))
    ifft2d_machine = pyfftw.FFTW(ifft2d_in, ifft2d_out, direction='FFTW_BACKWARD', flags=('FFTW_MEASURE',))

    if debug:
        print('setup took', time.time()-start)

    if write_wisdom:
        wisdom = pyfftw.export_wisdom()
        import pickle
        pickle.dump( wisdom, open( "wisdom.wis", "wb" ) )


    code_samples = np.arange(samples_per_chunk)
    # Generate a local signal sampled at the right sampling rate, and no phase change.
    start = time.time()
    code = gold(samples_per_ms, PRN, 1)
    if debug:
        print('gold code gen took', time.time()-start)

    start = time.time()
    codefreq = np.conj(np.fft.fft(code))

    phase_const = 2.0j*np.pi*sampling_period
    phasepoints = code_samples * phase_const
    if debug:
        print('gold code numpyfft took', time.time()-start)


    xcorr = np.zeros((numberOfFrqBins, samples_per_chunk)) # arguments need to be type int.
    best_epoch = -1
    best_xcorr = -1
    best_sn0 = -1
    start = time.time()

    if debug:
        print('num epochs:', epochs_available)
    for epoch in range(0, int(epochs_available)):
        # Get the correspi
        start_index = int(epoch*samples_per_ms) # int!    Changed by max
        stop_index = start_index + int(samples_per_ms) # int! Changed by max
        epoch_data = x[start_index:stop_index]
        epoch_xcorr = acquire_aux(epoch_data, sampling_freq, fc, numberOfFrqBins, PRN, code_samples, codefreq, phasepoints)
        # epoch_xcorr2 = acquire_aux2D(epoch_data, sampling_freq, fc, numberOfFrqBins, PRN, code_samples, codefreq, phasepoints)
        # epoch_xcorr = epoch_xcorr2
        #print '1d vs 2D', (epoch_xcorr2-epoch_xcorr).__eq__(0).all()

        sn0 = epoch_xcorr.max()#/epoch_xcorr.std()
        if sn0 > best_sn0:
        # if (epoch_xcorr.max() > best_xcorr):
            best_sn0 = sn0
            if debug:
                print('better sn0 for: ',PRN , sn0, '>', best_sn0)
            # print epoch_xcorr.max(), '>', best_xcorr
            best_epoch_xcorr = epoch_xcorr
            # best_xcorr = epoch_xcorr.max()
            # best_epoch = epoch
            best_data = epoch_data
        # xcorr = xcorr + epoch_xcorr

    if debug:
        print('corr took', time.time()-start)

    
    signal_strength = best_sn0 #best_xcorr/best_epoch_xcorr.std()
    # Calculate the average correlation
    result = best_epoch_xcorr

    frequencies = result.max(1)
    freq = frequencies.argmax()

    phases = result.max(0)
    codephase = phases.argmax() % int(samples_per_ms)
    codephase_frac = codephase / samples_per_ms

    # frequency = fc[freq] - center_freq
    # print 'PRN %02i, freqmax %3.3e phase %3.3e' % (PRN, fc[freq], codephase)


    #frequencies = best_epoch_xcorr.max(1)
    #freq = frequencies.argmax()
    frequency = fc[freq] - center_freq     # This is the doppler
    
    #phases = best_epoch_xcorr.max(0)
    #codephase = phases.argmax() % int(samples_per_ms)
    #codephase_frac = codephase / samples_per_ms


    #autocorrelation = result[freq,:]

    peak = result.max()
    #sd = result.std()
    ## signal_strength_sd = ((peak - result.mean()) / sd)
    #signal_strength_sd = (peak/ sd)

    offfreq = (freq + numberOfFrqBins//2) % numberOfFrqBins
    other_peak = result[offfreq].max()
    signal_strength_off = peak / other_peak

    #blank_start = max(int(codephase - 3*samples_per_chip), 0)
    #blank_end = min(int(codephase + 3*samples_per_chip), int(samples_per_ms))
    #autocorrelation[blank_start:blank_end] = 0
    #secondPeakSize = autocorrelation.max()
    #signal_strength_ratio = (peak / secondPeakSize)

    #signal_strength = signal_strength_sd

    start = time.time()
    try:
        #print("[%02d] best_xcorr=%f peak=%f sd=%f off=%f ratio=%f" % (PRN, best_xcorr, peak, signal_strength_sd, signal_strength_off, signal_strength_ratio))
        if ((signal_strength_off / 1.3) > 1.0):
            [codephase_frac_o, frequency_o, peak_o] = optimize_fit(PRN, best_data, 1, samples_per_ms, sampling_period, fc, freq, frequency, center_freq, codephase)
            #[codephase_frac_o, frequency_o, peak_o] = optimize_fit(PRN, x, epochs_available, samples_per_ms, sampling_period, fc, freq, frequency,center_freq, codephase)
            if (peak_o > best_xcorr):
                codephase_frac = codephase_frac_o
                frequency = frequency_o
                peak = peak_o
    except Exception as e:
        print(f"Optimization Failed {e.stacktrace()}")
    if debug:
        print(f"Optimization took {time.time()-start}")

    # assert type(frequency) is np.float64, "frequency is not an float: %r" % frequency
    return [PRN, signal_strength, codephase_frac, frequency]


def optimize_fit(PRN, x, epochs_available, samples_per_ms, sampling_period, fc, freq, frequency, center_freq, codephase):
    num_samples = int(epochs_available*samples_per_ms) # int! changed by max
    phase_const = 2.0j*np.pi*sampling_period;
    code_samples = np.linspace(0,num_samples-1,num_samples)
    signal = x[0:num_samples]
    # Generate a local signal sampled at the right sampling rate, and no phase change.
    code = gold(samples_per_ms, PRN, epochs_available)

    codefreq = np.conj(pyfftw.interfaces.numpy_fft.fft(code))

    phasepoints = code_samples * phase_const
    #print "Num Samples %d" % num_samples
    #print "Phasepoints %d" % np.size(phasepoints)
    #print "code %d" % np.size(code)
    #print "codefreq %d" % np.size(codefreq)
    #print "signal %d" % np.size(signal)

    func = lambda f: (1.0 / correlate_aux(f, signal, phasepoints, codefreq).max())
    f_max = scipy.optimize.fminbound(func,fc[max(freq-1,0)],fc[min(freq+1, len(fc)-1)], xtol=1e-2)
    cp_peak = 1.0 / func(f_max)

    #print("        [%02d] Freq Optim %d -> %f (peak %f)" % (PRN, frequency, f_max - center_freq, cp_peak))
    frequency = f_max - center_freq
    # Find the new codephase at this frequency
    autocorrelation = np.abs(pyfftw.interfaces.numpy_fft.ifft(pyfftw.interfaces.numpy_fft.fft(np.exp(phasepoints*f_max) * signal) * codefreq)) / np.sqrt(num_samples)
    cp_max = autocorrelation.argmax()
    cp_peak = autocorrelation.max()
    #print("        [%02d] Codephase Optim %d -> %d (peak %f)" % (PRN, codephase, cp_max % samples_per_ms, cp_peak))
    codephase_frac = cp_max / samples_per_ms;
    # Now fit a gaussian to the codephase.
    cw = 2
    phases = np.arange(codephase-cw,codephase+cw)
    values = autocorrelation[codephase-cw:codephase+cw]
    popt, success = peak_fit(phases, values, [values.max(), 1e-1, codephase])
    #print(popt, success)
    if (success <= 4.0):
        #print("        [%02d] Codephase Optim %d -> %f" % (PRN, codephase, popt[2]))
        codephase_frac = popt[2] / samples_per_ms

    return [codephase_frac, frequency, cp_peak]

def correlate_aux(frequency, signal, phasepoints, codefreq):
    expfreq = np.exp(phasepoints*frequency)

    global fft_in, fft_machine, ifft_in, ifft_machine
    fft_in[:] = expfreq * signal
    IQfreq1 = fft_machine()

    ifft_in[:] = IQfreq1 * codefreq
    ifft_result = ifft_machine()

    corr = np.abs(ifft_result) / np.sqrt(len(signal))
    return corr

def acquire_aux(x, sampling_freq, fc, numberOfFrqBins, PRN, code_samples, codefreq, phasepoints):
    samples_per_ms = sampling_freq/1000.0    # Samples per millisecond
    samples_per_chip = samples_per_ms / 1023
    epochs_available = np.floor(np.size(x)/samples_per_ms)
    signal1 = x[code_samples]
    samples_per_chunk = int(samples_per_ms)

    result1 = np.empty((numberOfFrqBins, samples_per_chunk))
    for i in range(0, numberOfFrqBins):
        acqRes1 = correlate_aux(fc[i], signal1, phasepoints, codefreq)
        result1[i,:] = acqRes1
    return result1


def acquire_aux2D(x, sampling_freq, fc, numberOfFrqBins, PRN, code_samples, codefreq, phasepoints):
    global fft2d_in, fft2d_machine, ifft2d_in, ifft2d_machine
    fft2d_in[:] = np.exp(np.outer(fc, phasepoints)) * x
    # generating array takes twice as long as rest of this function... (fc,phasepoints,x)
    IQfreq1 = fft2d_machine()
    ifft2d_in[:] = IQfreq1 * codefreq
    ifft2d_result = ifft2d_machine()
    corr = np.abs(ifft2d_result) / np.sqrt(len(x))
    return corr
