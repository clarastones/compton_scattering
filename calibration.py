import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np

def process_data(file_pattern, PMT_num, prom, width):
    peaks_all = []; fwhm_all = []
    for i in range(len(sources)):
        data = pd.read_csv(file_pattern.format(sources[i]), skiprows=1, header=None)

        # Extract voltage and counts columns
        voltage = data.iloc[:, 0]
        counts = data.iloc[:, PMT_num]
        
        #find peak(s)
        peaks, _ = find_peaks(counts, prominence=prom, width = width)
        
        if PMT_num == 2:
            if sources[i] in ('Cs', 'Am'):
                peaks = peaks[1:]
            if sources[i] == 'Ba':
                peaks = peaks[1:-1]
        
        half_max = []; x_left_interp = []; x_right_interp = []; fwhm = []
        for j in range(len(peaks)):
            peak_index = peaks[j]
            HM = counts[peak_index]/2

            # Find nearest points where y crosses half max
            x_left_interp.append(voltage[np.where(counts[:peak_index] < HM)[0][-1]])
            x_right_interp.append(voltage[np.where(counts[peak_index:] < HM)[0][0]+peak_index])
        
            fwhm.append(x_right_interp[j] - x_left_interp[j]); half_max.append(HM)
            
        #counts_rev = counts[::-1]
        #peaks_rev, _ = find_peaks(counts_rev, prominence=125, width = 2)
        #peaks_rev = [len(counts) - p for p in peaks_rev]
        #peaks_avg = [(peaks[i]+peaks_rev[i])/2 for i in range(len(peaks))]
        
        # Plot histogram
        plt.figure(figsize=(8, 6))
        plt.bar(voltage, counts, width=0.05, color='b', alpha=0.65)
        for j in range(len(peaks)):
            plt.axvline(x=voltage[peaks[j]], color='k', linestyle='--', label='Forward Peaks')
            plt.hlines(half_max[j], xmin=x_left_interp[j], xmax=x_right_interp[j], color='r', linestyle='dashed')
                
        #for k, p in enumerate(voltage[peaks_rev]):
            #if k == 0:
                #plt.axvline(x=p, color='r', linestyle='--', label='Forward Peaks')
            #else:
                #plt.axvline(x=p, color='r', linestyle='--')
                
        #plt.axvline(counts[peaks], "r", linestyle='--')
        plt.xlabel("Voltage (V)")
        plt.ylabel("Counts")
        plt.title(f"Voltage vs Count for PMT {PMT_num}: {sources[i]}")
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.savefig(f"calibration_histogram_PMT{PMT_num}_{sources[i]}", dpi=300)
        #plt.show()
        plt.close()
        peaks_all.append(list(voltage[peaks])); fwhm_all.append(list(fwhm))
        
    return peaks_all, fwhm_all
        
def calibration (file_pattern, coarse_gain, PMT_num, prom, width):
    peaks, fwhm = process_data(file_pattern, PMT_num, prom, width)

    #error in E is slope*fwhm
    if PMT_num == 2:
        exp_peaks[0] = exp_peaks[0][1:-1]

    v_peaks = []; E_peaks = []; FWHM = []
    #scale data by gain
    for i in range(len(peaks)):
        for j in range(len(peaks[i])):
            v_peaks.append(peaks[i][j]*10/coarse_gain[i])
            E_peaks.append(exp_peaks[i][j])
            FWHM.append(fwhm[i][j])
    

    v_peaks = np.array(v_peaks); E_peaks = np.array(E_peaks); FWHM = np.array(FWHM)

    p,cov = np.polyfit(v_peaks,E_peaks,1, cov=True)
    a=p[0]; b=p[1]
    sig_a = np.sqrt(cov[0][0]); sig_b = np.sqrt(cov[1][1])
    
    E_err = FWHM*a

    E_est = a * v_peaks + b

    plt.figure(figsize=(8, 6))
    plt.errorbar(v_peaks, E_peaks, xerr=FWHM, yerr=E_err, fmt='.', label="Experimental Data")
    plt.plot(v_peaks, E_est, 'r', label= 'E = ({:.2e})V + {:.2e}'.format(a, b))
    plt.title(f"Voltage-Energy Conversion for PMT {PMT_num}")
    plt.xlabel("Voltage (V)"); plt.ylabel("Energy (keV)")
    plt.legend()
    plt.savefig(f"calibration_PMT{PMT_num}", dpi=300)
    plt.show()
    
    sig_E = abs(b)*sig_a
    
    return a, b, sig_a, sig_b, FWHM, E_err, sig_E
    


sources = ["Ba", "Cs", "Na", "Am"]
# average Ba 2nd peak: 289.624
#first peak: 80.998 [80.998, 160.59, 223.22, 276.397]
exp_peaks = [[53.12, 80.998, 289.624], [661.638], [511.006], [59.537]]
coarse_gain_PMT1 = [40,10,10,10]
coarse_gain_PMT2 = [10,10,10,10]

file_patterns = ["calibration_{}.csv", "calibration_{}_PMT2.csv"]; gains = [coarse_gain_PMT1, coarse_gain_PMT2]
proms = [125, 400]; widths = [5, 1]
for i in range(2):
    slope, intercept, slope_err, intercept_err, fwhm, yerr, E_err = calibration(file_patterns[i], gains[i], i+1, proms[i], widths[i])
    print(f"PMT{i+1} conversion: ({slope} +/- {slope_err}) keV/V + ({intercept} +/- {intercept_err}) keV")
    print(f"FWHM: {fwhm} and error in Y: {yerr}")
    print(f"error in E: {E_err}")