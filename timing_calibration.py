import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np


#sum all counts from each PMT; collect 9 different data sets

#3 diff plots corresponding to 3 diff delay values on one sca while holding other constant

#do this three times on first sca

#plot delay vs counts with counts as sum of all counts of both PMTs

file_pattern = "timing_pmt1-0{}_pmt2-0{}.csv"
PMT2_delay = [0.1, 0.3, 0.5, 0.7]
PMT1_delay = [0.1, 0.2, 0.3]

#process data for each pmt1
def process_data(file_pattern, PMT1_num):
    counts_tot = []
    for i in range(4):
        #extract time
        data = pd.read_csv(file_pattern.format(PMT1_num, 2*i+1))
        time = int(str(data.iloc[0]).split(' ')[2])
        
        #extract counts columns and sum
        data = pd.read_csv(file_pattern.format(PMT1_num, 2*i+1), skiprows=1, header=None)
        counts_PMT1 = data.iloc[:, 1]; counts_PMT2 = data.iloc[:, 2]
        counts_tot.append((sum(counts_PMT1) + sum(counts_PMT2))/time)
        
    return counts_tot, max(counts_tot)
        
        
def plot_counts_delay (PMT1_delay, PMT2_delay, PMT1_num):
    
    counts, max_counts = process_data(file_pattern, PMT1_num)
    
    p = np.polyfit(PMT2_delay,counts,2)
    a=p[0]; b=p[1]; c=p[2]
    
    def parabola (x): return a*x**2 + b*x + c
    
    delay_range = np.linspace(PMT2_delay[0], PMT2_delay[-1], 100)
    
    plt.figure(figsize=(8, 6))
    plt.plot(PMT2_delay, counts, 'k.', label="Experimental Data")
    plt.plot(delay_range, parabola(delay_range),'b', label='C = {:.2f}d$^2$ + {:.2f}d + {:.2f}'.format(a,b,c))
    plt.title(f"Total Counts for PMT 1 delay: {PMT1_delay}")
    plt.xlabel("PMT 2 Delay (\u03BCs)"); plt.ylabel(r"Counts normalized by time (s$^{-1}$)")
    plt.legend()
    plt.savefig(f"timing_calibration_PMT1_delay_{PMT1_delay}.jpg", dpi=300)
    plt.show()
    
    print(f"PMT1 delay of {PMT1_delay}:")
    print(f"Maximum counts per second: {max_counts:.2f}")
    print(f"Polinomial fit: \nC = = {a:.2f}d$^2$ + {b:.2f}d + {c:.2f} \n")
    

#for i in range(3):
#    plot_counts_delay (PMT1_delay[i], PMT2_delay, i+1)
    
    

def plot_histogram(file_name):
    data = pd.read_csv(file_name, skiprows=1, header=None)

    # Extract voltage and counts columns
    voltage = data.iloc[:, 0]
    counts_PMT1 = data.iloc[:, 1]
    counts_PMT2 = data.iloc[:, 2]
    
    def PMT1_conv(v): return 129.0612079187347*v + 14.36301691501804
    def PMT2_conv(v): return 165.7468924263398*v + 13.375654389954574
    
    # Plot voltage histogram
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    axs[0].bar(voltage, counts_PMT1, width=0.05, color="blue", alpha=0.65)
    axs[0].set_ylabel("Counts")
    axs[0].set_title("PMT 1 Counts")
    
    axs[1].bar(voltage, counts_PMT2, width=0.05, color="blue", alpha=0.65)
    axs[1].set_ylabel("Counts")
    axs[1].set_title("PMT 2 Counts")
    
    axs[1].set_xlabel("Voltage (V)")

    plt.tight_layout()
    plt.savefig(f"Na_histogram_voltage", dpi=300)
    plt.show()
    #plt.close()
    
    # Plot energy histogram
    #voltage = np.array(voltage)
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    axs[0].bar(PMT1_conv(voltage), counts_PMT1, width=3, color="blue", alpha=0.65)
    axs[0].set_ylabel("Counts")
    axs[0].set_title("PMT 1 Counts")
    
    axs[1].bar(PMT2_conv(voltage), counts_PMT2, width=5, color="blue", alpha=0.65)
    axs[1].set_ylabel("Counts")
    axs[1].set_title("PMT 2 Counts")
    
    axs[1].set_xlabel("Energy (keV)")

    plt.tight_layout()
    plt.savefig(f"Na_histogram_energy", dpi=300)
    plt.show()
    #plt.close()
    
plot_histogram('10min_pmt1-01_pmt2-03_data.csv')
    
