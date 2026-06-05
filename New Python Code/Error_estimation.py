#Error estimation
from Read_Fork_Data import ForkData
import numpy as np


R_file = "Data/R_CalFile_20-28kHz_2000points_100mVpk.txt"

air_files = [
    "Data/Y_error estimation/Air error/TF1_airsweep1.txt",
    "Data/Y_error estimation/Air error/TF1_airsweep2.txt",
    "Data/Y_error estimation/Air error/TF1_airsweep3.txt",
    "Data/Y_error estimation/Air error/TF1_airsweep4.txt",
    "Data/Y_error estimation/Air error/TF1_airsweep5.txt",
    "Data/Y_error estimation/Air error/TF1_airsweep6.txt",
    "Data/Y_error estimation/Air error/TF1_airsweep7.txt",
    "Data/Y_error estimation/Air error/TF1_airsweep8.txt",
    "Data/Y_error estimation/Air error/TF1_airsweep9.txt",
    "Data/Y_error estimation/Air error/TF1_airsweep10.txt"
]

fluid_files = [
    "Data/Y_error estimation/W1G1 error/TF1_fluidsweep1.txt",
    "Data/Y_error estimation/W1G1 error/TF1_fluidsweep2.txt",
    "Data/Y_error estimation/W1G1 error/TF1_fluidsweep3.txt",
    "Data/Y_error estimation/W1G1 error/TF1_fluidsweep4.txt",
    #"Data/Y_error estimation/W1G1 error/TF1_fluidsweep5.txt",
    "Data/Y_error estimation/W1G1 error/TF1_fluidsweep6.txt",
    "Data/Y_error estimation/W1G1 error/TF1_fluidsweep7.txt",
    "Data/Y_error estimation/W1G1 error/TF1_fluidsweep8.txt",
    "Data/Y_error estimation/W1G1 error/TF1_fluidsweep9.txt",
    "Data/Y_error estimation/W1G1 error/TF1_fluidsweep10.txt"
]
realY_air_array = []
imagY_air_array = []
C0_air_array = []

f0_air_array = []
fp_air_array = []
Q_air_array = []
for i in range(len(air_files)):
    TF1 = ForkData(air_files,[0.2 for j in range(len(air_files))],100.6e3,R_file,0.1)
    FreqData, X_data, Y_data = TF1.Read_XY_Voltage(air_files[i])
    Z_fork = TF1.CalibrationResistorValue * (TF1.DrivingVoltages[i]/TF1.CalibrationResistorVoltage) * \
            TF1.CalibrationVoltage / (X_data + 1j*Y_data)
    Y_fork = 1 / Z_fork
    #realY_air_array.append(np.array(np.real(Y_fork)))
    #imagY_air_array.append(np.array(np.imag(Y_fork)))
    C0 = TF1.parallel_capacitance_estimation(Y_fork,FreqData)
    C0_air_array.append(C0)
    Y_m = Y_fork - 1j*2*np.pi*np.array(FreqData)*C0
    realY_air_array.append(np.array(np.real(Y_m)))
    imagY_air_array.append(np.array(np.imag(Y_m)))
    
    res = TF1.LorentzianFit(FreqData,abs(Y_m))
    f0_air_array.append(res[0])
    Q_air_array.append(res[2])
    res = TF1.LorentzianFit(FreqData,abs(Z_fork))
    fp_air_array.append(res[0])

def stats(Y):
    mean_out = np.mean(Y, axis=0)
    std_out = np.std(Y, axis=0, ddof=1)
    stderr_out = std_out / np.sqrt(np.array(Y).shape[0])
    return mean_out, std_out, stderr_out

Re_air_mean, Re_air_std, Re_air_stderr = stats(realY_air_array)
Im_air_mean, Im_air_std, Im_air_stderr = stats(imagY_air_array)
C0_air_mean, C0_air_std, C0_air_stderr = stats(C0_air_array)

print("C0 in air:",C0_air_mean,"+-",C0_air_stderr)
f0_air_mean,f0_air_std,f0_air_stderr = stats(f0_air_array)
print("Resonant Frequency in air:",f0_air_mean,"+-",f0_air_stderr)
fp_air_mean,fp_air_std,fp_air_stderr = stats(fp_air_array)
print("Parallel-Resonant Frequency in air:",fp_air_mean,"+-",fp_air_stderr)
Q_air_mean,Q_air_std,Q_air_stderr = stats(Q_air_array)
print("Q-factor in air",Q_air_mean,"+-",Q_air_stderr)

imagY_air_array = imagY_air_array - 2*np.pi*FreqData*C0_air_mean
Y_mag_array = [np.abs(r + 1j*i) for r, i in zip(realY_air_array, imagY_air_array)]
Y_mag_mean, Y_mag_std, Y_mag_stderr = stats(Y_mag_array)

from matplotlib import pyplot as plt
"""plt.plot(FreqData,(Y_mag_std/Y_mag_mean)*1e6)
plt.axvline(f0_air_mean,label='Series-Resonant Frequency',color='firebrick',alpha=0.5,ls='--')
plt.axvline(fp_air_mean,label='Parallel-Resonant Frequency',color='firebrick',alpha=0.5,ls='--')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Relative Standard Deviation of Admittance')
plt.legend();plt.grid(alpha=0.3)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()"""
#plt.plot(FreqData,(Re_air_std/Re_air_mean)*100,label="Real Uncertainty")
#plt.plot(FreqData,(Im_air_std/Im_air_mean)*100,label="Imaginary Uncertainty")
plt.plot(np.array(FreqData)*1e-3,(Y_mag_std/Y_mag_mean)*100,label="Magnitude Uncertainty")
plt.axvline(f0_air_mean*1e-3,label='Series-Resonant Frequency',color='firebrick',alpha=0.5,ls='--')
plt.axvline(fp_air_mean*1e-3,label='Parallel-Resonant Frequency',color='firebrick',alpha=0.5,ls='--')
plt.xlabel('Frequency (kHz)',fontsize=15)
plt.ylabel('Relative Standard Deviation of Admittance',fontsize=15)
plt.legend();plt.grid(alpha=0.3)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

realY_fluid_array = []
imagY_fluid_array = []
C0_fluid_array = []

f0_fluid_array = []
fp_fluid_array = []
Q_fluid_array = []
for i in range(len(fluid_files)):
    TF1 = ForkData(fluid_files,[0.2 for j in range(len(fluid_files))],100.6e3,R_file,0.1)
    FreqData, X_data, Y_data = TF1.Read_XY_Voltage(fluid_files[i])
    Z_fork = TF1.CalibrationResistorValue * (TF1.DrivingVoltages[i]/TF1.CalibrationResistorVoltage) * \
            TF1.CalibrationVoltage / (X_data + 1j*Y_data)
    Y_fork = 1 / Z_fork
    #realY_fluid_array.append(np.array(np.real(Y_fork)))
    #imagY_fluid_array.append(np.array(np.imag(Y_fork)))
    C0 = TF1.parallel_capacitance_estimation(Y_fork,FreqData)
    C0_fluid_array.append(C0)
    
    Y_m = Y_fork - 1j*2*np.pi*C0*np.array(FreqData)
    realY_fluid_array.append(np.array(np.real(Y_m)))
    imagY_fluid_array.append(np.array(np.imag(Y_m)))

    res = TF1.LorentzianFit(FreqData,abs(Y_m))
    f0_fluid_array.append(res[0])
    Q_fluid_array.append(res[2])
    res = TF1.LorentzianFit(FreqData,abs(Z_fork))
    fp_fluid_array.append(res[0])
    """plt.plot(np.abs(Y_fork))

plt.xlabel('Frequency (Hz)')
plt.ylabel('Admittance')
plt.yscale('log')
plt.legend();plt.grid(alpha=0.3)
plt.show()
plt.show()"""
C0_fluid_mean, C0_fluid_std, C0_fluid_stderr = stats(C0_fluid_array)
Re_fluid_mean, Re_fluid_std, Re_fluid_stderr = stats(realY_fluid_array)
#print(fluid_mean)
#print("",fluid_stderr[argmax(Y_fork)])
imagY_fluid_array = imagY_fluid_array - 2*np.pi*FreqData*C0_fluid_mean
Im_fluid_mean, Im_fluid_std, Im_fluid_stderr = stats(imagY_fluid_array)
print("C0 in fluid:",C0_fluid_mean,"+-",C0_fluid_stderr)
f0_fluid_mean,f0_fluid_std,f0_fluid_stderr = stats(f0_fluid_array)
print("Resonant Frequency in fluid:",f0_fluid_mean,"+-",f0_fluid_stderr)
fp_fluid_mean,fp_fluid_std,fp_fluid_stderr = stats(fp_fluid_array)
print("Parallel-Resonant Frequency in fluid:",fp_fluid_mean,"+-",fp_fluid_stderr)
Q_fluid_mean,Q_fluid_std,Q_fluid_stderr = stats(Q_fluid_array)
print("Q-factor in fluid",Q_fluid_mean,"+-",Q_fluid_stderr)

Y_mag_array = [np.abs(1/(r + 1j*i)) for r, i in zip(realY_fluid_array, imagY_fluid_array)]
Y_mag_mean, Y_mag_std, Y_mag_stderr = stats(Y_mag_array)

#plt.plot(FreqData,(Re_air_std/Re_air_mean)*100,label="Real Uncertainty")
#plt.plot(FreqData,(Im_air_std/Im_air_mean)*100,label="Imaginary Uncertainty")
plt.plot(np.array(FreqData)*1e-3,(Y_mag_std / Y_mag_mean)*100,label="Magnitude Uncertainty")
plt.axvline(f0_fluid_mean*1e-3,label='Series-Resonant Frequency',color='firebrick',alpha=0.5,ls='--')
plt.axvline(fp_fluid_mean*1e-3,label='Parallel-Resonant Frequency',color='firebrick',alpha=0.5,ls='--')
plt.xlabel('Frequency (kHz)',fontsize=15)
plt.ylabel('Relative Standard Deviation of Motional Impedance',fontsize=15)
plt.legend();plt.grid(alpha=0.3)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

plt.plot(np.array(FreqData)*1e-3,(Y_mag_std/Y_mag_mean)*100)
plt.axvline(f0_fluid_mean*1e-3,label='Series-Resonant Frequency',color='firebrick',alpha=0.5,ls='--')
plt.axvline(fp_fluid_mean*1e-3,label='Parallel-Resonant Frequency',color='firebrick',alpha=0.5,ls='--')
plt.xlabel('Frequency (kHz)',fontsize=15)
plt.ylabel('Relative Standard Deviation of Admittance',fontsize=15)
plt.legend();plt.grid(alpha=0.3)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()