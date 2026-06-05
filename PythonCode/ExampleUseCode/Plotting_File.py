from matplotlib import pyplot as plt
from ResonatorClasses.Read_Fork_Data import ForkData
import numpy as np

plt.rcParams.update({
        "font.size": 10,
        "font.family": "sans-serif",
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
})

ResistorFilename = "Data/R_CalFile_20-28kHz_2000points_100mVpk.txt"
File1 = "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_22.30-25.51kHz_W1G0.5_200mVpk_1100points_35C.txt"
File2 = "Data/Water and glycerol data/2 electrodes/TF1_A&B_22.51-25.57kHz_W1G0.5_200mVpk_35C.txt"
airFile = "Data/Initial data collected (pre week 5)/TF1_airsweep_25-28kHz_20251024.txt"
TF_Cal = ForkData(airFile,0.2,100.6e3,ResistorFilename,0.1)
TF_Cal2 = ForkData(airFile,0.2,100.6e3,ResistorFilename,0.1)
R_FreqData, R_X_Data, R_Y_Data = TF_Cal.Read_XY_Voltage(ResistorFilename)
FreqData, X_Data, Y_Data = TF_Cal.Read_XY_Voltage(airFile)
FreqData2, X_Data2, Y_Data2 = TF_Cal.Read_XY_Voltage(File2)
plt.plot(np.array(R_FreqData)*1e-3,R_X_Data,color='firebrick')
plt.plot(np.array(R_FreqData)*1e-3,R_Y_Data,color='mediumblue')

plt.grid(alpha=0.4)
        
# Ticks (major+minor, inward, on all sides)
plt.minorticks_on()
plt.tick_params(
    axis="both", 
    which="major",
    direction="in",
    length=6,
    width=1,
    bottom=True,
    left=True
    )
plt.tick_params(
    axis="both",
    which="minor",
    direction="in",
    length=4,
    width=1,
    bottom=True,
    left=True,
)

plt.xlabel('Frequency (kHz)',fontsize=15)
plt.ylabel('Voltage output (V)',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

plt.plot(np.array(FreqData)*1e-3,X_Data,color='firebrick')
plt.plot(np.array(FreqData)*1e-3,Y_Data,color='mediumblue')

plt.grid(alpha=0.4)
        
# Ticks (major+minor, inward, on all sides)
plt.minorticks_on()
plt.tick_params(
    axis="both", 
    which="major",
    direction="in",
    length=6,
    width=1,
    bottom=True,
    left=True
    )
plt.tick_params(
    axis="both",
    which="minor",
    direction="in",
    length=4,
    width=1,
    bottom=True,
    left=True,
)
plt.xlabel('Frequency (kHz)',fontsize=15)
plt.ylabel('Voltage output (V)',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

Z_air = 100.6e3 * (0.2/0.1) * \
            (np.average(R_X_Data) + 1j*np.average(R_Y_Data)) / (X_Data + 1j*Y_Data)
Y_air = 1/Z_air
C0 = TF_Cal.parallel_capacitance_estimation(Y_air,FreqData)

Y_motional = np.array(Y_air) - 1j*np.array(FreqData)*2*np.pi*C0
Z_air2 = 100.6e3 * (0.2/0.1) * \
            (np.average(R_X_Data) + 1j*np.average(R_Y_Data)) / (X_Data2 + 1j*Y_Data2)
Y_air2 = 1/Z_air2

C02 = TF_Cal.parallel_capacitance_estimation(Y_air2,FreqData2)
Y_motional2 = np.array(Y_air2) - 1j*np.array(FreqData2)*2*np.pi*C02
print(C02)
plt.plot(np.array(FreqData)*1e-3,np.real(Y_air)*1e6,color='firebrick')
plt.plot(np.array(FreqData)*1e-3,np.imag(Y_air)*1e6,color='mediumblue')
#plt.plot(FreqData,2*np.pi*np.array(FreqData) * C0 * 1e6, color='mediumblue')
#plt.plot(FreqData,np.imag(Y_air)*1e6,color='mediumblue')
#plt.yscale('log')
plt.grid(alpha=0.4)
        
# Ticks (major+minor, inward, on all sides)
plt.minorticks_on()
plt.tick_params(
    axis="both", 
    which="major",
    direction="in",
    length=6,
    width=1,
    bottom=True,
    left=True
    )
plt.tick_params(
    axis="both",
    which="minor",
    direction="in",
    length=4,
    width=1,
    bottom=True,
    left=True,
)
plt.xlabel('Frequency (kHz)',fontsize=15)
plt.ylabel(r'Fork Admittance ($\mu$S)',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.show()
