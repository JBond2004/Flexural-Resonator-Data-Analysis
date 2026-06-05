from Fork_Calibration import ForkCalibration
from Fork_Sensor import ForkSensor
from AqGlyc_Calculator import AqGly_rho_eta
import numpy as np
print('yo')
ResistorFilename = "Data/R_CalFile_20-28kHz_2000points_100mVpk.txt"
CalibrationResistorDrivingVoltage = 0.1
AirFilename = "Data/Initial data collected (pre week 5)/TF1_airsweep_25-28kHz_20251024.txt"
AirDrivingVoltage = 0.2

file_list = [
    #"Data/Water and glycerol data/3 electrodes/TF1_3electrode_23.24-27.08kHz_water_200mVpk_35C_take2.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_22.62-26.04kHz_W1G0.1_200mVpk_1100points_35C.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_22.39-25.12kHz_W1G0.2_200mVpk_1100points_35C.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_22.30-25.51kHz_W1G0.3_200mVpk_1100points_35C.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_22.41-25.36kHz_W1G0.4_200mVpk_1100points_35C.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_22.30-25.51kHz_W1G0.5_200mVpk_1100points_35C.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_22.98-25.49kHz_W1G0.6_200mVpk_1000points_35C.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_22.5-25.5kHz_W1G0.7_200mVpk_1000points_35C.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_21.72-26.38kHz_W1G0.8_200mVpk_1000points_35C.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_22.43-26.03kHz_W1G0.9_200mVpk_1000points_35C.txt",
    "Data/Water and glycerol data/3 electrodes/TF1_3electrodes_20.91-26.77kHz_W1G1_200mVpk_1000points_35C.txt"
]

CalibrationFilenames = file_list[::2]

DrivingVoltages = [
    0.2 for i in range(len(file_list))
]
CalibrationDrivingVoltages = DrivingVoltages[::2]

rho_list = [AqGly_rho_eta(30,0.1+0.1*i,1)[0]/1000 for i in range(len(file_list))]
eta_list = [AqGly_rho_eta(30,0.1+0.1*i,1)[1] for i in range(len(file_list))]
print(str(rho_list[0]),"to",str(rho_list[-1]))
print(str(eta_list[0]),"to",str(eta_list[-1]))

rho_uncertainty_list = [
    abs(AqGly_rho_eta(32, 0.1+0.1*i, 1)[0] - AqGly_rho_eta(30, 0.1+0.1*i, 1)[0]) / 2000
    for i in range(len(file_list))
]
eta_uncertainty_list = [
    abs(AqGly_rho_eta(32, 0.1+0.1*i, 1)[1] - AqGly_rho_eta(30, 0.1+0.1*i, 1)[1]) / 2
    for i in range(len(file_list))
]
calibration_rho_list = rho_list[::2]
calibration_eta_list = eta_list[::2]

TF1_Calibrate = ForkCalibration(file_list,DrivingVoltages,100.6e3,ResistorFilename,
                      0.1,CalibrationFilenames,CalibrationDrivingVoltages,
                      None,None,None,None,None,
                      AirFilename,0.2,None,None,None)
A_Re, A_std, B_Re, B_std, omega_s_vac, omega_45_vac = TF1_Calibrate.calibrate_Reittinger(calibration_rho_list,calibration_eta_list)
print(A_Re, "+-", A_std)
print(B_Re, "+-", B_std)
TF1 = ForkSensor(calibration=TF1_Calibrate,
    Filenames=file_list,
    DrivingVoltages=DrivingVoltages,
    CalibrationResistorValue=100.6e3,
    CalibrationResistorFile=ResistorFilename,
    CalibrationResistorVoltage=0.1)

TF1.plot_spectra(admittance=False)

rho_fit_list_Reittinger, eta_fit_list_Reittinger = TF1.sensor_Reittinger(A_Re,B_Re,omega_s_vac,omega_45_vac)

f0, Q0 = TF1_Calibrate.calibrate_fluid_dynamics(calibration_rho_list,calibration_eta_list)

hydro_rho_list, hydro_eta_list = TF1.sensor_Hydro(f0)

Rm, Lm, Cm, C0 = TF1_Calibrate.calibrate_air()

print(f"Rm = {Rm}")
print(f"Lm = {Lm}")
print(f"Cm = {Cm}")
print(f"C0 = {C0}")

A, B, sigma_A, sigma_B = TF1_Calibrate.calibrate_A_R_B_R_LeastSquares(calibration_rho_list,calibration_eta_list,rho_uncertainty_list[::2],eta_uncertainty_list[::2])

print(f"A = {A} +- {sigma_A}")
print(f"B = {B} +- {sigma_B}")

rho_fit_list, eta_fit_list, C0_list = TF1.sensor_AB()

print(rho_fit_list)
print(eta_fit_list)
print(C0_list)
from matplotlib import pyplot as plt

plt.errorbar(rho_list,100*(rho_fit_list-np.array(rho_list))/np.array(rho_list),
100 * np.array(rho_uncertainty_list) / np.array(rho_list),
fmt='o',capsize=3,label='Least Squares',color='firebrick')
plt.errorbar(rho_list,100*(hydro_rho_list-np.array(rho_list))/np.array(rho_list),
100 * np.array(rho_uncertainty_list) / np.array(rho_list),
fmt='o',capsize=3,label='Hydrodynamic Model',color='mediumblue')
plt.errorbar(rho_list,100*(rho_fit_list_Reittinger-np.array(rho_list))/np.array(rho_list),
100 * np.array(rho_uncertainty_list) / np.array(rho_list),
fmt='o',capsize=3,label="Reittinger's Model",color='seagreen')
plt.xlabel(r'Density ($gcm^{-3}$)',fontsize=15)
plt.ylabel('Percentage Scatter (%)',fontsize=15)
plt.grid(alpha=0.4)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=15)
plt.show()

plt.errorbar(eta_list,100*(eta_fit_list-np.array(eta_list))/np.array(eta_list),
100 * np.array(eta_uncertainty_list) / np.array(eta_list),
fmt='o',capsize=3,label='Least Squares',color='firebrick')
plt.errorbar(eta_list,100*(hydro_eta_list-np.array(eta_list))/np.array(eta_list),
100 * np.array(eta_uncertainty_list) / np.array(eta_list),
fmt='o',capsize=3,label='Hydrodynamic Model',color='mediumblue')
plt.errorbar(eta_list,100*(eta_fit_list_Reittinger-np.array(eta_list))/np.array(eta_list),
100 * np.array(eta_uncertainty_list) / np.array(eta_list),
fmt='o',capsize=3,label="Reittinger's Model",color='seagreen')
plt.xlabel('Viscosity (cP)',fontsize=15)
plt.ylabel('Percentage Scatter (%)',fontsize=15)
plt.grid(alpha=0.4)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=15)
plt.show()
plt.rcParams.update({
        "font.size": 10,
        "font.family": "sans-serif",
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    })
plt.minorticks_on()
plt.tick_params(
    axis="both", 
    which="major",
    direction="in",
    length=6,
    width=1,
# top=True,    # enable top ticks
    bottom=True,
    left=True,
    #right=True,  # enable right ticks
)
plt.tick_params(
    axis="both",
    which="minor",
    direction="in",
    length=4,
    width=1,
    #top=True,
    bottom=True,
    left=True,
    #right=True,
)

file = ["Data/Water and glycerol data/3 electrodes/TF1_3electrodes_21.62-23.48kHz_P9900_21.4C_100mVpk_900points.txt"]

TF1_P99 = ForkSensor(
    calibration=TF1_Calibrate,
    Filenames=file,
    DrivingVoltages=[0.1],
    CalibrationResistorValue=100.6e3,
    CalibrationResistorFile=ResistorFilename,
    CalibrationResistorVoltage=0.1
)

rho_ls, eta_ls, C0_list = TF1_P99.sensor_AB()
rho_Re, eta_Re = TF1_P99.sensor_Reittinger(A_Re,B_Re,omega_s_vac,omega_45_vac)
rho_Z, eta_Z = TF1_P99.sensor_Hydro(f0)

print(rho_ls,rho_Re,rho_Z)
print(eta_ls,eta_Re,eta_Z)
plt.errorbar(rho_list,rho_fit_list,rho_uncertainty_list,
fmt='o',capsize=3,label='Least Squares',color='firebrick')
plt.errorbar(rho_list,rho_fit_list_Reittinger,rho_uncertainty_list,
fmt='o',capsize=3,label="Reittinger's model",color='mediumblue')
plt.errorbar(rho_list,hydro_rho_list,rho_uncertainty_list,
fmt='o',capsize=3,label='Hydrodynamic Model',color='seagreen')
plt.plot([min(rho_list), max(rho_list)], [min(rho_list), max(rho_list)], 'k--',color='black')
plt.xlabel(r'Density ($gcm^{-3}$)',fontsize=14)
plt.ylabel(r'Fitted Density ($gcm^{-3}$)',fontsize=14)
plt.grid(alpha=0.4)
plt.legend()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()

plt.errorbar(eta_list,eta_fit_list,eta_uncertainty_list,
fmt='o',capsize=3,label='Least Squares',color='firebrick')
plt.errorbar(eta_list,eta_fit_list_Reittinger,eta_uncertainty_list,
fmt='o',capsize=3,label="Reittinger's model",color='mediumblue')
plt.errorbar(eta_list,hydro_eta_list,eta_uncertainty_list,
fmt='o',capsize=3,label='Hydrodynamic Model',color='seagreen')
plt.plot([min(eta_list), max(eta_list)], [min(eta_list), max(eta_list)], 'k--',color='black')
plt.xlabel('Viscosity (cP)',fontsize=14)
plt.ylabel('Fitted Viscosity (cP)',fontsize=14)
plt.grid(alpha=0.4)
plt.legend()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()

file_list = [
    "Data/ZnSO4/0.5M/TF1_3electrodes_23.15-25.43kHz_ZnSO40.5M_35C_data1.txt",
    "Data/ZnSO4/0.5M/TF1_3electrodes_23.15-25.43kHz_ZnSO40.5M_35C_data2.txt",
    "Data/ZnSO4/0.5M/TF1_3electrodes_23.25-25.37kHz_100mVpk_ZnSO40.5M_35C_data3.txt",

    "Data/ZnSO4/1.0M/TF1_23.14-25.34kHz_100mVpk_35C_ZnSO4_1.0M_data1.txt",
    "Data/ZnSO4/1.0M/TF1_23.14-25.34kHz_100mVpk_35C_ZnSO4_1.0M_data2.txt",
    "Data/ZnSO4/1.0M/TF1_23.14-25.34kHz_100mVpk_35C_ZnSO4_1.0M_data3.txt",

    "Data/ZnSO4/1.5M/35C/TF1_22.98-25.05kHz_100mVpk_35C_ZnSO4_1.5M_data1.txt",
    "Data/ZnSO4/1.5M/35C/TF1_22.98-25.05kHz_100mVpk_35C_ZnSO4_1.5M_data2.txt",
    "Data/ZnSO4/1.5M/35C/TF1_22.98-25.05kHz_100mVpk_35C_ZnSO4_1.5M_data3.txt",

    "Data/ZnSO4/2.0M/35C/TF1_22.85-25.00kHz_ZnSO4_2.0M_35C_100mVpk_data2.txt",
    "Data/ZnSO4/2.0M/35C/TF1_22.98-25.05kHz_ZnSO4_2.0M_35C_100mVpk_data1.txt",
    "Data/ZnSO4/2.0M/35C/TF1_22.98-25.05kHz_ZnSO4_2.0M_35C_100mVpk_data3.txt",

    "Data/ZnSO4/2.5M/35C/TF1_22.48-24.88kHz_ZnSO4_2.5M_35C_100mVpk_data1.txt",
    "Data/ZnSO4/2.5M/35C/TF1_22.48-24.88kHz_ZnSO4_2.5M_35C_100mVpk_data2.txt",
    "Data/ZnSO4/2.5M/35C/TF1_22.48-24.88kHz_ZnSO4_2.5M_35C_100mVpk_data3.txt"
]

DrivingVoltages = [0.1 for i in range(len(file_list))]

TF1_Zn = ForkSensor(
    calibration=TF1_Calibrate,
    Filenames=file_list,
    DrivingVoltages=DrivingVoltages,
    CalibrationResistorValue=100.6e3,
    CalibrationResistorFile=ResistorFilename,
    CalibrationResistorVoltage=0.1
)

rho_list, eta_list, C0_list = TF1_Zn.sensor_AB()
rho_raw = rho_list
eta_raw = eta_list

rho_list = [np.mean(rho_raw[0:3]),
np.mean(rho_raw[3:6]),np.mean(rho_raw[6:7]),
np.mean(rho_raw[9:12]),np.mean(rho_raw[12:15])]

eta_list = [np.mean(eta_raw[0:3]),
np.mean(eta_raw[3:6]),np.mean(eta_raw[6:9]),
np.mean(eta_raw[9:12]),np.mean(eta_raw[12:15])]

rho_unc_list = [np.std(rho_raw[0:3], ddof=1)/np.sqrt(3),
           np.std(rho_raw[3:6], ddof=1)/np.sqrt(3),
           np.std(rho_raw[6:9], ddof=1)/np.sqrt(3),
           np.std(rho_raw[9:12], ddof=1)/np.sqrt(3),
           np.std(rho_raw[12:15], ddof=1)/np.sqrt(3)]

eta_unc_list = [np.std(eta_raw[0:3])/np.sqrt(3),
np.std(eta_raw[3:6])/np.sqrt(3),np.std(eta_raw[6:9])/np.sqrt(3),
np.std(eta_raw[9:12])/np.sqrt(3),np.std(eta_raw[12:15])/np.sqrt(3)]

rho_list_R, eta_list_R = TF1_Zn.sensor_Reittinger(A_Re,B_Re,omega_s_vac,omega_45_vac)
rho_raw_R = rho_list_R
eta_raw_R = eta_list_R

rho_list_R = [np.mean(rho_raw_R[0:3]),
np.mean(rho_raw_R[3:6]),np.mean(rho_raw_R[6:7]),
np.mean(rho_raw_R[9:12]),np.mean(rho_raw_R[12:15])]

eta_list_R = [np.mean(eta_raw_R[0:3]),
np.mean(eta_raw_R[3:6]),np.mean(eta_raw_R[6:9]),
np.mean(eta_raw_R[9:12]),np.mean(eta_raw_R[12:15])]

rho_unc_list_R = [np.std(rho_raw_R[0:3], ddof=1)/np.sqrt(3),
           np.std(rho_raw_R[3:6], ddof=1)/np.sqrt(3),
           np.std(rho_raw_R[6:9], ddof=1)/np.sqrt(3),
           np.std(rho_raw_R[9:12], ddof=1)/np.sqrt(3),
           np.std(rho_raw_R[12:15], ddof=1)/np.sqrt(3)]

eta_unc_list_R = [np.std(eta_raw_R[0:3])/np.sqrt(3),
np.std(eta_raw_R[3:6])/np.sqrt(3),np.std(eta_raw_R[6:9])/np.sqrt(3),
np.std(eta_raw_R[9:12])/np.sqrt(3),np.std(eta_raw_R[12:15])/np.sqrt(3)]

rho_list_H, eta_list_H = TF1_Zn.sensor_Hydro(f0)
rho_raw_H = rho_list_H
eta_raw_H = eta_list_H

rho_list_H = [np.mean(rho_raw_H[0:3]),
np.mean(rho_raw_H[3:6]),np.mean(rho_raw_H[6:7]),
np.mean(rho_raw_H[9:12]),np.mean(rho_raw_H[12:15])]

eta_list_H = [np.mean(eta_raw_H[0:3]),
np.mean(eta_raw_H[3:6]),np.mean(eta_raw_H[6:9]),
np.mean(eta_raw_H[9:12]),np.mean(eta_raw_H[12:15])]

rho_unc_list_H = [np.std(rho_raw_H[0:3], ddof=1)/np.sqrt(3),
           np.std(rho_raw_H[3:6], ddof=1)/np.sqrt(3),
           np.std(rho_raw_H[6:9], ddof=1)/np.sqrt(3),
           np.std(rho_raw_H[9:12], ddof=1)/np.sqrt(3),
           np.std(rho_raw_H[12:15], ddof=1)/np.sqrt(3)]

eta_unc_list_H = [np.std(eta_raw_H[0:3])/np.sqrt(3),
np.std(eta_raw_H[3:6])/np.sqrt(3),np.std(eta_raw_H[6:9])/np.sqrt(3),
np.std(eta_raw_H[9:12])/np.sqrt(3),np.std(eta_raw_H[12:15])/np.sqrt(3)]

plt.errorbar([0.5,1.0,1.5,2.0,2.5],rho_list,rho_unc_list,
fmt='o',markersize=4,label='Least Squares',color='firebrick')
plt.errorbar([0.5,1.0,1.5,2.0,2.5],rho_list_R,rho_unc_list_R,
fmt='o',markersize=4,label="Reittinger's model",color='mediumblue')
plt.errorbar([0.5,1.0,1.5,2.0,2.5],rho_list_H,rho_unc_list_H,
fmt='o',markersize=4,label='Hydrodynamic Model',color='seagreen')
plt.xlabel(r'$ZnSO_4$ concentration (M)',fontsize=14)
plt.ylabel(r'Fitted Density ($gcm^{-3}$)',fontsize=14)
plt.grid(alpha=0.4)
plt.legend()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()
plt.errorbar([0.5,1.0,1.5,2.0,2.5],eta_list,eta_unc_list,
fmt='o',markersize=4,label='Least Squares',color='firebrick')
plt.errorbar([0.5,1.0,1.5,2.0,2.5],eta_list_R,eta_unc_list_R,
fmt='o',markersize=4,label="Reittinger's model",color='mediumblue')
plt.errorbar([0.5,1.0,1.5,2.0,2.5],eta_list_H,eta_unc_list_H,
fmt='o',markersize=4,label='Hydrodynamic Model',color='seagreen')
plt.xlabel(r'$ZnSO_4$ concentration (M)',fontsize=14)
plt.ylabel('Fitted Viscosity (cP)',fontsize=14)
plt.grid(alpha=0.4)
plt.legend()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()
