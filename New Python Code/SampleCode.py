#Sample code to run density and viscosity measurements with TF
from Read_Fork_Data import ForkData
from Fork_Calibration import ForkCalibration
from Fork_Sensor import ForkSensor
from AqGlyc_Calculator import AqGly_rho_eta
import numpy as np

#Boolean variables for each fitting method, set to false if you 
#do not want to use the corresponding method
Least_Squares = True
Reittinger = True
Hydrodynamic = True

#All files to be inputted can be inputted as strings of the path to the file
#The ResistorFilename includes an example of this
ResistorFilename = "Data/R_CalFile_20-28kHz_2000points_100mVpk.txt"
CalibrationResistorDrivingVoltage = 0.1 #mVpk
CalibrationResistorValue = 100.6e3 #Ohms
AirFilename = ""
AirDrivingVoltage = 0.2 #for example

#Files that you wish to test
file_list = [

]
#Driving voltages for the files above, in mVpk
DrivingVoltages = [

]

#Filenames of the Fork in 3-5 solutions of AqGlyc/other solution
CalibrationFilenames = [

]
#Driving voltages for the calibration files in mVpk
CalibrationDrivingVoltages = [

]
#Fluid parameter values for calibration
calibration_rho_list = [

]
calibration_eta_list = [

]
#If uncertainties are available, use them, if not, set each uncertainty to 0.
calibration_rho_uncertainty_list = [

]
calibration_eta_uncertainty_list = [

]

#Initialise Calibration Class
TF1_Calibrate = ForkCalibration(file_list,DrivingVoltages,
                      CalibrationResistorValue,
                      ResistorFilename,
                      CalibrationResistorDrivingVoltage,
                      CalibrationFilenames,
                      CalibrationDrivingVoltages,
                      None,None,None,None,None,
                      AirFilename,AirDrivingVoltage,
                      None,None,None)

if Least_Squares:
    Rm, Lm, Cm, C0 = TF1_Calibrate.calibrate_air()

    print(f"Rm = {Rm}")
    print(f"Lm = {Lm}")
    print(f"Cm = {Cm}")
    print(f"C0 = {C0}")

    A, B, sigma_A, sigma_B = TF1_Calibrate.calibrate_A_R_B_R_LeastSquares(
    calibration_rho_list,calibration_eta_list,
    calibration_rho_uncertainty_list,calibration_eta_uncertainty_list)

    print(f"A = {A} +- {sigma_A}")
    print(f"B = {B} +- {sigma_B}")
if Reittinger:
    A_Re, A_std, B_Re, B_std, omega_s_vac, omega_45_vac = \
    TF1_Calibrate.calibrate_Reittinger(calibration_rho_list,calibration_eta_list)
if Hydrodynamic:
    #Hydrodynamic parameters are saved to class variables
    f0, Q0 = TF1_Calibrate.calibrate_fluid_dynamics(calibration_rho_list,calibration_eta_list)

TF1 = ForkSensor(calibration=TF1_Calibrate,
    Filenames=file_list,
    DrivingVoltages=DrivingVoltages,
    CalibrationResistorValue=CalibrationResistorValue,
    CalibrationResistorFile=ResistorFilename,
    CalibrationResistorVoltage=CalibrationResistorDrivingVoltage)
#Perform measurements based off of calibration
if Least_Squares:
    rho_fit_list, eta_fit_list, C0_list = TF1.sensor_AB()
    for index, item in enumerate(rho_fit_list):
        print(f"{index}: {item}")
    for index, item in enumerate(eta_fit_list):
        print(f"{index}: {item}")
if Reittinger:
    rho_fit_list_Reittinger, eta_fit_list_Reittinger = \
    TF1.sensor_Reittinger(A_Re,B_Re,omega_s_vac,omega_45_vac)
    for index, item in enumerate(rho_fit_list_Reittinger):
        print(f"{index}: {item}")
    for index, item in enumerate(eta_fit_list_Reittinger):
        print(f"{index}: {item}")
if Hydrodynamic:
    hydro_rho_list, hydro_eta_list = TF1.sensor_Hydro(f0)
    for index, item in enumerate(hydro_rho_list):
        print(f"{index}: {item}")
    for index, item in enumerate(hydro_eta_list):
        print(f"{index}: {item}")

#For an example of plotting see Running_File.py or Plotting_File.py