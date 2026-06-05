from Read_Fork_Data import ForkData
from Fork_Calibration import ForkCalibration
import numpy as np

class ForkSensor(ForkData):
    def __init__(self, calibration: ForkCalibration,Filenames,DrivingVoltages,
                 CalibrationResistorValue,CalibrationResistorFile,
                 CalibrationResistorVoltage,):
        self.calibration = calibration
        super().__init__(Filenames,DrivingVoltages,CalibrationResistorValue,CalibrationResistorFile,
                 CalibrationResistorVoltage)
        
        # pull calibrated parameters
        self.Rm = calibration.Rm
        self.Lm = calibration.Lm
        self.Cm = calibration.Cm
        self.A = calibration.A
        self.B = calibration.B
    
    def sensor_AB(self):
        #Define Impedance model
        def impedance(f,A_R,B_R):
            w = 2*np.pi*f
            Z_air = self.calibration.Rm + \
            1j * w * self.calibration.Lm - 1j/(w*self.calibration.Cm)
            Z_fluid = 1j * A_R * w + \
            (1+1j)*B_R*np.sqrt(w)
            return (Z_air + Z_fluid)
        
        from scipy.optimize import least_squares
        #Define residual function for least_squares
        def residuals(params, f, Z_meas):
            A_R, B_R = params
            Z_model = impedance(f, A_R, B_R)

            #Estimate uncertainty
            Z_abs = np.abs(Z_meas)
            Z_norm = Z_abs / np.max(Z_abs)   
            Z_inv = 1 - Z_norm 

            rel_sigma = 0.025 + 0.04 * Z_inv  
            sigma = rel_sigma * np.maximum(np.abs(Z_meas), 1e-12)

            return np.concatenate([
                (Z_model.real - Z_meas.real)/sigma,
                (Z_model.imag - Z_meas.imag)/sigma
            ])
        #Reduced Chi^2 parameter to estimate strength of fit
        def reduced_chi_squared(f, Z_meas, params):
            A_R, B_R = params
            Z_model = impedance(f, A_R, B_R)
            #Estimate uncertainty
            Z_abs = np.abs(Z_meas)
            Z_norm = Z_abs / np.max(Z_abs)
            Z_inv = 1 - Z_norm     

            rel_sigma = 0.025 + 0.04 * Z_inv  
            sigma = rel_sigma * np.maximum(np.abs(Z_meas), 1e-12)

            # χ²: include real + imaginary parts
            chi2 = np.sum(
                ((Z_meas.real - Z_model.real)/sigma)**2 +
                ((Z_meas.imag - Z_model.imag)/sigma)**2
            )

            N = len(f)
            p = 2  # number of parameters 
            dof = 2*N - p  # factor 2 because complex → real+imag
            return chi2 / dof
        #Iterate for all input files
        A_R_list = []
        B_R_list = []
        C0_list = []
        for i, file in enumerate(self.Filenames):
            #Calculate admittance
            FreqData,X_data,Y_data = self.Read_XY_Voltage(file)
            Z_fluid = self.CalibrationResistorValue * (self.DrivingVoltages[i]/self.CalibrationResistorVoltage) * \
            self.CalibrationVoltage / (X_data + 1j*Y_data)
            #Admittance and impedance
            Y_fluid = 1 / Z_fluid
            C0 = self.parallel_capacitance_estimation(Y_fluid,FreqData)
            Y_motional = np.array(Y_fluid) - 1j*2*np.pi*C0*np.array(FreqData)
            Z_motional = 1 / Y_motional
            
            #Initial guesses of rho and eta are 1g/cm^3 and 1cP
            A_R0 = self.calibration.A * 1 
            B_R0 = self.calibration.B * np.sqrt(1 * 1)

            x0 = [A_R0, B_R0]
            #Perform least_squares fit
            res = least_squares(
                residuals,
                x0,
                bounds=([0, 0], [np.inf, np.inf]),
                args=(FreqData, Z_motional)
            )

            A_R, B_R = res.x
            A_R_list.append(A_R)
            B_R_list.append(B_R)
            C0_list.append(C0)

            print(f"Reduced Chi^2 for fit {i} = {reduced_chi_squared(FreqData,Z_motional,(A_R, B_R))}")
        
        rho_list = [A_R_list[i] / self.calibration.A for i in range(len(A_R_list))]

        sqrt_rho_eta_list = [B_R_list[i] / self.calibration.B for i in range(len(B_R_list))]

        eta_list = [
            (sqrt_rho_eta_list[i]**2) / rho_list[i]
            for i in range(len(rho_list))
        ]
        print(self.Filenames[-1])
        import matplotlib.pyplot as plt
        plt.plot(np.array(FreqData)*1e-3,np.array(np.abs(Z_fluid))*1e-6,color='firebrick')
        plt.plot(np.array(FreqData)*1e-3,
        np.array(np.abs(1/(1/impedance(FreqData,A_R,B_R) + 1j*2*np.pi*FreqData*C0)))*1e-6,
        color='mediumblue',ls='-.')
        plt.grid(alpha=0.4)
        plt.xlabel('Frequency (kHz)',fontsize=15)#;plt.ylabel(r'Admittance ($\mu$S)',fontsize=15)
        plt.ylabel(r'Impedance (M$\Omega$)',fontsize=15)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.yscale('log')
        plt.show()
        return rho_list, eta_list, C0_list
    
    def sensor_Hydro(self,f0):
        from scipy.optimize import least_squares
        def solve_rho_eta(A, B, C, omega0, omega1, Q1):

            y1 = (1 / Q1) * (omega0 / omega1)**2
            y2 = ((omega0 / omega1)**2) - 1

            def residuals(vars):
                rho, eta = vars

                u = np.sqrt(eta * rho / omega1)
                return [
                    A * u + B * (eta / omega1) - y1,
                    A * u + C * rho - y2
                ]

            result = least_squares(
                residuals,
                x0=[1, 1],
                bounds=([1e-3, 1e-6], [20, 20])
            )

            return result.x
        rho_list = []
        eta_list = []
        omega0 = 2*np.pi*f0
        for i, file in enumerate(self.Filenames):
            #Calculate admittance
            FreqData,X_data,Y_data = self.Read_XY_Voltage(file)
            Z_fluid = self.CalibrationResistorValue * (self.DrivingVoltages[i]/self.CalibrationResistorVoltage) * \
            self.CalibrationVoltage / (X_data + 1j*Y_data)
            #Admittance and impedance magnitudes
            Y_fluid = 1 / Z_fluid
            C0 = self.parallel_capacitance_estimation(Y_fluid,FreqData)
            Y_motional = np.array(Y_fluid) - 1j*2*np.pi*C0*np.array(FreqData)

            res = ForkData.LorentzianFit(FreqData,np.real(Y_motional))
            f1 = res[0]
            Q1 = res[2]
            omega1 = 2*np.pi*f1
            rho, eta = solve_rho_eta(self.calibration.A_Hydro,
                                     self.calibration.B_Hydro,
                                     self.calibration.C_Hydro,
                                     omega0,omega1,Q1)
            rho_list.append(rho)
            eta_list.append(eta)
        return rho_list, eta_list
    
    def sensor_Reittinger(self,A_Rei,B_Rei,omega_s_vac,omega_45_vac):
        rho_list = []
        rho_eta_list = []
        for i, file in enumerate(self.Filenames):
            FreqData,X_data,Y_data = self.Read_XY_Voltage(file)
            Z_fluid = self.CalibrationResistorValue * (self.DrivingVoltages[i]/self.CalibrationResistorVoltage) * \
            self.CalibrationVoltage / (X_data + 1j*Y_data) 
            #Admittance and impedance magnitudes
            Y_fluid = 1 / Z_fluid
            
            C0 = self.parallel_capacitance_estimation(Y_fluid,FreqData)
            Y_motional = np.array(Y_fluid) - 1j*2*np.pi*C0*np.array(FreqData)
            res = self.LorentzianFit(FreqData,np.real(Y_motional))
            omega_s = res[0] * 2 * np.pi
            omega_45 = self.find_first_intersection(FreqData,Y_motional) * 2 * np.pi

            num1 = ((omega_s_vac)/(omega_45))**2 - ((omega_s_vac)/(omega_s))**2 \
                - 2*((omega_s_vac-omega_45_vac)/(omega_45))
            den1 = B_Rei*(2/(np.sqrt(omega_45)) - 1/np.sqrt(omega_s))
            rho_eta_list.append((num1/den1)**2)

            num2 = ((omega_s_vac)/(omega_s))**2 - B_Rei*(np.sqrt(((num1/den1)**2)/omega_s)) - 1
            den2 = A_Rei

            rho_list.append(num2/den2)
        eta_list = [rho_eta_list[i] / rho_list[i] for i in range(len(rho_list))]
        return rho_list, eta_list

    def plot_spectra(self,admittance=True,logy=True,motional=False):
        from matplotlib import pyplot as plt
        plotcol = ["firebrick", "seagreen", 
                "mediumblue", "darkgoldenrod", 
                "saddlebrown", "slategray", "teal", 
                "peru", "steelblue","olive",
                "yellowgreen", "lightcoral", 
                "limegreen", "thistle", "indigo"]
        plt.rcParams.update({
        "font.size": 10,
        "font.family": "sans-serif",
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    })
        fig, axs = plt.subplots(2,1,sharex=True)
        for i, file in enumerate(self.Filenames[::3]):
            print(file)
            FreqData,X_data,Y_data = self.Read_XY_Voltage(file)
            Z_fluid = self.CalibrationResistorValue * (self.DrivingVoltages[i]/self.CalibrationResistorVoltage) * \
            self.CalibrationVoltage / (X_data + 1j*Y_data)
            #Admittance and impedance magnitudes
            Y_fluid = 1 / Z_fluid
            C0 = self.parallel_capacitance_estimation(Y_fluid,FreqData)
            Y_motional = np.array(Y_fluid) - 1j*2*np.pi*C0*np.array(FreqData)
            Z_motional = 1 / Y_motional

            if admittance and motional: 
                data = np.abs(Y_motional) * 1e6
                phase = np.angle(Y_motional,deg=True)
            elif admittance and not motional: 
                data = np.abs(Y_fluid) * 1e6
                phase = np.angle(Y_fluid,deg=True)
            elif not admittance and motional: 
                data = np.abs(Z_motional) * 1e-6
                phase = np.angle(Z_motional,deg=True)
            else: 
                data = np.abs(Z_fluid) * 1e-6
                phase = np.angle(Z_fluid,deg=True)
            
            axs[0].plot(np.array(FreqData)*1e-3, data, 
            color = plotcol[i] if i < len(plotcol) else None)
            axs[1].plot(np.array(FreqData)*1e-3,phase,
            color = plotcol[i] if i < len(plotcol) else None)

        
        axs[1].set_xlabel('Frequency (kHz)',fontsize=15)
        if admittance:
            axs[0].set_ylabel(r'Admittance ($\mu S$)',fontsize=15)
        else:
            axs[0].set_ylabel(r'Impedance ($M \Omega$)',fontsize=15)
        axs[1].set_ylabel('Phase (°)',fontsize=15)
        for ax in axs:
            ax.grid(True, which='both', alpha=0.3)
        
        # Ticks (major+minor, inward, on all sides)
        axs[0].minorticks_on()
        axs[1].minorticks_on()
        """fig.tick_params(
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
        fig.tick_params(
            axis="both",
            which="minor",
            direction="in",
            length=4,
            width=1,
            #top=True,
            bottom=True,
            left=True,
            #right=True,
        )"""
        if logy: axs[0].set_yscale('log')
        plt.xticks(fontsize=15)
        axs[0].tick_params(axis='y', labelsize=15)
        axs[1].tick_params(axis='y', labelsize=15)
        plt.show()