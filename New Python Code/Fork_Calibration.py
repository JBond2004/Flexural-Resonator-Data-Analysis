from Read_Fork_Data import ForkData
import numpy as np

class ForkCalibration(ForkData):
    def __init__(self,Filenames,DrivingVoltages,CalibrationResistorValue,CalibrationResistorFile,
                 CalibrationResistorVoltage,CalibrationFilenames,
                 CalibrationDrivingVoltages,A,B,
                 Rm,Lm,Cm,AirCalibrationFile,AirCalibrationDrivingVoltage,
                 A_Hydro,B_Hydro,C_Hydro):
        super().__init__(Filenames,DrivingVoltages,CalibrationResistorValue,CalibrationResistorFile,
                 CalibrationResistorVoltage)
        self.CalibrationFilenames = CalibrationFilenames
        self.CalibrationDrivingVoltages = CalibrationDrivingVoltages
        self.A = A
        self.B = B 
        self.Rm = Rm
        self.Lm = Lm
        self.Cm = Cm
        self.AirCalibrationFile = AirCalibrationFile
        self.AirCalibrationDrivingVoltage = AirCalibrationDrivingVoltage
        self.A_Hydro = A_Hydro
        self.B_Hydro = B_Hydro
        self.C_Hydro = C_Hydro

    def calibrate_air(self):
        #Calculate admittance
        FreqData,X_data,Y_data = self.Read_XY_Voltage(self.AirCalibrationFile)
        Z_air = self.CalibrationResistorValue * (self.AirCalibrationDrivingVoltage/self.CalibrationResistorVoltage) * self.CalibrationVoltage / (X_data + 1j*Y_data) 
        #Admittance and impedance magnitudes
        Y_air = 1 / Z_air
        
        C0 = self.parallel_capacitance_estimation(Y_air,FreqData)
        print("C0 = ",C0)
        Y_motional = np.array(Y_air) - 1j*2*np.pi*C0*np.array(FreqData)

        series_ResFreq, Sigma_ResFreq, Q, Sigma_Q, ResPeak, Sigma_ResPeak = \
        self.LorentzianFit(FreqData,np.real(Y_air))
        parallel_ResFreq, _, _, _, _, _ = \
        self.LorentzianFit(FreqData,np.real(1/Y_air))
        Cm = C0 * ((parallel_ResFreq/series_ResFreq)**2 - 1)
        Lm = 1 / (((2*np.pi*series_ResFreq) ** 2) * Cm)
        Rm = 1/ResPeak

        self.Rm = Rm
        self.Lm = Lm
        self.Cm = Cm

        def reduced_chi_squared(f, Y_meas, params, sigma_rel=0.02):
            Rm, Lm, Cm, C0 = params

            def admittance(f, Rm, Lm, Cm, C0):
                w = 2*np.pi*f
                return 1j*w*C0 + 1/(Rm + 1j*w*Lm - 1j/(w*Cm))

            Y_model = admittance(f, Rm, Lm, Cm, C0)
            Y_abs = np.abs(Y_meas)
            Y_norm = Y_abs / np.max(Y_abs)               

            rel_sigma = 0.03 + 0.07 * Y_norm  
            sigma = rel_sigma * np.maximum(np.abs(Y_meas), 1e-12)

            # Avoid divide-by-zero
            sigma[sigma == 0] = np.min(sigma[sigma > 0])

            # χ²: include real + imaginary parts
            chi2 = np.sum(
                ((Y_meas.real - Y_model.real)/sigma)**2 +
                ((Y_meas.imag - Y_model.imag)/sigma)**2
            )

            N = len(f)
            p = 4  # number of parameters (Rm, Lm, Cm, C0)

            dof = 2*N - p  # factor 2 because complex → real+imag

            return chi2 / dof
        print(f"Reduced Chi^2 = {reduced_chi_squared(FreqData,Y_air,(Rm,Lm,Cm,C0))}")
        def admittance(f, Rm, Lm, Cm, C0):
                w = 2*np.pi*f
                return 1j*w*C0 + 1/(Rm + 1j*w*Lm - 1j/(w*Cm))
        
        from matplotlib import pyplot as plt 
        plt.plot(np.array(FreqData)*1e-3,np.abs(1/Y_motional)*1e-6,label="Data")
        plt.plot(np.array(FreqData)*1e-3,
        np.abs(1/(admittance(FreqData,Rm,Lm,Cm,C0) - 1j*2*np.pi*C0*FreqData))*1e-6,label="Fit")
        plt.legend()
        plt.xlabel('Frequency (kHz)',fontsize=15)
        plt.ylabel(r'Impedance (M$\Omega$)',fontsize=15)
        #plt.ylabel(r'Admittance ($\mu$S)',fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.grid(alpha=0.4)
        plt.yscale('log')
        plt.show()
        return Rm, Lm, Cm, C0
    
    def calibrate_AB_LeastSquares(self,rho,eta):
        #Only use this function if you're calibrating with 1 file
        #Calculate admittance
        FreqData,X_data,Y_data = self.Read_XY_Voltage(self.CalibrationFilenames)
        Z_air = self.CalibrationResistorValue * (self.CalibrationDrivingVoltages/self.CalibrationResistorVoltage) * \
        self.CalibrationVoltage / (X_data + 1j*Y_data) 
        #Admittance and impedance magnitudes
        Y_fluid = 1 / Z_air
        C0 = self.parallel_capacitance_estimation(Y_fluid,FreqData)
        Y_motional = np.array(Y_fluid) - 1j*2*np.pi*C0*np.array(FreqData)

        def motional_fork_admittance(f,A,B):
            w = 2 * np.pi * f
            Z_air = self.Rm + 1j * w * self.Lm - 1j/(w * self.Cm)
            Z_fluid = 1j * A * w * rho + (1+1j)*B*np.sqrt(w * rho * eta)
            return 1/(Z_air + Z_fluid)
        from scipy.optimize import least_squares

        def residuals(params, f, Y_meas):
            A, B = params
            Y_model = motional_fork_admittance(f, A, B)
            return np.concatenate([
                (Y_model.real - Y_meas.real),
                (Y_model.imag - Y_meas.imag)
            ])

        res = least_squares(
            residuals,
            x0=[0.1, 5],
            args=(FreqData, Y_motional)
        )

        A_fit, B_fit = res.x
        self.A = A_fit
        self.B = B_fit

        def reduced_chi_squared(f, Y_meas, params, sigma_rel=0.02):
            """
            f         : frequency array
            Y_meas    : measured complex admittance (same length as f)
            params    : tuple/list (Rm, Lm, Cm, C0)
            sigma_rel : relative uncertainty (default 2%)
            """
            A, B = params

            Y_model = motional_fork_admittance(f, A, B)

            # Estimate uncertainties (2% of magnitude)
            sigma = sigma_rel * np.abs(Y_meas)

            # Avoid divide-by-zero
            sigma[sigma == 0] = np.min(sigma[sigma > 0])

            # χ²: include real + imaginary parts
            chi2 = np.sum(
                ((Y_meas.real - Y_model.real)/sigma)**2 +
                ((Y_meas.imag - Y_model.imag)/sigma)**2
            )

            N = len(f)
            p = 2  # number of parameters (Rm, Lm, Cm, C0)

            dof = 2*N - p  # factor 2 because complex → real+imag

            return chi2 / dof
        print(f"Reduced Chi^2 = {reduced_chi_squared(FreqData,Y_motional,(A_fit, B_fit))}")
        return A_fit, B_fit
    
    def calibrate_A_R_B_R_LeastSquares(self,rho_list,eta_list,rho_err_list,eta_err_list):

        def fork_impedance(f,A,B):
            w = 2 * np.pi * f
            Z_air = self.Rm + 1j * w * self.Lm - 1j/(w * self.Cm)
            Z_fluid = 1j * A * w + (1+1j) * B * np.sqrt(w)
            return (Z_air + Z_fluid)
        from scipy.optimize import least_squares

        def residuals(params, f, Z_meas):
            A, B = params
            Z_model = fork_impedance(f, A, B)

            Z_abs = np.abs(Z_meas)
            Z_norm = Z_abs / np.max(Z_abs)   
            Z_inv = 1 - Z_norm
            rel_sigma = 0.025 + 0.04 * Z_inv 
            sigma = rel_sigma * np.abs(Z_meas)

            return np.concatenate([
                (Z_model.real - Z_meas.real)/sigma,
                (Z_model.imag - Z_meas.imag)/sigma
            ])
        def reduced_chi_squared(f, Z_meas, params, sigma_rel=0.02):
            A, B = params

            Z_model = fork_impedance(f, A, B)

            Z_abs = np.abs(Z_meas)
            Z_norm = Z_abs / np.max(Z_abs)   
            Z_inv = 1 - Z_norm               

            rel_sigma = 0.025 + 0.04 * Z_inv 
            sigma = rel_sigma * np.abs(Z_meas)

            # χ²: include real + imaginary parts
            chi2 = np.sum(
                ((Z_meas.real - Z_model.real)/sigma)**2 +
                ((Z_meas.imag - Z_model.imag)/sigma)**2
            )

            N = len(f)
            p = 2  # number of parameters (Rm, Lm, Cm, C0)
            dof = 2*N - p  # factor 2 because complex → real+imag
            return chi2 / dof
        
        A_fit_list = []
        B_fit_list = []
        A_sigma_list = []
        B_sigma_list = []

        # A_R = A * rho
        # B_R = B * sqrt(rho * eta)
        sqrt_product_list = [np.sqrt(rho_list[i] * eta_list[i]) for i in range(len(rho_list))]

        for i, file in enumerate(self.CalibrationFilenames):
            #Calculate admittance
            FreqData,X_data,Y_data = self.Read_XY_Voltage(file)
            Z_fluid = self.CalibrationResistorValue * \
            (self.CalibrationDrivingVoltages[i]/self.CalibrationResistorVoltage) * \
            self.CalibrationVoltage / (X_data + 1j*Y_data) 
            #Admittance and impedance magnitudes
            Y_fluid = 1 / Z_fluid
            C0 = self.parallel_capacitance_estimation(Y_fluid,FreqData)
            #print("C0 = ", C0)
            Y_motional = np.array(Y_fluid) - 1j*2*np.pi*C0*np.array(FreqData)
            Z_motional = 1/Y_motional
            res = least_squares(
                residuals,
                x0=[0.1, 5*np.sqrt(1)],
                args=(FreqData, Z_motional)
            )

            A_R_fit, B_R_fit = res.x
            A_fit_list.append(A_R_fit)
            B_fit_list.append(B_R_fit)

            # Jacobian
            J = res.jac

            # Degrees of freedom
            N = len(FreqData)
            p = 2
            dof = 2*N - p

            # Residual variance estimate
            residual_var = np.sum(res.fun**2) / dof

            # Covariance matrix
            cov = residual_var * np.linalg.inv(J.T @ J)

            # Parameter uncertainties (standard deviations)
            A_sigma_list.append(np.sqrt(cov[0, 0]))
            B_sigma_list.append(np.sqrt(cov[1, 1]))

            print(f"Reduced Chi^2 = {reduced_chi_squared(FreqData,Z_motional,(A_R_fit, B_R_fit))}")
        A_fit = np.sum(np.array(rho_list) * np.array(A_fit_list)) / np.sum(np.array(rho_list)**2)
        B_fit = np.sum(np.array(sqrt_product_list) * np.array(B_fit_list)) / np.sum(np.array(sqrt_product_list)**2)

        self.A = A_fit
        self.B = B_fit 
        
        rho = np.array(rho_list)
        eta = np.array(eta_list)
        A_R = np.array(A_fit_list)
        B_R = np.array(B_fit_list)

        sigma_rho = np.array(rho_err_list)
        sigma_eta = np.array(eta_err_list)
        
        dx_drho = 0.5 * np.sqrt(eta / rho)
        dx_deta = 0.5 * np.sqrt(rho / eta)

        sigma_x = np.sqrt((dx_drho * sigma_rho)**2 +
                        (dx_deta * sigma_eta)**2)
        
        sigma_A = np.sqrt((np.sum((A_R - A_fit*rho)**2))/((len(self.CalibrationFilenames) - 1)*np.sum(rho**2)))

        sigma_B = np.sqrt((np.sum((B_R - B_fit*np.array(sqrt_product_list))**2))/((len(self.CalibrationFilenames) - 1)*np.sum(np.array(sqrt_product_list)**2)))
        from matplotlib import pyplot as plt 
        plt.errorbar(rho_list,A_fit_list,A_sigma_list,rho_err_list,fmt='o',markersize='3',color='firebrick')
        plt.plot(rho_list, np.array(rho_list)*A_fit,color='black',ls='--')
        plt.xlabel(r'Density ($gcm^{-3}$)',fontsize=15)
        plt.ylabel(r'Reduced $A_R$ ($A = A_R \rho$) (H)',fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.grid(alpha=0.4)
        plt.show()
        plt.errorbar(sqrt_product_list,B_fit_list,B_sigma_list,sigma_x,fmt='o',markersize='3',color='firebrick')
        plt.plot(sqrt_product_list, np.array(sqrt_product_list)*B_fit,color='black',ls='--')
        plt.xlabel(r'$\sqrt{\rho\eta}$ ($cPgcm^{-3}$)',fontsize=15)
        plt.ylabel(r'Reduced $B_R$ ($B = B_R \sqrt{\rho\eta}$) ($\Omega$)',fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.grid(alpha=0.4)
        plt.show()
        return A_fit, B_fit, sigma_A, sigma_B
    
    def calibrate_fluid_dynamics(self, rho_list, eta_list):
        """
        Robust calibration of A, B, C using Zhang hydrodynamic model
        Handles scaling, conditioning, and noise properly
        """
        import numpy as np
        from scipy.optimize import least_squares

        #Air calibration
        FreqData, X_data, Y_data = self.Read_XY_Voltage(self.AirCalibrationFile)

        Z_air = self.CalibrationResistorValue * (self.AirCalibrationDrivingVoltage/self.CalibrationResistorVoltage) \
                * self.CalibrationVoltage / (X_data + 1j*Y_data)

        Y_air = 1 / Z_air
        C0_air = self.parallel_capacitance_estimation(Y_air, FreqData)

        Y_motional = Y_air - 1j * 2*np.pi*C0_air*FreqData

        res = self.LorentzianFit(FreqData, np.real(Y_motional))
        f0 = res[0]
        Q0 = res[2]
        omega0 = 2*np.pi*f0

        omega1_list = []
        Q1_list = []

        y1_list = []
        y2_list = []
        for i, file in enumerate(self.CalibrationFilenames):

            FreqData, X_data, Y_data = self.Read_XY_Voltage(file)

            Z = self.CalibrationResistorValue * (self.CalibrationDrivingVoltages[i]/self.CalibrationResistorVoltage) \
                * self.CalibrationVoltage / (X_data + 1j*Y_data)

            Y = 1 / Z
            C0 = self.parallel_capacitance_estimation(Y, FreqData)

            Y_motional = Y - 1j * 2*np.pi*C0*FreqData

            res = self.LorentzianFit(FreqData, np.real(Y_motional))

            omega1 = 2*np.pi*res[0]
            omega1_list.append(omega1)
            Q1_list.append(res[2])

            #Zhang targets
            y1_list.append((1 / res[2]) * (omega0**2 / omega1**2))
            y2_list.append((omega0**2 / omega1**2) - 1)

        def residuals(params):
            A, B, C = params
            res_list = []
            for i in range(len(omega1_list)):
                u = np.sqrt(eta_list[i] * rho_list[i] / omega1_list[i])

                eq1 = A * u + B * (eta_list[i] / omega1_list[i]) - y1_list[i]
                eq2 = A * u + C * rho_list[i] - y2_list[i]

                res_list.append(eq1)
                res_list.append(eq2)
            return res_list
        
        initial_guess = [1, -30, 1]
        result = least_squares(
            residuals,
            x0=initial_guess,
            method='trf'
        )

        if not result.success:
            raise RuntimeError("Calibration failed: " + result.message)

        A, B, C = result.x
        J = result.jac  # Jacobian matrix
        residual = result.fun
        n = len(residual)          # number of data points
        p = len(result.x)          # number of parameters (3 here)

        # Residual variance (reduced chi-square estimate)
        sigma2 = np.sum(residual**2) / (n - p)

        # Covariance matrix
        cov = np.linalg.inv(J.T @ J) * sigma2

        # Standard errors (1-sigma uncertainties)
        errors = np.sqrt(np.diag(cov))

        A_err, B_err, C_err = errors
        print("\n--- Calibration Results ---")
        print("A =", A, "+-", A_err)
        print("B =", B, "+-", B_err)
        print("C =", C, "+-", C_err)

        #Save results
        self.A_Hydro = A
        self.B_Hydro = B
        self.C_Hydro = C
        return f0, Q0
    
    def calibrate_Reittinger(self,rho_list,eta_list):
        
        
        FreqData,X_data,Y_data = self.Read_XY_Voltage(self.AirCalibrationFile)
        Z_air = self.CalibrationResistorValue * (self.AirCalibrationDrivingVoltage/self.CalibrationResistorVoltage) * self.CalibrationVoltage / (X_data + 1j*Y_data) 
        #Admittance and impedance magnitudes
        Y_air = 1 / Z_air
        
        C0 = self.parallel_capacitance_estimation(Y_air,FreqData)
        Y_motional = np.array(Y_air) - 1j*2*np.pi*C0*np.array(FreqData)
        
        res = self.LorentzianFit(FreqData,np.real(Y_motional))
        omega_s_vac = res[0] * 2 * np.pi
        omega_45_vac = self.find_first_intersection(FreqData,Y_motional) * 2 * np.pi

        A_list = []
        B_list = []
        
        for i, file in enumerate(self.CalibrationFilenames):
            FreqData,X_data,Y_data = self.Read_XY_Voltage(file)
            Z_fluid = self.CalibrationResistorValue * (self.CalibrationDrivingVoltages[i]/self.CalibrationResistorVoltage) * \
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
            den1 = np.sqrt(rho_list[i]*eta_list[i])*(2/(np.sqrt(omega_45)) - 1/np.sqrt(omega_s))

            B_list.append(num1/den1)

            num2 = ((omega_s_vac)/(omega_s))**2 - (num1/den1)*(np.sqrt(rho_list[i]*eta_list[i]/omega_s)) - 1
            den2 = rho_list[i]

            A_list.append(num2/den2)
        
        A = np.mean(A_list)
        B = np.mean(B_list)
        A_std = np.std(A_list)
        B_std = np.std(B_list)
        from matplotlib import pyplot as plt
        plt.plot(FreqData*1e-3,np.real(Y_motional)*1e6,color='firebrick')
        plt.plot(FreqData*1e-3,np.imag(Y_motional)*1e6,color='mediumblue')
        plt.axvline(1e-3*omega_s/(2*np.pi),ls='-.')
        plt.axvline(1e-3*omega_45/(2*np.pi),ls='-.')
        plt.xlabel('Frequency (kHz)')
        plt.ylabel(r'Motional Admittance ($\mu$S)')
        plt.show()
        return A, A_std, B, B_std, omega_s_vac, omega_45_vac
        
            