import numpy as np
import scipy.optimize as scp

class ForkData():
    def __init__(self,Filenames,DrivingVoltages,CalibrationResistorValue,CalibrationResistorFile,
                 CalibrationResistorDrivingVoltage):
        self.Filenames = Filenames
        self.DrivingVoltages = DrivingVoltages
        R_FreqData, R_X_data, R_Y_data = self.Read_XY_Voltage(CalibrationResistorFile)
        self.CalibrationVoltage = np.mean(R_X_data) + 1j*np.mean(R_Y_data)
        self.CalibrationResistorVoltage = CalibrationResistorDrivingVoltage
        self.CalibrationResistorValue = CalibrationResistorValue
    
    @staticmethod
    def Read_XY_Voltage(file):
        """
        Function to read the X and Y data of a file generated 
        by Zurich instruments during a frequency sweep.
        This function is specific to the file type generated, though it
        likely can be adapted for other files of similar nature.

        Note to self - Potentially make a function for a more general file
        type with different variables or sweeps.
        """
        import numpy as np
        X_data = []
        Y_data = []
        freq_data = []
        #Open the file 
        with open(file, "r", encoding="utf-8") as f:
            for line in f:
                #Remove spaces from the file
                line = line.strip()
                # Skip comments and empty lines
                if not line or line.startswith("%"):
                    continue

                # Split by semicolon
                parts = [p.strip() for p in line.split(";")]
                if len(parts) != 3:
                    continue  # skip malformed lines

                try:
                    freq = float(parts[0])
                    x_val = float(parts[1])
                    y_val = float(parts[2])

                    freq_data.append(freq)
                    X_data.append(x_val)
                    Y_data.append(y_val)
                except ValueError:
                    print("ValueError")
                    # Skip lines that can't be converted to floats
                    continue
        X_data = np.array(X_data)
        Y_data = np.array(Y_data)
        freq_data = np.array(freq_data)
        return freq_data, X_data, Y_data
    
    @staticmethod
    def LorentzianFit(FreqData, ReYData, plot=False):
        """
        Fits a Lorentzian to frequency vs. magnitude data.

        Returns:
        ResFreq, Sigma_ResFreq, Q, Sigma_Q, ResPeak, Sigma_ResPeak
        """

        # --- Lorentzian function ---
        def lorentzian(f, y0, A, f0, Gamma):
            return y0 + A * ((Gamma/2)**2 / ((f - f0)**2 + (Gamma/2)**2))

        # --- Initial guesses ---
        guess_y0 = np.median(ReYData)
        guess_A = max(ReYData) - guess_y0
        guess_f0 = FreqData[np.argmax(ReYData)]

        # Half-max width estimation for Gamma
        half_max = guess_y0 + guess_A/2
        indices_above_half = np.where(ReYData > half_max)[0]
        if len(indices_above_half) >= 2:
            f_left = FreqData[indices_above_half[0]]
            f_right = FreqData[indices_above_half[-1]]
            guess_Gamma = f_right - f_left
        else:
            guess_Gamma = (max(FreqData) - min(FreqData)) / 50  # fallback
        p0 = [guess_y0, guess_A, guess_f0, guess_Gamma]

        # --- Bounds ---
        lower_bounds = [0, 0, min(FreqData), 0]
        upper_bounds = [np.inf, np.inf, max(FreqData), np.inf]

        # --- Weight near peak to emphasize it ---
        weights = 1 + 5 * np.exp(-((FreqData - guess_f0)**2)/(2*(guess_Gamma/2)**2))

        # --- Fit ---
        popt, pcov = scp.curve_fit(
            lorentzian, FreqData, ReYData, p0=p0,
            sigma=1/weights, absolute_sigma=False,
            bounds=(lower_bounds, upper_bounds), maxfev=50000
        )
        perr = np.sqrt(np.diag(pcov))

        # --- Extract parameters ---
        y0, A, ResFreq, Gamma = popt
        sigma_y0, sigma_A, Sigma_ResFreq, sigma_Gamma = perr

        Q = ResFreq / Gamma
        Sigma_Q = Q * np.sqrt((Sigma_ResFreq/ResFreq)**2 + (sigma_Gamma/Gamma)**2)

        ResPeak = y0 + A
        Sigma_ResPeak = np.sqrt(sigma_y0**2 + sigma_A**2)

        # --- Optional plot ---
        from matplotlib import pyplot as plt
        if plot:
            plt.figure(figsize=(6,4))
            plt.plot(FreqData, ReYData, 'b.', label='Data')
            plt.plot(FreqData, lorentzian(FreqData, y0, A, ResFreq, Gamma), 'r-', label='Lorentzian Fit')
            plt.xlabel('Frequency [Hz]')
            plt.ylabel('Re(Admittance) [S]')
            plt.title(f'QTF Fit: f0={ResFreq:.3f} Hz, Q={Q:.1f}')
            plt.grid(alpha=0.4)
            plt.legend()
            plt.show()

        return ResFreq, Sigma_ResFreq, Q, Sigma_Q, ResPeak, Sigma_ResPeak
    
    @staticmethod
    def parallel_capacitance_estimation(Y_data,freq_data):
        """
        Takes Y_data as the complex admittance spectra against frequency sweep of
        freq_data
        """
        return sum(np.imag(np.array(Y_data)) / (2*np.pi*np.array(freq_data))) / len(freq_data)
    
    @staticmethod
    def find_first_intersection(freq, Y):
        freq = np.array(freq)
        Y = np.array(Y)
        
        diff_val = np.real(Y) - np.imag(Y)
        crossings = []

        # Loop through consecutive points
        for i in range(len(diff_val)-1):
            # Exact zero at a point
            if diff_val[i] == 0:
                crossings.append(freq[i])
            # Sign change indicates crossing between points
            elif diff_val[i] * diff_val[i+1] < 0:
                # Linear interpolation for zero crossing
                f_cross = freq[i] + (freq[i+1] - freq[i]) * (-diff_val[i]) / (diff_val[i+1] - diff_val[i])
                crossings.append(f_cross)

        if len(crossings) == 0:
            raise RuntimeError("No intersection found where Re(Y) = Im(Y)")
    
        return crossings[0]  # Return the first crossing