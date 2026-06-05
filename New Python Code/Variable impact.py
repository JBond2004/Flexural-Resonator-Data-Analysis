#Test to see impact of different variables
import matplotlib.pyplot as plt
import numpy as np
from Read_Fork_Data import ForkData


def total_impedance(f, 
                    Rm, Cm, Lm,
                    C0,
                    A,B,
                    rho, eta):
            w = 2 * np.pi * f

            Zm = Rm + 1j*w*Lm - 1j/(w*Cm)
            Z_hyd = 1j*A*w*rho + (1+1j)*B*np.sqrt(w*rho*eta)

            Y_total = 1/(Zm + Z_hyd) + 1j*w*C0
            return (1 / Y_total)

#Test C0
freq_data = np.linspace(22e3,26e3,int(1e4),dtype=np.float64)

Rm = 11720
Lm = 141.2
Cm = 0.2601e-12
C0 = 4.0e-12
A = 26.05
B = 244.3
rho_real = 1
eta_real = 1
real_density_list = [(900 + 5*i)*1e-3 for i in range(50)]
real_viscosity_list = [0.1 + 0.1*i for i in range(50)]

res_freq_C0 = []
sigma_res_freq_C0 = []
antires_freq_C0 = []
sigma_antires_freq_C0 = []
res_imp_C0 = []
sigma_res_imp_C0 = []
antires_imp_C0 = []
sigma_antires_imp_C0 = []
Q_C0 = []
sigma_Q_C0 = []

res_freq_rho = []
sigma_res_freq_rho = []
antires_freq_rho = []
sigma_antires_freq_rho = []
res_imp_rho = []
sigma_res_imp_rho = []
antires_imp_rho = []
sigma_antires_imp_rho = []
Q_rho = []
sigma_Q_rho = []

res_freq_eta = []
sigma_res_freq_eta = []
antires_freq_eta = []
sigma_antires_freq_eta = []
res_imp_eta = []
sigma_res_imp_eta = []
antires_imp_eta = []
sigma_antires_imp_eta = []
Q_eta = []
sigma_Q_eta = []

C0_list = []


for i in range(10,21):
    C0 = i * 2.5e-13
    C0_list.append(C0 * 1e12)
    Z_data = [total_impedance(f,Rm,Cm,Lm,C0,A,B,rho_real,eta_real) for f in freq_data]
    Y_data = 1/np.array(Z_data)
    Y_m = Y_data - 1j*np.array(freq_data)*2*np.pi*C0
    Z_m = 1/Y_m

    res_series = ForkData.LorentzianFit(freq_data,np.real(Y_m))
    res_parallel = ForkData.LorentzianFit(freq_data,np.abs(Z_data))

    res_freq_C0.append(res_series[0])
    sigma_res_freq_C0.append(res_series[1])
    antires_freq_C0.append(res_parallel[0])
    sigma_antires_freq_C0.append(res_parallel[1])
    res_imp_C0.append(1/res_series[4])
    sigma_res_imp_C0.append(res_series[5]/(res_series[4]**2))
    antires_imp_C0.append(res_parallel[4])
    sigma_antires_imp_C0.append(res_parallel[5])
    Q_C0.append(res_series[2])
    sigma_Q_C0.append(res_series[3])


"""plt.errorbar(C0_list,res_freq_C0,sigma_res_freq_C0,fmt='o',
             label='Series Resonant Frequency',markersize=3,color='firebrick')
plt.errorbar(C0_list,antires_freq_C0,sigma_antires_freq_C0,
             fmt='o',label='Parallel Resonant Frequency',markersize=3,color='mediumblue')
plt.legend(fontsize=14)
plt.grid(alpha=0.4)
plt.xlabel('Parallel capacitance (pF)',fontsize=14)
plt.ylabel('Frequency (Hz)',fontsize=14)
plt.tight_layout()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.show()
plt.errorbar(C0_list,np.abs(np.array(res_imp_C0) + np.array(res_freq_C0)*2*np.pi*np.array(C0_list))*1e-6,np.array(sigma_res_imp_C0)*1e-6,fmt='o',
                label='Resonant Impedance',markersize=3,color='firebrick')
    
plt.errorbar(C0_list,np.array(antires_imp_C0)*1e-6,np.array(sigma_antires_imp_C0)*1e-6,fmt='o',
             label='Parallel Resonant Impedance',markersize=3,color='mediumblue')
plt.legend(fontsize=14)
plt.grid(alpha=0.4)
plt.xlabel('Parallel capacitance (pF)',fontsize=15)
plt.ylabel(r'Impedance (M$\Omega$)',fontsize=15)
plt.tight_layout()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()

plt.errorbar(C0_list,Q_C0,sigma_Q_C0,fmt='o',
            markersize=3,color='firebrick')
plt.grid(alpha=0.4)
plt.xlabel('Parallel capacitance (pF)',fontsize=15)
plt.ylabel('Q-factor',fontsize=15)
plt.ylim((Q_C0[0]-10,Q_C0[0]+10))
plt.ticklabel_format(style='plain', axis='y', useOffset=False)
plt.tight_layout()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()"""
#Density
C0 = 4.0e-12

for i in range(len(real_density_list)):
    rho_real = real_density_list[i]
    Z_data = [total_impedance(f,Rm,Cm,Lm,C0,A,B,rho_real,eta_real) for f in freq_data]
    Y_data = 1/np.array(Z_data)
    Y_m = Y_data - 1j*np.array(freq_data)*2*np.pi*C0
    Z_m = 1/Y_m

    res_series = ForkData.LorentzianFit(freq_data,np.real(Y_m))
    res_parallel = ForkData.LorentzianFit(freq_data,np.abs(Z_data))

    res_freq_rho.append(res_series[0])
    sigma_res_freq_rho.append(res_series[1])
    antires_freq_rho.append(res_parallel[0])
    sigma_antires_freq_rho.append(res_parallel[1])
    res_imp_rho.append(1/res_series[4])
    sigma_res_imp_rho.append(res_series[5]/(res_series[4]**2))
    antires_imp_rho.append(res_parallel[4])
    sigma_antires_imp_rho.append(res_parallel[5])
    Q_rho.append(res_series[2])
    sigma_Q_rho.append(res_series[3])

deltaf = [26251.746 - f for f in res_freq_rho]
"""
plt.errorbar(real_density_list,res_freq_rho,sigma_res_freq_rho,fmt='o',label='Series Resonant Frequency')
plt.errorbar(real_density_list,antires_freq_rho,sigma_antires_freq_rho,fmt='o',label='Parallel Resonant Frequency')
plt.legend(fontsize=15)
plt.grid(alpha=0.4)
plt.xlabel(r'Mass Density (gcm$^{-3}$)',fontsize=15)
plt.ylabel('Frequency (Hz)',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()"""

"""fig, axs = plt.subplots(2, 1, sharex=True)

axs[1].errorbar(real_density_list, np.abs(np.array(res_imp_rho) + np.array(res_freq_rho)*2*np.pi*C0)*1e-6,
                np.array(sigma_res_imp_rho)*1e-6, fmt='o',
                label='Resonant Impedance', markersize=3, color='firebrick')

axs[0].errorbar(real_density_list, np.array(antires_imp_rho)*1e-6,
                np.array(sigma_antires_imp_rho)*1e-6, fmt='o',
                label='Parallel Resonant Impedance', markersize=3, color='mediumblue')

axs[0].grid(alpha=0.4)
axs[1].grid(alpha=0.4)
axs[0].ticklabel_format(style='plain', axis='y')
axs[1].ticklabel_format(style='plain', axis='y')
axs[1].set_xlabel(r'Mass Density (gcm$^{-3}$)', fontsize=15)
axs[0].set_ylabel(r'Impedance (M$\Omega$)', fontsize=15)
axs[1].set_ylabel(r'Impedance (M$\Omega$)', fontsize=15)

# Control tick sizes
axs[0].tick_params(axis='both', labelsize=12)
axs[1].tick_params(axis='both', labelsize=15)

# Remove top x-axis labels
axs[0].tick_params(axis='x', labelbottom=False)

plt.show()

# Convert to arrays (just in case)
rho = np.array(real_density_list)
Q = np.array(Q_rho)
Q_err = np.array(sigma_Q_rho)

# --- Fit in log-log space ---
log_rho = np.log10(rho)
log_Q = np.log10(Q)

# Linear fit: log(Q) = m*log(rho) + c
m, c = np.polyfit(log_rho, log_Q, 1)

# Generate smooth fit line
rho_fit = np.linspace(min(rho), max(rho), 200)
Q_fit = 10**(m * np.log10(rho_fit) + c)

plt.errorbar(real_density_list,Q_rho,sigma_Q_rho,fmt='o',
            markersize=3,color='firebrick')
plt.plot(rho_fit, Q_fit,
         color='black',
         linestyle='--',
         label=f'Fit: slope = {m:.2f}, intercept = {c:.2f}')
plt.grid(alpha=0.4)
plt.xlabel(r'Mass Density (gcm$^{-3}$)',fontsize=15)
plt.ylabel('Q-factor',fontsize=15)
plt.ticklabel_format(style='plain', axis='y', useOffset=False)
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()"""
#Viscosity
rho_real = 1

rho_eta_list=[]
for i in range(len(real_viscosity_list)):
    eta_real = real_viscosity_list[i]
    rho_eta_list.append(1/np.sqrt(rho_real*real_viscosity_list[i]))
    Z_data = [total_impedance(f,Rm,Cm,Lm,C0,A,B,rho_real,eta_real) for f in freq_data]
    Y_data = 1/np.array(Z_data)
    Y_m = Y_data - 1j*np.array(freq_data)*2*np.pi*C0
    Z_m = 1/Y_m

    res_series = ForkData.LorentzianFit(freq_data,np.real(Y_m))
    res_parallel = ForkData.LorentzianFit(freq_data,np.abs(Z_data))

    res_freq_eta.append(res_series[0])
    sigma_res_freq_eta.append(res_series[1])
    antires_freq_eta.append(res_parallel[0])
    sigma_antires_freq_eta.append(res_parallel[1])
    res_imp_eta.append(1/res_series[4])
    sigma_res_imp_eta.append(res_series[5]/(res_series[4]**2))
    antires_imp_eta.append(res_parallel[4])
    sigma_antires_imp_eta.append(res_parallel[5])
    Q_eta.append(res_series[2])
    sigma_Q_eta.append(res_series[3])

s_p_frequency = [abs(res_freq_eta[i] - antires_freq_eta[i]) for i in range(len(res_freq_eta))]
deltaf = [26251.746 - f for f in res_freq_eta]

"""fig,axs = plt.subplots(2,1)

axs[1].errorbar(real_viscosity_list,np.array(res_freq_eta)*1e-3,np.array(sigma_res_freq_eta)*1e-3,fmt='o',color='firebrick')
axs[0].errorbar(real_viscosity_list,np.array(antires_freq_eta)*1e-3,np.array(sigma_antires_freq_eta)*1e-3,fmt='o',color='mediumblue')
axs[0].grid(alpha=0.4)
axs[1].grid(alpha=0.4)
axs[1].set_xlabel('Viscosity (cP)',fontsize=15)
axs[0].set_ylabel('Frequency (kHz)',fontsize=15)
axs[1].set_ylabel('Frequency (kHz)',fontsize=15)
axs[0].ticklabel_format(useOffset=False)
axs[1].ticklabel_format(useOffset=False)

# Remove top x-axis labels
axs[0].tick_params(axis='x', labelbottom=False)

plt.tight_layout()
plt.show()"""

"""fig, axs = plt.subplots(2,1)

axs[1].errorbar(real_viscosity_list,np.array(res_imp_eta)*1e-6,np.array(sigma_res_imp_eta)*1e-6,fmt='o',
             label='Resonant Impedance',markersize=3,color='firebrick')
axs[0].errorbar(real_viscosity_list,np.array(antires_imp_eta)*1e-6,np.array(sigma_antires_imp_eta)*1e-6,fmt='o',
             label='Parallel Resonant Impedance',markersize=3,color='mediumblue')

axs[0].grid(alpha=0.4)
axs[1].grid(alpha=0.4)
axs[1].set_xlabel('Viscosity (cP)',fontsize=15)

axs[0].set_ylabel(r'Impedance (M$\Omega$)',fontsize=15)
axs[1].set_ylabel(r'Impedance (M$\Omega$)',fontsize=15)

# Remove top x-axis labels
axs[0].tick_params(axis='x', labelbottom=False)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
plt.show()

plt.errorbar(real_viscosity_list,Q_eta,sigma_Q_eta,fmt='o',
            markersize=3,color='firebrick')
plt.grid(alpha=0.4)
plt.xlabel('Viscosity (cP)',fontsize=15)
plt.ylabel('Q-factor',fontsize=15)
plt.ticklabel_format(style='plain', axis='y', useOffset=False)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#plt.yscale('log')
#plt.xscale('log')
plt.show()"""

rho_vals = np.array(real_density_list)
eta_vals = np.array(real_viscosity_list)

RHO, ETA = np.meshgrid(rho_vals, eta_vals)

res_freq_grid = np.zeros_like(RHO)
Q_grid = np.zeros_like(RHO)
for i in range(len(eta_vals)):
    print(i)
    for j in range(len(rho_vals)):
        rho_real = rho_vals[j]
        eta_real = eta_vals[i]

        Z_data = [total_impedance(f, Rm, Cm, Lm, C0, A, B, rho_real, eta_real)
                  for f in freq_data]

        Y_data = 1 / np.array(Z_data)
        Y_m = Y_data - 1j * np.array(freq_data) * 2 * np.pi * C0

        res_series = ForkData.LorentzianFit(freq_data, np.real(Y_m))

        # Store results
        res_freq_grid[i, j] = res_series[0]
        Q_grid[i, j] = res_series[2]

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(RHO, ETA, res_freq_grid,cmap='viridis')

ax.set_xlabel('Density (g cm$^{-3}$)', fontsize=15, labelpad=15)
ax.set_ylabel('Viscosity (cP)', fontsize=15, labelpad=15)
ax.set_zlabel('Resonant Frequency (Hz)', fontsize=15, labelpad=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(RHO, ETA, Q_grid,cmap='viridis')

ax.set_xlabel('Density (g cm$^{-3}$)', fontsize=15, labelpad=15)
ax.set_ylabel('Viscosity (cP)', fontsize=15, labelpad=15)
ax.set_zlabel('Q-factor', fontsize=15, labelpad=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.show()