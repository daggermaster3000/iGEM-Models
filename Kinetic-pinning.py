"""
Created on Thu Aug 19 17:05:31 2021

@author: quillan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#define model equations
def model_eq_14(concentration,chi):
    #Use model equation 14 for proteins with c>a and c>b
    b = 69.171
    c = 178.764
    p = 1   #probability of binding to ice plane
    l0 = 1  #constant (not sure about the unit)
    lambd = 12000 #not sure about the units but fits the paper's graph
    log_term_1 = ((b+c)/(2*b)+1)
    return((chi*b)*(lambd/c)**(1/4)*(np.sqrt(2*np.pi*(np.log10(log_term_1)*l0*p*concentration))))

def model_eq_13(concentration,chi):
    a = 15
    b = 12
    p = 0.55
    lambd = 12000                    #not sure about the units but fits the paper's graph
    beta = np.log10((a+b)/2*b)
    l0 = 1
    r1 = (chi*b)*(lambd/a)**(1/4)
    r2 = np.sqrt(2*np.pi*(beta+1)*l0*p*concentration)
    right_term = r1 * r2
    
    return(right_term)

#experimental data
x_data = [0.5, 1, 2]
y_data = [1.33, 2.54, 3.29]

#concentrations
Conc = np.arange(0,3,0.1)   
 

#get best fit
best_fit_ab, covar = curve_fit(model_eq_14, x_data, y_data)
sigma_ab = np.sqrt(np.diagonal(covar))
chi = str(round(best_fit_ab[0],3))
#print(f'best fit param: {chi}')

#get TH for each concentration

TH = model_eq_14(Conc,*best_fit_ab)

#get confidence interval for fitting curve
bound_upper = model_eq_14(Conc, *(best_fit_ab + sigma_ab)) 
bound_lower = model_eq_14(Conc, *(best_fit_ab - sigma_ab)) 




# Plot results.
fig, ax = plt.subplots()
ax.text(1.5, 1, f'best \u03C7: {chi}', style='italic', bbox={'facecolor': 'grey', 'alpha': 0.1, 'pad': 10})
ax.plot(Conc,TH , label="ffibp-HIS")        #Model
ax.scatter(x=x_data, y=y_data,facecolor = 'silver', edgecolor = 'k', s = 10, alpha = 1, label='experimental values')    #experimental values
plt.fill_between(Conc, bound_lower, bound_upper,color = 'black', alpha = 0.15)  #confidence interval
ax.set_xlabel("C (mM)")
ax.set_ylabel("TH")
ax.legend()

ax.set(title='TH as function of concentration')
plt.show()