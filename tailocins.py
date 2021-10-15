import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Function implementing the differential equations rules.
def equations(t, y):
    # Get values from the tuple.
    P, T = y
    # Compute differential values.
    
    if t==100:
        T=1
    dP = (alpha * (1 - P) - beta * T - delta_P) * P
    dT = -delta_T *T
    # Return differential values.
    return (dP, dT)

# Define parameters.
alpha   = 0.1
beta    = 1
delta_P = 0.001
delta_T = 0.1

# Define time span and evaluation points.S
t_span = (0., 500.)
t_eval = np.arange(t_span[0], t_span[1], 0.01)

# Define initial conditions: (P: pseudomonas, T: tailocins).
y0 = (1., 1.)

# Solve differential equation.
# See documentation at: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
sol = solve_ivp(equations, t_span, y0, t_eval=t_eval)

# Plot results.
fig, ax = plt.subplots()
ax.plot(sol.t, sol.y[0], label="pseudomonas (P)")
ax.plot(sol.t, sol.y[1], label="tailocins (T)")

ax.set_xlabel("time (a.u)")
ax.set_ylabel("P, T")
ax.legend()
plt.annotate(f'\u03B1={alpha}, \u03B2={alpha}, \u0394P={delta_P}, \u0394T={delta_T}', xy=(0.25, 1.1), xytext=(12, -12), va='top', xycoords='axes fraction', textcoords='offset points')
plt.show()