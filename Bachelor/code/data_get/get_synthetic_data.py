import numpy as np
from scipy.integrate import odeint

'''
direction_rules = {
    (N1, N2): 1,
    (N1, N3): 1,
    (N2, N4): 1,
    (N3, N4): -1,
    (N4, N5): 1,
    (N5, N6): 1,
    (N6, N1): -1,
}
'''

'''Equation
dx/dt = alpha_{i} * (Sigma_{j = 1}^{6} W_{ji} * x_{j}^{n} / K_{j}^{n} + x_{j}^{n}) - beta_{i} * x_{i} + sigma * epsilon_{t}
'''

# Define weights for the edges
W = np.array([
    [0, 0.8, 0.6, 0, 0, -0.3],
    [0, 0, 0, -0.9, 0, 0],
    [0, 0, 0, 0, 0.7, 0],
    [0, 0, 0, 0, 0, 0.5],
    [0, 0, 0, 0, 0, 0]
])

# Define parameters
alpha = np.array([1.0, 1.2, 0.9, 1.0, 0.8, 0.7])
beta = np.array([0.1, 0.3, 0.4, 0.2, 0.3, 0.2])
K = 0.5
n = 2
# Define original x values
x0 = np.array([0.8, 0.2, 0.3, 0.1, 0.0, 0.0])

# Define ODE system
def model(x, t):
    dxdt = np.zeros(6)
    for i in range(6):

        # define Hill function
        Hill = sum(W[j, i] * (x[j] ** n / (K ** n + x[j] ** n)) for j in range(5))
        # dynamic parameters
        if 30  <= t <= 50 and i == 0:
            alpha_i = 2.0
        else:
            alpha_i = alpha[i]
        if t >= 70 and i == 3:
            Hill += (-0.9 - (-0.2)) * (x[2] ** n / (K ** n + x[2] ** n))
        
        dxdt[i] = alpha_i * Hill - beta[i] * x[i] + 0.1 * np.random.normal()
    
    return dxdt

t = np.linspace(0, 100, 100)
X = odeint(model, x0, t)

print(X)

