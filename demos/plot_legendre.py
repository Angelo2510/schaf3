import numpy as np
import matplotlib.pyplot as plt

# Load data (skip the header row)
data = np.loadtxt("legendre_data.txt", skiprows=1)

x = data[:, 0]
P = data[:, 1:]  # columns P0..P5

degrees = P.shape[1]  # should be 6 (0..5)

plt.figure()

for n in range(degrees):
    plt.plot(x, P[:, n], label=f"P{n}(x)")

plt.xlabel("x")
plt.ylabel("P_n(x)")
plt.title("Legendre Polynomials P0..P5")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()