# schaf3 — ASC-ODE: Mechanical Systems + ODE Solver Library

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)]()

_A C++ library for solving ODEs with automatic differentiation, plus a 3D mass-spring-constraint system simulator with Python bindings and interactive visualization._

### Sebastian Wagner, 12307039

### Victor Pawlek, 12320223 (already has done SciComp)

### Angelo Stjepanovic, 12302549 (already has done SciComp)

### Carlos Aleksandar Duenas Lopez, 12308486 (already has done SciComp)

### Fabian Senftner, 12307538 (already has done SciComp)

---

## 🚀 What is schaf3?

**schaf3** (ASC-ODE) is a C++ library designed to help you:

- **Simulate mechanical systems** — masses, springs, and rigid distance constraints in 2D/3D with Python bindings
- **Solve ODEs** — multiple time-stepping schemes (Euler, Crank-Nicolson, Runge-Kutta, Newmark, generalized-α)
- **Automatic differentiation** — compute derivatives and Jacobians without manual calculus
- **Visualize in real-time** — Python integration with `pythreejs` for interactive 3D rendering in Jupyter notebooks

Under the hood, it uses implicit time-stepping methods to handle stiff ODEs and rigid constraints, with automatic differentiation for computing Jacobians in the Newton solver. Perfect for scientific computing coursework, mechanical system prototyping, or experimenting with physics simulations.

---

## 📦 Quick Example

```python
from mass_spring import *

# Create a pendulum
mss = MassSpringSystem3d()
mss.gravity = (0, 0, -9.81)

fix = mss.add(Fix((0, 0, 1.0)))                    # anchor at top
mass = mss.add(Mass(1.0, (0.1, 0, 0.0)))           # hanging mass
mss.add(Spring(1.0, 100.0, (fix, mass)))           # spring connection

# Simulate
mss.simulate(dt=0.01, substeps=100)
print(mass.pos, mass.vel)
```

### Springs vs. Constraints

```python
# Elastic spring (can stretch)
mss.add(Spring(length, stiffness, (conn1, conn2)))

# Rigid constraint (perfectly fixed distance)
mss.add(DistanceConstraint(length, (conn1, conn2)))
```

### Automatic Differentiation (C++)

```cpp
AutoDiff<2> x = Variable<0>(2.0);  // x = 2, dx/dx = 1
AutoDiff<2> y = Variable<1>(3.0);  // y = 3, dy/dy = 1
auto f = x * sin(y);               // compute f and all derivatives simultaneously
```

---

## 🛠️ Getting Started

### Requirements

- **C++ Compiler:** C++17 or newer (GCC, Clang, MSVC)
- **CMake:** ≥ 3.20
- **For Python bindings:** Python ≥ 3.8, `pybind11`
- **For visualization:** `numpy`, `pythreejs` (for Jupyter notebooks)

### Installation

```bash
# Clone the repository
git clone https://github.com/Angelo2510/schaf3
cd schaf3

# Build with CMake
mkdir build && cd build
cmake ..
cmake --build .
```

### Try It Out

```bash
# Run C++ ODE demos
./test_ode
./demo_autodiff

# Python mass-spring simulator
python mechsystem/test_mass_spring.py

# Jupyter notebooks with visualization
jupyter notebook mechsystem/finale.ipynb
```

---

## 📁 What's Inside

```
src/               ODE solvers, time steppers, automatic differentiation
mechsystem/        Mass-spring simulator with Python bindings
  ├── mass_spring.hpp/cpp      C++ implementation
  ├── bind_mass_spring.cpp     Python bindings (pybind11)
  ├── finale.ipynb             Demo: pendulums (springs vs constraints)
  └── kreisel.ipynb            Demo: spinning top
nanoblas/          Lightweight linear algebra
demos/             C++ example programs
```

---

## 🎯 Features

**Mechanical Systems:**
- 2D/3D masses, fixed anchors, springs, distance constraints
- Gravity and external forces
- Access to positions, velocities, masses
- Python bindings for easy scripting

**ODE Solvers (C++):**
- Explicit/Implicit Euler, Crank-Nicolson, Runge-Kutta
- Newmark and generalized-α for second-order ODEs
- Template-based automatic differentiation

**Visualization:**
- Real-time 3D rendering in Jupyter with `pythreejs`
- Interactive orbit controls
- See examples in `mechsystem/*.ipynb`

---

## 📚 Example Notebooks

- **finale.ipynb** — Compare springs vs constraints with single and double pendulums
- **kreisel.ipynb** — Spinning top with rigid body constraints
- **mass_spring.ipynb** — Basic mass-spring chain

---

## 🤔 Why not just use NGSolve?

Because sometimes you want to simulate a pendulum without compiling a finite element library the size of a small planet.

_(But seriously, if you need industrial-grade FEM: [NGSolve](https://ngsolve.org). This is for learning and fun.)_

---

## 📄 License

MIT License — see [LICENSE](LICENSE)
