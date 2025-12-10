# schaf3 — ASC-ODE: A C++ Package for Solving Ordinary Differential Equations

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)]()

_A lightweight, header-only (mostly) C++ library for ODE integration and automatic differentiation, plus sample usage / plotting tools._

---

## 🚀 What is schaf3?

schaf3 (ASC-ODE) is a C++ library designed to help you:

-   Define and solve ordinary differential equations (ODEs) in a flexible way.
-   Use automatic differentiation (AD) to compute derivatives — enabling e.g. sensitivity analysis or Jacobians without manual calculus.
-   Easily switch between plain numeric types (e.g. `double`) and AD-enabled types, thanks to templated code.
-   Run demos that export data and visualize results using Python + Matplotlib (e.g. for polynomial examples or ODE trajectories).

In short: **schaf3** aims to make solving ODEs + computing derivatives easier and more expressive, without sacrificing performance or readability.

---

## 📦 What’s inside

Features:

-   Templated `AutoDiff<N, T>` types (for arbitrary dimension `N`)
-   Basic operations: addition, multiplication, division, plus standard math (`sin`, `cos`, `exp`, …) for AD types
-   Simple `Variable<I, T>` helper for initializing AD variables
-   Sample code to export data and use Python + Matplotlib to plot values and derivatives

---

## 🛠️ Getting Started

### Requirements

-   A C++17 (or newer) compiler
-   CMake (if you want to build demos)
-   For demos using plotting: a working Python (≥ 3.7), plus `numpy` and `matplotlib` installed

### Build and Run (C++ only)

```bash
git clone https://github.com/YourUsername/schaf3.git
cd schaf3
mkdir build && cd build
cmake ..
make
# Then run any demo binary, e.g.:
./demos/legendre_demo   # or whatever demo names you have
```
