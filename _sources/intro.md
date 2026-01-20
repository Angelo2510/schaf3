# Mass-Spring System Simulator

Simulate 3D mechanical systems with masses, springs, and rigid distance constraints. Perfect for pendulums, spinning tops, and other mechanical systems.

Under the hood, the library uses implicit time-stepping methods (Newmark/generalized-α) to handle stiff ODEs and rigid constraints, with automatic differentiation for computing Jacobians in the Newton solver. But you don't need to worry about that — just build your system and hit simulate.

## What can you do?

Build mechanical systems with:
- **Masses** — objects with mass that move under forces
- **Fixed points** — anchors that don't move
- **Springs** — elastic connections with stiffness and rest length
- **Distance constraints** — perfectly rigid connections (no stretch)

Simulate physics in real-time and visualize in Jupyter notebooks with interactive 3D rendering.

## Quick Example

```python
from mass_spring import *

# Create a system with gravity
mss = MassSpringSystem3d()
mss.gravity = (0, 0, -9.81)

# Add a fixed anchor point
fix = mss.add(Fix((0, 0, 1.0)))

# Add a mass
mass = mss.add(Mass(1.0, (0.1, 0, 0.0)))  # mass, (x, y, z)

# Connect with a spring
mss.add(Spring(1.0, 100.0, (fix, mass)))  # length, stiffness, connectors

# Simulate
mss.simulate(dt=0.01, substeps=100)

# Access current state
print(mass.pos)  # position
print(mass.vel)  # velocity
```

## Key Concepts

### Creating a System

```python
# 2D or 3D
mss = MassSpringSystem2d()
mss = MassSpringSystem3d()

# Set gravity
mss.gravity = (0, 0, -9.81)  # x, y, z components
```

### Adding Components

**Fixed points** (anchors):
```python
fix = mss.add(Fix((x, y, z)))
```

**Masses** (position set at creation, cannot be changed later):
```python
mass = mss.add(Mass(mass_value, (x, y, z)))
```

**Springs** (elastic, can stretch):
```python
spring = mss.add(Spring(rest_length, stiffness, (connector1, connector2)))
```

**Distance constraints** (rigid, no stretch):
```python
constraint = mss.add(DistanceConstraint(distance, (connector1, connector2)))
```

### Accessing State

Components are stored in lists:
```python
for m in mss.masses:
    print(m.pos)  # read-only
    print(m.vel)  # can set: m.vel = (vx, vy, vz)
    print(m.mass)  # can modify: m.mass = 2.0

for f in mss.fixes:
    print(f.pos)  # fixed position

for s in mss.springs:
    print(s.length, s.stiffness)
    print(s.connectors)  # tuple of two connectors
```

Get position of any connector:
```python
pos = mss[connector].pos  # works for both masses and fixes
```

### Simulation

```python
# dt: time step
# substeps: solver iterations per step (higher = more stable for constraints)
mss.simulate(dt=0.01, substeps=100)
```

## Visualization with pythreejs

```python
from pythreejs import *

# Create sphere for each mass
masses_vis = []
for m in mss.masses:
    sphere = Mesh(
        SphereBufferGeometry(0.15, 16, 16),
        MeshStandardMaterial(color='red'),
        position=m.pos
    )
    masses_vis.append(sphere)

# Create lines for springs/constraints
positions = []
for s in mss.springs:
    pA = mss[s.connectors[0]].pos
    pB = mss[s.connectors[1]].pos
    positions.append([pA, pB])

lines = LineSegments2(
    LineSegmentsGeometry(positions=positions),
    LineMaterial(linewidth=2, color='cyan')
)

# Setup scene
camera = PerspectiveCamera(position=[3, 3, 3], aspect=16/9)
scene = Scene(children=[*masses_vis, lines, camera, DirectionalLight()])
renderer = Renderer(camera=camera, scene=scene, width=800, height=600)
renderer
```

Update in simulation loop:
```python
for i in range(1000):
    mss.simulate(0.01, 100)
    
    # Update positions
    for m, vis in zip(mss.masses, masses_vis):
        vis.position = m.pos
    
    # Update lines
    positions = [[ mss[s.connectors[0]].pos, mss[s.connectors[1]].pos ] 
                 for s in mss.springs]
    lines.geometry = LineSegmentsGeometry(positions=positions)
```

## Examples

See the Jupyter notebooks:
- `mechsystem/mass_spring.ipynb` — Basic chain of masses
- `mechsystem/finale.ipynb` — Pendulums with springs vs. constraints
- `mechsystem/kreisel.ipynb` — Spinning top with rigid body constraints

## Installation

Build from source:
```bash
git clone <repository-url>
cd schaf3
mkdir build && cd build
cmake ..
cmake --build .
```

Requirements: Python ≥3.8, `pybind11`, `numpy`, `pythreejs`

---

**Note:** This is a teaching tool for understanding mechanical systems and constraints. If you need industrial-grade multibody dynamics, there are... other options. But where's the fun in that? 🎢 (Looking at you, NGSolve.)




   
