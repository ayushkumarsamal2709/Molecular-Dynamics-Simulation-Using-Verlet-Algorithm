# Molecular Dynamics Simulation Using Verlet Algorithm (Scilab)

## 📌 Overview

This project presents a **Molecular Dynamics (MD) simulation implemented in Scilab**, where particles evolve over time under inter-particle forces. The simulation demonstrates how an initially ordered configuration evolves into a **distorted/disordered state** due to interactions.

---

## 🧠 Theory

### Newton’s Equation of Motion

The motion of each particle is governed by:

F = m a

For N-particle system:
m d²r/dt² = F

---

### Interatomic Interaction

Particles interact via a model potential (e.g., Lennard-Jones type), which includes:

* Short-range repulsion
* Long-range attraction

---

### Verlet Algorithm

The position update is:

r(t+Δt) = 2r(t) - r(t-Δt) + a(t)Δt²

✔ Stable
✔ Energy conserving
✔ Suitable for MD simulations

---

## ⚙️ Methodology

1. Initialize particle positions (ordered structure)
2. Assign initial velocities
3. Compute forces between particles
4. Update positions using Verlet algorithm
5. Repeat for multiple time steps
6. Observe structural evolution

---

## 📂 Project Structure

```
molecular-dynamics-simulation/
│── src/              # Scilab code (.sce / .sci)
│── input/            # initial configuration
│── output/           # simulation data
│── plots/            # images (initial & distorted states)
│── README.md
```

---

## 💻 Code Implementation (Scilab)

```scilab
clc;
clear;

N = 100;        // number of particles
dt = 0.01;      // time step

// initialize positions and velocities
// compute forces
// update using Verlet algorithm
```

---

## 📊 Results

The simulation shows **structural evolution of particles**:

### 🔹 Initial Configuration

![Initial State 1](plots/initial1.png)
![Initial State 2](plots/initial2.png)

👉 Particles are arranged in an ordered configuration.

---

### 🔹 After Distortion (Time Evolution)

![Distorted State 1](plots/final1.png)
![Distorted State 2](plots/final2.png)

👉 Due to inter-particle forces, the system evolves into a **distorted/disordered state**.

---

## 📈 Observations

* Initially ordered system becomes unstable over time
* Particle interactions lead to distortion
* Demonstrates fundamental MD behavior
* Sensitive to time step and initial conditions

---

## 🚀 How to Run

1. Open Scilab
2. Load the script:

```
exec('main.sce');
```

3. Run simulation

---

## 📉 Error & Stability

* Smaller time step (dt) → better accuracy
* Large dt → instability
* Energy conservation used for validation

---

## 🔬 Applications

* Condensed matter physics
* Materials science
* Statistical mechanics
* Particle system simulations

---

## 🔧 Future Improvements

* Add energy vs time plots
* Implement Velocity-Verlet
* Extend to 3D systems
* Add temperature control (thermostats)

---

## ✍️ Author

**Ayush Kumar Samal**
