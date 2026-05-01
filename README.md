# Molecular Dynamics Simulation Using Verlet Algorithm (Scilab)

## 📌 Overview

This project implements a **classical Molecular Dynamics (MD) simulation** using **Scilab**. The simulation studies how an initially ordered system of particles evolves over time into a **distorted/disordered configuration** due to inter-particle interactions.

---

## 🧠 Theory

### Newton’s Equation of Motion

The motion of particles is governed by:

F = m a

For an N-particle system:

m d²r/dt² = F

---

### Interatomic Interaction

Particles interact through a model potential (Lennard-Jones type), which includes:

* Short-range repulsion
* Long-range attraction

---

### Verlet Algorithm

The position update equation is:

r(t+Δt) = 2r(t) - r(t-Δt) + a(t)Δt²

✔ Simple and stable
✔ Good energy conservation
✔ Widely used in MD simulations

---

## ⚙️ Methodology

1. Initialize particle positions (ordered structure)
2. Assign initial velocities
3. Compute inter-particle forces
4. Update positions using Verlet algorithm
5. Repeat for multiple time steps
6. Observe system evolution

---

## 📂 Project Structure

```
molecular-dynamics-simulation/
│── src/              # Scilab code
│── input/            # initial conditions
│── output/           # simulation data
│── initial1.png
│── initial2.png
│── final1.png
│── final2.png
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

## 📊 Simulation Results

### 🔹 Initial Configuration

<p align="center">
  <img src="initial1.png" width="45%">
  <img src="initial2.png" width="45%">
</p>

👉 Particles are arranged in an ordered configuration at the beginning.

---

### 🔹 After Distortion (Time Evolution)

<p align="center">
  <img src="final1.png" width="45%">
  <img src="final2.png" width="45%">
</p>

👉 Due to inter-particle interactions, the system evolves into a distorted/disordered state.

---

## 📈 Observations

* Ordered structure becomes unstable over time
* Inter-particle forces lead to distortion
* System evolution depends on initial conditions
* Time step affects stability and accuracy

---

## 🚀 How to Run

1. Open **Scilab**
2. Load your script:

```
exec('main.sce');
```

3. Run the simulation

---

## 📉 Error & Stability

* Smaller time step (dt) → better accuracy
* Larger dt → numerical instability
* Energy conservation used for validation

---

## 🔬 Applications

* Condensed matter physics
* Statistical mechanics
* Materials science
* Particle simulations


---

## ✍️ Author

**Ayush Kumar Samal**
