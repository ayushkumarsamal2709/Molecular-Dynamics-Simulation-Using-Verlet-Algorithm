# Molecular Dynamics Simulation Using Verlet Algorithm (Scilab)

## 📌 Overview

This project demonstrates a **Molecular Dynamics (MD) simulation** using Scilab. An initially ordered system of particles evolves into a **distorted/disordered configuration** due to inter-particle interactions.

---

## 🧠 Theory

### Newton’s Equation of Motion

F = m a

m d²r/dt² = F

---

### Verlet Algorithm

r(t+Δt) = 2r(t) - r(t-Δt) + a(t)Δt²

---

## ⚙️ Methodology

1. Initialize particle positions
2. Assign velocities
3. Compute forces
4. Update positions using Verlet algorithm
5. Repeat over time

---

## 📊 Simulation Results

### 🔹 Initial Configuration

![Initial 1](./initial1.png)

[Click if not visible](./initial1.png)

![Initial 2](./initial2.png)

[Click if not visible](./initial2.png)

---

### 🔹 After Distortion

![Final 1](./final1.png)

[Click if not visible](./final1.png)

![Final 2](./final2.png)

[Click if not visible](./final2.png)

---

## 📈 Observations

* Ordered structure becomes distorted over time
* Particle interactions drive system evolution
* Results depend on time step and initial conditions

---

## 🚀 How to Run

1. Open Scilab
2. Run:

```
exec('main.sce');
```

---

## ✍️ Author

Ayush Kumar Samal
