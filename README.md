# Molecular Dynamics Simulation Using Verlet Algorithm (Scilab)

## 📌 Overview

This project presents a **classical Molecular Dynamics (MD) simulation** implemented in Scilab. The simulation models a system of interacting particles confined in a three-dimensional box and evolves their motion over time using the **Verlet integration algorithm**.

The system begins with an ordered lattice configuration and evolves dynamically under inter-particle forces, demonstrating the transition toward a more complex state governed by physical interactions.

---

## 🧠 Theoretical Background

### Newton’s Equation of Motion

The motion of each particle is governed by Newton’s second law:

F = m a

For an N-particle system:

m d²r/dt² = F

---

### Lennard-Jones Potential

The interaction between particles is modeled using the Lennard-Jones potential:

V(r) = 4ε [ (σ/r)¹² − (σ/r)⁶ ]

This potential represents:

* Strong **repulsion** at short distances
* Weak **attraction** at larger distances

---

### Periodic Boundary Conditions

To simulate an infinite system using a finite box, **periodic boundary conditions** are applied. When particles cross the simulation box boundary, they re-enter from the opposite side.

---

### Verlet Integration Algorithm

The time evolution of particle positions is computed using:

r(t + Δt) = 2r(t) − r(t − Δt) + a(t)Δt²

This method is:

* Numerically stable
* Efficient for large systems
* Widely used in MD simulations

---

## ⚙️ Simulation Methodology

1. Generate initial particle positions using a cubic lattice
2. Assign random initial velocities
3. Scale velocities to match desired temperature
4. Compute inter-particle forces using Lennard-Jones potential
5. Apply periodic boundary conditions
6. Integrate equations of motion using Verlet algorithm
7. Calculate potential energy over time
8. Perform averaging after equilibration

---

## 📂 Project Structure

```id="s1"
Molecular-Dynamics-Simulation/
│── main.sce
│── initial1.png
│── initial2.png
│── final1.png
│── final2.png
│── README.md
```

---

## 💻 Scilab Implementation

### 🔹 Lattice Generator

```scilab id="s2"
function [a1,a2,a3] = lattice_pos(index,Lx,Ly,Lz,a)
n = ceil(index^(1/3));
count = 1;

for i=0:n-1
    for j=0:n-1
        for k=0:n-1
            if count == index then
                a1 = i*a;
                a2 = j*a;
                a3 = k*a;
                return
            end
            count = count + 1;
        end
    end
end
endfunction
```

---

### 🔹 Initialization

```scilab id="s3"
function [x,y,z,vx,vy,vz,xm,ym,zm] = init(npart,Lx,Ly,Lz,a,temp,dt)

sumv2 = 0;

for i=1:npart
    [x(i),y(i),z(i)] = lattice_pos(i,Lx,Ly,Lz,a);

    vx(i)=rand()-0.5;
    vy(i)=rand()-0.5;
    vz(i)=rand()-0.5;

    sumv2 = sumv2 + vx(i)^2 + vy(i)^2 + vz(i)^2;
end

fs = sqrt((3*temp*npart)/sumv2);

for i=1:npart
    vx(i)=vx(i)*fs;
    vy(i)=vy(i)*fs;
    vz(i)=vz(i)*fs;

    xm(i)=x(i)-vx(i)*dt;
    ym(i)=y(i)-vy(i)*dt;
    zm(i)=z(i)-vz(i)*dt;
end

endfunction
```

---

### 🔹 Force Calculation (Lennard-Jones)

```scilab id="s4"
function [fx,fy,fz,en] = force(npart,x,y,z,L,rc)

rc2 = rc*rc;
en = 0;

fx = zeros(1,npart);
fy = zeros(1,npart);
fz = zeros(1,npart);

for i=1:npart-1
    for j=i+1:npart

        xr = x(i)-x(j);
        yr = y(i)-y(j);
        zr = z(i)-z(j);

        xr = xr - L*round(xr/L);
        yr = yr - L*round(yr/L);
        zr = zr - L*round(zr/L);

        r2 = xr^2 + yr^2 + zr^2;

        if r2 < rc2 then
            r2i = 1/r2;
            r6i = r2i^3;

            ff = 48*r2i*r6i*(r6i - 0.5);

            fx(i)=fx(i)+ff*xr;
            fx(j)=fx(j)-ff*xr;

            fy(i)=fy(i)+ff*yr;
            fy(j)=fy(j)-ff*yr;

            fz(i)=fz(i)+ff*zr;
            fz(j)=fz(j)-ff*zr;

            en = en + 4*r6i*(r6i - 1);
        end
    end
end

endfunction
```

---

### 🔹 Verlet Integration

```scilab id="s5"
function [x,y,z,xm,ym,zm] = integ(fx,fy,fz,npart,x,y,z,xm,ym,zm,dt)

for i=1:npart
    xx = 2*x(i) - xm(i) + dt^2*fx(i);
    yy = 2*y(i) - ym(i) + dt^2*fy(i);
    zz = 2*z(i) - zm(i) + dt^2*fz(i);

    xm(i)=x(i);
    ym(i)=y(i);
    zm(i)=z(i);

    x(i)=xx;
    y(i)=yy;
    z(i)=zz;
end

endfunction
```

---

### 🔹 Main Simulation

```scilab id="s6"
clear;
clc;

Nlist = [100 200 300 400 500 600 700];

L = 15;
a = 1.1;
temp = 1.0;
dt = 0.001;
rc = 2.5;

steps = 20000;
equil = 5000;

Volume = L^3;

avgE = zeros(1,length(Nlist));
errE = zeros(1,length(Nlist));
rho = zeros(1,length(Nlist));

scf(1);
clf();
legends = [];

for nindex = 1:length(Nlist)

    npart = Nlist(nindex);
    rho(nindex) = npart/Volume;

    [x,y,z,vx,vy,vz,xm,ym,zm] = init(npart,L,L,L,a,temp,dt);

    PE = zeros(1,steps);

    for t=1:steps
        [fx,fy,fz,en] = force(npart,x,y,z,L,rc);
        [x,y,z,xm,ym,zm] = integ(fx,fy,fz,npart,x,y,z,xm,ym,zm,dt);

        PE(t) = en/npart;
    end

    plot(1:steps,PE);
    legends(nindex) = string(npart);

    data = PE(equil:steps);
    avgE(nindex) = mean(data);
    errE(nindex) = stdev(data);

end

xlabel("Time");
ylabel("Potential Energy");
title("Potential Energy vs Time");
legend(legends);

scf(2);
clf();

plot(rho,avgE,'-o');

xlabel("Number Density");
ylabel("Average Energy");
title("Number Density vs Average Energy");

for i=1:length(rho)
    xpoly([rho(i) rho(i)], [avgE(i)-errE(i) avgE(i)+errE(i)]);
end
```

---

## 📊 Results and Observations

* The system evolves dynamically under inter-particle forces
* Potential energy varies with time and stabilizes after equilibration
* Average energy depends on particle density
* Error bars indicate statistical fluctuations in energy

---

## 🚀 How to Run

1. Open Scilab
2. Load the script
3. Execute:

```id="run1"
exec('main.sce');
```

---

## ✍️ Author

**Ayush Kumar Samal**
