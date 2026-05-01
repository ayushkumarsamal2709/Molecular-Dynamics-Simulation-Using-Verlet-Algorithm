# Molecular Dynamics Simulation Using Verlet Algorithm (Scilab)

## 📌 Overview

This project presents a **Molecular Dynamics (MD) simulation** using Scilab. It models the evolution of a system of interacting particles using the **Verlet integration method**, starting from an ordered configuration and evolving over time due to inter-particle forces.

---

## 🧠 Theory

### Newton’s Equation of Motion

F = m a

For a system of particles:

m d²r/dt² = F

---

### Lennard-Jones Potential

The interaction between particles is given by:

V(r) = 4ε [ (σ/r)¹² − (σ/r)⁶ ]

* Short-range repulsion
* Long-range attraction

---

### Verlet Algorithm

The position update equation is:

r(t + Δt) = 2r(t) − r(t − Δt) + a(t)Δt²

---

## ⚙️ Methodology

1. Initialize particles in a lattice
2. Assign random velocities
3. Normalize velocities using temperature
4. Compute forces using Lennard-Jones potential
5. Apply periodic boundary conditions
6. Update positions using Verlet algorithm
7. Calculate potential energy over time

---

## 📊 Simulation Results

### 🔹 Initial Configuration

![Initial 1](initial1.png)

![Initial 2](initial2.png)

---

### 🔹 After Distortion

![Final 1](final1.png)

![Final 2](final2.png)

---

## 📂 Project Structure

```
Molecular-Dynamics-Simulation/
│── main.sce
│── initial1.png
│── initial2.png
│── final1.png
│── final2.png
│── README.md
```

---

## 💻 Scilab Code (Main Parts)

### Lattice Generator

```scilab
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

### Initialization

```scilab
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

### Force Calculation

```scilab
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

## 🚀 How to Run

Open Scilab and run:

```
exec('main.sce');
```

---

## ✍️ Author

Ayush Kumar Samal
