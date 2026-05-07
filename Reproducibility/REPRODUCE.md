# Reproducibility Instructions

## GitHub Repository

Repository Link:

https://github.com/adityanimbalkar245/Docking-Adapter-Latch-Optimization

---

# Step-by-Step Reproduction Guide

## 1. Clone Repository

git clone https://github.com/adityanimbalkar245/Docking-Adapter-Latch-Optimization.git

---

## 2. Navigate to Code Directory

cd Docking-Adapter-Latch-Optimization/3_Code/src

---

## 3. Install Required Dependencies

pip install -r ../requirements.txt

---

## 4. Run Main Optimization Script

python main.py

---

# Expected Outputs

The script performs:
- RRRP linkage optimization,
- trajectory generation,
- mechanism visualization,
- and docking performance evaluation.

The following outputs are expected:

## Console Outputs
- Optimized linkage parameters
- Docking performance metrics
- Axial tolerance values
- Reachability margins

## Generated Files

Outputs are automatically saved in:

../outputs/

Generated files include:
- trajectory.png
- optimization_results.txt

---

# Notebook Reproduction

A notebook implementation is also provided:

../notebooks/mechanism_optimization.ipynb

The notebook can be executed sequentially to reproduce:
- optimization results,
- trajectory plots,
- and metric calculations.

---

# Expected Optimization Behavior

The optimization should converge to:
- a physically feasible RRRP linkage,
- monotonic inward/downward latch motion,
- and successful seating near the target docking position.

Minor numerical variation may occur depending on:
- operating system,
- SciPy version,
- and floating-point implementation.