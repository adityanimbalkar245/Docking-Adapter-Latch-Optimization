# Design Decisions

## 1. Selection of RRRP Mechanism

An RRRP linkage configuration was selected due to:
- compact packaging capability,
- controllable latch-tip trajectory,
- and compatibility with the available docking envelope.

The mechanism also provided sufficient design flexibility through variation of:
- link lengths,
- slider position,
- and coupler offset coordinates.

---

## 2. Trajectory-Based Optimization

Instead of designing the mechanism using only discrete precision points, the complete latch trajectory over the slider stroke was evaluated during optimization.

This allowed:
- monotonic inward/downward motion enforcement,
- trajectory smoothness evaluation,
- and axial tolerance estimation.

---

## 3. Objective Function Selection

The optimization objective was formulated as a weighted combination of:
- seating accuracy,
- capture capability,
- axial tolerance,
- monotonicity,
- smoothness,
- and motion transmission metrics.

Higher weights were assigned to quantities directly affecting successful docking.

---

## 4. Smooth Penalty Formulation

Most geometric requirements were implemented using smooth quadratic penalty functions instead of discontinuous logical constraints.

This decision was made to:
- improve convergence stability,
- maintain differentiability,
- and avoid abrupt optimization behavior.

---

## 5. Choice of Optimization Algorithm

Sequential Least Squares Programming (SLSQP) was selected because:
- the problem contains nonlinear objectives,
- bounded design variables,
- and nonlinear inequality constraints.

SLSQP also integrates well with smooth penalty formulations.

---

## 6. Planar Modeling Assumption

The mechanism was modeled as a planar system to simplify the synthesis and optimization process.

Dynamic effects, structural flexibility, friction, and contact mechanics were not included in the current formulation.
