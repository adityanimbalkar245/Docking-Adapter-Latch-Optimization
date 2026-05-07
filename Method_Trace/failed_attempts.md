# Failed Attempts and Observations

## 1. Direct Precision-Point Design

Initial attempts focused on manually selecting linkage dimensions to satisfy only:
- open position,
- and closed position requirements.

However, these configurations frequently produced:
- undesirable intermediate trajectories,
- non-monotonic motion,
- and geometric infeasibility during the slider stroke.

---

## 2. Hard Constraint Formulation

Several trajectory conditions were initially implemented using strict logical constraints.

This resulted in:
- unstable optimization behavior,
- discontinuous objective variation,
- and poor convergence performance.

The approach was later replaced with smooth quadratic penalties.

---

## 3. Aggressive Weighting of Single Objectives

Early optimization trials assigned excessively large weights to individual metrics such as capture margin.

This caused the optimizer to produce:
- unrealistic linkage geometries,
- poor seating accuracy,
- or unstable trajectories.

The final weights were selected empirically after multiple optimization iterations.

---

## 4. Invalid Linkage Configurations

Without explicit reachability constraints, the optimizer occasionally generated mechanisms where:
- the linkage could not physically assemble,
- or links separated during the slider stroke.

Additional geometric feasibility constraints were introduced to prevent this issue.

---

## 5. Non-Monotonic Trajectories

Several early solutions achieved acceptable end positions but contained:
- radial reversal,
- upward motion,
- or zig-zag trajectory segments.

Monotonicity penalties were introduced to discourage such behavior.