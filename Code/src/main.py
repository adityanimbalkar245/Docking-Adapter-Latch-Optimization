"""
Optimization of an RRRP-Based Latching Mechanism
for CubeSat Docking Adapter

Main executable script.
"""

import os
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# ==========================================
# OUTPUT DIRECTORY
# ==========================================
OUTPUT_DIR = "../Data_Results"
OUTPUT_DIR1 = "../Data_Results/figures"
OUTPUT_DIR2 = "../Data_Results/outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR1, exist_ok=True)
os.makedirs(OUTPUT_DIR2, exist_ok=True)

# ==========================================
# 1. FIXED SYSTEM PARAMETERS
# ==========================================
Hboss_t = 30.00
Hboss_c = 20.00
H_dock_c = 30.00
H_dock_t = H_dock_c + Hboss_c + Hboss_t
R_hex = 27.28

# Latch Seat
H_seat = H_dock_c + Hboss_c + 5.00
R_seat = 20.00

# Constraints & Strokes
a3_min = 10.00
a3_max = 25.00
Rmax = 50.00
X_MAX = 60.00

# ==========================================
# 2. EXACT KINEMATICS
# ==========================================
def rrrp_kinematics(a1, a2, R, s, u, v, a3):
    O = np.array([R, 0.0])
    B = np.array([s, a3])

    d = float(np.linalg.norm(B - O))

    # Feasibility Check
    if d > (a1 + a2) or d < abs(a1 - a2) or d < 1e-6:
        return None, None, None

    # Circle Intersection
    l = (a1**2 - a2**2 + d**2) / (2 * d)
    h_sq = a1**2 - l**2
    h = np.sqrt(max(0.0, h_sq))

    dir_u = (B - O) / d
    dir_u_x, dir_u_y = dir_u

    dir_v = np.array([-dir_u_y, dir_u_x])

    A1 = O + l * dir_u + h * dir_v
    A2 = O + l * dir_u - h * dir_v

    A1_x, A1_y = A1
    A2_x, A2_y = A2

    # Upward elbow selection
    A = A1 if A1_y > A2_y else A2

    # Coupler frame
    AB = B - A
    norm_AB = float(np.linalg.norm(AB))

    t_vec = AB / norm_AB

    t_vec_x, t_vec_y = t_vec

    n_vec = np.array([-t_vec_y, t_vec_x])

    # Tip Position
    tip = A + u * t_vec + v * n_vec

    return tip, A, B

# ==========================================
# 3. TRAJECTORY GENERATION
# ==========================================
def compute_trajectory_full(params, N=25):

    a1, a2, R, s, u, v = map(float, params)

    a3_vals = np.linspace(a3_min, a3_max, N)

    traj = []
    A_list = []
    B_list = []

    for a3 in a3_vals:

        tip, A, B = rrrp_kinematics(a1, a2, R, s, u, v, a3)

        if tip is None:
            return None, None, None, None

        traj.append(tip)
        A_list.append(A)
        B_list.append(B)

    return (
        np.array(traj),
        np.array(A_list),
        np.array(B_list),
        a3_vals
    )

# ==========================================
# 4. COST FUNCTION
# ==========================================
def cost_function(params):

    a1, a2, R, s, u, v = map(float, params)

    traj, A_list, B_list, a3_vals = compute_trajectory_full(params)

    if traj is None:

        O = np.array([R, 0.0])
        B_mid = np.array([s, (a3_min + a3_max)/2])

        d_mid = float(np.linalg.norm(B_mid - O))

        fail_cost = 1e6 + (d_mid - (a1 + a2))**2 * 1000

        return float(fail_cost)

    Z = 0
    L = -1

    x = traj[:, Z]
    y = traj[:, 1]

    penalty = 0.0

    # --------------------------------------
    # 1. Closed Position
    # --------------------------------------
    J_closed = (
        (x[Z] - R_seat)**2 +
        (y[Z] - H_seat)**2
    )

    # --------------------------------------
    # 2. Capture
    # --------------------------------------
    x_open = x[L]
    y_open = y[L]

    J_capture = -(x_open - R_hex)

    # --------------------------------------
    # 3. Axial Tolerance
    # --------------------------------------
    axial_tol = -1e3

    for i in range(len(x) - 1):

        crossed = (
            (x[i] >= R_hex >= x[i+1]) or
            (x[i] <= R_hex <= x[i+1])
        )

        if crossed:

            t = (
                (R_hex - x[i]) /
                (x[i+1] - x[i] + 1e-9)
            )

            y_interp = (
                y[i] +
                t * (y[i+1] - y[i])
            )

            axial_tol = y_interp - H_seat

            break

    J_axial = -axial_tol

    # --------------------------------------
    # 4. Monotonicity
    # --------------------------------------
    dx = np.diff(x)
    dy = np.diff(y)

    J_mono = (
        np.sum(np.maximum(0, dx)**2) +
        np.sum(np.maximum(0, dy)**2)
    )

    # --------------------------------------
    # 5. Smoothness
    # --------------------------------------
    J_smooth = np.sum(
        np.diff(traj, n=2, axis=0)**2
    )

    # --------------------------------------
    # 6. Force Transmission
    # --------------------------------------
    dx_g = np.gradient(x)
    dy_g = np.gradient(y)

    vel = np.sqrt(dx_g**2 + dy_g**2) + 1e-6

    alignment = -dy_g / vel

    J_force = -np.mean(alignment / vel)

    # --------------------------------------
    # Smooth Penalties
    # --------------------------------------
    penalty += np.sum(np.maximum(0, R_seat - x)**2)

    penalty += np.maximum(0, R_hex - x_open)**2

    penalty += np.maximum(0, H_seat - y_open)**2

    penalty += np.maximum(0, y_open - H_dock_t)**2

    penalty += np.maximum(0, a1 - (Rmax - R))**2

    penalty += np.maximum(0, 5 - abs(R - s))**2

    penalty += np.sum(
        np.maximum(0, -(A_list[:,Z] - R))**2
    )

    penalty += np.sum(
        np.maximum(0, -A_list[:,1])**2
    )

    penalty += np.sum(
        np.maximum(0, B_list[:,Z] - A_list[:,Z])**2
    )

    penalty += np.sum(
        np.maximum(0, -(B_list[:,1] - A_list[:,1]))**2
    )

    # --------------------------------------
    # Final Objective
    # --------------------------------------
    J_total = (
        20 * J_closed +
        22 * J_capture +
        15 * J_axial +
        15 * J_mono +
        2 * J_smooth +
        8 * J_force +
        50 * penalty
    )

    return float(np.sum(J_total))

# ==========================================
# 5. PLOTTING FUNCTION
# ==========================================
def plot_docking_system(params):

    result = compute_trajectory_full(params, N=50)

    fig, ax = plt.subplots(figsize=(10, 8))

    # Geometry
    ax.axhline(
        H_dock_c,
        color='blue',
        linestyle='--',
        alpha=0.5,
        label='H_dock_c'
    )

    ax.axhline(
        H_seat,
        color='orange',
        linestyle='--',
        alpha=0.7,
        label='H_seat'
    )

    ax.axhline(
        H_dock_t,
        color='red',
        linestyle='--',
        alpha=0.5,
        label='H_dock_t'
    )

    ax.axvline(
        R_hex,
        color='purple',
        linestyle=':',
        linewidth=2,
        label='R_hex'
    )

    ax.axvline(
        R_seat,
        color='green',
        linestyle=':',
        linewidth=2,
        label='R_seat'
    )

    ax.plot(
        R_seat,
        H_seat,
        'ko',
        markersize=8,
        label='Seat Target'
    )

    R_val = params[2]

    ax.plot(
        R_val,
        0,
        'ks',
        markersize=8,
        label='Origin'
    )

    Z = 0
    L = -1

    if result[Z] is not None:

        traj, A_list, B_list, _ = result

        O_x, O_y = R_val, 0.0

        # Trajectory
        traj_x = traj[:, Z]
        traj_y = traj[:, 1]

        ax.plot(
            traj_x,
            traj_y,
            'k-',
            linewidth=2,
            label='Tip Trajectory'
        )

        # Closed Position
        A_c_x, A_c_y = A_list[Z]
        B_c_x, B_c_y = B_list[Z]
        tip_c_x, tip_c_y = traj[Z]

        ax.plot(
            [O_x, A_c_x],
            [O_y, A_c_y],
            'b-',
            linewidth=3
        )

        ax.plot(
            [A_c_x, B_c_x],
            [A_c_y, B_c_y],
            'g-',
            linewidth=3
        )

        ax.plot(
            [A_c_x, tip_c_x],
            [A_c_y, tip_c_y],
            'r-',
            linewidth=2
        )

        ax.plot(
            tip_c_x,
            tip_c_y,
            'ro',
            markersize=6
        )

        ax.text(
            tip_c_x + 2,
            tip_c_y,
            'Closed',
            fontweight='bold'
        )

        # Open Position
        A_o_x, A_o_y = A_list[L]
        B_o_x, B_o_y = B_list[L]
        tip_o_x, tip_o_y = traj[L]

        ax.plot(
            [O_x, A_o_x],
            [O_y, A_o_y],
            'b--',
            linewidth=2,
            alpha=0.6
        )

        ax.plot(
            [A_o_x, B_o_x],
            [A_o_y, B_o_y],
            'g--',
            linewidth=2,
            alpha=0.6
        )

        ax.plot(
            [A_o_x, tip_o_x],
            [A_o_y, tip_o_y],
            'r--',
            linewidth=2,
            alpha=0.6
        )

        ax.plot(
            tip_o_x,
            tip_o_y,
            'r^',
            markersize=6
        )

        ax.text(
            tip_o_x + 2,
            tip_o_y,
            'Open',
            fontweight='bold'
        )

    ax.set_aspect('equal')

    ax.set_xlim(-10, 70)
    ax.set_ylim(-10, 90)

    ax.set_xlabel('Radial distance x (mm)')
    ax.set_ylabel('Docking axis y (mm)')

    ax.set_title(
        'Docking Geometry & Optimized RRRP Mechanism'
    )

    ax.grid(True, linestyle=':', alpha=0.6)

    ax.legend(
        loc='upper right',
        bbox_to_anchor=(1.35, 1)
    )

    plt.tight_layout()

    # Save Figure
    fig.savefig(
        os.path.join(
            OUTPUT_DIR1,
            "trajectory.png"
        ),
        dpi=300,
        bbox_inches='tight'
    )

    plt.show()

# ==========================================
# 6. METRICS FUNCTION
# ==========================================
def print_comprehensive_metrics(params):

    idx_0 = 0
    idx_1 = 1
    last = -1

    a1, a2, R, s, u, v = params

    traj, A_list, B_list, a3_vals = compute_trajectory_full(
        params,
        N=50
    )

    if traj is None:
        print("Trajectory is invalid.")
        return

    results_text = []

    def log(line):
        print(line)
        results_text.append(line)

    # --------------------------------------
    # FIXED GEOMETRY
    # --------------------------------------
    log("=== FIXED SYSTEM GEOMETRY ===")

    log(f"R_seat (Target X)    : {R_seat:.2f} mm")
    log(f"H_seat (Target Y)    : {H_seat:.2f} mm")
    log(f"R_hex (Entry X limit): {R_hex:.2f} mm")
    log(f"H_dock_top (Max Y)   : {H_dock_t:.2f} mm")
    log(f"Slider Stroke (a3)   : {a3_min:.2f} to {a3_max:.2f} mm\n")

    # --------------------------------------
    # PARAMETERS
    # --------------------------------------
    log("=== OPTIMIZED MECHANISM PARAMETERS ===")

    log(f"Link 1 (a1)          : {a1:.3f} mm")
    log(f"Link 2 (a2)          : {a2:.3f} mm")
    log(f"Origin (R)           : {R:.3f} mm")
    log(f"Slider Pos (s)       : {s:.3f} mm")
    log(f"Coupler offset (u,v) : ({u:.3f}, {v:.3f}) mm\n")

    # --------------------------------------
    # CLOSED POSITION
    # --------------------------------------
    A_c = A_list[idx_0]
    B_c = B_list[idx_0]
    T_c = traj[idx_0]

    log("=== CLOSED POSITION COORDINATES ===")

    log(f"Origin (O) : ({R:.2f}, 0.00)")
    log(f"Slider (B) : ({B_c[idx_0]:.2f}, {B_c[idx_1]:.2f})")
    log(f"Elbow (A)  : ({A_c[idx_0]:.2f}, {A_c[idx_1]:.2f})")
    log(f"Tip (T)    : ({T_c[idx_0]:.2f}, {T_c[idx_1]:.2f})\n")

    # --------------------------------------
    # OPEN POSITION
    # --------------------------------------
    A_o = A_list[last]
    B_o = B_list[last]
    T_o = traj[last]

    log("=== OPEN POSITION COORDINATES ===")

    log(f"Origin (O) : ({R:.2f}, 0.00)")
    log(f"Slider (B) : ({B_o[idx_0]:.2f}, {B_o[idx_1]:.2f})")
    log(f"Elbow (A)  : ({A_o[idx_0]:.2f}, {A_o[idx_1]:.2f})")
    log(f"Tip (T)    : ({T_o[idx_0]:.2f}, {T_o[idx_1]:.2f})\n")

    # --------------------------------------
    # PERFORMANCE
    # --------------------------------------
    log("=== PERFORMANCE METRICS ===")

    err_x = abs(T_c[idx_0] - R_seat)
    err_y = abs(T_c[idx_1] - H_seat)

    total_err = np.sqrt(err_x**2 + err_y**2)

    log(f"Seat Target Error     : {total_err:.3f} mm")

    capture_margin = T_o[idx_0] - R_hex

    log(f"Lateral Capture Margin: {capture_margin:.3f} mm")

    log(f"Open Tip Height       : {T_o[idx_1]:.3f} mm")

    axial_tol = None

    traj_x = traj[:, idx_0]
    traj_y = traj[:, idx_1]

    for i in range(len(traj_x) - 1):

        crossed = (
            (traj_x[i] >= R_hex >= traj_x[i+1]) or
            (traj_x[i] <= R_hex <= traj_x[i+1])
        )

        if crossed:

            t = (
                (R_hex - traj_x[i]) /
                (traj_x[i+1] - traj_x[i] + 1e-9)
            )

            y_interp = (
                traj_y[i] +
                t * (traj_y[i+1] - traj_y[i])
            )

            axial_tol = y_interp - H_seat

            break

    if axial_tol is not None:
        log(f"Axial Tolerance       : {axial_tol:.3f} mm")

    max_reach_req = np.sqrt(
        (R - s)**2 + a3_max**2
    )

    log(f"Max Reach Required    : {max_reach_req:.3f} mm")
    log(f"Link Reach Available  : {(a1 + a2):.3f} mm")

    reach_margin = (
        (a1 + a2) - max_reach_req
    )

    log(f"Reach Safety Margin   : {reach_margin:.3f} mm")

    # Save Results
    with open(
        os.path.join(
            OUTPUT_DIR2,
            "optimization_results.txt"
        ),
        "w"
    ) as f:

        for line in results_text:
            f.write(line + "\n")

# ==========================================
# 7. OPTIMIZATION EXECUTION
# ==========================================
if __name__ == "__main__":

    p0 = [20.0, 25.0, 35.0, 20.0, 15.0, -5.0]

    bounds = [
        (10.0, 30.0),
        (10.0, 40.0),
        (20.0, Rmax),
        (20.0, 50.0),
        (-80.0, 80.0),
        (-80.0, 80.0)
    ]

    def reach_constraint(params):

        a1, a2, R, s, u, v = params

        max_dist = np.sqrt(
            (R - s)**2 + a3_max**2
        )

        min_dist = np.sqrt(
            (R - s)**2 + a3_min**2
        )

        return np.array([
            (a1 + a2) - max_dist - 0.5,
            min_dist - abs(a1 - a2) - 0.5
        ])

    cons = {
        'type': 'ineq',
        'fun': reach_constraint
    }

    print("Running Optimization...\n")

    result = minimize(
        cost_function,
        p0,
        method='SLSQP',
        bounds=bounds,
        constraints=cons,
        options={
            'maxiter': 500,
            'disp': True
        }
    )

    print("\nFinal Parameters:\n")

    labels = ['a1', 'a2', 'R', 's', 'u', 'v']

    for name, val in zip(labels, result.x):
        print(f"{name}: {val:.3f}")

    print("\nGenerating Trajectory Plot...\n")

    plot_docking_system(result.x)

    print("\nGenerating Comprehensive Metrics...\n")

    print_comprehensive_metrics(result.x)

    print("\nOutputs saved in:")
    print(os.path.abspath(OUTPUT_DIR))