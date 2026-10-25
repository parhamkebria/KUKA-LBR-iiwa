# KUKA LBR iiwa kinematics

Forward and inverse kinematics for the 7-DoF KUKA LBR iiwa, as an importable
Python module (NumPy only). 

The arm is redundant, so a pose has infinitely many IK solutions. They are
selected with:

- **`param`** – the arm angle ψ (rad): rotation of the elbow around the
  shoulder–wrist axis, measured from a reference plane (the arm with joint 3 = 0).
- **`conf`** – configuration number: `bit0` = 1 if joint 2 < 0 (arm),
  `bit1` = 1 if joint 4 < 0 (elbow), `bit2` = 1 if joint 6 < 0 (wrist).

Units: metres and radians. Poses are 4×4 homogeneous matrices (flange in base frame).
Link lengths: 0.34, 0.4, 0.4, 0.126 m.

## Usage

```python
import numpy as np
from KUKA import forward_kinematics, inverse_kinematics, within_joint_limits

q = np.deg2rad([10, 30, -20, 60, 15, 40, 5])

pose, param, conf, _ = forward_kinematics(q)
joints, _, _ = inverse_kinematics(pose, param, conf)   # recovers q

print(within_joint_limits(joints))
```

Sweeping the elbow (arm angle) for a fixed end-effector pose:

```python
for psi in np.linspace(-1, 1, 5):
    q_psi, _, _ = inverse_kinematics(pose, psi, conf)
```

## API

| Function | Description |
|---|---|
| `forward_kinematics(joints)` | → `(pose, param, conf, jout)`; `joints` has 7 values. |
| `inverse_kinematics(pose, param=0.0, conf=0)` | → `(joints, s_mat, w_mat)`; `s_mat`/`w_mat` are the shoulder/wrist matrices `[A, B, C]` (3×3×3) used for joint-limit handling. |
| `configuration(rconf)` | Decode to `(arm, elbow, wrist)` ∈ {±1}. |
| `dh_matrix(a, alpha, d, theta)` | Single DH transform. |
| `within_joint_limits(joints)` | Checks ±170, 120, 170, 120, 170, 120, 175°. |

`inverse_kinematics` and `forward_kinematics` use an internal reference-plane
computation (virtual arm with joint 3 = 0) that defines the zero of `param`.
It is also where the reachability and elbow-singularity checks happen, so
those errors surface from both functions.

`inverse_kinematics` raises `KukaError` (a `ValueError`) when the wrist is
outside the reachable workspace or at the elbow singularity. Joint limits are
**not** enforced by IK; check with `within_joint_limits`.

## Self-test

```
python KUKA.py
```

Runs random FK → IK → FK round trips (max pose error ≈ 1e-12).

# KUKA LBR Med & LBR iiwa DH Parameters

This repository contains the standard Denavit-Hartenberg (DH) parameters for the **KUKA LBR Med** and **KUKA LBR iiwa** 7-DOF collaborative robot arms. Because the LBR Med is built directly on the LBR iiwa chassis, models sharing the same payload/reach designation utilize identical geometric forward and inverse kinematics models.

## Kinematic Structure
Both variants feature a 7-axis anthropomorphic design with alternating roll-pitch-roll configurations, introducing a kinematically redundant null space. In the zero position ($q_i = 0$), the arm stands fully upright along the vertical Z-axis.

### Parameter Convention
These parameters adhere to the **Standard DH Convention**. The coordinate frame transformation matrix from frame $i-1$ to frame $i$ is defined as:

$$T_i^{i-1} = Rot_z(\theta_i) \cdot Trans_z(d_i) \cdot Trans_x(a_i) \cdot Rot_x(\alpha_i)$$

---

## 1. KUKA LBR 7 R800 (Med / iiwa)
* **Maximum Reach:** 800 mm
* **Payload Capacity:** 7 kg

| Joint ($i$) | $\alpha_{i-1}$ (rad) | $a_{i-1}$ (mm) | $d_i$ (mm) | $\theta_i$ (rad) | Joint Range |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **1** | $0$ | $0$ | $310$ | $q_1$ | $\pm 170^\circ$ |
| **2** | $-\pi/2$ | $0$ | $0$ | $q_2$ | $\pm 120^\circ$ |
| **3** | $\pi/2$ | $0$ | $400$ | $q_3$ | $\pm 170^\circ$ |
| **4** | $-\pi/2$ | $0$ | $0$ | $q_4$ | $\pm 120^\circ$ |
| **5** | $\pi/2$ | $0$ | $400$ | $q_5$ | $\pm 170^\circ$ |
| **6** | $-\pi/2$ | $0$ | $0$ | $q_6$ | $\pm 120^\circ$ |
| **7** | $\pi/2$ | $0$ | $126$ | $q_7$ | $\pm 175^\circ$ |

---

## 2. KUKA LBR 14 R820 (Med / iiwa)
* **Maximum Reach:** 820 mm
* **Payload Capacity:** 14 kg

| Joint ($i$) | $\alpha_{i-1}$ (rad) | $a_{i-1}$ (mm) | $d_i$ (mm) | $\theta_i$ (rad) | Joint Range |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **1** | $0$ | $0$ | $360$ | $q_1$ | $\pm 170^\circ$ |
| **2** | $-\pi/2$ | $0$ | $0$ | $q_2$ | $\pm 120^\circ$ |
| **3** | $\pi/2$ | $0$ | $420$ | $q_3$ | $\pm 170^\circ$ |
| **4** | $-\pi/2$ | $0$ | $0$ | $q_4$ | $\pm 120^\circ$ |
| **5** | $\pi/2$ | $0$ | $400$ | $q_5$ | $\pm 170^\circ$ |
| **6** | $-\pi/2$ | $0$ | $0$ | $q_6$ | $\pm 120^\circ$ |
| **7** | $\pi/2$ | $0$ | $126$ | $q_7$ | $\pm 175^\circ$ |

---

## Technical Notes

* **Flange Face Transformation:** The $d_7$ parameter ($126\text{ mm}$) transforms the mathematical model directly to the native mechanical tool flange face. When integrating custom end-effectors, add your tool's center-of-mass/tool-center-point (TCP) parameters directly onto $d_7$ or append it as an explicit rigid-body extension matrix.
* **Modified DH Conversion:** If your inverse kinematics algorithm or environment (such as certain ROS configurations or legacy solvers) requires **Modified DH (Craig's Convention)**, shift the row index alignments down by one step. The absolute physical geometry metrics remain invariant.
