# KUKA LBR iiwa kinematics (`KUKA.py`)

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
