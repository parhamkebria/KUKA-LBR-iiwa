# KUKA LBR iiwa / LBR Med kinematics

Forward and inverse kinematics for the 7-DoF KUKA LBR iiwa and LBR Med (7 R800
and 14 R820), as an importable Python module ([`python/KUKA.py`](python/KUKA.py), NumPy only) and an
equivalent dependency-free C++17 library ([`cpp/`](cpp/)).

The arm is redundant, so a pose has infinitely many IK solutions. They are
selected with:

- **`param`** – the arm angle ψ (rad): rotation of the elbow around the
  shoulder–wrist axis, measured from a reference plane (the arm with joint 3 = 0).
- **`conf`** – configuration number: `bit0` = 1 if joint 2 < 0 (arm),
  `bit1` = 1 if joint 4 < 0 (elbow), `bit2` = 1 if joint 6 < 0 (wrist).

Units: metres and radians. Poses are 4×4 homogeneous matrices (flange in base frame).
The robot model is chosen with the `payload` argument (kg), default `7`:

| `payload` | Model | d1, d3, d5, d7 (mm) |
|:---:|---|---|
| `7` | LBR 7 R800 | 310, 400, 400, 126 |
| `14` | LBR 14 R820 | 360, 420, 400, 126 |

Any other value raises `ValueError`. Joint angles follow the DH tables in the [DH Parameters](#dh-parameters) section below.

## Python usage

Run from `python/` (or add it to `PYTHONPATH`).

```python
import numpy as np
from KUKA import forward_kinematics, inverse_kinematics, within_joint_limits

q = np.deg2rad([10, 30, -20, 60, 15, 40, 5])

payload = 14   # 7 -> LBR 7 R800, 14 -> LBR 14 R820

pose, param, conf, _ = forward_kinematics(q, payload)
joints, _, _ = inverse_kinematics(pose, param, conf, payload)   # recovers q

print(within_joint_limits(joints, payload))
```

Sweeping the elbow (arm angle) for a fixed end-effector pose:

```python
for psi in np.linspace(-1, 1, 5):
    q_psi, _, _ = inverse_kinematics(pose, psi, conf, payload)
```

## API

| Function | Description |
|---|---|
| `forward_kinematics(joints, payload=7)` | → `(pose, param, conf, jout)`; `joints` has 7 values. |
| `inverse_kinematics(pose, param=0.0, conf=0, payload=7)` | → `(joints, s_mat, w_mat)`; `s_mat`/`w_mat` are the shoulder/wrist matrices `[A, B, C]` (3×3×3) used for joint-limit handling. |
| `configuration(conf)` | Decode to `(arm, elbow, wrist)` ∈ {±1}. |
| `dh_table(payload=7)` | DH table `[alpha_{i-1}, a_{i-1}, d_i]` of the selected model. |
| `dh_matrix(alpha, a, d, theta)` | Single DH transform (README convention). |
| `within_joint_limits(joints, payload=7)` | Checks ±170, 120, 170, 120, 170, 120, 175°. |

`inverse_kinematics` and `forward_kinematics` use an internal reference-plane
computation (virtual arm with joint 3 = 0) that defines the zero of `param`.
It is also where the reachability and elbow-singularity checks happen, so
those errors surface from both functions.

Both functions raise `KukaError` (a `ValueError`) when the wrist is outside the
reachable workspace or the elbow is singular (fully stretched or folded). Joint limits are
**not** enforced by IK; check with `within_joint_limits`.

## Self-test

```
python python/KUKA.py
```

Runs random FK → IK → FK round trips for both models (max pose error ≈ 1e-12).

## C++ usage

[`cpp/KUKA.h`](cpp/KUKA.h) / [`cpp/KUKA.cpp`](cpp/KUKA.cpp) provide the same functions in the
`kuka` namespace (C++17, standard library only). Names, argument order, defaults, units, conventions
and `payload` handling are identical to the Python module; matrices are `std::array`s
(`Mat4`, `Mat3`, `Joints`) and results are structs.

```cpp
#include "KUKA.h"

kuka::Joints q = {0.17, 0.52, -0.35, 1.05, 0.26, 0.70, 0.09};  // rad
int payload = 14;                                              // 7 or 14

kuka::FKResult fk = kuka::forward_kinematics(q, payload);      // .pose .param .conf .jout
kuka::IKResult ik = kuka::inverse_kinematics(fk.pose, fk.param, fk.conf, payload);  // .joints .s_mat .w_mat
bool ok = kuka::within_joint_limits(ik.joints, payload);
```

| Python | C++ |
|---|---|
| `forward_kinematics(joints, payload=7)` | `kuka::forward_kinematics(const Joints&, int payload = 7)` → `FKResult` |
| `inverse_kinematics(pose, param=0.0, conf=0, payload=7)` | `kuka::inverse_kinematics(const Mat4&, double param = 0, int conf = 0, int payload = 7)` → `IKResult` |
| `configuration`, `dh_table`, `dh_matrix`, `within_joint_limits` | same names in `kuka::` |
| `KukaError` (`ValueError`) | `kuka::KukaError` (`std::runtime_error`) for unreachable/singular poses; `std::invalid_argument` for a bad `payload` |

Build and run the self-test (also cross-checks fixed values against the Python module):

```
cd cpp
cmake -S . -B build && cmake --build build && ./build/test_kuka
# or: c++ -std=c++17 -O2 KUKA.cpp test_kuka.cpp -o test_kuka && ./test_kuka
```

## DH Parameters

This section lists the modified (Craig) Denavit-Hartenberg parameters for the **KUKA LBR Med** and **KUKA LBR iiwa** 7-DOF collaborative robot arms. Because the LBR Med is built directly on the LBR iiwa chassis, models sharing the same payload/reach designation utilize identical geometric forward and inverse kinematics models.

### Kinematic Structure
Both variants feature a 7-axis anthropomorphic design with alternating roll-pitch-roll configurations, introducing a kinematically redundant null space. In the zero position ($q_i = 0$), the arm stands fully upright along the vertical Z-axis.

### Parameter Convention
These parameters adhere to the **Modified (Craig) DH Convention**, which is what `dh_matrix` implements. The transformation from frame $i-1$ to frame $i$ is:

$$T_{i-1}^{i} = Rot_x(\alpha_{i-1}) \cdot Trans_x(a_{i-1}) \cdot Rot_z(\theta_i) \cdot Trans_z(d_i)$$

The analytic solver internally uses an equivalent standard-DH chain; the two agree on the full flange pose once joint 4 is sign-converted, which both implementations do automatically, so all inputs and outputs follow the tables below.

---

### 1. KUKA LBR 7 R800 (Med / iiwa)
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

### 2. KUKA LBR 14 R820 (Med / iiwa)
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

### Technical Notes

* **Flange Face Transformation:** The $d_7$ parameter ($126\text{ mm}$) transforms the mathematical model directly to the native mechanical tool flange face. When integrating custom end-effectors, add your tool's center-of-mass/tool-center-point (TCP) parameters directly onto $d_7$ or append it as an explicit rigid-body extension matrix.
* **Standard DH Conversion:** If your solver requires the standard DH convention, shift the row alignments by one step ($\alpha_{i-1}, a_{i-1}$ move to row $i-1$). The physical geometry is unchanged.
