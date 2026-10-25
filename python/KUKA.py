"""Forward and inverse kinematics for the 7-DoF KUKA LBR iiwa / LBR Med.

Supported models, selected with the ``payload`` argument (kg):

    payload = 7   -> LBR 7 R800   (d = 310, 400, 400, 126 mm)
    payload = 14  -> LBR 14 R820  (d = 360, 420, 400, 126 mm)

The redundancy is parameterised by the arm angle ``param`` (psi) and the robot
configuration ``conf`` (bit0 = arm, bit1 = elbow, bit2 = wrist; a bit is set
when joint 2 / 4 / 6 is negative).

Units: metres and radians. Poses are 4x4 homogeneous matrices (flange in base
frame). Joint angles follow the DH tables in README.md.

    import numpy as np
    from KUKA import forward_kinematics, inverse_kinematics

    q = np.deg2rad([10, 30, -20, 60, 15, 40, 5])
    pose, param, conf, _ = forward_kinematics(q, payload=7)
    joints, _, _ = inverse_kinematics(pose, param, conf, payload=7)
"""

import numpy as np

__all__ = [
    "MODELS",
    "JOINT_LIMITS",
    "KukaError",
    "configuration",
    "dh_matrix",
    "dh_table",
    "forward_kinematics",
    "inverse_kinematics",
    "within_joint_limits",
]

# Link lengths d1, d3, d5, d7 in metres (README DH tables)
MODELS = {
    7: {"name": "LBR 7 R800", "l": (0.31, 0.40, 0.40, 0.126)},
    14: {"name": "LBR 14 R820", "l": (0.36, 0.42, 0.40, 0.126)},
}

# +/- limits of joints 1..7 in radians (identical for both models)
JOINT_LIMITS = np.deg2rad([170, 120, 170, 120, 170, 120, 175])

# The README tables define joint 4 with the opposite sense of rotation to the
# convention used by the analytic solver below; this vector converts between them.
_SIGN = np.array([1, 1, 1, -1, 1, 1, 1.0])

_TOL = 1e-8


class KukaError(ValueError):
    """Raised for unreachable poses and singular configurations."""


def _lengths(payload):
    try:
        return MODELS[int(payload)]["l"]
    except (KeyError, TypeError, ValueError):
        raise ValueError(
            f"payload must be one of {sorted(MODELS)} (7 -> LBR 7 R800, 14 -> LBR 14 R820), "
            f"got {payload!r}") from None


def dh_table(payload=7):
    """DH table of the selected model, rows [alpha_{i-1}, a_{i-1}, d_i] (README convention)."""
    l = _lengths(payload)
    h = np.pi / 2
    return np.array([
        [0, 0, l[0]],
        [-h, 0, 0],
        [h, 0, l[1]],
        [-h, 0, 0],
        [h, 0, l[2]],
        [-h, 0, 0],
        [h, 0, l[3]],
    ])


def _solver_table(l):
    # Rows [a, alpha, d] of the equivalent chain used by the analytic solver
    h = np.pi / 2
    return np.array([
        [0, -h, l[0]],
        [0, h, 0],
        [0, h, l[1]],
        [0, -h, 0],
        [0, -h, l[2]],
        [0, h, 0],
        [0, 0, l[3]],
    ])


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def _skew(v):
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])


def dh_matrix(alpha, a, d, theta):
    """Transform of one README DH row:
    Rot_x(alpha_{i-1}) * Trans_x(a_{i-1}) * Rot_z(theta_i) * Trans_z(d_i)."""
    ct, st = np.cos(theta), np.sin(theta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([[ct, -st, 0, a],
                     [st * ca, ct * ca, -sa, -sa * d],
                     [st * sa, ct * sa, ca, ca * d],
                     [0, 0, 0, 1.0]])


def _std_dh(a, alpha, d, theta):
    ct, st = np.cos(theta), np.sin(theta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([[ct, -st * ca, st * sa, a * ct],
                     [st, ct * ca, -ct * sa, a * st],
                     [0, sa, ca, d],
                     [0, 0, 0, 1.0]])


def configuration(conf):
    """Decode ``conf`` into (arm, elbow, wrist), each +1 or -1."""
    conf = int(conf)
    return (-1 if conf & 1 else 1,
            -1 if conf & 2 else 1,
            -1 if conf & 4 else 1)


def _shoulder_wrist(pose, l):
    xs = np.array([0.0, 0.0, l[0]])
    xw = pose[:3, 3] - pose[:3, :3] @ np.array([0.0, 0.0, l[3]])
    return xs, xw


def _elbow_angle(dsw, l, elbow, tol):
    """Joint 4 (solver convention) from the shoulder-wrist distance.

    Raises KukaError if the wrist is unreachable or the elbow is singular.
    """
    lse, lew = l[1], l[2]
    if not (abs(lse - lew) < dsw < lse + lew):
        raise KukaError("Specified pose outside reachable workspace")
    c4 = (dsw**2 - lse**2 - lew**2) / (2 * lse * lew)
    if 1 - abs(c4) <= tol:
        raise KukaError("Elbow singularity (arm fully stretched or folded)")
    return elbow * np.arccos(c4)


def _reference_plane(pose, elbow, l):
    """Reference plane of the virtual arm with joint 3 locked at 0.

    ``elbow`` uses the solver convention. Returns (normal, R03, joints).
    The reachability and singularity checks live here so FK and IK share them.
    """
    dh = _solver_table(l)
    lse, lew = l[1], l[2]
    xs, xw = _shoulder_wrist(pose, l)
    xsw = xw - xs
    dsw = np.linalg.norm(xsw)

    joints = np.zeros(7)
    joints[3] = _elbow_angle(dsw, l, elbow, 1e-6)

    # Joint 1 is undefined when the wrist is on the vertical axis
    if np.linalg.norm(np.cross(xsw, [0, 0, 1])) > 1e-6:
        joints[0] = np.arctan2(xsw[1], xsw[0])

    phi = np.arccos(np.clip((lse**2 + dsw**2 - lew**2) / (2 * lse * dsw), -1, 1))
    joints[1] = np.arctan2(np.hypot(xsw[0], xsw[1]), xsw[2]) + elbow * phi

    T = [_std_dh(*dh[i], joints[i]) for i in range(4)]
    T03 = T[0] @ T[1] @ T[2]
    T04 = T03 @ T[3]

    normal = np.cross(_unit(T04[:3, 3] - xs), _unit(xw - xs))
    return normal, T03[:3, :3], joints


def forward_kinematics(joints, payload=7):
    """Forward kinematics.

    Parameters
    ----------
    joints : array-like, shape (7,), radians
    payload : 7 (LBR 7 R800) or 14 (LBR 14 R820)

    Returns
    -------
    pose : (4, 4) flange transform in the base frame
    param : arm angle psi in radians
    conf : int, robot configuration (bit0 arm, bit1 elbow, bit2 wrist)
    jout : (7,) joint values of the virtual arm (joint 3 = 0) used as reference
    """
    l = _lengths(payload)
    joints = np.asarray(joints, dtype=float).ravel()
    if joints.shape != (7,):
        raise ValueError("joints must contain exactly 7 values")

    conf = int(joints[1] < 0) + 2 * int(joints[3] < 0) + 4 * int(joints[5] < 0)
    elbow = -configuration(conf)[1]  # solver convention

    dh = dh_table(payload)
    tr = []
    T = np.eye(4)
    for i in range(7):
        T = T @ dh_matrix(dh[i, 0], dh[i, 1], dh[i, 2], joints[i])
        tr.append(T)

    xs, xe, xw = tr[0][:3, 3], tr[3][:3, 3], tr[5][:3, 3]
    xsw = xw - xs
    pose = tr[-1]

    vv, _, jout = _reference_plane(pose, elbow, l)
    vc = np.cross(_unit(xe - xs), _unit(xsw))

    cos_ns = np.clip(np.dot(_unit(vv), _unit(vc)), -1, 1)
    v3 = np.cross(_unit(vv), _unit(vc))
    if np.linalg.norm(v3) > _TOL:
        param = np.sign(np.dot(v3, xsw)) * np.arccos(cos_ns)
    elif np.linalg.norm(vv - vc) < _TOL:
        param = 0.0
    else:
        param = np.pi
    return pose, float(param), conf, jout * _SIGN


def inverse_kinematics(pose, param=0.0, conf=0, payload=7):
    """Inverse kinematics for a given arm angle and configuration.

    Parameters
    ----------
    pose : (4, 4) desired flange transform in the base frame
    param : arm angle psi in radians
    conf : int, robot configuration (bit0 arm, bit1 elbow, bit2 wrist)
    payload : 7 (LBR 7 R800) or 14 (LBR 14 R820)

    Returns
    -------
    joints : (7,) joint values in radians (limits are NOT enforced; see
        ``within_joint_limits``)
    s_mat : (3, 3, 3) shoulder matrices As, Bs, Cs
    w_mat : (3, 3, 3) wrist matrices Aw, Bw, Cw

    Raises
    ------
    KukaError if the pose is unreachable or at an elbow singularity.
    """
    l = _lengths(payload)
    pose = np.asarray(pose, dtype=float)
    if pose.shape != (4, 4):
        raise ValueError("pose must be a 4x4 homogeneous matrix")

    arm, elbow_user, wrist = configuration(conf)
    elbow = -elbow_user  # solver convention
    dh = _solver_table(l)

    xs, xw = _shoulder_wrist(pose, l)
    xsw = xw - xs
    usw = _unit(xsw)

    joints = np.zeros(7)
    joints[3] = _elbow_angle(np.linalg.norm(xsw), l, elbow, _TOL)
    R34 = _std_dh(*dh[3], joints[3])[:3, :3]

    _, R03_o, _ = _reference_plane(pose, elbow, l)

    skew_usw = _skew(usw)
    As = skew_usw @ R03_o
    Bs = -skew_usw @ skew_usw @ R03_o
    Cs = np.outer(usw, usw) @ R03_o

    psi = param
    R03 = As * np.sin(psi) + Bs * np.cos(psi) + Cs

    joints[0] = np.arctan2(arm * R03[1, 1], arm * R03[0, 1])
    joints[1] = arm * np.arccos(np.clip(R03[2, 1], -1, 1))
    joints[2] = np.arctan2(arm * -R03[2, 2], arm * -R03[2, 0])

    R07 = pose[:3, :3]
    Aw = R34.T @ As.T @ R07
    Bw = R34.T @ Bs.T @ R07
    Cw = R34.T @ Cs.T @ R07
    R47 = Aw * np.sin(psi) + Bw * np.cos(psi) + Cw

    joints[4] = np.arctan2(wrist * R47[1, 2], wrist * R47[0, 2])
    joints[5] = wrist * np.arccos(np.clip(R47[2, 2], -1, 1))
    joints[6] = np.arctan2(wrist * R47[2, 1], wrist * -R47[2, 0])

    return joints * _SIGN, np.stack([As, Bs, Cs], axis=2), np.stack([Aw, Bw, Cw], axis=2)


def within_joint_limits(joints, payload=7):
    """True if every joint is within its +/- limit (radians)."""
    _lengths(payload)
    return bool(np.all(np.abs(np.asarray(joints, dtype=float)) <= JOINT_LIMITS))


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    for payload, model in MODELS.items():
        worst, n_ok = 0.0, 0
        for _ in range(1000):
            q = rng.uniform(-JOINT_LIMITS, JOINT_LIMITS)
            try:
                pose, psi, conf, _ = forward_kinematics(q, payload)
                q_ik, _, _ = inverse_kinematics(pose, psi, conf, payload)
                pose2, *_ = forward_kinematics(q_ik, payload)
            except KukaError:
                continue
            worst = max(worst, np.abs(pose - pose2).max())
            n_ok += 1
        print(f"{model['name']}: {n_ok} FK->IK->FK round trips, max pose error {worst:.2e}")
