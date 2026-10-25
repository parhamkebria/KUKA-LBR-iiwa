"""Forward and inverse kinematics for the 7-DoF KUKA LBR iiwa.

Units: metres and radians. Poses are 4x4 homogeneous matrices (base frame).

    import numpy as np
    from KUKA import forward_kinematics, inverse_kinematics

    pose, param, conf, _ = forward_kinematics(np.deg2rad([10, 30, -20, 60, 15, 40, 5]))
    joints, _, _ = inverse_kinematics(pose, param, conf)
"""

import numpy as np

__all__ = [
    "LINK_LENGTHS",
    "JOINT_LIMITS",
    "KukaError",
    "configuration",
    "dh_matrix",
    "forward_kinematics",
    "inverse_kinematics",
    "within_joint_limits",
]

# base-shoulder, shoulder-elbow, elbow-wrist, wrist-flange
LINK_LENGTHS = (0.34, 0.4, 0.4, 0.126)

# +/- limits of A1..A7 in radians (170, 120, 170, 120, 170, 120, 175 deg)
JOINT_LIMITS = np.deg2rad([170, 120, 170, 120, 170, 120, 175])

_TOL = 1e-8


class KukaError(ValueError):
    """Raised for unreachable poses and singular configurations."""


def _dh_table():
    l = LINK_LENGTHS
    # columns: a, alpha, d
    return np.array([
        [0, -np.pi / 2, l[0]],
        [0, np.pi / 2, 0],
        [0, np.pi / 2, l[1]],
        [0, -np.pi / 2, 0],
        [0, -np.pi / 2, l[2]],
        [0, np.pi / 2, 0],
        [0, 0, l[3]],
    ])


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def _skew(v):
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])


def dh_matrix(a, alpha, d, theta):
    """Homogeneous transform for one Denavit-Hartenberg row."""
    ct, st = np.cos(theta), np.sin(theta)
    ca, sa = np.cos(alpha), np.sin(alpha)
    return np.array([[ct, -st * ca, st * sa, a * ct],
                     [st, ct * ca, -ct * sa, a * st],
                     [0, sa, ca, d],
                     [0, 0, 0, 1.0]])


def configuration(rconf):
    """Decode ``rconf`` into (arm, elbow, wrist), each +1 or -1."""
    rconf = int(rconf)
    return (-1 if rconf & 1 else 1,
            -1 if rconf & 2 else 1,
            -1 if rconf & 4 else 1)


def _wrist_position(pose):
    pose = np.asarray(pose, dtype=float)
    xs = np.array([0.0, 0.0, LINK_LENGTHS[0]])
    xw = pose[:3, 3] - pose[:3, :3] @ np.array([0.0, 0.0, LINK_LENGTHS[3]])
    return xs, xw


def _reference_plane(pose, elbow):
    """Reference plane of the virtual arm with joint 3 locked at 0.

    Returns (normal_vector, R03_reference, virtual_joints).
    """
    pose = np.asarray(pose, dtype=float)
    dh = _dh_table()
    _, lse, lew, _ = LINK_LENGTHS
    xs, xw = _wrist_position(pose)
    xsw = xw - xs
    dsw = np.linalg.norm(xsw)

    if not (lse - lew < dsw < lse + lew):
        raise KukaError("Specified pose outside reachable workspace")
    if abs(dsw**2 - lse**2 - lew**2 - 2 * lse * lew) <= 1e-6:
        raise KukaError("Elbow singularity. Tip at reach limit.")

    joints = np.zeros(7)
    joints[3] = elbow * np.arccos((dsw**2 - lse**2 - lew**2) / (2 * lse * lew))

    # Joint 1 is undefined when the wrist is on the vertical axis
    if np.linalg.norm(np.cross(xsw, [0, 0, 1])) > 1e-6:
        joints[0] = np.arctan2(xsw[1], xsw[0])

    r = np.hypot(xsw[0], xsw[1])
    phi = np.arccos(np.clip((lse**2 + dsw**2 - lew**2) / (2 * lse * dsw), -1, 1))
    joints[1] = np.arctan2(r, xsw[2]) + elbow * phi

    T = [dh_matrix(*dh[i], joints[i]) for i in range(4)]
    T04 = T[0] @ T[1] @ T[2] @ T[3]
    R03 = (T[0] @ T[1] @ T[2])[:3, :3]

    normal = np.cross(_unit(T04[:3, 3] - xs), _unit(xw - xs))
    return normal, R03, joints


def forward_kinematics(joints):
    """Forward kinematics.

    Parameters
    ----------
    joints : array-like, shape (7,), radians

    Returns
    -------
    pose : (4, 4) end-effector (flange) transform in the base frame
    param : arm angle in radians
    conf : int, robot configuration (bit0 arm, bit1 elbow, bit2 wrist)
    jout : joint values of the virtual (joint 3 = 0) arm
    """
    joints = np.asarray(joints, dtype=float).ravel()
    if joints.shape != (7,):
        raise ValueError("joints must contain exactly 7 values")

    rconf = int(joints[1] < 0) + 2 * int(joints[3] < 0) + 4 * int(joints[5] < 0)
    _, elbow, _ = configuration(rconf)

    dh = _dh_table()
    tr = []
    T = np.eye(4)
    for i in range(7):
        T = T @ dh_matrix(*dh[i], joints[i])
        tr.append(T)

    xs = tr[0][:3, 3]
    xe = tr[3][:3, 3]
    xw = tr[5][:3, 3]
    xsw = xw - xs
    pose = tr[-1]

    vv, _, jout = _reference_plane(pose, elbow)
    vc = np.cross(_unit(xe - xs), _unit(xsw))

    cos_ns = np.clip(np.dot(_unit(vv), _unit(vc)), -1, 1)
    v3 = np.cross(_unit(vv), _unit(vc))
    if np.linalg.norm(v3) > _TOL:
        nsparam = np.sign(np.dot(v3, xsw)) * np.arccos(cos_ns)
    elif np.linalg.norm(vv - vc) < _TOL:
        nsparam = 0.0
    else:
        nsparam = np.pi
    return pose, float(nsparam), rconf, jout


def inverse_kinematics(pose, param=0.0, conf=0):
    """Inverse kinematics for a given arm angle and configuration.

    Parameters
    ----------
    pose : (4, 4) desired flange transform in the base frame
    param : arm angle in radians
    conf : int, robot configuration (bit0 arm, bit1 elbow, bit2 wrist)

    Returns
    -------
    joints : (7,) joint values in radians (joint limits are NOT enforced;
        see ``within_joint_limits``)
    s_mat : (3, 3, 3) shoulder matrices As, Bs, Cs
    w_mat : (3, 3, 3) wrist matrices Aw, Bw, Cw

    Raises
    ------
    KukaError if the pose is unreachable or at an elbow singularity.
    """
    pose = np.asarray(pose, dtype=float)
    if pose.shape != (4, 4):
        raise ValueError("pose must be a 4x4 homogeneous matrix")

    arm, elbow, wrist = configuration(conf)
    dh = _dh_table()
    _, lse, lew, _ = LINK_LENGTHS

    xs, xw = _wrist_position(pose)
    xsw = xw - xs
    dsw = np.linalg.norm(xsw)
    usw = _unit(xsw)

    if not (lse - lew < dsw < lse + lew):
        raise KukaError("Specified pose outside reachable workspace")
    if abs(dsw**2 - lse**2 - lew**2 - 2 * lse * lew) <= _TOL:
        raise KukaError("Elbow singularity. Tip at reach limit.")

    joints = np.zeros(7)
    joints[3] = elbow * np.arccos((dsw**2 - lse**2 - lew**2) / (2 * lse * lew))
    R34 = dh_matrix(*dh[3], joints[3])[:3, :3]

    _, R03_o, _ = _reference_plane(pose, elbow)

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

    return joints, np.stack([As, Bs, Cs], axis=2), np.stack([Aw, Bw, Cw], axis=2)


def within_joint_limits(joints, limits=JOINT_LIMITS):
    """True if every joint is within its +/- limit (radians)."""
    return bool(np.all(np.abs(np.asarray(joints, dtype=float)) <= limits))


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    worst = 0.0
    n_ok = 0
    for _ in range(1000):
        q = rng.uniform(-JOINT_LIMITS, JOINT_LIMITS)
        try:
            pose, param, conf, _ = forward_kinematics(q)
            q_ik, _, _ = inverse_kinematics(pose, param, conf)
            pose2, *_ = forward_kinematics(q_ik)
        except KukaError:
            continue
        worst = max(worst, np.abs(pose - pose2).max())
        n_ok += 1
    print(f"FK->IK->FK round trip on {n_ok} random configs, max pose error: {worst:.2e}")
