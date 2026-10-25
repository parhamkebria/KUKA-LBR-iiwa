// Forward and inverse kinematics for the 7-DoF KUKA LBR iiwa / LBR Med.
// C++17 port of KUKA.py (no external dependencies).
//
// Supported models, selected with the `payload` argument (kg):
//   payload = 7   -> LBR 7 R800   (d = 310, 400, 400, 126 mm)
//   payload = 14  -> LBR 14 R820  (d = 360, 420, 400, 126 mm)
//
// The redundancy is parameterised by the arm angle `param` (psi) and the robot
// configuration `conf` (bit0 = arm, bit1 = elbow, bit2 = wrist; a bit is set
// when joint 2 / 4 / 6 is negative).
//
// Units: metres and radians. Poses are 4x4 homogeneous matrices (flange in
// base frame). Joint angles follow the DH tables in README.md.
//
//   #include "KUKA.h"
//   kuka::Joints q = {0.17, 0.52, -0.35, 1.05, 0.26, 0.70, 0.09};
//   auto fk = kuka::forward_kinematics(q, 7);
//   auto ik = kuka::inverse_kinematics(fk.pose, fk.param, fk.conf, 7);

#ifndef KUKA_H
#define KUKA_H

#include <array>
#include <map>
#include <stdexcept>
#include <string>

namespace kuka {

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;
using Mat4 = std::array<std::array<double, 4>, 4>;
using Joints = std::array<double, 7>;

// Thrown for unreachable poses and singular configurations.
class KukaError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct Model {
    std::string name;
    std::array<double, 4> l;  // d1, d3, d5, d7 in metres
};

struct FKResult {
    Mat4 pose;     // flange transform in the base frame
    double param;  // arm angle psi (rad)
    int conf;      // robot configuration (bit0 arm, bit1 elbow, bit2 wrist)
    Joints jout;   // joints of the virtual arm (joint 3 = 0) used as reference
};

struct IKResult {
    Joints joints;                // limits are NOT enforced; see within_joint_limits
    std::array<Mat3, 3> s_mat;    // shoulder matrices As, Bs, Cs
    std::array<Mat3, 3> w_mat;    // wrist matrices Aw, Bw, Cw
};

struct Configuration {
    int arm, elbow, wrist;  // each +1 or -1
};

// Model table keyed by payload (7, 14).
const std::map<int, Model>& models();

// +/- limits of joints 1..7 in radians (identical for both models).
extern const Joints JOINT_LIMITS;

// DH table of the selected model, rows [alpha_{i-1}, a_{i-1}, d_i] (README convention).
// Throws std::invalid_argument for an unsupported payload.
std::array<std::array<double, 3>, 7> dh_table(int payload = 7);

// Transform of one README DH row:
// Rot_x(alpha_{i-1}) * Trans_x(a_{i-1}) * Rot_z(theta_i) * Trans_z(d_i).
Mat4 dh_matrix(double alpha, double a, double d, double theta);

// Decode `conf` into (arm, elbow, wrist), each +1 or -1.
Configuration configuration(int conf);

// Throws std::invalid_argument (bad payload) or KukaError (unreachable / singular elbow).
FKResult forward_kinematics(const Joints& joints, int payload = 7);
IKResult inverse_kinematics(const Mat4& pose, double param = 0.0, int conf = 0, int payload = 7);

// True if every joint is within its +/- limit (radians).
bool within_joint_limits(const Joints& joints, int payload = 7);

}  // namespace kuka

#endif  // KUKA_H
