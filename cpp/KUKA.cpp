#include "KUKA.h"

#include <algorithm>
#include <cmath>

namespace kuka {

namespace {

const double PI = 3.14159265358979323846;
const double TOL = 1e-8;

using DhRow = std::array<double, 3>;
using DhTable = std::array<DhRow, 7>;

// The README tables define joint 4 with the opposite sense of rotation to the
// convention used by the analytic solver; this vector converts between them.
const Joints SIGN = {1, 1, 1, -1, 1, 1, 1};

const std::array<double, 4>& lengths(int payload) {
    const auto& m = models();
    auto it = m.find(payload);
    if (it == m.end()) {
        throw std::invalid_argument("payload must be one of [7, 14] (7 -> LBR 7 R800, 14 -> LBR 14 R820), got " +
                                    std::to_string(payload));
    }
    return it->second.l;
}

// ---- small linear algebra helpers ----------------------------------------

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Vec3 unit(const Vec3& v) {
    double n = norm(v);
    return n > 0 ? Vec3{v[0] / n, v[1] / n, v[2] / n} : v;
}
double sign(double x) { return (x > 0) - (x < 0); }
double clip1(double x) { return std::max(-1.0, std::min(1.0, x)); }

Mat3 mul(const Mat3& A, const Mat3& B) {
    Mat3 C{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) C[i][j] += A[i][k] * B[k][j];
    return C;
}
Mat4 mul(const Mat4& A, const Mat4& B) {
    Mat4 C{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k) C[i][j] += A[i][k] * B[k][j];
    return C;
}
Mat3 transpose(const Mat3& A) {
    Mat3 T{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) T[i][j] = A[j][i];
    return T;
}
Vec3 mul(const Mat3& A, const Vec3& v) {
    return {dot(A[0], v), dot(A[1], v), dot(A[2], v)};
}
Mat3 rotation(const Mat4& T) {
    Mat3 R{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) R[i][j] = T[i][j];
    return R;
}
Vec3 translation(const Mat4& T) { return {T[0][3], T[1][3], T[2][3]}; }
Mat4 identity4() {
    Mat4 I{};
    for (int i = 0; i < 4; ++i) I[i][i] = 1.0;
    return I;
}
Mat3 skew(const Vec3& v) {
    return {{{0, -v[2], v[1]}, {v[2], 0, -v[0]}, {-v[1], v[0], 0}}};
}
Mat3 scale(const Mat3& A, double s) {
    Mat3 C = A;
    for (auto& r : C)
        for (auto& x : r) x *= s;
    return C;
}
Mat3 add(const Mat3& A, const Mat3& B) {
    Mat3 C{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) C[i][j] = A[i][j] + B[i][j];
    return C;
}
Mat3 outer(const Vec3& a, const Vec3& b) {
    Mat3 C{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) C[i][j] = a[i] * b[j];
    return C;
}

// ---- kinematic chain ------------------------------------------------------

// Rows [a, alpha, d] of the equivalent chain used by the analytic solver.
DhTable solver_table(const std::array<double, 4>& l) {
    const double h = PI / 2;
    return {{{0, -h, l[0]}, {0, h, 0}, {0, h, l[1]}, {0, -h, 0}, {0, -h, l[2]}, {0, h, 0}, {0, 0, l[3]}}};
}

Mat4 std_dh(double a, double alpha, double d, double theta) {
    double ct = std::cos(theta), st = std::sin(theta), ca = std::cos(alpha), sa = std::sin(alpha);
    return {{{ct, -st * ca, st * sa, a * ct}, {st, ct * ca, -ct * sa, a * st}, {0, sa, ca, d}, {0, 0, 0, 1}}};
}

void shoulder_wrist(const Mat4& pose, const std::array<double, 4>& l, Vec3& xs, Vec3& xw) {
    xs = {0.0, 0.0, l[0]};
    Vec3 off = mul(rotation(pose), Vec3{0.0, 0.0, l[3]});
    xw = sub(translation(pose), off);
}

// Joint 4 (solver convention) from the shoulder-wrist distance.
double elbow_angle(double dsw, const std::array<double, 4>& l, double elbow, double tol) {
    double lse = l[1], lew = l[2];
    if (!(std::fabs(lse - lew) < dsw && dsw < lse + lew)) throw KukaError("Specified pose outside reachable workspace");
    double c4 = (dsw * dsw - lse * lse - lew * lew) / (2 * lse * lew);
    if (1 - std::fabs(c4) <= tol) throw KukaError("Elbow singularity (arm fully stretched or folded)");
    return elbow * std::acos(c4);
}

struct RefPlane {
    Vec3 normal;
    Mat3 R03;
    Joints joints;
};

// Reference plane of the virtual arm with joint 3 locked at 0. `elbow` uses the
// solver convention. Reachability and singularity checks live here so FK and IK share them.
RefPlane reference_plane(const Mat4& pose, double elbow, const std::array<double, 4>& l) {
    DhTable dh = solver_table(l);
    double lse = l[1], lew = l[2];
    Vec3 xs, xw;
    shoulder_wrist(pose, l, xs, xw);
    Vec3 xsw = sub(xw, xs);
    double dsw = norm(xsw);

    Joints j{};
    j[3] = elbow_angle(dsw, l, elbow, 1e-6);

    // Joint 1 is undefined when the wrist is on the vertical axis
    if (norm(cross(xsw, Vec3{0, 0, 1})) > 1e-6) j[0] = std::atan2(xsw[1], xsw[0]);

    double phi = std::acos(clip1((lse * lse + dsw * dsw - lew * lew) / (2 * lse * dsw)));
    j[1] = std::atan2(std::hypot(xsw[0], xsw[1]), xsw[2]) + elbow * phi;

    Mat4 T03 = identity4();
    for (int i = 0; i < 3; ++i) T03 = mul(T03, std_dh(dh[i][0], dh[i][1], dh[i][2], j[i]));
    Mat4 T04 = mul(T03, std_dh(dh[3][0], dh[3][1], dh[3][2], j[3]));

    RefPlane rp;
    rp.normal = cross(unit(sub(translation(T04), xs)), unit(sub(xw, xs)));
    rp.R03 = rotation(T03);
    rp.joints = j;
    return rp;
}

}  // namespace

const std::map<int, Model>& models() {
    static const std::map<int, Model> m = {
        {7, {"LBR 7 R800", {0.31, 0.40, 0.40, 0.126}}},
        {14, {"LBR 14 R820", {0.36, 0.42, 0.40, 0.126}}},
    };
    return m;
}

const Joints JOINT_LIMITS = {170 * PI / 180, 120 * PI / 180, 170 * PI / 180, 120 * PI / 180,
                             170 * PI / 180, 120 * PI / 180, 175 * PI / 180};

std::array<std::array<double, 3>, 7> dh_table(int payload) {
    const auto& l = lengths(payload);
    const double h = PI / 2;
    return {{{0, 0, l[0]}, {-h, 0, 0}, {h, 0, l[1]}, {-h, 0, 0}, {h, 0, l[2]}, {-h, 0, 0}, {h, 0, l[3]}}};
}

Mat4 dh_matrix(double alpha, double a, double d, double theta) {
    double ct = std::cos(theta), st = std::sin(theta), ca = std::cos(alpha), sa = std::sin(alpha);
    return {{{ct, -st, 0, a},
             {st * ca, ct * ca, -sa, -sa * d},
             {st * sa, ct * sa, ca, ca * d},
             {0, 0, 0, 1}}};
}

Configuration configuration(int conf) {
    return {(conf & 1) ? -1 : 1, (conf & 2) ? -1 : 1, (conf & 4) ? -1 : 1};
}

FKResult forward_kinematics(const Joints& joints, int payload) {
    const auto& l = lengths(payload);

    int conf = (joints[1] < 0) + 2 * (joints[3] < 0) + 4 * (joints[5] < 0);
    double elbow = -configuration(conf).elbow;  // solver convention

    auto dh = dh_table(payload);
    std::array<Mat4, 7> tr;
    Mat4 T = identity4();
    for (int i = 0; i < 7; ++i) {
        T = mul(T, dh_matrix(dh[i][0], dh[i][1], dh[i][2], joints[i]));
        tr[i] = T;
    }

    Vec3 xs = translation(tr[0]), xe = translation(tr[3]), xw = translation(tr[5]);
    Vec3 xsw = sub(xw, xs);
    Mat4 pose = tr[6];

    RefPlane rp = reference_plane(pose, elbow, l);
    Vec3 vv = rp.normal;
    Vec3 vc = cross(unit(sub(xe, xs)), unit(xsw));

    double cos_ns = clip1(dot(unit(vv), unit(vc)));
    Vec3 v3 = cross(unit(vv), unit(vc));
    double param;
    if (norm(v3) > TOL) {
        param = sign(dot(v3, xsw)) * std::acos(cos_ns);
    } else if (norm(sub(vv, vc)) < TOL) {
        param = 0.0;
    } else {
        param = PI;
    }

    Joints jout;
    for (int i = 0; i < 7; ++i) jout[i] = rp.joints[i] * SIGN[i];
    return {pose, param, conf, jout};
}

IKResult inverse_kinematics(const Mat4& pose, double param, int conf, int payload) {
    const auto& l = lengths(payload);

    Configuration c = configuration(conf);
    double arm = c.arm, wrist = c.wrist;
    double elbow = -c.elbow;  // solver convention
    DhTable dh = solver_table(l);

    Vec3 xs, xw;
    shoulder_wrist(pose, l, xs, xw);
    Vec3 xsw = sub(xw, xs);
    Vec3 usw = unit(xsw);

    Joints j{};
    j[3] = elbow_angle(norm(xsw), l, elbow, TOL);
    Mat3 R34 = rotation(std_dh(dh[3][0], dh[3][1], dh[3][2], j[3]));

    RefPlane rp = reference_plane(pose, elbow, l);
    const Mat3& R03_o = rp.R03;

    Mat3 S = skew(usw);
    Mat3 As = mul(S, R03_o);
    Mat3 Bs = scale(mul(mul(S, S), R03_o), -1.0);
    Mat3 Cs = mul(outer(usw, usw), R03_o);

    double s = std::sin(param), co = std::cos(param);
    Mat3 R03 = add(add(scale(As, s), scale(Bs, co)), Cs);

    j[0] = std::atan2(arm * R03[1][1], arm * R03[0][1]);
    j[1] = arm * std::acos(clip1(R03[2][1]));
    j[2] = std::atan2(arm * -R03[2][2], arm * -R03[2][0]);

    Mat3 R07 = rotation(pose);
    Mat3 R34t = transpose(R34);
    Mat3 Aw = mul(mul(R34t, transpose(As)), R07);
    Mat3 Bw = mul(mul(R34t, transpose(Bs)), R07);
    Mat3 Cw = mul(mul(R34t, transpose(Cs)), R07);
    Mat3 R47 = add(add(scale(Aw, s), scale(Bw, co)), Cw);

    j[4] = std::atan2(wrist * R47[1][2], wrist * R47[0][2]);
    j[5] = wrist * std::acos(clip1(R47[2][2]));
    j[6] = std::atan2(wrist * R47[2][1], wrist * -R47[2][0]);

    IKResult out;
    for (int i = 0; i < 7; ++i) out.joints[i] = j[i] * SIGN[i];
    out.s_mat = {As, Bs, Cs};
    out.w_mat = {Aw, Bw, Cw};
    return out;
}

bool within_joint_limits(const Joints& joints, int payload) {
    lengths(payload);  // validates payload
    for (int i = 0; i < 7; ++i)
        if (std::fabs(joints[i]) > JOINT_LIMITS[i]) return false;
    return true;
}

}  // namespace kuka
