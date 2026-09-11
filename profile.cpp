// profile - rest-to-rest motion profiles for a single axis.
//
// Two generators:
//
//   trapezoidal   ramp up at amax, cruise at vmax, ramp down. Acceleration
//                 steps instantly, so jerk is infinite at four points. Every
//                 hobby servo library does this, and it is why cheap gantries
//                 ring after every move.
//
//   double-S      jerk-limited. Acceleration itself ramps, so jerk stays
//                 bounded. The move takes longer on paper and settles far
//                 faster in practice, because you stop exciting the first
//                 structural mode of the machine.
//
// Both handle the cases that break naive implementations: a move too short
// to reach vmax, a move too short to even reach amax, and moves in the
// negative direction.
//
// Build:  g++ -std=c++17 -O2 -o profile profile.cpp
// Test:   ./profile --test
// CSV:    ./profile --csv --distance 400 --vmax 300 --amax 1200 --jmax 9000

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace mp {

struct Limits {
    double vmax = 1.0;   // units per second
    double amax = 1.0;   // units per second squared
    double jmax = 1.0;   // units per second cubed, double-S only
};

struct State {
    double t = 0, q = 0, v = 0, a = 0, j = 0;
};

// A stretch of time over which jerk is constant.
struct Segment {
    double dt = 0;
    double j = 0;
    const char* name = "";
};

// -------------------------------------------------------------------------
// A profile is just a list of constant-jerk segments plus a start state.
// Integrating a cubic is exact, so sampling never accumulates error.
// -------------------------------------------------------------------------

class Profile {
public:
    Profile(double q0, std::vector<Segment> segs, double a_start = 0.0)
        : q0_(q0), seg_(std::move(segs)), a_start_(a_start) {
        double t = 0;
        for (const auto& s : seg_) t += s.dt;
        duration_ = t;
    }

    double duration() const { return duration_; }
    const std::vector<Segment>& segments() const { return seg_; }

    State sample(double t) const {
        State s{0, q0_, 0, a_start_, 0};
        double elapsed = 0;
        for (const auto& seg : seg_) {
            const double dt = std::min(std::max(0.0, t - elapsed), seg.dt);
            advance(s, seg.j, dt);
            elapsed += seg.dt;
            if (t <= elapsed) {
                s.t = t;
                // Jerk is right-continuous inside a segment, zero past the end.
                s.j = (t < elapsed) ? seg.j : 0.0;
                return s;
            }
        }
        s.t = t;
        s.j = 0;
        s.a = 0;
        s.v = 0;
        return s;
    }

private:
    static void advance(State& s, double j, double dt) {
        s.q += s.v * dt + 0.5 * s.a * dt * dt + j * dt * dt * dt / 6.0;
        s.v += s.a * dt + 0.5 * j * dt * dt;
        s.a += j * dt;
    }

    double q0_;
    std::vector<Segment> seg_;
    double a_start_;
    double duration_ = 0;
};

// -------------------------------------------------------------------------
// trapezoidal
// -------------------------------------------------------------------------
// Modelled as constant-jerk segments of zero jerk with the acceleration set
// by an impulse at each boundary. To keep one evaluator, the impulses are
// folded in as very short segments would be; instead the acceleration steps
// are expressed by starting each phase from the right state, which the
// builder below does by handing over a matching a_start per phase. Simpler:
// build it as three zero-jerk phases and carry acceleration through a small
// wrapper.

class Trapezoid {
public:
    Trapezoid(double q0, double q1, Limits lim) : q0_(q0), lim_(lim) {
        const double D = std::fabs(q1 - q0);
        dir_ = (q1 >= q0) ? 1.0 : -1.0;

        double t_ramp = lim.vmax / lim.amax;
        double d_ramp = 0.5 * lim.amax * t_ramp * t_ramp;

        if (2 * d_ramp <= D) {                       // vmax is reached
            v_peak_ = lim.vmax;
            t_ramp_ = t_ramp;
            t_flat_ = (D - 2 * d_ramp) / lim.vmax;
        } else {                                     // triangular
            t_ramp_ = std::sqrt(D / lim.amax);
            v_peak_ = lim.amax * t_ramp_;
            t_flat_ = 0;
        }
        duration_ = 2 * t_ramp_ + t_flat_;
        reached_vmax_ = (t_flat_ > 0);
    }

    double duration() const { return duration_; }
    double peak_velocity() const { return v_peak_; }
    bool reached_vmax() const { return reached_vmax_; }

    State sample(double t) const {
        State s;
        s.t = t;
        const double a = lim_.amax;
        double q, v, acc;
        if (t <= 0) {
            q = 0; v = 0; acc = 0;
        } else if (t < t_ramp_) {
            acc = a;
            v = a * t;
            q = 0.5 * a * t * t;
        } else if (t < t_ramp_ + t_flat_) {
            const double td = t - t_ramp_;
            acc = 0;
            v = v_peak_;
            q = 0.5 * a * t_ramp_ * t_ramp_ + v_peak_ * td;
        } else if (t < duration_) {
            const double td = t - t_ramp_ - t_flat_;
            acc = -a;
            v = v_peak_ - a * td;
            q = 0.5 * a * t_ramp_ * t_ramp_ + v_peak_ * t_flat_ +
                v_peak_ * td - 0.5 * a * td * td;
        } else {
            q = 0.5 * a * t_ramp_ * t_ramp_ + v_peak_ * t_flat_ +
                0.5 * a * t_ramp_ * t_ramp_;
            v = 0;
            acc = 0;
        }
        s.q = q0_ + dir_ * q;
        s.v = dir_ * v;
        s.a = dir_ * acc;
        s.j = 0;                    // undefined at the corners, zero elsewhere
        return s;
    }

private:
    double q0_, dir_, v_peak_ = 0, t_ramp_ = 0, t_flat_ = 0, duration_ = 0;
    bool reached_vmax_ = false;
    Limits lim_;
};

// -------------------------------------------------------------------------
// double-S, jerk limited, rest to rest
// -------------------------------------------------------------------------

struct DoubleSPlan {
    double Tj = 0;     // duration of each constant-jerk ramp
    double Ta = 0;     // total acceleration phase
    double Tv = 0;     // cruise phase
    double v_peak = 0;
    double a_peak = 0;
    bool reached_vmax = false;
    bool reached_amax = false;
    double duration = 0;
};

DoubleSPlan plan_double_s(double D, Limits lim) {
    DoubleSPlan p;
    D = std::fabs(D);
    if (D <= 0) return p;

    // Step 1: what does the acceleration phase look like if we do reach vmax?
    double Tj, Ta;
    bool amax_hit;
    if (lim.vmax * lim.jmax >= lim.amax * lim.amax) {
        amax_hit = true;
        Tj = lim.amax / lim.jmax;
        Ta = Tj + lim.vmax / lim.amax;
    } else {
        amax_hit = false;
        Tj = std::sqrt(lim.vmax / lim.jmax);
        Ta = 2 * Tj;
    }

    // A symmetric ramp from rest to v covers exactly v*Ta/2, whatever shape
    // the acceleration takes, because velocity is antisymmetric about v/2.
    const double d_ramp = lim.vmax * Ta / 2.0;

    if (D >= 2 * d_ramp) {                         // cruise exists
        p.reached_vmax = true;
        p.reached_amax = amax_hit;
        p.Tj = Tj;
        p.Ta = Ta;
        p.v_peak = lim.vmax;
        p.Tv = (D - 2 * d_ramp) / lim.vmax;
    } else {
        p.reached_vmax = false;
        // Solve v * Ta(v) = D. First assume amax is still reached, so
        // Ta = amax/jmax + v/amax, giving a quadratic in v.
        const double k = lim.amax / lim.jmax;
        const double v = lim.amax *
                         (-k + std::sqrt(k * k + 4 * D / lim.amax)) / 2.0;
        if (v >= lim.amax * lim.amax / lim.jmax) {
            p.reached_amax = true;
            p.v_peak = v;
            p.Tj = lim.amax / lim.jmax;
            p.Ta = p.Tj + v / lim.amax;
        } else {
            // amax is never reached either: the acceleration is a pure
            // triangle, so Ta = 2*sqrt(v/jmax) and v^(3/2) = D*sqrt(j)/2.
            p.reached_amax = false;
            p.v_peak = std::cbrt(D * D * lim.jmax / 4.0);
            p.Tj = std::sqrt(p.v_peak / lim.jmax);
            p.Ta = 2 * p.Tj;
        }
        p.Tv = 0;
    }

    p.a_peak = lim.jmax * p.Tj;
    p.duration = 2 * p.Ta + p.Tv;
    return p;
}

Profile build_double_s(double q0, double q1, Limits lim, DoubleSPlan* out) {
    const double D = q1 - q0;
    const double dir = (D >= 0) ? 1.0 : -1.0;
    DoubleSPlan p = plan_double_s(D, lim);
    if (out) *out = p;

    const double j = dir * lim.jmax;
    const double flat = std::max(0.0, p.Ta - 2 * p.Tj);
    std::vector<Segment> segs{
        {p.Tj, j, "jerk up"},
        {flat, 0, "hold amax"},
        {p.Tj, -j, "jerk down"},
        {p.Tv, 0, "cruise"},
        {p.Tj, -j, "jerk down"},
        {flat, 0, "hold -amax"},
        {p.Tj, j, "jerk up"},
    };
    return Profile(q0, segs);
}

}  // namespace mp

// -------------------------------------------------------------------------
// tests
// -------------------------------------------------------------------------

namespace {

int g_fail = 0;

void check(bool ok, const char* what, const std::string& note = "") {
    if (!ok) ++g_fail;
    std::printf("  [%s] %-50s %s\n", ok ? "pass" : "FAIL", what, note.c_str());
}

template <typename P>
void extremes(const P& prof, double T, double& vmax, double& amax,
              double& jmax, double& qend) {
    vmax = amax = jmax = 0;
    const int N = 20000;
    for (int i = 0; i <= N; ++i) {
        auto s = prof.sample(T * i / N);
        vmax = std::max(vmax, std::fabs(s.v));
        amax = std::max(amax, std::fabs(s.a));
        jmax = std::max(jmax, std::fabs(s.j));
    }
    qend = prof.sample(T).q;
}

int run_tests() {
    using namespace mp;
    char buf[160];
    Limits lim{300.0, 1200.0, 9000.0};

    // --- trapezoid, long move --------------------------------------------
    {
        Trapezoid tp(0, 100, lim);
        double v, a, j, qend;
        extremes(tp, tp.duration(), v, a, j, qend);
        std::snprintf(buf, sizeof buf, "T %.5f s, end %.6f", tp.duration(), qend);
        check(std::fabs(qend - 100) < 1e-9, "trapezoid lands exactly on target", buf);
        std::snprintf(buf, sizeof buf, "v %.4f <= %.1f, a %.4f <= %.1f", v,
                      lim.vmax, a, lim.amax);
        check(v <= lim.vmax + 1e-6 && a <= lim.amax + 1e-6,
              "trapezoid respects vmax and amax", buf);
        check(tp.reached_vmax(), "long move reaches cruise velocity");
    }

    // --- trapezoid, short move becomes triangular ------------------------
    {
        Trapezoid tp(0, 2.0, lim);
        double v, a, j, qend;
        extremes(tp, tp.duration(), v, a, j, qend);
        std::snprintf(buf, sizeof buf, "peak v %.4f of vmax %.1f",
                      tp.peak_velocity(), lim.vmax);
        check(!tp.reached_vmax() && std::fabs(qend - 2.0) < 1e-9,
              "short move collapses to a triangle", buf);
    }

    // --- double-S, long move ---------------------------------------------
    {
        DoubleSPlan p;
        auto pr = build_double_s(0, 400, lim, &p);
        double v, a, j, qend;
        extremes(pr, pr.duration(), v, a, j, qend);
        std::snprintf(buf, sizeof buf, "T %.5f s, end %.9f", pr.duration(), qend);
        check(std::fabs(qend - 400) < 1e-7, "double-S lands exactly on target", buf);
        std::snprintf(buf, sizeof buf, "v %.3f  a %.3f  j %.1f", v, a, j);
        check(v <= lim.vmax + 1e-6 && a <= lim.amax + 1e-6 &&
                  j <= lim.jmax + 1e-6,
              "double-S respects vmax, amax and jmax", buf);
        check(p.reached_vmax && p.reached_amax,
              "long move reaches both plateaus");
        std::snprintf(buf, sizeof buf, "final v %.3e, final a %.3e",
                      pr.sample(pr.duration()).v, pr.sample(pr.duration()).a);
        check(std::fabs(pr.sample(pr.duration()).v) < 1e-9 &&
                  std::fabs(pr.sample(pr.duration()).a) < 1e-9,
              "double-S ends at rest with zero acceleration", buf);
    }

    // --- double-S, vmax never reached ------------------------------------
    {
        DoubleSPlan p;
        auto pr = build_double_s(0, 5.0, lim, &p);
        double v, a, j, qend;
        extremes(pr, pr.duration(), v, a, j, qend);
        std::snprintf(buf, sizeof buf, "peak v %.4f, amax reached %s",
                      p.v_peak, p.reached_amax ? "yes" : "no");
        check(!p.reached_vmax && std::fabs(qend - 5.0) < 1e-8 &&
                  v <= lim.vmax + 1e-6,
              "short move drops the cruise phase", buf);
    }

    // --- double-S, amax never reached either ------------------------------
    {
        DoubleSPlan p;
        auto pr = build_double_s(0, 0.05, lim, &p);
        double v, a, j, qend;
        extremes(pr, pr.duration(), v, a, j, qend);
        std::snprintf(buf, sizeof buf, "peak a %.4f of amax %.1f", p.a_peak,
                      lim.amax);
        check(!p.reached_amax && std::fabs(qend - 0.05) < 1e-9 &&
                  a <= lim.amax + 1e-6,
              "very short move keeps acceleration triangular", buf);
    }

    // --- negative direction ----------------------------------------------
    {
        DoubleSPlan p;
        auto pr = build_double_s(50, -30, lim, &p);
        double v, a, j, qend;
        extremes(pr, pr.duration(), v, a, j, qend);
        std::snprintf(buf, sizeof buf, "end %.9f", qend);
        check(std::fabs(qend + 30) < 1e-7 && v <= lim.vmax + 1e-6,
              "negative moves mirror correctly", buf);
    }

    // --- the whole point: same move, less jerk ---------------------------
    {
        Trapezoid tp(0, 400, lim);
        DoubleSPlan p;
        auto ds = build_double_s(0, 400, lim, &p);
        const double penalty =
            (ds.duration() - tp.duration()) / tp.duration() * 100.0;
        std::snprintf(buf, sizeof buf,
                      "trap %.5f s, double-S %.5f s, %.1f%% slower",
                      tp.duration(), ds.duration(), penalty);
        check(ds.duration() > tp.duration() && penalty < 40.0,
              "jerk limiting costs time, but not much", buf);
    }

    // --- position is the integral of velocity ----------------------------
    {
        DoubleSPlan p;
        auto pr = build_double_s(0, 100, lim, &p);
        const int N = 200000;
        const double h = pr.duration() / N;
        double q = 0;
        for (int i = 0; i < N; ++i)                 // trapezoid rule
            q += 0.5 * (pr.sample(i * h).v + pr.sample((i + 1) * h).v) * h;
        std::snprintf(buf, sizeof buf, "integrated %.6f vs sampled %.6f", q,
                      pr.sample(pr.duration()).q);
        check(std::fabs(q - 100) < 1e-4,
              "numerically integrating v reproduces q", buf);
    }

    std::printf("\n  %s\n",
                g_fail == 0 ? "all checks passed" : "SOME CHECKS FAILED");
    return g_fail == 0 ? 0 : 1;
}

// -------------------------------------------------------------------------

int run_csv(double distance, mp::Limits lim, int samples) {
    using namespace mp;
    Trapezoid tp(0, distance, lim);
    DoubleSPlan p;
    auto ds = build_double_s(0, distance, lim, &p);
    const double T = std::max(tp.duration(), ds.duration());

    std::fprintf(stderr,
                 "trapezoid %.5f s   double-S %.5f s   penalty %.1f%%\n",
                 tp.duration(), ds.duration(),
                 (ds.duration() - tp.duration()) / tp.duration() * 100);
    std::fprintf(stderr, "double-S peak v %.3f  peak a %.3f  vmax hit %s  amax hit %s\n",
                 p.v_peak, p.a_peak, p.reached_vmax ? "yes" : "no",
                 p.reached_amax ? "yes" : "no");

    std::printf("t,trap_q,trap_v,trap_a,ds_q,ds_v,ds_a,ds_j\n");
    for (int i = 0; i <= samples; ++i) {
        const double t = T * 1.05 * i / samples;
        auto A = tp.sample(std::min(t, tp.duration()));
        auto B = ds.sample(std::min(t, ds.duration()));
        std::printf("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", t, A.q, A.v,
                    A.a, B.q, B.v, B.a, B.j);
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    mp::Limits lim{300.0, 1200.0, 9000.0};
    double distance = 400.0;
    int samples = 2000;
    bool csv = false;

    for (int i = 1; i < argc; ++i) {
        auto next = [&](double& dst) {
            if (i + 1 < argc) dst = std::stod(argv[++i]);
        };
        if (!std::strcmp(argv[i], "--test")) return run_tests();
        else if (!std::strcmp(argv[i], "--csv")) csv = true;
        else if (!std::strcmp(argv[i], "--distance")) next(distance);
        else if (!std::strcmp(argv[i], "--vmax")) next(lim.vmax);
        else if (!std::strcmp(argv[i], "--amax")) next(lim.amax);
        else if (!std::strcmp(argv[i], "--jmax")) next(lim.jmax);
        else if (!std::strcmp(argv[i], "--samples")) {
            double s = samples; next(s); samples = static_cast<int>(s);
        }
    }

    if (!csv) {
        std::printf(
            "usage:\n"
            "  profile --test\n"
            "  profile --csv --distance 400 --vmax 300 --amax 1200 "
            "--jmax 9000 > move.csv\n"
            "\nDistances are in whatever unit you like; velocity is unit/s,\n"
            "acceleration unit/s^2, jerk unit/s^3.\n");
        return 1;
    }
    return run_csv(distance, lim, samples);
}
