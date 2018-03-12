//
// Created by kevin on 3/12/18.
//
#include "r_s_planner/dubins.h"

#include <queue>
#include <boost/math/constants/constants.hpp>

namespace {
    const double twopi = 2. * boost::math::constants::pi<double>();
    const double DUBINS_EPS = 1e-6;
    const double DUBINS_ZERO = -1e-7;

    inline double mod2pi(double x) {
        if (x < 0 && x > DUBINS_ZERO)
            return 0;
        double xm = x - twopi * floor(x / twopi);
        if (twopi - xm < .5 * DUBINS_EPS) xm = 0.;
        return xm;
    }

    DubinsStateSpace::DubinsPath dubinsLSL(double d, double alpha, double beta) {
        double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
        double tmp = 2. + d * d - 2. * (ca * cb + sa * sb - d * (sa - sb));
        if (tmp >= DUBINS_ZERO) {
            double theta = atan2(cb - ca, d + sa - sb);
            double t = mod2pi(-alpha + theta);
            double p = sqrt(std::max(tmp, 0.));
            double q = mod2pi(beta - theta);
            assert(fabs(p * cos(alpha + t) - sa + sb - d) < 2 * DUBINS_EPS);
            assert(fabs(p * sin(alpha + t) + ca - cb) < 2 * DUBINS_EPS);
            assert(mod2pi(alpha + t + q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
            return DubinsStateSpace::DubinsPath(DubinsStateSpace::dubinsPathType[0], t, p, q);
        }
        return DubinsStateSpace::DubinsPath();
    }

    DubinsStateSpace::DubinsPath dubinsRSR(double d, double alpha, double beta) {
        double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
        double tmp = 2. + d * d - 2. * (ca * cb + sa * sb - d * (sb - sa));
        if (tmp >= DUBINS_ZERO) {
            double theta = atan2(ca - cb, d - sa + sb);
            double t = mod2pi(alpha - theta);
            double p = sqrt(std::max(tmp, 0.));
            double q = mod2pi(-beta + theta);
            assert(fabs(p * cos(alpha - t) + sa - sb - d) < 2 * DUBINS_EPS);
            assert(fabs(p * sin(alpha - t) - ca + cb) < 2 * DUBINS_EPS);
            assert(mod2pi(alpha - t - q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
            return DubinsStateSpace::DubinsPath(DubinsStateSpace::dubinsPathType[1], t, p, q);
        }
        return DubinsStateSpace::DubinsPath();
    }

    DubinsStateSpace::DubinsPath dubinsRSL(double d, double alpha, double beta) {
        double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
        double tmp = d * d - 2. + 2. * (ca * cb + sa * sb - d * (sa + sb));
        if (tmp >= DUBINS_ZERO) {
            double p = sqrt(std::max(tmp, 0.));
            double theta = atan2(ca + cb, d - sa - sb) - atan2(2., p);
            double t = mod2pi(alpha - theta);
            double q = mod2pi(beta - theta);
            assert(fabs(p * cos(alpha - t) - 2. * sin(alpha - t) + sa + sb - d) < 2 * DUBINS_EPS);
            assert(fabs(p * sin(alpha - t) + 2. * cos(alpha - t) - ca - cb) < 2 * DUBINS_EPS);
            assert(mod2pi(alpha - t + q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
            return DubinsStateSpace::DubinsPath(DubinsStateSpace::dubinsPathType[2], t, p, q);
        }
        return DubinsStateSpace::DubinsPath();
    }

    DubinsStateSpace::DubinsPath dubinsLSR(double d, double alpha, double beta) {
        double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
        double tmp = -2. + d * d + 2. * (ca * cb + sa * sb + d * (sa + sb));
        if (tmp >= DUBINS_ZERO) {
            double p = sqrt(std::max(tmp, 0.));
            double theta = atan2(-ca - cb, d + sa + sb) - atan2(-2., p);
            double t = mod2pi(-alpha + theta);
            double q = mod2pi(-beta + theta);
            assert(fabs(p * cos(alpha + t) + 2. * sin(alpha + t) - sa - sb - d) < 2 * DUBINS_EPS);
            assert(fabs(p * sin(alpha + t) - 2. * cos(alpha + t) + ca + cb) < 2 * DUBINS_EPS);
            assert(mod2pi(alpha + t - q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
            return DubinsStateSpace::DubinsPath(DubinsStateSpace::dubinsPathType[3], t, p, q);
        }
        return DubinsStateSpace::DubinsPath();
    }

    DubinsStateSpace::DubinsPath dubinsRLR(double d, double alpha, double beta) {
        double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
        double tmp = .125 * (6. - d * d + 2. * (ca * cb + sa * sb + d * (sa - sb)));
        if (fabs(tmp) < 1.) {
            double p = twopi - acos(tmp);
            double theta = atan2(ca - cb, d - sa + sb);
            double t = mod2pi(alpha - theta + .5 * p);
            double q = mod2pi(alpha - beta - t + p);
            assert(fabs(2. * sin(alpha - t + p) - 2. * sin(alpha - t) - d + sa - sb) < 2 * DUBINS_EPS);
            assert(fabs(-2. * cos(alpha - t + p) + 2. * cos(alpha - t) - ca + cb) < 2 * DUBINS_EPS);
            assert(mod2pi(alpha - t + p - q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
            return DubinsStateSpace::DubinsPath(DubinsStateSpace::dubinsPathType[4], t, p, q);
        }
        return DubinsStateSpace::DubinsPath();
    }

    DubinsStateSpace::DubinsPath dubinsLRL(double d, double alpha, double beta) {
        double ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta);
        double tmp = .125 * (6. - d * d + 2. * (ca * cb + sa * sb - d * (sa - sb)));
        if (fabs(tmp) < 1.) {
            double p = twopi - acos(tmp);
            double theta = atan2(-ca + cb, d + sa - sb);
            double t = mod2pi(-alpha + theta + .5 * p);
            double q = mod2pi(beta - alpha - t + p);
            assert(fabs(-2. * sin(alpha + t - p) + 2. * sin(alpha + t) - d - sa + sb) < 2 * DUBINS_EPS);
            assert(fabs(2. * cos(alpha + t - p) - 2. * cos(alpha + t) + ca - cb) < 2 * DUBINS_EPS);
            assert(mod2pi(alpha + t - p + q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
            return DubinsStateSpace::DubinsPath(DubinsStateSpace::dubinsPathType[5], t, p, q);
        }
        return DubinsStateSpace::DubinsPath();
    }

    DubinsStateSpace::DubinsPath dubins(double d, double alpha, double beta) {
        if (d < DUBINS_EPS && fabs(alpha - beta) < DUBINS_EPS)
            return DubinsStateSpace::DubinsPath(DubinsStateSpace::dubinsPathType[0], 0, d, 0);

        DubinsStateSpace::DubinsPath path(dubinsLSL(d, alpha, beta)), tmp(dubinsRSR(d, alpha, beta));
        double len, minLength = path.length();

        if ((len = tmp.length()) < minLength) {
            minLength = len;
            path = tmp;
        }
        tmp = dubinsRSL(d, alpha, beta);
        if ((len = tmp.length()) < minLength) {
            minLength = len;
            path = tmp;
        }
        tmp = dubinsLSR(d, alpha, beta);
        if ((len = tmp.length()) < minLength) {
            minLength = len;
            path = tmp;
        }
        tmp = dubinsRLR(d, alpha, beta);
        if ((len = tmp.length()) < minLength) {
            minLength = len;
            path = tmp;
        }
        tmp = dubinsLRL(d, alpha, beta);
        if ((len = tmp.length()) < minLength)
            path = tmp;
        return path;
    }
}

const DubinsStateSpace::DubinsPathSegmentType DubinsStateSpace::dubinsPathType[6][3] = {
        {DUBINS_LEFT, DUBINS_STRAIGHT, DUBINS_LEFT},
        {DUBINS_RIGHT, DUBINS_STRAIGHT, DUBINS_RIGHT},
        {DUBINS_RIGHT, DUBINS_STRAIGHT, DUBINS_LEFT},
        {DUBINS_LEFT, DUBINS_STRAIGHT, DUBINS_RIGHT},
        {DUBINS_RIGHT, DUBINS_LEFT, DUBINS_RIGHT},
        {DUBINS_LEFT, DUBINS_RIGHT, DUBINS_LEFT}
};

double DubinsStateSpace::distance(double q0[3], double q1[3]) {
    return rho_ * dubins(q0, q1).length();
}

DubinsStateSpace::DubinsPath DubinsStateSpace::dubins(double q0[3], double q1[3]) {
    double x1 = q0[0], y1 = q0[1], th1 = q0[2];
    double x2 = q1[0], y2 = q1[1], th2 = q1[2];
    double dx = x2 - x1, dy = y2 - y1, d = sqrt(dx * dx + dy * dy) / rho_, th = atan2(dy, dx);
    double alpha = mod2pi(th1 - th), beta = mod2pi(th2 - th);
    return ::dubins(d, alpha, beta);
}

void DubinsStateSpace::interpolate(double q0[3], DubinsPath &path, double seg, double s[3]) {

    if (seg < 0.0) seg = 0.0;
    if (seg > path.length()) seg = path.length();

    double phi, v;

    s[0] = s[1] = 0.0;
    s[2] = q0[2];

    for (unsigned int i = 0; i < 3 && seg > 0; ++i) {
        v = std::min(seg, path.length_[i]);
        seg -= v;
        phi = s[2];
        switch (path.type_[i]) {
            case DUBINS_LEFT:
                s[0] += ( sin(phi+v) - sin(phi));
                s[1] += (-cos(phi+v) + cos(phi));
                s[2] = phi + v;
                break;
            case DUBINS_RIGHT:
                s[0] += (-sin(phi-v) + sin(phi));
                s[1] += ( cos(phi-v) - cos(phi));
                s[2] = phi - v;
                break;
            case DUBINS_STRAIGHT:
                s[0] += (v * cos(phi));
                s[1] += (v * sin(phi));
                break;
        }
    }

    s[0] = s[0] * rho_ + q0[0];
    s[1] = s[1] * rho_ + q0[1];

}

void DubinsStateSpace::sample(double q0[3], double q1[3], double step_size, double &length, std::vector<std::vector<double> > &points) {
    DubinsPath path = dubins(q0, q1);
    length = rho_ * path.length();

    for (double seg=0.0; seg<=length; seg+=step_size){
        double qnew[3] = {};
        interpolate(q0, path, seg/rho_, qnew);
        std::vector<double> v(qnew, qnew + sizeof qnew / sizeof qnew[0]);
        points.push_back(v);
    }
    return;
}










