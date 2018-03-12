/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2010, Rice University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Rice University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Mark Moll */

#ifndef R_S_PLANNER_DUBINS_H
#define R_S_PLANNER_DUBINS_H

#include <boost/math/constants/constants.hpp>
#include <cassert>

class DubinsStateSpace {
public:
    /** \brief The Dubins path segment type */
    enum DubinsPathSegmentType {
        DUBINS_LEFT = 0, DUBINS_STRAIGHT = 1, DUBINS_RIGHT = 2
    };
    /** \brief Dubins path types */
    static const DubinsPathSegmentType dubinsPathType[6][3];

    /** \brief Complete description of a Dubins path */
    class DubinsPath {
    public:
        DubinsPath(const DubinsPathSegmentType *type = dubinsPathType[0], double t = 0.,
                   double p = std::numeric_limits<double>::max(), double q = 0.) : type_(type) {
            length_[0] = t;
            length_[1] = p;
            length_[2] = q;
            assert(t >= 0.);
            assert(p >= 0.);
            assert(q >= 0.);
        }

        double length() const {
            return length_[0] + length_[1] + length_[2];
        }

        /** Path segment types */
        const DubinsPathSegmentType *type_;
        /** Path segment lengths */
        double length_[3];
    };

    DubinsStateSpace(double turningRadius = 1.0) : rho_(turningRadius) {}

    void sample(double q0[3], double q1[3], double step_size, double &length, std::vector<std::vector<double> > &points);

    double distance(double q0[3], double q1[3]);


    /** \brief Return the shortest Dubins path from SE(2) state state1 to SE(2) state state2 */
    DubinsPath dubins(double q0[3], double q1[3]);

protected:

    /** \brief Turning radius */
    double rho_;

    void interpolate(double q0[3], DubinsPath &path, double seg, double s[3]);

};


#endif //R_S_PLANNER_DUBINS_H
