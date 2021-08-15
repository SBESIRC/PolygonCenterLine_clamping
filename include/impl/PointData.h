#pragma once
#include "CenterLineSolver.h"
namespace CenterLineSolver {
    template <typename K, class Poly_with_holes, class Poly>
    struct CenterLineSolver<K, Poly_with_holes, Poly>::PointData {
        std::vector<Location> &locations;
        int start_loc, end_loc; // start time & location; end time & location;
        int point_id;
        Vector_2 src_vector, dest_vector; // src_vector; dest_vector; (unit vector)
        //PointData *branches[2]; // branches[2]: PointData. Points created when splitting
        std::vector<std::pair<PointData *, bool>> branches[2]; // 0: start; 1: end
        PointData *prev, *next;
        bool is_ans;

        PointData(std::vector<Location> &loc, int _start_loc,
                  const Vector_2 &src, const Vector_2 &dest)
            : locations(loc), start_loc(_start_loc), end_loc(-1), src_vector(src), dest_vector(dest), prev(this), next(this)
        {
            std::cout << "start_loc = " << _start_loc << std::endl;
            locations[_start_loc].branches.emplace_back(this, 0);
            FT tmp_inner_product = inner_product(src_vector, dest_vector);
            FT absolute_scale = src_vector.squared_length() * dest_vector.squared_length();
            is_ans = (outer_product(src_vector, dest_vector) >= 0 && tmp_inner_product < 0 &&
                      tmp_inner_product * tmp_inner_product > FT(0.25 - eps) * absolute_scale);
        }
        PointData(std::vector<Location> &loc, int _start_loc, int _end_loc, int id)
            : locations(loc), start_loc(_start_loc), end_loc(_end_loc), is_ans(true), prev(this), next(this), point_id(id) {}
        void set_end(int loc)
        {
            end_loc = loc;
            locations[loc].branches.emplace_back(this, 1);
        }
        Point_2 source() const { return locations[start_loc].point; }
        Point_2 target() const
        {
            if (end_loc != -1)
                return locations[end_loc].point;
            else
                return Point_2(FT(0) / FT(0), FT(0) / FT(0));
        }
        Point_2 point(bool port) { return port ? target() : source(); }
        int location(bool port) { return port ? end_loc : start_loc; }
        FT start_time() const { return locations[start_loc].time; }
        FT end_time() const {
            if(end_loc != -1) return locations[end_loc].time;
            else return FT(1) / FT(0);
        }
        Vector_2 speed() const
        {
            // TODO: 处理极端情况(夹角接近+-180°)
            return calc_speed(src_vector, dest_vector);
        }
        Ray_2 path() const
        {
            return Ray_2(source(), speed());
        }
        Point_2 new_loc(FT time) const
        {
            return source() + calc_real_speed(src_vector, dest_vector) * (CGAL::sqrt(time) - CGAL::sqrt(start_time()));
        }
    };
}