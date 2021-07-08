#pragma once
#include "CenterLineSolver.h"
namespace CenterLineSolver{
    template <typename K, class Poly_with_holes, class Poly>
    struct CenterLineSolver<K, Poly_with_holes, Poly>::EdgeData {
        int corr_line;
        PointData *src_point, *dest_point;
        EdgeData *prev, *next;
        bool del;

        EdgeData(int line, PointData *src, PointData *dest)
            : corr_line(line), src_point(src), dest_point(dest), prev(this), next(this), del(false) {}
        bool is_active() const
        {
            return !del && src_point->end_loc == -1 && dest_point->end_loc == -1;
        }
        Vector_2 to_vector() const { return src_point->dest_vector; }
        Segment_2 new_loc(FT new_time)
        {
            return Segment_2(src_point->new_loc(new_time),
                             dest_point->new_loc(new_time));
        }
        Segment_2 new_loc(const std::vector<Line_2> &line_loc, const FT &new_time) const
        {
            if (is_debugging) {
                std::cout << src_point->path() << std::endl;
                std::cout << dest_point->path() << std::endl;
                std::cout << line_loc[corr_line] << std::endl;
            }
            auto a = CGAL::intersection(line_loc[corr_line], src_point->path().supporting_line());
            auto b = CGAL::intersection(line_loc[corr_line], dest_point->path().supporting_line());
            // assert(a  && b);
            if (!a || !b) {
                std::cerr << "no intersection :\n " << a << std::endl
                          << b << std::endl;
                throw("no intersection");
            }
            Point_2 *p = boost::get<Point_2>(&*a), *q = boost::get<Point_2>(&*b), l, r;
            if (!p || !q) {
                std::cout << line_loc[corr_line].to_vector() << " " << line_loc[prev->corr_line].to_vector() << " "
                          << line_loc[next->corr_line].to_vector() << std::endl;
                std::cerr << "intersection not Point_2 :\n"
                          << a << std::endl
                          << b << std::endl;
                throw("intersection not Point_2");
            }
            if (new_time == src_point->start_time())
                l = src_point->source();
            else
                l = *p;
            if (new_time == dest_point->start_time())
                r = dest_point->source();
            else
                r = *q;
            return Segment_2(l, r);
        }
    };
}