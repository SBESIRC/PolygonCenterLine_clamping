#pragma once
#include "CenterLineSolver.h"
namespace CenterLineSolver{
    template <typename K, class Poly_with_holes, class Poly>
    struct CenterLineSolver<K, Poly_with_holes, Poly>::Interval {
        FT begin, end;
        Point_2 begin_point, end_point;
        bool is_inf, is_empty;
        // -1/-2: always endpoint; 0: parallel or intersect; 1/2: collision starts at
        // source/target
        int is_endpoint;
        std::string to_string()
        {
            std::stringstream ss;
            CGAL::set_pretty_mode(ss);
            ss << "------------------" << std::endl;
            ss << "begin = " << begin << std::endl;
            ss << "end = " << end << std::endl;
            ss << "begin_point = " << begin_point << std::endl;
            ss << "end_point = " << end_point << std::endl;
            ss << "is_inf=" << is_inf << std::endl;
            ss << "is_empty=" << is_empty << std::endl;
            ss << "is_endpoint=" << is_endpoint << std::endl;
            ss << "------------------" << std::endl;
            return ss.str();
        }
        Interval()
            : is_inf(false), is_empty(true), is_endpoint(0)
        {
        }
        void init() { is_empty = true; }
        void set(FT &l, Point_2 &L, int end = 0)
        {
            init();
            begin = l;
            begin_point = L;
            is_inf = true;
            is_empty = false;
            is_endpoint = end;
        }
        void set(FT &l, FT &r, Point_2 &L, Point_2 &R, int is_end = 0)
        {
            init();
            begin = l, end = r;
            begin_point = L, end_point = R;
            is_empty = !less(l, r); //(l >= r);
            is_endpoint = is_end;
        }
        Interval operator+(const Interval &b) const
        {
            Interval ans;
            if (is_empty || b.is_empty) {
                ans.is_empty = true;
                return ans;
            }
            ans.is_empty = false;

            if (is_endpoint < 0 && b.is_endpoint < 0) {
                assert(false); // 相邻边，不考虑
            }
            else if (is_endpoint < 0 || b.is_endpoint < 0) {
                // 一方有一端恒在平分线上，另一方也是相同端与之相遇
                if (is_endpoint && is_endpoint + b.is_endpoint == 0) {
                    ans.is_empty = true;
                    return ans;
                }
            }
            // 点与点相撞，且头对头/脚对脚
            else if (is_endpoint != 0 && is_endpoint == b.is_endpoint && equal(begin, b.begin)) {
                ans.is_empty = true;
                return ans;
            }

            if (begin < b.begin)
                ans.begin = b.begin, ans.begin_point = b.begin_point;
            else
                ans.begin = begin, ans.begin_point = begin_point;

            if (is_inf && b.is_inf)
                ans.is_inf = true;
            else if (is_inf) {
                // if(b.end <= ans.begin)
                if (!less(ans.begin, b.end))
                    ans.is_empty = true;
                else
                    ans.end = b.end, ans.end_point = b.end_point;
            }
            else if (b.is_inf) {
                // if(end <= ans.begin)
                if (!less(ans.begin, end))
                    ans.is_empty = true;
                else
                    ans.end = end, ans.end_point = end_point;
            }
            else {
                // if(end <= ans.begin || b.end <= ans.begin)
                if (!less(ans.begin, end) || !less(ans.begin, b.end))
                    ans.is_empty = true;
                else if (end < b.end)
                    ans.end = end, ans.end_point = end_point;
                else
                    ans.end = b.end, ans.end_point = b.end_point;
            }
            // if(ans.is_inf || ans.begin < ans.end)
            if (!ans.is_inf && !less(ans.begin, ans.end))
                ans.is_empty = true;
            return ans;
        }
    };
}