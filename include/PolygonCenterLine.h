/*****************************************************************
To use this, you should define "PolygonCenterLine_Implementation"
before including this file in exactly one source file.
****************************************************************/
#ifndef POLYGON_CENTERLINE_H
#define POLYGON_CENTERLINE_H
#include <iostream>
#include <string>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_with_holes_2.h>

#include "CenterLineSolver.h"
#include "convertKernel.h"
#include "parse.h"

namespace CenterLine {
    class PolygonCenterLine {
        using K = CGAL::Exact_predicates_exact_constructions_kernel;
        //https://github.com/CGAL/cgal/issues/1873
        using Point_2 = typename K::Point_2;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Segment_2 = CGAL::Segment_2<K>;
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

        using InnerK = CGAL::Simple_cartesian<CGAL::Gmpfr>;
        // Gmpq precisions
        static const int _input_precision = 20; // 小数点后20二进制，约1e-6

        KernelConverter::KernelConverter<K, InnerK, KernelConverter::NumberConverter<typename K::FT, InnerK::FT, _input_precision>> kernel_converter;
        //KernelConverter::KernelConverter<K, InnerK, KernelConverter::NumberConverter<typename K::FT, InnerK::FT>> kernel_converter;

        KernelConverter::KernelConverter<InnerK, K, KernelConverter::NumberConverter<InnerK::FT, typename K::FT>> result_converter;
        // results:
        std::vector<Segment_2> _segments, _sub_segments;

        void init()
        {
            _segments.clear();
            _sub_segments.clear();
        }

    public:
        // input_precision: 从to_double(K)转化为Gmpfr时的精度; inner_precision: 内部运算时的默认精度
        PolygonCenterLine()
        {
            CGAL::Gmpfr::set_default_precision(256);
            CGAL::set_pretty_mode(std::cout);
            CGAL::set_pretty_mode(std::cerr);
        }
        void setDefaultPrecision(int p) { CGAL::Gmpfr::set_default_precision(p); }
        const std::vector<Segment_2> &centerline() const { return _segments; }
        const std::vector<Segment_2> &sub_centerline() const { return _sub_segments; }
        bool calcCenterLine(const Polygon_with_holes_2 &space);
    };
} // namespace CenterLine
#endif // POLYGON_CENTERLINE_H

#ifdef PolygonCenterLine_Implementation
namespace CenterLine {
    bool PolygonCenterLine::calcCenterLine(const Polygon_with_holes_2 &space)
    {
        init();

        CGAL::Polygon_with_holes_2<InnerK> new_space = kernel_converter.convert(space);
        //std::vector<CGAL::Segment_2<InnerK>> segments, sub_segments;
        //std::vector<CGAL::Polygon_2<InnerK>> new_PolyParts;

        // calculate relative coordinates ( relative_poly )
        CGAL::Polygon_2<InnerK> outer = new_space.outer_boundary();
        std::vector<CGAL::Polygon_2<InnerK>> holes(new_space.holes_begin(), new_space.holes_end());
        CGAL::Vector_2<InnerK> polygon_offset(outer.left_vertex()->x(), outer.bottom_vertex()->y());
        for (auto it = outer.vertices_begin(); it != outer.vertices_end(); ++it) {
            std::cout << "p=" << (*it - polygon_offset) << std::endl;
            outer.set(it, *it - polygon_offset);
        }
        for (auto &poly : holes)
            for (auto i = poly.vertices_begin(); i != poly.vertices_end(); ++i) {
                std::cout << "p=" << (*i - polygon_offset) << std::endl;
                poly.set(i, *i - polygon_offset);
            }
        CGAL::Polygon_with_holes_2<InnerK> relative_poly(outer, holes.begin(), holes.end());

        // calculate with relative_poly
        CenterLineSolver::CenterLineSolver<InnerK, CGAL::Polygon_with_holes_2<InnerK>, CGAL::Polygon_2<InnerK>> solver;
        bool success = false;
        try {
            solver(relative_poly);
            success = true;
        }
        catch (const char *str){
            std::cerr << "error = " << str << std::endl;
        }

        const auto &segments = solver.res_segments;
        const auto &sub_segments = solver.sub_segments;

        for (size_t i = 0; i < segments.size(); ++i) {
            auto &seg = segments[i];
            CGAL::Segment_2<InnerK> new_seg(seg.source() + polygon_offset, seg.target() + polygon_offset);
            _segments.push_back(result_converter(new_seg));
        }
        for (size_t i = 0; i < sub_segments.size(); ++i) {
            auto &seg = sub_segments[i];
            CGAL::Segment_2<InnerK> new_seg(seg.source() + polygon_offset, seg.target() + polygon_offset);
            _sub_segments.push_back(result_converter(new_seg));
        }
        return success;
    } // void PolygonCenterLine::showCenterLine(double interval)

} // namespace CenterLine
#undef PolygonCenterLine_Implementation
#endif // PolygonCenterLine_Implementation
