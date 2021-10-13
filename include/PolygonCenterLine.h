/*****************************************************************
To use this, you should define "PolygonCenterLine_Implementation"
before including this file in exactly one source file.
****************************************************************/
#ifndef POLYGON_CENTERLINE_H
#define POLYGON_CENTERLINE_H
#include <iostream>
#include <string>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>
#include <CGAL/intersections.h>
#include <CGAL/partition_2.h>
#include <boost/optional/optional_io.hpp>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Straight_skeleton_converter_2.h>

#include "CenterLineSolver.h"
#include "UcsPartitionSolver.h"
#include "PartitionSolver.h"
#include "convertKernel.h"
#include "CenterLineGeoJSON.h"

namespace CenterLine {
    struct PolygonCenterLine {
        using K = CGAL::Exact_predicates_exact_constructions_kernel;
        using FT = K::FT;
        //https://github.com/CGAL/cgal/issues/1873
        using Point_2 = K::Point_2;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Segment_2 = CGAL::Segment_2<K>;
        using Vector_2 = CGAL::Vector_2<K>;
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

		using Gps_traits_2 = CGAL::Gps_circle_segment_traits_2<K>;
		using Offset_polygon_2 = Gps_traits_2::Polygon_2;
		using Offset_polygon_with_holes_2 = Gps_traits_2::Polygon_with_holes_2;

        //using NType = CGAL::Gmpfr;
        //class InnerK : public CGAL::Filtered_kernel_adaptor<CGAL::Type_equality_wrapper<CGAL::Simple_cartesian<NType>::Base<InnerK>::Type, InnerK>,
        //#ifdef CGAL_NO_STATIC_FILTERS
        //            false >
        //#else
        //            true >
        //#endif
        //{};
        //using InnerK = CGAL::Filtered_kernel_adaptor<CGAL::Simple_cartesian<NType>>;
        using InnerK = CGAL::Epick;

        using InnerSs = CGAL::Straight_skeleton_2<InnerK>;
        using Ss = CGAL::Straight_skeleton_2<K>;

        // Gmpq precisions
        static const int _input_precision = 20; // С�����20�����ƣ�Լ1e-6

        //KernelConverter::NumberConverter<K::FT, InnerK::FT, _input_precision> nt_converter;
        KernelConverter::NumberConverter<K::FT, InnerK::FT> nt_converter;
        KernelConverter::NumberConverter<InnerK::FT, K::FT> res_nt_converter;
        KernelConverter::KernelConverter<InnerK, K, KernelConverter::NumberConverter<InnerK::FT, K::FT>> result_converter;
        //KernelConverter::KernelConverter<K, InnerK, KernelConverter::NumberConverter<K::FT, InnerK::FT, _input_precision>> kernel_converter;
        KernelConverter::KernelConverter<K, InnerK, KernelConverter::NumberConverter<K::FT, InnerK::FT>> kernel_converter;

        Polygon_with_holes_2 relative_poly;
        Vector_2 polygon_offset;
        std::vector<FT> point_contour_distance;
        // results:
        size_t seg_cnt_before_connect;
        boost::shared_ptr<Ss> skeleton;
        std::vector<Segment_2> relative_segments, relative_sub_segments;
        std::vector<Segment_2> _segments, _sub_segments;
        std::vector<std::pair<FT, FT>> segment_dis;
        std::vector<Offset_polygon_with_holes_2> rel_remainders;
        std::vector<Polygon_2> rel_convex_remainders;
        std::vector<Polygon_with_holes_2> rel_parts, centerline_parts, rel_ucs_parts;
        std::vector<Vector_2> corr_ucs;

        void init()
        {
            _segments.clear();
            _sub_segments.clear();
        }

    public:
        // input_precision: ��to_double(K)ת��ΪGmpfrʱ�ľ���; inner_precision: �ڲ�����ʱ��Ĭ�Ͼ���
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
        bool calcCenterLine(std::string geojson) { return calcCenterLine(geojson_to_poly(geojson)); }
        // calc with _segments
        bool calcPartition(FT R);
        // 弧度值差距不超过eps被认为是同一个ucs
        bool calcUcsPartition(double eps);

        std::string centerline_geojson() const { return segments_to_geojson(_segments); }
        std::string sub_centerline_geojson() const { return segments_to_geojson(_sub_segments); }
        std::string parts_geojson() const { return multipoly_to_geojson(rel_parts, centerline_parts, polygon_offset); }
        std::string ucs_parts_geojson() const { return multipoly_to_geojson(rel_ucs_parts, corr_ucs, polygon_offset); }

        static std::string segments_to_geojson(const std::vector<Segment_2> &segs){
            std::vector<std::pair<Point, Point>> res;
            for(auto &seg : segs){
                Point_2 src = seg.source(), dest = seg.target();
                Point p0(CGAL::to_double(src.x()), CGAL::to_double(src.y()));
                Point p1(CGAL::to_double(dest.x()), CGAL::to_double(dest.y()));
                res.push_back(std::make_pair(p0, p1));
            }
            return out2str(res);
        }
        static std::string multipoly_to_geojson(const std::vector<Polygon_with_holes_2> &polygons, const std::vector<Vector_2> &corr_ucs, Vector_2 offset = Vector_2(0, 0)){
            std::vector<Block> ucs_blocks;
            std::vector<Point> ucs;
            for(auto &polygon : polygons){
                Block block;
                block.coords.push_back(convert_points(polygon.outer_boundary(), offset));
                for(auto h_it = polygon.holes_begin(); h_it != polygon.holes_end(); ++h_it){
                    block.coords.push_back(convert_points(*h_it, offset));
                }
                ucs_blocks.push_back(block);
            }
            for(auto &co_ucs : corr_ucs){
                Point dir(CGAL::to_double(co_ucs.x().exact()), CGAL::to_double(co_ucs.y().exact()));
                ucs.push_back(dir);
            }
            return out2str(ucs_blocks, ucs);
        }
        static std::string multipoly_to_geojson(const std::vector<Polygon_with_holes_2> &polygons, const std::vector<Polygon_with_holes_2> &centerline_parts, Vector_2 offset = Vector_2(0, 0)){
            std::vector<Block> rect_blocks, centerline_blocks;
            for(auto &polygon : polygons){
                Block block;
                block.coords.push_back(convert_points(polygon.outer_boundary(), offset));
                for(auto h_it = polygon.holes_begin(); h_it != polygon.holes_end(); ++h_it){
                    block.coords.push_back(convert_points(*h_it, offset));
                }
                rect_blocks.push_back(block);
            }
            for(auto &polygon : centerline_parts){
                Block block;
                block.coords.push_back(convert_points(polygon.outer_boundary(), offset));
                for(auto h_it = polygon.holes_begin(); h_it != polygon.holes_end(); ++h_it){
                    block.coords.push_back(convert_points(*h_it, offset));
                }
                centerline_blocks.push_back(block);
            }
            return out2str(rect_blocks, centerline_blocks);
        }
        static Polygon_2 convert_poly(std::vector<Point> &points){
            std::vector<Point_2> pts;
            Point_2 tmp;
            for(auto &p : points){
                Point_2 pt(p.x, p.y);
                if(!pts.empty() && tmp == pt) continue;
                pts.push_back(pt);
                tmp = pt;
            }
            if(pts.front() == pts.back()) pts.pop_back();
            return Polygon_2(pts.begin(), pts.end());
        }
        static std::vector<Point> convert_points(const Polygon_2 &poly, Vector_2 offset){
            std::vector<Point> ans;
            auto it = poly.vertices_begin();
            if(it == poly.vertices_end()){
                throw("poly is empty");
            }
            for(;it != poly.vertices_end();++it)
                ans.emplace_back(CGAL::to_double(it->x() + offset.x()), CGAL::to_double(it->y() + offset.y()));
            ans.emplace_back(CGAL::to_double(poly.vertices_begin()->x() + offset.x()), CGAL::to_double(poly.vertices_begin()->y() + offset.y()));
            return ans;
        }
        static Polygon_with_holes_2 geojson_to_poly(std::string geojson){
            parseout res;
            parse_geojson(geojson, res);
            auto &data = res.data[0];
            Polygon_2 poly = convert_poly(data.coords[0]);
            if(poly.is_clockwise_oriented()) poly.reverse_orientation();
            std::vector<Polygon_2> holes;
            for(int i = 1;i < data.coords.size();++i){
                Polygon_2 hole = convert_poly(data.coords[i]);
                if(hole.is_counterclockwise_oriented()) hole.reverse_orientation();
                holes.push_back(hole);
            }
            if(holes.empty()) return Polygon_with_holes_2(poly);
            else return Polygon_with_holes_2(poly, holes.begin(), holes.end());
        }
    };
} // namespace CenterLine
#endif // POLYGON_CENTERLINE_H

#ifdef PolygonCenterLine_Implementation
namespace CenterLine {
    bool PolygonCenterLine::calcCenterLine(const Polygon_with_holes_2 &space)
    {
        init();

        //std::vector<CGAL::Segment_2<InnerK>> segments, sub_segments;
        //std::vector<CGAL::Polygon_2<InnerK>> new_PolyParts;

        // calculate relative coordinates ( relative_poly )
        Polygon_2 outer = space.outer_boundary();
        std::vector<Polygon_2> holes(space.holes_begin(), space.holes_end());
        polygon_offset = Vector_2(outer.left_vertex()->x(), outer.bottom_vertex()->y());
        for (auto it = outer.vertices_begin(); it != outer.vertices_end(); ++it) {
            std::cout << "p=" << (*it - polygon_offset) << std::endl;
            outer.set(it, *it - polygon_offset);
        }
        for (auto it = holes.begin(); it != holes.end();++it)
            for (auto i = it->vertices_begin(); i != it->vertices_end(); ++i) {
                std::cout << "p=" << (*i - polygon_offset) << std::endl;
                it->set(i, *i - polygon_offset);
            }

        relative_poly = Polygon_with_holes_2(outer, holes.begin(), holes.end());

        CGAL::Polygon_with_holes_2<InnerK> new_poly = kernel_converter.convert(relative_poly);
        // calculate with relative_poly
        CenterLineSolver::CenterLineSolver<InnerK, CGAL::Polygon_with_holes_2<InnerK>, CGAL::Polygon_2<InnerK>> solver;
        bool success = false;
        try {
            solver(new_poly);
            success = true;
        }
        catch (const char *str){
            std::cerr << "error = " << str << std::endl;
        }

        relative_segments.clear(); relative_sub_segments.clear(); segment_dis.clear();
        this->seg_cnt_before_connect = solver.seg_cnt_before_connect;
        std::transform(solver.res_segments.begin(), solver.res_segments.end(), std::back_inserter(relative_segments), result_converter);
        std::transform(solver.sub_segments.begin(), solver.sub_segments.end(), std::back_inserter(relative_sub_segments), result_converter);
        for(auto p : solver.res_seg_dis){
            segment_dis.emplace_back(res_nt_converter(p.first), res_nt_converter(p.second));
        }
        // for auto selecting R for PartitionSolver Test
        for(auto p_it : solver.point_data_pool) if(p_it->is_ans) {
            FT x = res_nt_converter(p_it->start_time());
            x = int(CGAL::sqrt(CGAL::to_double(x)) / 20 + 0.5) * 20;
            point_contour_distance.push_back(x);

            x = res_nt_converter(p_it->end_time());
            x = int(CGAL::sqrt(CGAL::to_double(x)) / 20 + 0.5) * 20;
            point_contour_distance.push_back(x);
        }
        std::sort(point_contour_distance.begin(), point_contour_distance.end());
        std::vector<FT>::iterator it = std::unique(point_contour_distance.begin(), point_contour_distance.end());
        point_contour_distance.erase(it, point_contour_distance.end());

        //relative_segments = solver.res_segments;
        //relative_sub_segments = solver.sub_segments;

        using ItemsCvt = CGAL::Straight_skeleton_items_converter_2<InnerSs, Ss>;
        CGAL::Straight_skeleton_converter_2<InnerSs, Ss, ItemsCvt> Ss_converter;
        skeleton = Ss_converter(*solver.skeleton);

        for (size_t i = 0; i < relative_segments.size(); ++i) {
            auto &seg = relative_segments[i];
            //CGAL::Segment_2<InnerK> new_seg(seg.source() + polygon_offset, seg.target() + polygon_offset);
            //_segments.push_back(result_converter(new_seg));
            _segments.emplace_back(seg.source() + polygon_offset, seg.target() + polygon_offset);
        }
        for (size_t i = 0; i < relative_sub_segments.size(); ++i) {
            auto &seg = relative_sub_segments[i];
            //CGAL::Segment_2<InnerK> new_seg(seg.source() + polygon_offset, seg.target() + polygon_offset);
            //_sub_segments.push_back(result_converter(new_seg));
            _sub_segments.emplace_back(seg.source() + polygon_offset, seg.target() + polygon_offset);
        }
        return success;
    } // bool PolygonCenterLine::calcCenterLine(const Polygon_with_holes_2 &space)

    bool PolygonCenterLine::calcUcsPartition(double eps){
        bool success = true;

        UcsPartitionSolver<K> solver(skeleton);
        success = solver(eps);
        if(success){
            for(size_t i = 0;i < solver.ucs_parts.size();++i){
                rel_ucs_parts.push_back(solver.ucs_parts[i]);
                corr_ucs.push_back(solver.ucs_direction[i]);
            }
        }
        return success;
    }

    bool PolygonCenterLine::calcPartition(FT R) {
        if(relative_segments.empty()) return false;
        PartitionSolver<K> solver;
        solver.space = relative_poly;
        solver.segments = relative_segments;
        solver.seg_dis = segment_dis;
        bool res = false;
        try {
            res = solver.solve(R);
        }
        catch (const char *str){
            std::cerr << "error = " << str << std::endl;
        }
        if(!res) return false;
        rel_remainders = solver.remainders;
        rel_convex_remainders = solver.convex_remainders;
        rel_parts = solver.parts;
        centerline_parts = solver.cencerline_parts;
        
        return true;
    } // bool PolygonCentrerLine::calcPartition(const Polygon_with_holes_2 &space);
} // namespace CenterLine
#undef PolygonCenterLine_Implementation
#endif // PolygonCenterLine_Implementation
