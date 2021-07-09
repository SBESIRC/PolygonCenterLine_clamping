#ifndef GET_CENTERLINE_H
#define GET_CENTERLINE_H
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <sstream>
#include <unordered_map>

#include <CGAL/Partition_traits_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/intersections.h>
#include <CGAL/partition_2.h>
#include <boost/optional/optional_io.hpp>

namespace CenterLineSolver {
    static bool is_debugging = false;
    struct PairHash {
        std::size_t operator()(const std::pair<int, int> &p) const
        {
            return std::hash<std::uint64_t>{}(std::uint64_t(p.first) << 32 | p.second);
        }
    };
    // 确保<a, b> = <b, a>
    inline std::pair<int, int> make_mono_pair(int a, int b)
    {
        return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
    }

    const double eps = 1. / (1 << 20), eps2 = eps * eps,
                 eps_limit = 0.25 / (1ll << 62);

    template <typename K, class Poly_with_holes, class Poly>
    struct CenterLineSolver {
        using FT = typename K::FT;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Point_2 = CGAL::Point_2<K>;
        using Vector_2 = CGAL::Vector_2<K>;
        using Segment_2 = CGAL::Segment_2<K>;
        using Line_2 = CGAL::Line_2<K>;
        using Ray_2 = CGAL::Ray_2<K>;
        struct PointData;
        struct EdgeData;
        struct Interval;
        struct Event;
        struct Location;
        enum class EventType : std::uint8_t { PARALLEL,
                                              NOT_PARALLEL,
                                              TRIANGLE };

        const double CornerCosineThreshold = -0.5;
        std::vector<Line_2> line_data_pool, line_loc;
        std::vector<Vector_2> line_speed;
        std::vector<PointData *> point_data_pool;
        std::vector<EdgeData *> edge_data_pool;
        std::priority_queue<Event> events;
        Point_2 top_right; // 右边界与上边界，用于判断交点是否落在区域外
        FT cur_time;       // events最新pop出的事件的发生时间. 初始为0
        //bool is_debugging;
        std::unordered_map<std::pair<int, int>, Line_2, PairHash> bisector_map;
        // results
        std::vector<Segment_2> res_segments, sub_segments;
        std::vector<Location> locations;

        static bool equal(const FT &a, const FT &b, FT epsilon = eps)
        {
            return a + epsilon > b && a < b + epsilon;
        }
        static bool less(const FT &a, const FT &b, FT epsilon = eps)
        {
            return a + epsilon <= b;
        }
        static FT Cosine(const Segment_2 &a, const Segment_2 &b)
        {
            auto x = a.to_vector(), y = b.to_vector();
            return (x * y) / CGAL::sqrt(x * x) / CGAL::sqrt(y * y);
        }
        static FT inner_product(const Vector_2 &a, const Vector_2 &b)
        {
            return a.x() * b.x() + a.y() * b.y();
        }
        static FT outer_product(const Vector_2 &a, const Vector_2 &b)
        {
            return a.x() * b.y() - a.y() * b.x();
        }
        static FT Cosine(const Vector_2 &a, const Vector_2 &b)
        {
            return inner_product(a, b) / CGAL::sqrt(a.squared_length() * b.squared_length());
        }
        static Vector_2 calc_real_speed(const Vector_2 &src_vector, const Vector_2 &dest_vector)
        {
            Vector_2 src = src_vector / CGAL::sqrt(src_vector.squared_length()),
                     dest = dest_vector / CGAL::sqrt(dest_vector.squared_length());
            Vector_2 v_bisector, &v_norm = src.perpendicular(CGAL::LEFT_TURN);
            if (inner_product(src_vector, dest_vector) > 0) {
                v_bisector = v_norm + dest.perpendicular(CGAL::LEFT_TURN);
            }
            else {
                v_bisector = dest - src;
                if (outer_product(src_vector, v_bisector) < 0)
                    v_bisector = -v_bisector;
            }
            FT tmp_inner_product = inner_product(v_norm, v_bisector);
            // speed的参考值为squared_length=1
            if (tmp_inner_product * tmp_inner_product < v_bisector.squared_length() * FT(eps_limit)) {
                std::cerr << "speed too large" << std::endl;
                assert(false); //"speed too large"
            }
            return v_bisector / tmp_inner_product;
        }
        static Vector_2 calc_speed(const Vector_2 &src_vector, const Vector_2 &dest_vector)
        {
            Vector_2 v_bisector;
            if (inner_product(src_vector, dest_vector) > 0) {
                Line_2 src_norm(CGAL::ORIGIN, src_vector.perpendicular(CGAL::LEFT_TURN));
                Line_2 dest_norm(CGAL::ORIGIN,
                                 dest_vector.perpendicular(CGAL::LEFT_TURN));
                // v_bisector = v_norm + dest.perpendicular(CGAL::LEFT_TURN);
                v_bisector = CGAL::bisector(src_norm, dest_norm).to_vector();
            }
            else {
                // v_bisector = dest - src;
                Line_2 src(CGAL::ORIGIN, -src_vector), dest(CGAL::ORIGIN, dest_vector);
                v_bisector = CGAL::bisector(src, dest).to_vector();
                if (outer_product(src_vector, v_bisector) < 0)
                    v_bisector = -v_bisector;
            }
            return v_bisector;

            // FT tmp_inner_product = inner_product(v_norm, v_bisector);
            //// speed的参考值为squared_length=1
            // if(tmp_inner_product * tmp_inner_product < v_bisector.squared_length() *
            // FT(eps_limit)) { 	assert(false); //"speed too large"
            //}
            // return v_bisector / tmp_inner_product;
        }
        static bool is_parallel(Vector_2 a, Vector_2 b)
        {
            FT absolute_eps = a.squared_length() * b.squared_length() * FT(eps2);
            FT tmp_outer_product = outer_product(a, b);
            return (tmp_outer_product * tmp_outer_product < absolute_eps);
        }
        static boost::optional<boost::variant<Point_2, Ray_2>> intersection(const Line_2 &line, const Ray_2 &ray, const Line_2 &bound)
        {
            Line_2 ray_line = ray.supporting_line();
            auto x = CGAL::intersection(line, ray_line);
            if (!x)
                return boost::none;
            Line_2 *inter_line;
            Point_2 *point;
            if ((inter_line = boost::get<Line_2>(&*x)))
                return ray;
            point = boost::get<Point_2>(&*x);
            FT point_loc = calc_loc_on_line(ray_line, *point);
            FT side = (bound.a() * point->x() + bound.b() * point->y() + bound.c()) / calc_line_scale(bound);
            if (is_debugging) {
                std::cout << ray << std::endl;
                std::cout << "ray = " << ray.to_vector() << std::endl;
                std::cout << "line = " << line.to_vector() << std::endl;
                std::cout << "x = " << x << std::endl;
                std::cout << "side = " << side << std::endl;
                std::cout << point_loc << std::endl;
            }
            if (less(side, 0))
                return boost::none;
            return *point;
        }
        static FT calc_line_scale(const Line_2 &line)
        {
            return CGAL::sqrt(line.a() * line.a() + line.b() * line.b());
        }
        // 沿着line的方向递增
        static FT calc_loc_on_line(const Line_2 &line, const Point_2 &point)
        {
            // std::cout << line << std::endl << point << std::endl;
            return (line.b() * point.x() - line.a() * point.y()) / calc_line_scale(line);
        }
        // 在bisector上的位置区间(按先后顺序从先到后)
        bool calc_interval_on_line(const Line_2 &bisector, const EdgeData &edge, const EdgeData &op_edge,
                                   Interval &interval, bool parallel = false);
        void append_collision_events(EdgeData *edge);
        void start_debug() { is_debugging = true; }
        void end_debug() { is_debugging = false; }
        void append_line_data(const Polygon_2 &poly, size_t &line_it);
        // type=0: base的source至edge的source; type=1: edge的target至base的target
        // TODO: 点的终止 (对于另一端相交于同一点的情况，考虑两边相邻的情况)
        void cut_half_edge(EdgeData *base, EdgeData *edge, bool type, int location,
                           std::vector<PointData *> &point_pool,
                           std::vector<EdgeData *> &edge_pool);
        void cut_half_point(EdgeData *prev, EdgeData *next, int location,
                            std::vector<PointData *> &point_pool,
                            std::vector<EdgeData *> &edge_pool);
        void cut_edge(EdgeData *base, EdgeData *prev, EdgeData *next, int location,
                      std::vector<PointData *> &point_pool,
                      std::vector<EdgeData *> &edge_pool);
        bool in_box(Point_2 point) { return point.x() >= 0 && point.x() <= top_right.x() && point.y() >= 0 && point.y() <= top_right.y(); }
        void process_event(const Event &event);
        void connect_segments();
        void operator()(const Poly_with_holes &polygon);
        CenterLineSolver() { is_debugging = false; }
    }; // CenterLineSolver

    // PolygonCenterLine_Implementation

    template <typename K, class Poly_with_holes, class Poly>
    struct CenterLineSolver<K, Poly_with_holes, Poly>::Location {
        Point_2 point;
        FT time;
        std::vector<std::pair<PointData *, bool>> branches;
        Location(Point_2 p, FT t) : point(p), time(t) {}
    };

    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::cut_edge(
        EdgeData *base, EdgeData *prev, EdgeData *next, int location,
        std::vector<PointData *> &point_pool, std::vector<EdgeData *> &edge_pool)
    {
        if (is_debugging) {
            std::cout << "base= " << base->new_loc(line_loc, this->cur_time) << std::endl;
            std::cout << "prev= " << prev->new_loc(line_loc, this->cur_time) << std::endl;
            std::cout << "next= " << next->new_loc(line_loc, this->cur_time) << std::endl;
            std::cout << "location= " << location << std::endl;
        }
        cut_half_edge(base, prev, 1, location, point_pool, edge_pool);
        cut_half_edge(base, next, 0, location, point_pool, edge_pool);
    }

    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::operator()(
        const Poly_with_holes &polygon)
    {
        top_right = Point_2(polygon.outer_boundary().right_vertex()->x(), polygon.outer_boundary().top_vertex()->y());

        size_t tot_edges = polygon.outer_boundary().size(); // the size of line_data_pool.
        for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
            tot_edges += it->size();

        // get the corresbonding line of every edges;
        // build point_data_pool and edge_data_pool at the same time;
        line_data_pool.reserve(tot_edges);
        size_t line_it = 0;
        line_data_pool.clear();
        line_loc.clear();
        line_speed.clear();
        append_line_data(polygon.outer_boundary(), line_it);
        for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
            append_line_data(*it, line_it);

        // for each pair of non-adjacent edges, calculate the bisector of them and
        // when the endpoints of both edges intersect with the bisector. store the
        // edge pair indexed with the start time of the overlapping interval if it
        // exists.
        //				-> as a function: "for a new edge with each edge in
        //the Iterable list"
        // TODO: 预处理速度，减少求速度的消耗

        // set an event for each start time of the overlapping interval of an edge
        // pair.
        while (!events.empty())
            events.pop();
        cur_time = 0;
        for (auto e_it : edge_data_pool) {
            append_collision_events(e_it);
        }

        // fetch an event at a time in chronological order:
        // if(edge0.is_active() && edge1.is_active())
        // > case edge vs edge: split the edges (or maybe endpoints)
        //		(平行相遇:直接终止原来的端点轨迹，然后输出相交的片段，再把多余的片段截取出来，重建多边形)
        // > case edge vs point: split the edge; split the point
        //		(相遇位置位于某个线段内部:
        //把这个线段拆成两半，然后让点终止，分裂成两个点分别与两边的edge片段相，重建多边形)
        // > case point vs point: modify the linking relation of points and edges;
        //		(相遇位置为两个线段的顶点:
        //终止两个顶点轨迹，重建新的顶点，重新连接多边形)
        // change the speed of new points
        // for new edges: call "append_collision_events"

        const char *error = NULL;
        try {
            while (!events.empty()) {
                Event event = events.top();
                events.pop();
                if (!event.is_active())
                    continue;

                this->process_event(event);
                if(is_debugging){
                    std::cout << "segments:=====================================" << std::endl;
                    for (size_t i = 0; i < point_data_pool.size(); ++i) std::cout << i << ": " << point_data_pool[i]->source() << " => " << point_data_pool[i]->target() << std::endl;
                }
            }
        }
        catch (const char *str) {
            std::cerr << "unexpectedly exit: " << str << std::endl;
            error = str;
        }

        for(int i = 0;i < (int) point_data_pool.size(); ++i) point_data_pool[i]->point_id = i;

        std::cout << "segments:=====================================" << std::endl;
        for (size_t i = 0; i < point_data_pool.size(); ++i)
            std::cout << i << ": " << point_data_pool[i]->source() << " => " << point_data_pool[i]->target() << std::endl;

        std::unordered_map<PointData *, int> pd_id;
        pd_id[(PointData *)0] = -1;
        for (size_t i = 0; i < point_data_pool.size(); ++i)
            pd_id[point_data_pool[i]] = i;
        for (size_t i = 0; i < point_data_pool.size(); ++i) {
            std::cout << i << " : " << std::endl;
            PointData *pd = point_data_pool[i];
            std::cout << pd->src_vector << " " << pd->dest_vector << std::endl;
            std::cout << "<< " << locations[pd->start_loc].branches.size() << std::endl;
            //for (auto p : locations[pd->start_loc].branches)
            //    std::cout << "<" << pd_id[p.first] << ", " << p.second << ") ";
            //std::cout << std::endl
            //          << ">> ";
            //for (auto p : locations[pd->end_loc].branches)
            //    std::cout << "<" << pd_id[p.first] << ", " << p.second << ") ";
            //std::cout << std::endl;
        }

        // 补全连线
        try{
            this->connect_segments();
        }
        catch(const char *str){
            std::cerr << "connect segments error: " << str << std::endl;
        }

        int counter = 0;
        for (PointData *it : point_data_pool) {
            ++counter;
            if (it->end_loc == -1) {
                std::cout << "point " << counter - 1 << "not ended" << std::endl;
                continue;
            }
            if (it->is_ans) {
                this->res_segments.emplace_back(it->source(), it->target());
            }
            else
                this->sub_segments.emplace_back(it->source(), it->target());
        }
        if (error)
            throw(error);
    }
} // namespace CenterLineSolver

#include "impl/EdgeData.h"
#include "impl/Event.h"
#include "impl/Interval.h"
#include "impl/PointData.h"
#include "impl/append_collision_events.h"
#include "impl/append_line_data.h"
#include "impl/calc_interval_on_line.h"
#include "impl/connect_segments.h"
#include "impl/cut_half_edge.h"
#include "impl/cut_half_point.h"
#include "impl/process_event.h"

#endif // GET_CENTERLINE_H