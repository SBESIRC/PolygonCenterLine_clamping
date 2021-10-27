#ifndef GET_CENTERLINE_H
#define GET_CENTERLINE_H
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include <CGAL/Partition_traits_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/intersections.h>
#include <CGAL/partition_2.h>
#include <boost/optional/optional_io.hpp>
#include <CGAL/Straight_skeleton_builder_2.h>

#include "CenterLineContext.h"

namespace CenterLineSolver {
    static bool is_debugging = false;
    struct PairHash {
        std::size_t operator()(const std::pair<int, int> &p) const
        {
            return std::hash<std::uint64_t>{}(std::uint64_t(p.first) << 32 | p.second);
        }
    };
    // ȷ��<a, b> = <b, a>
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
        using Ss = CGAL::Straight_skeleton_2<K>;
        using Halfedge_iterator = typename Ss::Halfedge_iterator;
        using Halfedge_handle = typename Ss::Halfedge_handle;
        using Vertex_handle = typename Ss::Vertex_handle;
        using Vertex_iterator = typename Ss::Vertex_iterator;
        using SsBuilderTraits = CGAL::Straight_skeleton_builder_traits_2<K>;
        using SsBuilder = CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss>;

        struct PointData;
        struct EdgeData;
        struct Interval;
        struct Event;
        struct Location;
        enum class EventType : std::uint8_t { PARALLEL,
                                              NOT_PARALLEL,
                                              TRIANGLE };

        const CenterLine::Context _context;    
        const double CornerCosineThreshold = -0.5;
        const Poly_with_holes *origin_space;
        std::vector<PointData *> point_data_pool;
        std::priority_queue<Event> events;
        Point_2 top_right; // �ұ߽����ϱ߽磬�����жϽ����Ƿ�����������
        FT cur_time;       // events����pop�����¼��ķ���ʱ��. ��ʼΪ0
        //bool is_debugging;
        std::unordered_map<std::pair<int, int>, Line_2, PairHash> bisector_map;
        // results
        boost::shared_ptr<Ss> skeleton;
        size_t seg_cnt_before_connect;
        std::vector<Segment_2> res_segments, sub_segments;
        std::vector<std::pair<FT, FT>> res_seg_dis;
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
        static FT Sine(const Vector_2 &a, const Vector_2 &b)
        {
            return outer_product(a, b) / CGAL::sqrt(a.squared_length() * b.squared_length());
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
            // speed�Ĳο�ֵΪsquared_length=1
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
            //// speed�Ĳο�ֵΪsquared_length=1
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
        // ����line�ķ������
        static FT calc_loc_on_line(const Line_2 &line, const Point_2 &point)
        {
            // std::cout << line << std::endl << point << std::endl;
            return (line.b() * point.x() - line.a() * point.y()) / calc_line_scale(line);
        }
        // ��bisector�ϵ�λ������(���Ⱥ�˳����ȵ���)
        bool calc_interval_on_line(const Line_2 &bisector, const EdgeData &edge, const EdgeData &op_edge,
                                   Interval &interval, bool parallel = false);
        void append_collision_events(EdgeData *edge);
        void start_debug() { is_debugging = true; }
        void end_debug() { is_debugging = false; }
        void append_line_data(const Polygon_2 &poly, size_t &line_it);
        // type=0: base��source��edge��source; type=1: edge��target��base��target
        // TODO: �����ֹ (������һ���ཻ��ͬһ�������������������ڵ����)
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
        CenterLineSolver(const CenterLine::Context &context) : _context(context) { is_debugging = false; }
        ~CenterLineSolver(){
            for(auto ptr : point_data_pool) delete ptr;
        }
    }; // CenterLineSolver

    // PolygonCenterLine_Implementation

    template <typename K, class Poly_with_holes, class Poly>
    struct CenterLineSolver<K, Poly_with_holes, Poly>::Location {
        Point_2 point;
        FT time;
        std::vector<std::pair<PointData *, bool>> branches;
        Location(){}
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
    inline int ufs_find(std::vector<int> &f, size_t x){
        if(f[x] == x) return x;
        return f[x] = ufs_find(f, f[x]);
    }
    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::operator()(
        const Poly_with_holes &polygon)
    {
        origin_space = &polygon;
        top_right = Point_2(polygon.outer_boundary().right_vertex()->x(), polygon.outer_boundary().top_vertex()->y());

        size_t tot_edges = polygon.outer_boundary().size(); // the size of line_data_pool.
        for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
            tot_edges += it->size();

        // get the corresbonding line of every edges;
        // build point_data_pool and edge_data_pool at the same time;

        SsBuilder ssb;
        ssb.enter_contour(polygon.outer_boundary().vertices_begin(), polygon.outer_boundary().vertices_end());
        for(auto hole = polygon.holes_begin();hole != polygon.holes_end();++hole)
            ssb.enter_contour(hole->vertices_begin(), hole->vertices_end());
        try{
            skeleton = ssb.construct_skeleton();
        }
        catch(...){
            std::cerr << "construct_skeleton exception ...\n";
        }
        if(!skeleton) throw("skeleton was not correctly constructed");

        std::unordered_map<int, int> index;
        std::unordered_set<std::pair<int, int>, PairHash> edge_set;
        std::vector<int> f;
        std::cout << "size = " << skeleton->size_of_vertices() << std::endl;
        size_t cnt = 0;
        for(Vertex_iterator it = skeleton->vertices_begin();it != skeleton->vertices_end();++it){
            index[it->id()] = cnt;
            f.push_back(cnt);
            ++cnt;
        }
        for(Halfedge_iterator it = skeleton->halfedges_begin();it != skeleton->halfedges_end();++it) if(it->is_bisector()) {
            Vertex_handle from = it->opposite()->vertex(), to = it->vertex();
            if(equal(from->point().x(), to->point().x(), 1e-4) && equal(from->point().y(), to->point().y(), 1e-4)){
                int a = ufs_find(f, index[from->id()]), b = ufs_find(f, index[to->id()]);
                f[a] = b;
            }
        }
        std::vector<int> rk(cnt, -1);
        int rk_it = 0;
        for(int i = 0;i < cnt;++i){
            f[i] = ufs_find(f, i);
            if(rk[f[i]] == -1) rk[f[i]] = rk_it++;
            f[i] = rk[f[i]];
        }

        locations.resize(rk_it);
        for(Vertex_iterator it = skeleton->vertices_begin();it != skeleton->vertices_end();++it){
            int id = f[index[it->id()]];
            locations[id].point = it->point();
            //locations[cnt].time = it->time();
            locations[id].time = it->time() * it->time();
        }
        for(Halfedge_iterator it = skeleton->halfedges_begin();it != skeleton->halfedges_end();++it) if(it->is_bisector()) {
            Vertex_handle from = it->opposite()->vertex(), to = it->vertex();
            int a = f[index[from->id()]], b = f[index[to->id()]];
            Halfedge_handle l = it->defining_contour_edge(), r = it->opposite()->defining_contour_edge();
            //if(index[from->id()] < index[to->id()]){
            if(a < b){
                if(from->time() > to->time()){
                    swap(a, b);
                    swap(from, to);
                    swap(l, r);
                }
                if(edge_set.count(std::make_pair(a, b))) continue;
                edge_set.insert(std::make_pair(a, b));
                int e_id = point_data_pool.size();
                //PointData *edge = new PointData(locations, index[from->id()], index[to->id()], e_id);
                PointData *edge = new PointData(locations, a, b, e_id);
                point_data_pool.push_back(edge);
                //std::cout << "edge = " << from->point() << " " << to->point() << std::endl;
                Vector_2 v_l = l->vertex()->point() - l->opposite()->vertex()->point();
                Vector_2 v_r = r->vertex()->point() - r->opposite()->vertex()->point();
                std::cout << "l = " << l->opposite()->vertex()->point() << " " << l->vertex()->point() << std::endl;
                std::cout << "r = " << r->opposite()->vertex()->point() << " " << r->vertex()->point() << std::endl;
                std::cout << "v_l = " << v_l << "\nv_r = " << v_r << std::endl;
                FT inner = v_l.x() * v_r.x() + v_l.y() * v_r.y();
                FT outer = v_l.x() * v_r.y() - v_l.y() * v_r.x();
                FT absolute_scale = v_l.squared_length() * v_r.squared_length();
                edge->is_ans = (inner < 0 && outer >= 0 && inner * inner * 4 >= absolute_scale);
                std::cout << edge->is_ans << std::endl;
                //locations[index[from->id()]].branches.emplace_back(edge, 0);
                //locations[index[to->id()]].branches.emplace_back(edge, 1);
                locations[a].branches.emplace_back(edge, 0);
                locations[b].branches.emplace_back(edge, 1);
            }
        }

        size_t sep = point_data_pool.size();
        try{
            this->connect_segments();
        }
        catch(const char *str){
            std::cerr << "connect segments error: " << str << std::endl;
        }

        int counter = 0;
        for (PointData *it : point_data_pool) {
            ++counter;
            if(counter - 1 == sep) this->seg_cnt_before_connect = this->res_segments.size();
            if (it->end_loc == -1) {
                std::cout << "point " << counter - 1 << "not ended" << std::endl;
                continue;
            }
            if (it->is_ans) {
                if(it->end_loc != -1){
                    this->res_segments.emplace_back(it->source(), it->target());
                    this->res_seg_dis.emplace_back(it->start_time(), it->end_time());
                }
                else{
                    std::cerr << "point " << it->point_id << " not ended" << std::endl;
                }
            }
            else
                this->sub_segments.emplace_back(it->source(), it->target());
        }
        if(point_data_pool.size() == sep) this->seg_cnt_before_connect = this->res_segments.size();
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