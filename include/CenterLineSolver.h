#ifndef GET_CENTERLINE_H
#define GET_CENTERLINE_H
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
        enum class EventType : std::uint8_t { PARALLEL,
                                              NOT_PARALLEL,
                                              TRIANGLE };

        const double CornerCosineThreshold = -0.5;
        std::vector<Line_2> line_data_pool, line_loc;
        std::vector<Vector_2> line_speed;
        std::vector<PointData *> point_data_pool;
        std::vector<EdgeData *> edge_data_pool;
        std::priority_queue<Event> events;
        FT cur_time; // events最新pop出的事件的发生时间. 初始为0
        bool is_debugging;
        std::unordered_map<std::pair<int, int>, Line_2, PairHash> bisector_map;
        // results
        std::vector<Segment_2> res_segments, sub_segments;

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
        FT calc_line_scale(const Line_2 &line)
        {
            return CGAL::sqrt(line.a() * line.a() + line.b() * line.b());
        }
        // 沿着line的方向递增
        FT calc_loc_on_line(const Line_2 &line, const Point_2 &point)
        {
            // std::cout << line << std::endl << point << std::endl;
            return (line.b() * point.x() - line.a() * point.y()) / calc_line_scale(line);
        }
        // 在bisector上的位置区间(按先后顺序从先到后)
        bool calc_interval_on_line(const Line_2 &bisector, const EdgeData &edge,
                                   Interval &interval, bool parallel = false);
        void append_collision_events(EdgeData *edge);
        void start_debug() { is_debugging = true; }
        void end_debug() { is_debugging = false; }
        void append_line_data(const Polygon_2 &poly, size_t &line_it);
        // type=0: base的source至edge的source; type=1: edge的target至base的target
        // TODO: 点的终止 (对于另一端相交于同一点的情况，考虑两边相邻的情况)
        void CutHalfEdge(EdgeData *base, EdgeData *edge, bool type, Point_2 location,
                         std::vector<PointData *> &point_pool,
                         std::vector<EdgeData *> &edge_pool);
        void CutHalfPoint(EdgeData *prev, EdgeData *next, Point_2 location,
                          std::vector<PointData *> &point_pool,
                          std::vector<EdgeData *> &edge_pool);
        void CutEdge(EdgeData *base, EdgeData *prev, EdgeData *next, Point_2 location,
                     std::vector<PointData *> &point_pool,
                     std::vector<EdgeData *> &edge_pool);
        void operator()(const Poly_with_holes &polygon);
        CenterLineSolver() { is_debugging = false; }
    }; // CenterLineSolver

    // PolygonCenterLine_Implementation
    template <typename K, class Poly_with_holes, class Poly>
    struct CenterLineSolver<K, Poly_with_holes, Poly>::PointData {
        Point_2 start_loc, end_loc;
        FT start_time, end_time;          // start time & location; end time & location;
        Vector_2 src_vector, dest_vector; // src_vector; dest_vector; (unit vector)
        PointData *
            branches[2]; // branches[2]: PointData. Points created when splitting
        PointData *prev, *next;

        PointData(const Point_2 &_start_loc, const FT &_start_time,
                  const Vector_2 &src, const Vector_2 &dest)
            : start_loc(_start_loc), start_time(_start_time), end_time(0), src_vector(src), dest_vector(dest), branches{0}, prev(this), next(this) {}
        Vector_2 speed() const
        {
            // TODO: 处理极端情况(夹角接近+-180°)
            return calc_speed(src_vector, dest_vector);
        }
        Ray_2 path() const
        {
            return Ray_2(start_loc, speed());
        }
        Point_2 new_loc(FT time) const
        {
            return start_loc + calc_real_speed(src_vector, dest_vector) * (CGAL::sqrt(time) - CGAL::sqrt(start_time));
        }
    };
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
            return !del && src_point->end_time == FT(0) && dest_point->end_time == FT(0);
        }
        Segment_2 new_loc(FT new_time)
        {
            return Segment_2(src_point->new_loc(new_time),
                             dest_point->new_loc(new_time));
        }
        Segment_2 new_loc(const std::vector<Line_2> &line_loc) const
        {
            //std::cout << src_point->path() << std::endl;
            //std::cout << dest_point->path() << std::endl;
            //std::cout << line_loc[corr_line] << std::endl;
            auto a = CGAL::intersection(line_loc[corr_line], src_point->path());
            auto b = CGAL::intersection(line_loc[corr_line], dest_point->path());
            // assert(a  && b);
            if (!a || !b) {
                std::cerr << "no intersection :\n " << a << std::endl << b << std::endl;
                exit(1);
            }
            Point_2 *p = boost::get<Point_2>(&*a), *q = boost::get<Point_2>(&*b);
            if (!p || !q) {
                std::cerr << "intersection not Point_2 :\n" << p << std::endl << q << std::endl;
                exit(1);
            }
            return Segment_2(*p, *q);
        }
    };

    template <typename K, class Poly_with_holes, class Poly>
    struct CenterLineSolver<K, Poly_with_holes, Poly>::Event {
        FT time;
        EdgeData *a, *b;
        EventType type;
        Event(FT t, EdgeData *e0, EdgeData *e1, EventType tp)
            : time(t), a(e0), b(e1), type(tp)
        {
            assert(e0);
            assert(e1);
        }
        bool is_active() const { return a->is_active() && b->is_active(); }
        bool operator<(const Event &b) const { return time > b.time; }
    };

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

    template <typename K, class Poly_with_holes, class Poly>
    inline bool CenterLineSolver<K, Poly_with_holes, Poly>::calc_interval_on_line(
        const Line_2 &bisector, const EdgeData &edge, Interval &interval,
        bool parallel)
    {
        PointData *src = edge.src_point, *dest = edge.dest_point;
        Line_2 corr_line = line_data_pool[edge.corr_line];
        Vector_2 edge_vector = corr_line.to_vector(),
                 bi_vector = bisector.to_vector();
        FT tmp_inner_product = inner_product(edge_vector, bi_vector);
        FT tmp_outer_product = outer_product(edge_vector, bi_vector);
        // product(a, b) / sqrt((a*a) * (b*b)) ====== eps
        FT absolute_eps = edge_vector.squared_length() * bi_vector.squared_length() * FT(eps2);
        assert(tmp_outer_product >= 0);
        Point_2 *tmp_p = NULL, *tmp_q = NULL;
        Ray_2 *tmp_r = NULL;
        Segment_2 *tmp_s = NULL;
        int is_endpoint = 0;
        interval.init();
        std::cout << "bisector_v = " << bisector.to_vector() << std::endl;
        auto a = CGAL::intersection(bisector, src->path());
        auto b = CGAL::intersection(bisector, dest->path());
        std::cout << tmp_inner_product << " " << tmp_outer_product << std::endl;
        // parallel
        if (parallel) {
            if (!a || !b)
                return false;
            if (tmp_inner_product < 0)
                std::swap(a, b);
            if (!(tmp_p = boost::get<Point_2>(&*a)) || !(tmp_q = boost::get<Point_2>(&*b)))
                return false;
            std::cout << "tmp_p, tmp_q" << *tmp_p << " " << *tmp_q << std::endl;
            FT loc_p = calc_loc_on_line(bisector, *tmp_p),
               loc_q = calc_loc_on_line(bisector, *tmp_q);
            std::cout << "loc_p, loc_q = " << loc_p << " " << loc_q << std::endl;
            // if(loc_p >= loc_q)
            if (!less(loc_p, loc_q))
                return false;
            interval.set(loc_p, loc_q, *tmp_p, *tmp_q, is_endpoint);
            return true;
        }
        else { // tmp_outer_product > 0
            Segment_2 new_seg = edge.new_loc(line_loc);
            std::cout << "new_seg" << new_seg << std::endl;
            Point_2 new_src = new_seg.source(), new_dest = new_seg.target();
            Point_2 begin, end;
            auto o = CGAL::intersection(bisector, line_loc[edge.corr_line]);
            assert(bool(o));
            Point_2 origin = *boost::get<Point_2>(&*o);

            // std::cout << "bisector= " << bisector << std::endl;
            // std::cout << new_src << " " << new_dest << " " << origin << std::endl;
            FT src_loc = calc_loc_on_line(bisector, new_src);
            FT dest_loc = calc_loc_on_line(bisector, new_dest);
            FT origin_loc = calc_loc_on_line(bisector, origin);
            // （线段运动方向背对直线，如果依然能相交，说明运动方向偏转超过中线，说明相邻边偏转超过对位边，所以相邻边不会是对位边的相邻边，那么显然这个相邻边更能代表这一点的碰撞，所以不用考虑？）（若两边在中线同侧，则正对中线的一方一定是靠近另一边的一侧先接触中线，而对应的相邻边一定更接近背对中线的边，且角度更接近）
            //std::cout << src_loc << " " << dest_loc << " " << origin_loc << std::endl;
            // if(src_loc <= origin_loc && dest_loc <= origin_loc)
            if (equal(src_loc, dest_loc))
                return false; // 已经退化为一个点
            if (!less(origin_loc, src_loc) && !less(origin_loc, dest_loc))
                return false;

            auto x = CGAL::intersection(bisector, new_seg);
            if (!x) {
                if (tmp_inner_product < 0) {
                    std::swap(a, b);
                    std::swap(src, dest);
                    is_endpoint = 2;
                }
                else
                    is_endpoint = 1;
                if (!a)
                    return false;
                if ((tmp_r = boost::get<Ray_2>(&*a))) {
                    begin = tmp_r->source();
                    is_endpoint = -is_endpoint;
                }
                else if ((tmp_p = boost::get<Point_2>(&*a))) {
                    begin = *tmp_p;
                }

                std::cout << "begin = " << begin << std::endl;

                FT loc_p = calc_loc_on_line(bisector, begin), loc_q;

                if (!b)
                    interval.set(loc_p, begin, is_endpoint);
                else if ((tmp_q = boost::get<Point_2>(&*b))) {
                    end = *tmp_q;
                    std::cout << end << std::endl;
                    loc_q = calc_loc_on_line(bisector, end);
                    std::cout << "loc_q = " << loc_q << std::endl;
                    interval.set(loc_p, loc_q, begin, end, is_endpoint);
                }
                else
                    assert(false); // 线段和中线没有交点，end不可能在中线上
                return true;
            }
            // 线段一开始就与中线相交
            // 需要排除一开始从中线远离的情况
            else if ((tmp_p = boost::get<Point_2>(&*x))) {
                // std::cout << tmp_p << std::endl;

                // if(src_loc == origin_loc)
                if (equal(src_loc, origin_loc))
                    is_endpoint = 1;
                // else if(dest_loc == origin_loc)
                else if (equal(dest_loc, origin_loc))
                    is_endpoint = 2;

                FT loc_p = origin_loc, loc_q;
                begin = *tmp_p;
                if (!a && !b) {
                    interval.set(loc_p, begin, is_endpoint);
                }
                else if (a && b) {
                    // 此前is_endpoint只判断了origin_loc是否等于src_loc和dest_loc
                    if (is_endpoint == 1) {
                        if ((tmp_r = boost::get<Ray_2>(&*a))) {
                            is_endpoint = -1;
                        }
                        // 远离bisector
                        // else if(outer_product(src->speed(), bisector.to_vector()) > 0){
                        else if (less(0, outer_product(src->speed(), bisector.to_vector()),
                                      FT(eps) * calc_line_scale(bisector))) {
                            return false;
                        }
                        if ((tmp_q = boost::get<Point_2>(&*b))) {
                            end = *tmp_q;
                            loc_q = calc_loc_on_line(bisector, end);
                        }
                        else
                            assert(false); // dest不可能也在中线上

                        interval.set(loc_p, loc_q, begin, end, is_endpoint);
                    }
                    else if (is_endpoint == 2) {
                        if ((tmp_r = boost::get<Ray_2>(&*b))) {
                            is_endpoint = -2;
                        }
                        // 远离bisector
                        // else if(outer_product(dest->speed(), bisector.to_vector()) < 0){
                        else if (less(outer_product(dest->speed(), bisector.to_vector()), 0,
                                      FT(eps) * calc_line_scale(bisector))) {
                            return false;
                        }
                        if ((tmp_q = boost::get<Point_2>(&*a))) {
                            end = *tmp_q;
                            loc_q = calc_loc_on_line(bisector, end);
                        }
                        else
                            assert(false); // src 不可能也在中线上
                        interval.set(loc_p, loc_q, begin, end, is_endpoint);
                    }
                    else {
                        // 此时交点一定在内部
                        Point_2 *tmp_src = boost::get<Point_2>(&*a);
                        Point_2 *tmp_dest = boost::get<Point_2>(&*b);
                        assert(tmp_src && tmp_dest); // 线段与bisector交点在线段内部

                        //std::cout << bisector << std::endl
                        std::cout << "tmp_src, tmp_dest = " << *tmp_src << " " << *tmp_dest << std::endl;
                        FT tmp_src_loc = calc_loc_on_line(bisector, *tmp_src);
                        FT tmp_dest_loc = calc_loc_on_line(bisector, *tmp_dest);
                        if (tmp_src_loc > tmp_dest_loc) {
                            end = *tmp_dest;
                            loc_q = tmp_dest_loc;
                        }
                        else {
                            end = *tmp_src;
                            loc_q = tmp_src_loc;
                        }
                        is_endpoint = 0;
                        interval.set(loc_p, loc_q, begin, end, is_endpoint);
                    }
                }
                // a的路线有交点或b的路线有交点
                else {
                    if (a && (tmp_q = boost::get<Point_2>(&*a))) {
                        end = *tmp_q;
                    }
                    else if (b && (tmp_q = boost::get<Point_2>(&*b))) {
                        end = *tmp_q;
                    }
                    loc_q = calc_loc_on_line(bisector, end);
                    interval.set(loc_p, loc_q, begin, end, is_endpoint);
                }
            }
            else
                assert(false); // bisector 与segment不平行
        }
        return true;
    }

    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::append_collision_events(
        EdgeData *edge)
    {
        static int run_id = 0;
        std::cout << "run_id=" << run_id++ << std::endl;

        // 求出的bisector均沿着交点运动的方向(如果平行则沿着l0的方向)
        Line_2 bisector,
            edge_opposite_line = line_data_pool[edge->corr_line].opposite();

        // 退化为三角形的情况
        if (edge->next->next == edge->prev) {
            auto inter = CGAL::intersection(edge->src_point->path(), edge->dest_point->path());
            Point_2 *incenter;
            if (inter && (incenter = boost::get<Point_2>(&*inter))) {
                FT time = CGAL::squared_distance(*incenter, line_data_pool[edge->corr_line]);
                events.emplace(time, edge, edge, EventType::TRIANGLE);
            }
        }
        int e_it_id = 0;
        for (auto e_it : edge_data_pool) {
            ++e_it_id;
            if (e_it == edge)
                break; // 只处理edge_data_pool中当前edge与位于之前的边相交的事件
            if (!e_it->is_active())
                continue;
            // 相邻边不需要添加事件
            if (e_it->src_point == edge->dest_point || edge->src_point == e_it->dest_point)
                continue;
            const Line_2 &cur_line = line_data_pool[e_it->corr_line];
            Vector_2 it_vector = cur_line.to_vector(),
                     edge_vector = edge_opposite_line.to_vector();
            const FT &tmp_inner_product = inner_product(it_vector, edge_vector);
            const FT &tmp_outer_product = outer_product(it_vector, edge_vector);
            FT absolute_eps = it_vector.squared_length() * edge_vector.squared_length() * FT(eps2);
            Point_2 L0, R0, L1, R1, point_l, point_r, collision_point;
            FT l0, r0, l1, r1, loc_l, loc_r;

            auto line_pair = std::make_pair(e_it->corr_line, edge->corr_line);
            bool parallel = is_parallel(it_vector, edge_vector);
            // 平行且同向，不可能出现重叠
            // if(tmp_outer_product == 0 && tmp_inner_product < 0) continue;
            // is_parallel(it_vector, edge_vector)
            if (parallel && (tmp_inner_product < 0 && less(0, tmp_inner_product * tmp_inner_product, absolute_eps)))
                continue;
            std::cout << "run_id=" << run_id - 1 << " "
                      << "e_it_id=" << e_it_id - 1 << std::endl;
            std::cout << "edge_corr_line = " << edge->corr_line << std::endl;
            std::cout << "e_it_corr_line = " << e_it->corr_line << std::endl;
            std::cout << e_it->src_point->path() << std::endl;
            std::cout << e_it->dest_point->path() << std::endl;
            std::cout << edge->dest_point->path() << std::endl;
            std::cout << edge->src_point->path() << std::endl;
            std::cout << "cur_line= " << cur_line << std::endl;
            std::cout << "edge= " << edge_opposite_line << std::endl;
            bool bisector_not_found = !bisector_map.count(line_pair);
            if (bisector_not_found) {
                if (tmp_inner_product >= 0)
                    bisector = CGAL::bisector(cur_line, edge_opposite_line);
                else {
                    auto intersect = CGAL::intersection(cur_line, edge_opposite_line);
                    Point_2 *origin;
                    if (!intersect || !(origin = boost::get<Point_2>(&*intersect)))
                        continue;
                    bisector = CGAL::bisector(
                        cur_line.perpendicular(*origin),
                        edge_opposite_line.opposite().perpendicular(*origin));
                    // std::cout << *origin << std::endl;
                    // std::cout << cur_line.perpendicular(*origin) << std::endl;
                    // std::cout << edge_opposite_line.opposite().perpendicular(*origin) <<
                    // std::endl;
                }
                if (inner_product(bisector.to_vector(),
                                  cur_line.to_vector().perpendicular(CGAL::LEFT_TURN)) < 0)
                    bisector = bisector.opposite();
            }
            else
                bisector = bisector_map[line_pair];

            // 判断是否确定会相交
            // L{0,1}, R{0,1}: 两个线段的端点; point_l, point_r: 线段公共部分的端点
            // l{0,1}, r{0,1}: 两个线段端点与bisector相交的位置; loc_l, loc_r:
            // 相交部分在bisector上的位置
            Interval interval0, interval1;
            if ((run_id == 15 && e_it_id == 3)) {
                start_debug();
                std::cout << run_id - 1 << " " << e_it_id - 1 << ": " << bisector
                          << std::endl;
            }
            bool tmp1 = calc_interval_on_line(bisector, *e_it, interval0, parallel);
            bool tmp2 = calc_interval_on_line(bisector, *edge, interval1, parallel);
            std::cout << interval0.to_string() << interval1.to_string();
            if (!tmp1 || !tmp2)
                continue; // (!tmp1 && !tmp2)会漏掉情况
            // if(!calc_interval_on_line(bisector, *e_it, interval0) ||
            // !calc_interval_on_line(bisector, *edge, interval1)) continue;
            Interval interval = interval0 + interval1;
            std::cout << interval.to_string();
            if (interval.is_empty)
                continue;

            if (bisector_not_found)
                bisector_map[line_pair] = bisector;

            // if(tmp_outer_product < 0) collision_point = point_r;
            // else collision_point = point_l;
            collision_point = interval.begin_point;
            FT time = CGAL::squared_distance(collision_point, cur_line);
            std::cout << "time = " << time << std::endl;
            // parallel
            // TODO: 设置一个eps
            // is_parallel(it_vector, edge_vector)
            std::cout << "is_parallel = " << parallel << std::endl;
            if (parallel) {
                events.emplace(time, e_it, edge, EventType::PARALLEL);
            }
            else
                events.emplace(time, e_it, edge, EventType::NOT_PARALLEL);
        }
    }

    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::append_line_data(
        const Polygon_2 &poly, size_t &line_it)
    {
        std::vector<PointData *> tmp_point;
        std::vector<EdgeData *> tmp_edge;
        int len = poly.size();
        PointData *src_point = NULL;
        EdgeData *src_edge = NULL;
        Vector_2 src_vector, dest_vector = poly.edge(0).to_vector();
        // build all edge_data and point_data;
        for (int i = 0; i < len; ++i) {
            line_data_pool.emplace_back(poly.edge(i));
            std::cout << poly.edge(i) << std::endl;
            std::cout << line_data_pool[line_it] << std::endl;
            line_loc.emplace_back(line_data_pool[line_it]);
            Vector_2 speed = line_data_pool[line_it].to_vector().perpendicular(CGAL::LEFT_TURN);
            line_speed.emplace_back(speed / CGAL::sqrt(speed.squared_length()));
            src_vector = dest_vector;
            dest_vector = poly.edge((i + 1) % len).to_vector();
            //std::cout << poly.edge(i).to_vector() << std::endl;
            // PointData(const Point_2 &_start_loc, const FT &_start_time, const
            // Vector_2 &src, const Vector_2 &dest);
            PointData *new_point = new PointData(poly.edge(i).target(), FT(0), src_vector, dest_vector);
            src_vector = dest_vector;
            tmp_point.push_back(new_point);
            // EdgeData(int line, PointData<K> *src, PointData<K> *dest) :
            // corr_line(line), src_point(src), dest_point(dest), prev(this), next(this)
            // {}
            EdgeData *new_edge = new EdgeData(line_it, src_point, new_point);
            src_point = new_point;
            tmp_edge.push_back(new_edge);
            // std::cout << new_edge<<std::endl;
            ++line_it;
        }
        tmp_edge[0]->src_point = src_point;
        // link them together
        src_edge = tmp_edge[len - 1];
        for (int i = 0; i < len; ++i) {
            src_edge->next = tmp_edge[i];
            tmp_edge[i]->prev = src_edge;
            src_edge = tmp_edge[i];
            edge_data_pool.push_back(src_edge);
            src_point->next = tmp_point[i];
            tmp_point[i]->prev = src_point;
            src_point = tmp_point[i];
            point_data_pool.push_back(src_point);
        }
    }

    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::CutHalfEdge(
        EdgeData *base, EdgeData *edge, bool type, Point_2 location,
        std::vector<PointData *> &point_pool, std::vector<EdgeData *> &edge_pool)
    {
        std::cout << "base->prev=" << base->prev->corr_line << std::endl;
        std::cout << "base->next=" << base->next->corr_line << std::endl;
        std::cout << "edge->prev=" << edge->prev->corr_line << std::endl;
        std::cout << "edge->next=" << edge->next->corr_line << std::endl;
        base->del = edge->del = true;
        Line_2 base_line = line_data_pool[base->corr_line],
               edge_line = line_data_pool[edge->corr_line];
        Vector_2 base_v = base_line.to_vector(), edge_v = edge_line.to_vector();
        // parallel
        if (is_parallel(base_v, edge_v)) {
            assert(inner_product(base_v, edge_v) < 0);
            Line_2 bisector = CGAL::bisector(base_line, edge_line.opposite());
            PointData *inter_seg;
            // FT bisector_eps = FT(eps) * calc_line_scale(bisector);
            if (type == 0) {
                inter_seg = new PointData(location, this->cur_time, base->dest_point->src_vector,
                                          edge->src_point->dest_vector);
                point_pool.push_back(inter_seg);
                if (!edge->src_point->branches[0]) {
                    edge->src_point->end_time = this->cur_time;
                    edge->src_point->end_loc = location;
                    edge->src_point->branches[0] = inter_seg;
                }
                else
                    edge->src_point->branches[1] = inter_seg;

                auto inter0 = CGAL::intersection(bisector, base->src_point->path()),
                     inter1 = CGAL::intersection(bisector, edge->dest_point->path());
                assert(inter0 && inter1);
                Point_2 *p0 = boost::get<Point_2>(&*inter0),
                        *p1 = boost::get<Point_2>(&*inter1);
                assert(p0 && p1);
                FT l0 = calc_loc_on_line(bisector, *p0),
                   l1 = calc_loc_on_line(bisector, *p1);
                std::cout << "p0=" << *p0 << std::endl;
                std::cout << "p1=" << *p1 << std::endl;
                if (less(l0, l1)) {
                    inter_seg->end_time = this->cur_time;
                    inter_seg->end_loc = *p1;
                    if (edge->dest_point->end_time == FT(0))
                        CutHalfEdge(base, edge->next, 0, *p1, point_pool, edge_pool);
                    // TODO: 是否有可能循环回来?
                    inter_seg->branches[0] = edge->dest_point->branches[0];
                }
                else if (less(l1, l0)) {
                    inter_seg->end_time = this->cur_time;
                    inter_seg->end_loc = *p0;
                    if (base->src_point->end_time == FT(0))
                        CutHalfEdge(edge, base->prev, 1, *p0, point_pool, edge_pool);
                    inter_seg->branches[0] = base->src_point->branches[0];
                }
                else {
                    inter_seg->end_time = this->cur_time;
                    inter_seg->end_loc = *p0;
                    if (base->src_point->end_time == FT(0))
                        CutHalfPoint(base->prev, edge->next, *p0, point_pool, edge_pool);
                    inter_seg->branches[0] = base->src_point->branches[0];
                }
            }
            else { // type == 1
                inter_seg = new PointData(location, this->cur_time, edge->dest_point->src_vector,
                                          base->src_point->dest_vector);
                point_pool.push_back(inter_seg);
                if (!edge->dest_point->branches[0]) {
                    edge->dest_point->end_time = this->cur_time;
                    edge->dest_point->end_loc = location;
                    edge->dest_point->branches[0] = inter_seg;
                }
                else
                    edge->dest_point->branches[1] = inter_seg;

                auto inter0 = CGAL::intersection(bisector, base->dest_point->path()),
                     inter1 = CGAL::intersection(bisector, edge->src_point->path());
                assert(inter0 && inter1);
                Point_2 *p0 = boost::get<Point_2>(&*inter0),
                        *p1 = boost::get<Point_2>(&*inter1);
                assert(p0 && p1);
                FT r0 = calc_loc_on_line(bisector, *p0),
                   r1 = calc_loc_on_line(bisector, *p1);
                if (less(r0, r1)) {
                    inter_seg->end_time = this->cur_time;
                    inter_seg->end_loc = *p0;
                    if (base->dest_point->end_time == FT(0))
                        CutHalfEdge(edge, base->next, 0, *p0, point_pool, edge_pool);
                    inter_seg->branches[0] = base->dest_point->branches[0];
                }
                else if (less(r1, r0)) {
                    inter_seg->end_time = this->cur_time;
                    inter_seg->end_loc = *p1;
                    if (edge->src_point->end_time == FT(0))
                        CutHalfEdge(base, edge->prev, 1, *p1, point_pool, edge_pool);
                    inter_seg->branches[0] = edge->src_point->branches[0];
                }
                else {
                    inter_seg->end_time = this->cur_time;
                    inter_seg->end_loc = *p1;
                    if (edge->src_point->end_time == FT(0))
                        CutHalfPoint(edge->prev, base->next, *p1, point_pool, edge_pool);
                    inter_seg->branches[0] = edge->src_point->branches[0];
                }
            }
            return;
        }

        if (type == 0) {
            PointData *point = new PointData(location, this->cur_time, base->src_point->dest_vector,
                                             edge->src_point->dest_vector);
            point_pool.push_back(point);
            if (!edge->src_point->branches[0]) {
                edge->src_point->end_time = this->cur_time;
                edge->src_point->end_loc = location; // 终止一个点
                edge->src_point->branches[0] = point;
            }
            else
                edge->src_point->branches[1] = point;
            EdgeData *src = new EdgeData(base->corr_line, base->src_point, point);
            EdgeData *dest = new EdgeData(edge->corr_line, point, edge->dest_point);
            edge_pool.push_back(src);
            edge_pool.push_back(dest);
            src->prev = base->prev;
            src->prev->next = src;
            dest->next = edge->next;
            dest->next->prev = dest;
            src->next = dest;
            dest->prev = src;
        }
        else {
            PointData *point = new PointData(location, this->cur_time, edge->src_point->dest_vector,
                                             base->src_point->dest_vector);
            point_pool.push_back(point);
            if (!edge->dest_point->branches[0]) {
                edge->dest_point->end_time = this->cur_time;
                edge->dest_point->end_loc = location; // 终止一个点
                edge->dest_point->branches[0] = point;
            }
            else
                edge->dest_point->branches[1] = point;
            EdgeData *src = new EdgeData(edge->corr_line, edge->src_point, point);
            EdgeData *dest = new EdgeData(base->corr_line, point, base->dest_point);
            edge_pool.push_back(src);
            edge_pool.push_back(dest);
            src->prev = edge->prev;
            src->prev->next = src;
            dest->next = base->next;
            dest->next->prev = dest;
            src->next = dest;
            dest->prev = src;
        }
    }

    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::CutHalfPoint(
        EdgeData *prev, EdgeData *next, Point_2 location,
        std::vector<PointData *> &point_pool, std::vector<EdgeData *> &edge_pool)
    {
        if (prev == next) {
            prev->dest_point->end_time = next->src_point->end_time = this->cur_time;
            prev->dest_point->end_loc = next->src_point->end_loc = location;
            return;
        }
        Line_2 prev_line = line_data_pool[prev->corr_line],
               next_line = line_data_pool[next->corr_line];
        Vector_2 prev_v = prev_line.to_vector(), next_v = next_line.to_vector();
        PointData *point = new PointData(location, this->cur_time, prev->dest_point->src_vector,
                                         next->src_point->dest_vector);
        point_pool.push_back(point);
        if (!prev->dest_point->branches[0]) {
            prev->dest_point->end_time = this->cur_time;
            prev->dest_point->end_loc = location;
            prev->dest_point->branches[0] = point;
        }
        else
            prev->dest_point->branches[1] = point;

        if (!next->src_point->branches[0]) {
            next->src_point->end_time = this->cur_time;
            next->src_point->end_loc = location;
            next->src_point->branches[0] = point;
        }
        else
            next->src_point->branches[1] = point;
        // parallel and opposite
        if (is_parallel(prev_v, next_v) && inner_product(prev_v, next_v) < 0) {
            Line_2 bisector = CGAL::bisector(prev_line.opposite(), next_line);
            auto inter0 = CGAL::intersection(bisector, prev->src_point->path()),
                 inter1 = CGAL::intersection(bisector, next->dest_point->path());
            assert(inter0 && inter1);
            Point_2 *p0 = boost::get<Point_2>(&*inter0),
                    *p1 = boost::get<Point_2>(&*inter1);
            assert(p0 && p1);
            FT r0 = calc_loc_on_line(bisector, *p0),
               r1 = calc_loc_on_line(bisector, *p1);
            if (less(r0, r1)) {
                point->end_time = this->cur_time;
                point->end_loc = *p0;
                if (prev->src_point->end_time == FT(0))
                    CutHalfEdge(next, prev->prev, 1, *p0, point_pool, edge_pool);
                point->branches[0] = prev->src_point->branches[0];
            }
            else if (less(r1, r0)) {
                point->end_time = this->cur_time;
                point->end_loc = *p1;
                if (next->dest_point->end_time == FT(0))
                    CutHalfEdge(prev, next->next, 0, *p1, point_pool, edge_pool);
                point->branches[0] = next->dest_point->branches[0];
            }
            else {
                point->end_time = this->cur_time;
                point->end_loc = *p1;
                if (next->dest_point->end_time == FT(0))
                    CutHalfPoint(prev->prev, next->next, *p1, point_pool, edge_pool);
                point->branches[0] = next->dest_point->branches[0];
            }
            return;
        }

        prev->del = next->del = true;
        EdgeData *new_prev = new EdgeData(prev->corr_line, prev->src_point, point);
        EdgeData *new_next = new EdgeData(next->corr_line, point, next->dest_point);
        new_prev->prev = prev->prev;
        prev->prev->next = new_prev;
        new_next->next = next->next;
        next->next->prev = new_next;
        new_prev->next = new_next;
        new_next->prev = new_prev;
        edge_pool.push_back(new_prev);
        edge_pool.push_back(new_next);
    }

    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::CutEdge(
        EdgeData *base, EdgeData *prev, EdgeData *next, Point_2 location,
        std::vector<PointData *> &point_pool, std::vector<EdgeData *> &edge_pool)
    {
        std::cout << "base= " << base->new_loc(line_loc) << std::endl;
        std::cout << "prev= " << prev->new_loc(line_loc) << std::endl;
        std::cout << "next= " << next->new_loc(line_loc) << std::endl;
        std::cout << "location= " << location << std::endl;
        CutHalfEdge(base, prev, 1, location, point_pool, edge_pool);
        CutHalfEdge(base, next, 0, location, point_pool, edge_pool);
    }

    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::operator()(
        const Poly_with_holes &polygon)
    {
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

        FT last_time = 0;
        std::vector<PointData *> point_pool;
        std::vector<EdgeData *> edge_pool;
        while (!events.empty()) {
            Event event = events.top();
            events.pop();
            if (!event.is_active())
                continue;
            point_pool.clear();
            edge_pool.clear();
            this->cur_time = event.time;
            FT exact_time = CGAL::sqrt(event.time);
            for (size_t i = 0; i < line_data_pool.size(); ++i) {
                line_loc[i] = Line_2(line_data_pool[i].point(0) + line_speed[i] * exact_time,
                                     line_data_pool[i].to_vector());
            }
            std::cout << "time = sqrt(" << event.time << ") = " << exact_time
                      << std::endl;
            EdgeData *event_a = event.a, *event_b = event.b;
            Segment_2 seg0 = event_a->new_loc(line_loc),
                      seg1 = event_b->new_loc(line_loc);
            Line_2 bisector = bisector_map[std::make_pair(event_a->corr_line, event_b->corr_line)];
            FT l0 = calc_loc_on_line(bisector, seg0.source()),
               r0 = calc_loc_on_line(bisector, seg0.target());
            FT l1 = calc_loc_on_line(bisector, seg1.target()),
               r1 = calc_loc_on_line(bisector, seg1.source());
            std::cout << "event.a=" << event_a->corr_line
                      << "     event.b=" << event_b->corr_line << std::endl;
            std::cout << event_a << " " << event_b << std::endl;
            std::cout << "l0 r0 l1 r1= " << l0 << " " << r0 << std::endl
                      << l1 << " " << r1 << std::endl
                      << std::endl;
            std::cout << "seg0 seg1=" << seg0 << " " << seg1 << std::endl;
            std::cout << "is_parallel = " << (event.type == EventType::PARALLEL) << std::endl;
            if (event.type == EventType::TRIANGLE) {
                EdgeData *it = event_a;
                if (it->src_point->end_time == FT(0)) {
                    Line_2 line[3] = {line_loc[it->corr_line], line_loc[it->next->corr_line], line_loc[it->prev->corr_line]};
                    auto inter0 = CGAL::intersection(line[1], line[2]);
                    auto inter1 = CGAL::intersection(line[2], line[0]);
                    auto inter2 = CGAL::intersection(line[0], line[1]);
                    if (!inter0 || !inter1 || !inter2) {
                        std::cerr << "triangle: no intersection" << std::endl;
                        exit(1);
                    }
                    Point_2 *pts[3] = {boost::get<Point_2>(&*inter0), boost::get<Point_2>(&*inter1), boost::get<Point_2>(&*inter2)};
                    PointData *data[3] = {it->next->dest_point, it->src_point, it->dest_point};
                    if (!pts[0] || !pts[1] || !pts[2]) {
                        std::cerr << "triangle: intersection not Point_2" << std::endl;
                        exit(1);
                    }
                    Point_2 barycenter = CGAL::ORIGIN + (((*pts[0] - CGAL::ORIGIN) + (*pts[1] - CGAL::ORIGIN) + (*pts[2] - CGAL::ORIGIN)) / 3);
                    FT end_time = CGAL::squared_distance(barycenter, line[0]);
                    data[0]->end_time = data[1]->end_time = data[2]->end_time = end_time;
                    data[0]->end_loc = data[1]->end_loc = data[2]->end_loc = barycenter;
                }
            }
            else if (event.type == EventType::PARALLEL) {
                /* 输出公共部分 TODO: 改成创建点 */
                this->res_segments.emplace_back(
                    (l0 < l1) ? seg1.target() : seg0.source(),
                    (r0 < r1) ? seg0.target() : seg1.source());
                if (event_a->prev != event_b->next) {
                    // if(l0 < l1)
                    if (less(l0, l1))
                        CutHalfEdge(event_a, event_b->next, 0, seg1.target(), point_pool, edge_pool);
                    // else if(l0 > l1)
                    else if (less(l1, l0))
                        CutHalfEdge(event_b, event_a->prev, 1, seg0.source(), point_pool, edge_pool);
                    else {
                        CutHalfPoint(event_a->prev, event_b->next, seg0.source(), point_pool, edge_pool);
                    }
                    /* if(seg0.source()和seg1.target()相邻){ 直接跳过 }
else if(l0 < l1){ 用seg1.target()的位置建点，切割seg0，重连 }
else if(l0 > l1){ 用seg0.source()的位置建点，切割seg1，重连 }*/
                }
                else {
                    event_b->dest_point->end_time = event_a->src_point->end_time = event.time;
                    event_a->src_point->end_loc = seg0.source();
                    event_b->dest_point->end_loc = seg1.target();
                }
                if (event_a->next != event_b->prev) {
                    // if(r0 < r1)
                    if (less(r0, r1))
                        CutHalfEdge(event_b, event_a->next, 0, seg0.target(), point_pool, edge_pool);
                    // else if(r0 > r1)
                    else if (less(r1, r0))
                        CutHalfEdge(event_a, event_b->prev, 1, seg1.source(), point_pool, edge_pool);
                    else {
                        CutHalfPoint(event_b->prev, event_a->next, seg1.source(), point_pool, edge_pool);
                    }
                    /* if(seg0.source()和seg1.target()相邻){ 直接跳过 }
if(r0 < r1){ 用seg0.target()的位置建点，切割seg1，重连 }
else if(r0 > r1){ 用seg1.source()的位置建点，切割seg0，重连 } */
                }
                else {
                    event_b->src_point->end_time = event_a->dest_point->end_time = event.time;
                    event_a->dest_point->end_loc = seg0.target();
                    event_b->src_point->end_loc = seg1.source();
                }
                std::cout << "parallel: a=" << event_a->src_point->start_loc << " "
                          << event_a->src_point->end_loc << std::endl;
                std::cout << "parallel: a=" << event_a->dest_point->start_loc << " "
                          << event_a->dest_point->end_loc << std::endl;
                std::cout << "parallel: b=" << event_b->src_point->start_loc << " "
                          << event_b->src_point->end_loc << std::endl;
                std::cout << "parallel: b=" << event_b->dest_point->start_loc << " "
                          << event_b->dest_point->end_loc << std::endl;
            }
            else {
                if (l0 > r0) {
                    std::swap(seg0, seg1);
                    std::swap(l0, r1);
                    std::swap(l1, r0);
                    std::swap(event_a, event_b);
                }
                if (event_a->prev == event_b->next) {
                    CutHalfPoint(event_b, event_a, seg0.source(), point_pool, edge_pool);
                }
                // else if(l0 < l1){
                else if (less(l0, l1)) {
                    CutEdge(event_a, event_b, event_b->next, seg1.target(), point_pool, edge_pool);
                }
                // else if(l0 > l1){
                else if (less(l1, l0)) {
                    CutEdge(event_b, event_a->prev, event_a, seg0.source(), point_pool, edge_pool);
                }
                else {
                    CutHalfPoint(event_b, event_a, seg0.source(), point_pool, edge_pool);
                    CutHalfPoint(event_a->prev, event_b->next, seg1.target(), point_pool, edge_pool);
                }
            }

            for (EdgeData *it : edge_pool) {
                append_collision_events(it);
            }
            edge_data_pool.insert(edge_data_pool.end(), edge_pool.begin(), edge_pool.end());
            point_data_pool.insert(point_data_pool.end(), point_pool.begin(), point_pool.end());

            std::cout << "segments:=====================================" << std::endl;
            for (PointData *it : point_data_pool)
                std::cout << it->start_loc << " => " << it->end_loc << std::endl;
        }

        std::cout << "segments:=====================================" << std::endl;
        for (PointData *it : point_data_pool)
            std::cout << it->start_loc << " => " << it->end_loc << std::endl;

        int counter = 0;
        for (PointData *it : point_data_pool) {
            ++counter;
            if (it->end_time == FT(0)) {
                std::cout << "point " << counter - 1 << "not ended" << std::endl;
                continue;
            }
            FT tmp_inner_product = inner_product(it->src_vector, it->dest_vector);
            FT absolute_scale = it->src_vector.squared_length() * it->dest_vector.squared_length();
            if (outer_product(it->src_vector, it->dest_vector) >= 0 && tmp_inner_product < 0 && tmp_inner_product * tmp_inner_product > FT(0.25 - eps) * absolute_scale) {
                this->res_segments.emplace_back(it->start_loc, it->end_loc);
            }
            else
                this->sub_segments.emplace_back(it->start_loc, it->end_loc);
        }
    }
} // namespace CenterLineSolver
#endif // GET_CENTERLINE_H