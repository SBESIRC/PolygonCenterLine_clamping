#pragma once
#include <queue>
#include "CenterLineSolver.h"
#include "convertKernel.h"

namespace CenterLineSolver {
    using IK = CGAL::Epeck;
    using IFT = IK::FT;
    inline void get_intersection(const CGAL::Polygon_2<IK> &poly, const CGAL::Ray_2<IK> &ray, IK::Point_2 &new_p, IFT &cur_dis){
        CGAL::Point_2<IK> *tmp_p;
        CGAL::Segment_2<IK> *tmp_seg;
        for(auto it = poly.edges_begin(); it != poly.edges_end();++it){
            auto inter = CGAL::intersection(*it, ray);
            if(!inter) continue;
            if((tmp_p = boost::get<CGAL::Point_2<IK>>(&*inter))){
                IFT dis = CGAL::squared_distance(*tmp_p, ray.source());
                if(cur_dis < 0 || dis < cur_dis){
                    cur_dis = dis;
                    new_p = *tmp_p;
                }
            }
            else if((tmp_seg = boost::get<CGAL::Segment_2<IK>>(&*inter))){
                IFT dis_0 = CGAL::squared_distance(tmp_seg->source(), ray.source());
                IFT dis_1 = CGAL::squared_distance(tmp_seg->target(), ray.source());
                if(cur_dis < 0 || CGAL::min(dis_0, dis_1) < cur_dis){
                    if(dis_0 < dis_1){
                        cur_dis = dis_0;
                        new_p = tmp_seg->source();
                    }
                    else{
                        cur_dis = dis_1;
                        new_p = tmp_seg->target();
                    }
                }
            }
        }
    }
    template <typename K>
    inline typename K::Point_2 get_intersection(const CGAL::Polygon_with_holes_2<K> &poly, typename K::Point_2 point, typename K::Vector_2 dir){
        KernelConverter::KernelConverter<K, CGAL::Epeck, KernelConverter::NumberConverter<K::FT, IK::FT>> to_exact;
        KernelConverter::KernelConverter<CGAL::Epeck, K, KernelConverter::NumberConverter<CGAL::Epeck::FT, K::FT, 256>> to_Gmpfr;
        CGAL::Polygon_with_holes_2<IK> space = to_exact.convert(poly);
        CGAL::Point_2<IK> p = to_exact(point), new_p;
        IFT cur_dis = -1;
        CGAL::Vector_2<IK> v = to_exact(dir);
        CGAL::Ray_2<IK> ray(p, v);
        std::cout << "polygon = " << space << std::endl << "ray = " << ray << std::endl;
        get_intersection(space.outer_boundary(), ray, new_p, cur_dis);
        for(auto h_it = space.holes_begin(); h_it != space.holes_end();++h_it)
            get_intersection(*h_it, ray, new_p, cur_dis);
        if(cur_dis < 0){
            std::cerr << "get_intersection(poly, point, dir): no intersection" << std::endl << point << " " << dir << std::endl;
            throw("get_intersection(poly, point, dir): no intersection");
        }
        return to_Gmpfr(new_p);
    }
    template<typename K>
    inline bool evalDirection(CGAL::Vector_2<K> from, CGAL::Vector_2<K> base_v, CGAL::Vector_2<K> &used_vector, typename K::FT &value){
        using FT = typename K::FT;
        using Solver = CenterLineSolver<K, CGAL::Polygon_with_holes_2<K>, CGAL::Polygon_2<K>>;
        bool upd = false;
        FT Cos = Solver::Cosine(base_v, from);
        if(Solver::equal(Cos, 1)) {
            if(value < 1){
                used_vector = base_v;
                value = 1;
                upd = true;
            }
        }
        else if(Cos >= 0){
            if(value < Cos){
                used_vector = from;
                value = Cos;
                upd = true;
            }
        }
        else if(value < Cos) {
            used_vector = base_v;
            value = Cos;
            upd = true;
        }
        return upd;
    }
    template<typename K>
    inline typename K::FT evalDirection(CGAL::Vector_2<K> from, CGAL::Vector_2<K> base_v, CGAL::Vector_2<K> to){

    }
    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::connect_segments(){
        for(int i = 0;i < locations.size();++i){
            std::cout << "location " << i << ": <" << locations[i].point << ", " << locations[i].time << ">" << std::endl;
        }
        for(auto p : point_data_pool){
            std::cout << p->point_id << ":" << std::endl;
            std::cout << "start_loc = " << p->start_loc << std::endl;
            std::cout << "end_loc = " << p->end_loc << std::endl;
        }

        std::vector<bool> visited(locations.size(), 0), inQ(locations.size(), 0);
        std::vector<std::vector<Vector_2>> in_vectors(locations.size());
        std::priority_queue<std::pair<FT, int>> Q; // <-time, loc_id>

        for(auto p : point_data_pool) if(p->is_ans) {
            if(p->end_loc == -1){
                std::cerr << "point_data " << p->point_id << " not ended" << std::endl;
            }
            else{
                in_vectors[p->start_loc].push_back(p->point(0) - p->point(1));
                in_vectors[p->end_loc].push_back(p->point(1) - p->point(0));
            }
        }
        std::cout << "preprocessing finished" << std::endl;
        size_t locations_size = locations.size();
        for(int i = 0;i < locations.size();++i) if(locations[i].branches.size() == 1){ // leaves
            visited[i] = true;
            auto &e = locations[i].branches[0];
            if(e.first->end_loc == -1) continue;
            int new_loc = e.first->location(e.second ^ 1);
            std::cout << "i=" << i << "\te.first->point_id=" << e.first->point_id << std::endl;
            std::cout << "e.second = " << e.second << std::endl;
            std::cout << "new_loc = " << new_loc << std::endl;
            if(e.first->is_ans){
                Vector_2 v = e.first->target() - e.first->source();
                if(e.second == 1) v = -v;
                std::cout << v;
                in_vectors[new_loc].push_back(v);
            }
            if(!inQ[new_loc]){
                inQ[new_loc] = true;
                Q.emplace(-locations[new_loc].time, new_loc);
            }
        }
        while(!Q.empty()){
            int loc_id = Q.top().second; Q.pop();
            visited[loc_id] = true;
            
            for(auto &e : locations[loc_id].branches) if(e.second == false) {
                if(e.first->end_loc == -1){
                    std::cerr << "Point " << e.first->point_id << " not ended." << std::endl;
                    continue;
                }
                int new_loc = e.first->location(e.second ^ 1);
                std::cout << "\ncur = " << loc_id << "\tnew_loc="<< new_loc << std::endl;
                if(visited[new_loc]) continue;
                if(!e.first->is_ans && in_vectors[loc_id].size() != 0){
                    const Point_2 &src = e.first->point(e.second), &dest = e.first->point(e.second^1);
                    Vector_2 base_v = dest - src;
                    std::cout << "connecting PointData " << e.first->point_id << std::endl;
                    std::cout << "from: "; for(auto v : in_vectors[loc_id]) std::cout << v << " "; std::cout << std::endl;
                    std::cout << e.first->point(e.second) << " => " << e.first->point(e.second^1) << std::endl;
                    std::cout << "vect = " << e.first->point(e.second^1) - e.first->point(e.second) << std::endl;
                    std::cout << "to: "; for(auto v : in_vectors[new_loc]) std::cout << -v << " "; std::cout << std::endl;

                    if(equal(src.x(), dest.x()) && equal(src.y(), dest.y())){ // �˻��ĵ�
                        throw("repeated point");
                        //in_vectors[new_loc].insert(in_vectors[new_loc].end(), in_vectors[loc_id].begin(), in_vectors[loc_id].end());
                    }
                    else if(in_vectors[new_loc].size() == 0){ // ��ص�����PointData�����Ǵ�
                        Vector_2 used_vector;
                        FT value = -1;
                        for(auto v : in_vectors[loc_id]){
                            evalDirection(v, base_v, used_vector, value);
                        }
                        if(value == 1 || value < 0){
                            PointData *new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                            point_data_pool.push_back(new_point);
                            in_vectors[loc_id].push_back(-base_v);
                            in_vectors[new_loc].push_back(base_v);
                        }
                        else{
                            Point_2 foot = Line_2(src, used_vector).projection(dest);
                            int turning = locations.size();
                            locations.push_back(Location(foot, (locations[loc_id].time + locations[new_loc].time) / 2));
                            PointData *point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                            point_data_pool.push_back(point0);
                            in_vectors[loc_id].push_back(src - foot);

                            PointData *point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                            point_data_pool.push_back(point1);
                            in_vectors[new_loc].push_back(dest - foot);
                        }
                    }
                    else{ // �������in_vectors
                        Vector_2 used_vector0, used_vector1, used_v0, used_v1;
                        FT value = -2, cur_value;
                        for(auto v0 : in_vectors[loc_id]) {
                            FT Sin0 = outer_product(v0, base_v) / CGAL::sqrt(v0.squared_length() * base_v.squared_length());
                            FT Cos0 = Cosine(v0, base_v);
                            for(auto v1 : in_vectors[new_loc]){
                                FT Sin1 = outer_product(base_v, -v1) / CGAL::sqrt(v1.squared_length() * base_v.squared_length());
                                FT Cos1 = Cosine(-v1, base_v);
                                if(equal(Cos0, 1) || equal(Cos1, 1)){
                                    used_v0 = used_v1 = base_v;
                                    cur_value = 4;
                                }
                                else if(Cos0 < 0 || Cos1 < 0){
                                    used_v0 = used_v1 = base_v;
                                    cur_value = (Cos0 + Cos1) / 2;
                                }
                                else if(Sin0 * Sin1 < 0) {
                                    used_v0 = v0; used_v1 = -v1;
                                    //cur_value = -Sin0 * Sin1;
                                    cur_value = Cosine(v0, -v1);
                                }
                                else{
                                    used_v0 = v0; used_v1 = -v1;
                                    FT tmp = Sine(v0, -v1);
                                    cur_value = 2 * tmp * tmp;
                                }
                                if(cur_value > value){
                                    value = cur_value;
                                    used_vector0 = used_v0, used_vector1 = used_v1;
                                }
                            }
                        }
                        if(used_vector0 == base_v && used_vector1 == base_v){ // ֱ������
                            PointData *new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                            point_data_pool.push_back(new_point);
                            in_vectors[loc_id].push_back(-base_v);
                            in_vectors[new_loc].push_back(base_v);
                        }
                        else if(outer_product(used_vector0, base_v) * outer_product(base_v, used_vector1) < 0){
                            Vector_2 v0 = used_vector0 / CGAL::sqrt(used_vector0.squared_length());
                            Vector_2 v1 = used_vector1 / CGAL::sqrt(used_vector1.squared_length());
                            Point_2 mid = src + (dest - src) / 2, foot0, foot1;
                            if(inner_product(v0, base_v) < inner_product(v1, base_v)){
                                foot0 = Line_2(src, used_vector0).projection(mid);
                                Line_2 mid_line(mid, used_vector0.perpendicular(CGAL::LEFT_TURN));
                                auto inter = CGAL::intersection(mid_line, Line_2(dest, used_vector1));
                                Point_2 *tmp_p;
                                if(!inter || !(tmp_p = boost::get<Point_2>(&*inter))) throw("no intersection of mid_line and v1");
                                foot1 = *tmp_p;
                            }
                            else{
                                foot1 = Line_2(dest, used_vector1).projection(mid);
                                Line_2 mid_line(mid, used_vector1.perpendicular(CGAL::LEFT_TURN));
                                auto inter = CGAL::intersection(mid_line, Line_2(src, used_vector0));
                                Point_2 *tmp_p;
                                if(!inter || !(tmp_p = boost::get<Point_2>(&*inter))) throw("no intersection of mid_line and v0");
                                foot0 = *tmp_p;
                            }
                            int turning0 = locations.size();
                            locations.push_back(Location(foot0, (locations[loc_id].time * 3 + locations[new_loc].time) / 4));
                            int turning1 = locations.size();
                            locations.push_back(Location(foot1, (locations[loc_id].time + locations[new_loc].time * 3) / 4));
                            PointData *point0 = new PointData(locations, loc_id, turning0, point_data_pool.size());
                            point_data_pool.push_back(point0);
                            in_vectors[loc_id].push_back(src - foot0);

                            PointData *point1 = new PointData(locations, turning0, turning1, point_data_pool.size());
                            point_data_pool.push_back(point1);
                            PointData *point2 = new PointData(locations, turning1, new_loc, point_data_pool.size());
                            point_data_pool.push_back(point2);
                            in_vectors[new_loc].push_back(dest - foot1);
                        }
                        else{
                            auto inter = CGAL::intersection(Line_2(src, used_vector0), Line_2(dest, used_vector1));
                            Point_2 *tmp_p;
                            if(!inter || !(tmp_p = boost::get<Point_2>(&*inter))) throw("no intersection of v0 and v1");
                            int turning = locations.size();
                            locations.push_back(Location(*tmp_p, (locations[loc_id].time + locations[new_loc].time) / 2));
                            PointData *point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                            point_data_pool.push_back(point0);
                            in_vectors[loc_id].push_back(src - *tmp_p);

                            PointData *point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                            point_data_pool.push_back(point1);
                            in_vectors[new_loc].push_back(dest - *tmp_p);
                        }
                    }
                }
                else if(in_vectors[loc_id].size() == 0 && locations[loc_id].branches.size() > 3) {
                    const Point_2 &src = e.first->point(e.second), &dest = e.first->point(e.second^1);
                    Vector_2 base_v = dest - src;
                    if(in_vectors[new_loc].size()){
                        Vector_2 used_vector;
                        FT value = -1;
                        for(auto v : in_vectors[new_loc]){
                            evalDirection(v, -base_v, used_vector, value);
                        }
                        if(value == 1 || value < 0){ // simply link
                            PointData *new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                            point_data_pool.push_back(new_point);
                            in_vectors[new_loc].push_back(base_v);
                            in_vectors[loc_id].push_back(-base_v);
                        }
                        else{
                            Point_2 foot = Line_2(src, used_vector).projection(dest);
                            int turning = locations.size();
                            locations.push_back(Location(foot, (locations[loc_id].time + locations[new_loc].time) / 2));
                            PointData *point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                            point_data_pool.push_back(point0);
                            in_vectors[loc_id].push_back(src - foot);

                            PointData *point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                            point_data_pool.push_back(point1);
                            in_vectors[new_loc].push_back(dest - foot);
                        }
                    }
                    else{ // in_vectors of src and dest are both empty
                        if(equal(base_v.x(), 0) || equal(base_v.y(), 0)){
                            PointData *new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                            point_data_pool.push_back(new_point);
                            in_vectors[new_loc].push_back(base_v);
                            in_vectors[loc_id].push_back(-base_v);
                        }
                        else{
                            Point_2 foot;
                            if(base_v.x() < base_v.y()) foot = Point_2(src.x(), dest.y());
                            else foot = Point_2(dest.x(), src.y());
                            int turning = locations.size();
                            locations.push_back(Location(foot, (locations[loc_id].time + locations[new_loc].time) / 2));
                            PointData *point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                            point_data_pool.push_back(point0);
                            in_vectors[loc_id].push_back(src - foot);

                            PointData *point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                            point_data_pool.push_back(point1);
                            in_vectors[new_loc].push_back(dest - foot);
                        }
                    }
                }
                if(!inQ[new_loc]){
                    inQ[new_loc] = true;
                    Q.emplace(-locations[new_loc].time, new_loc);
                }
            }
        }
        // extend leaf segments
        std::vector<int> degree(locations.size(), 0);
        for(auto p : point_data_pool) if(p->is_ans) {
            ++degree[p->start_loc];
            std::cout << p->point_id << " is ans\n";
            std::cout << "ans=" << p->source() << " " << p->target()<< std::endl;
            if(p->end_loc != -1) ++degree[p->end_loc];
        }
        for(size_t loc_id = 0;loc_id != locations_size;++loc_id){
            //if(degree[loc_id] == 1){
            if(degree[loc_id] == 1){
                if(in_vectors[loc_id].size() == 0){
                    throw("connect_segments error: leaf node has no in_vector");
                }
                std::cout << "degree" << in_vectors[loc_id].size() << std::endl;
                std::cout << "leaf point " << loc_id << ": " << locations[loc_id].point << std::endl;
                std::cout << "vector=" << in_vectors[loc_id][0] << std::endl;
                try{
                    Point_2 new_p = get_intersection(*origin_space, locations[loc_id].point, in_vectors[loc_id][0]);
                    if(CGAL::squared_distance(new_p, locations[loc_id].point) < 1e-6){
                        continue;
                    }
                    std::cout << "intersection " << locations[loc_id].point << " " << new_p << std::endl;
                    size_t endpoint = locations.size();
                    locations.push_back(Location(new_p, locations[loc_id].time));
                    PointData *point0 = new PointData(locations, loc_id, endpoint, point_data_pool.size());
                    point_data_pool.push_back(point0);
                }
                catch (const char *str){
                    std::cerr << "connect_segments error = \n" << str << std::endl;
                }
            }
            else if(degree[loc_id] == 0 && locations[loc_id].branches.size() > 2){
                bool has_out = false;
                for(auto &branch : locations[loc_id].branches) if(branch.second == 0) {
                    has_out = true;
                    break;
                }
                if(!has_out){
                    std::cout << "isoloated point" << locations[loc_id].point << std::endl;
                    try{
                        Point_2 l = get_intersection(*origin_space, locations[loc_id].point, Vector_2(-1, 0));
                    std::cout << "intersection " << locations[loc_id].point << " " << l << std::endl;
                        size_t endpoint = locations.size();
                        locations.push_back(Location(l, locations[loc_id].time));
                        PointData *point0 = new PointData(locations, loc_id, endpoint, point_data_pool.size());
                        point_data_pool.push_back(point0);

                        Point_2 r = get_intersection(*origin_space, locations[loc_id].point, Vector_2(1, 0));
                    std::cout << "intersection " << locations[loc_id].point << " " << r << std::endl;
                        endpoint = locations.size();
                        locations.push_back(Location(r, locations[loc_id].time));
                        PointData *point1 = new PointData(locations, loc_id, endpoint, point_data_pool.size());
                        point_data_pool.push_back(point1);
                    }
                    catch (const char *str){
                        std::cerr << "connect_segments error = \n" << str << std::endl;
                    }
                }
            }
        }
    }
} // namespace CenterLineSolver