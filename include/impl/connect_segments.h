#pragma once
#include <queue>
#include "CenterLineSolver.h"

namespace CenterLineSolver {
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
            in_vectors[p->start_loc].push_back(p->point(0) - p->point(1));
            in_vectors[p->end_loc].push_back(p->point(1) - p->point(0));
        }
        std::cout << "preprocessing finished" << std::endl;
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
            inQ[new_loc] = true;
            Q.emplace(-locations[new_loc].time, new_loc);
        }
        while(!Q.empty()){
            int loc_id = Q.top().second; Q.pop();
            visited[loc_id] = true;
            for(auto &e : locations[loc_id].branches) {
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
                    std::cout << "to: "; for(auto v : in_vectors[new_loc]) std::cout << -v << " "; std::cout << std::endl;

                    if(equal(src.x(), dest.x()) && equal(src.y(), dest.y())){ // 退化的点
                        in_vectors[new_loc].insert(in_vectors[new_loc].end(), in_vectors[loc_id].begin(), in_vectors[loc_id].end());
                    }
                    else if(in_vectors[new_loc].size() == 0){ // 相关的所有PointData均不是答案
                        Vector_2 used_vector;
                        FT value = -1;
                        for(auto v : in_vectors[loc_id]){
                            FT Cos = Cosine(base_v, v);
                            if(equal(Cos, 1)) {
                                if(value < 1){
                                    used_vector = base_v;
                                    value = 1;
                                }
                            }
                            else if(Cos >= 0){
                                if(value < Cos){
                                    used_vector = v;
                                    value = Cos;
                                }
                            }
                            else if(value < Cos) {
                                used_vector = base_v;
                                value = Cos;
                            }
                        }
                        if(value == 1 || value < 0){
                            PointData *new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                            point_data_pool.push_back(new_point);
                            in_vectors[new_loc].push_back(base_v);
                        }
                        else{
                            Point_2 foot = Line_2(src, used_vector).projection(dest);
                            int turning = locations.size();
                            locations.emplace_back(foot, (locations[loc_id].time + locations[new_loc].time) / 2);
                            PointData *point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                            point_data_pool.push_back(point0);
                            PointData *point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                            point_data_pool.push_back(point1);
                            in_vectors[new_loc].push_back(dest - foot);
                        }
                    }
                    else{ // 两侧均有in_vectors
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
                                    cur_value = -Sin0 * Sin1;
                                }
                                else{
                                    used_v0 = v0; used_v1 = -v1;
                                    cur_value = Cosine(v0, -v1);
                                }
                                if(cur_value > value){
                                    value = cur_value;
                                    used_vector0 = used_v0, used_vector1 = used_v1;
                                }
                            }
                        }
                        if(used_vector0 == base_v && used_vector1 == base_v){ // 直接连线
                            PointData *new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                            point_data_pool.push_back(new_point);
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
                            locations.emplace_back(foot0, (locations[loc_id].time * 3 + locations[new_loc].time) / 4);
                            int turning1 = locations.size();
                            locations.emplace_back(foot1, (locations[loc_id].time + locations[new_loc].time * 3) / 4);
                            PointData *point0 = new PointData(locations, loc_id, turning0, point_data_pool.size());
                            point_data_pool.push_back(point0);
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
                            locations.emplace_back(*tmp_p, (locations[loc_id].time + locations[new_loc].time) / 2);
                            PointData *point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                            point_data_pool.push_back(point0);
                            PointData *point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                            point_data_pool.push_back(point1);
                            in_vectors[new_loc].push_back(dest - *tmp_p);
                        }
                    }
                }
                if(!inQ[new_loc]){
                    inQ[new_loc] = true;
                    Q.emplace(-locations[new_loc].time, new_loc);
                }
            }
        }
	}
} // namespace CenterLineSolver