#pragma once
#include "CenterLineSolver.h"
namespace CenterLineSolver{
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
        if (edge->next->next == edge->prev && outer_product(edge->src_point->src_vector, edge->src_point->dest_vector) > 0) {
            auto inter = CGAL::intersection(edge->src_point->path(), edge->dest_point->path());
            Point_2 *incenter;
            if (inter && (incenter = boost::get<Point_2>(&*inter))) {
                FT time = CGAL::squared_distance(*incenter, line_data_pool[edge->corr_line]);
                events.emplace(time, edge, edge, EventType::TRIANGLE);
            }
        }
        // 退化为平行四边形的情形
        else if (edge->next->next == edge->prev->prev &&
                 is_parallel(edge->src_point->src_vector, edge->dest_point->dest_vector) && is_parallel(edge->src_point->dest_vector, edge->prev->src_point->src_vector)) {
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
            if (is_debugging) {
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
            }
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
                    if (is_debugging) {
                        std::cout << *origin << std::endl;
                        std::cout << cur_line.perpendicular(*origin) << std::endl;
                        std::cout << edge_opposite_line.opposite().perpendicular(*origin) << std::endl;
                    }
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
            //if ((run_id-1 == 13 && e_it_id-1 == 1)) {
            //    start_debug();
            //    std::cout << run_id - 1 << " " << e_it_id - 1 << ": " << bisector << std::endl;
            //}
            bool tmp1 = calc_interval_on_line(bisector, *e_it, *edge, interval0, parallel);
            bool tmp2 = calc_interval_on_line(bisector, *edge, *e_it, interval1, parallel);
            if (is_debugging)
                std::cout << interval0.to_string() << interval1.to_string();
            if (!tmp1 || !tmp2)
                continue; // (!tmp1 && !tmp2)会漏掉情况
            // if(!calc_interval_on_line(bisector, *e_it, interval0) ||
            // !calc_interval_on_line(bisector, *edge, interval1)) continue;
            Interval interval = interval0 + interval1;
            if (is_debugging)
                std::cout << interval.to_string();
            if (interval.is_empty)
                continue;

            if (bisector_not_found)
                bisector_map[line_pair] = bisector;

            // if(tmp_outer_product < 0) collision_point = point_r;
            // else collision_point = point_l;
            collision_point = interval.begin_point;
            FT time = CGAL::squared_distance(collision_point, cur_line);
            //std::cout << "time = " << time << std::endl;
            // parallel
            // TODO: 设置一个eps
            // is_parallel(it_vector, edge_vector)
            //std::cout << "is_parallel = " << parallel << std::endl;
            if (parallel) {
                events.emplace(time, e_it, edge, EventType::PARALLEL);
            }
            else
                events.emplace(time, e_it, edge, EventType::NOT_PARALLEL);
        }
    }

}
