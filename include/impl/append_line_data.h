#pragma once
#include "CenterLineSolver.h"
namespace CenterLineSolver{
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
            std::cout << i << " " << poly.edge(i) << std::endl;
            std::cout << line_data_pool[line_it] << std::endl;
            line_loc.emplace_back(line_data_pool[line_it]);
            Vector_2 speed = line_data_pool[line_it].to_vector().perpendicular(CGAL::LEFT_TURN);
            line_speed.emplace_back(speed / CGAL::sqrt(speed.squared_length()));
            src_vector = dest_vector;
            dest_vector = poly.edge((i + 1) % len).to_vector();
            //std::cout << poly.edge(i).to_vector() << std::endl;
            // PointData(const Point_2 &_start_loc, const FT &_start_time, const
            // Vector_2 &src, const Vector_2 &dest);
            int loc_id = locations.size();
            locations.emplace_back(poly.edge(i).target(), 0);
            PointData *new_point = new PointData(locations, loc_id, src_vector, dest_vector);

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
}
