#include "CenterLineSolver.h"
namespace CenterLineSolver{
    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::cut_half_point(
        EdgeData *prev, EdgeData *next, int location,
        std::vector<PointData *> &point_pool, std::vector<EdgeData *> &edge_pool)
    {
        if (prev == next) {
            prev->dest_point->set_end(location);
            next->src_point->set_end(location);
            return;
        }
        Line_2 prev_line = line_data_pool[prev->corr_line],
               next_line = line_data_pool[next->corr_line];
        Vector_2 prev_v = prev_line.to_vector(), next_v = next_line.to_vector();
        PointData *point = new PointData(locations, location, prev->to_vector(),
                                         next->to_vector());
        point_pool.push_back(point);
        if (prev->dest_point->end_loc == -1) {
            prev->dest_point->set_end(location);
        }
        //prev->dest_point->branches[1].emplace_back(point, 0);
        //point->branches[0].emplace_back(prev->dest_point, 1);

        if (next->src_point->end_loc == -1) {
            next->src_point->set_end(location);
        }
        //next->src_point->branches[1].emplace_back(point, 0);
        //point->branches[0].emplace_back(next->src_point, 1);
        // parallel and opposite
        if (is_parallel(prev_v, next_v) && inner_product(prev_v, next_v) < 0) {
            Line_2 bisector = CGAL::bisector(prev_line.opposite(), next_line);
            //auto inter0 = CGAL::intersection(bisector, prev->src_point->path()),
            //     inter1 = CGAL::intersection(bisector, next->dest_point->path());
            auto inter0 = intersection(bisector, prev->src_point->path(), line_loc[prev->corr_line]),
                 inter1 = intersection(bisector, next->dest_point->path(), line_loc[next->corr_line]);
            assert(inter0 && inter1);
            Point_2 *p0 = boost::get<Point_2>(&*inter0),
                    *p1 = boost::get<Point_2>(&*inter1);
            assert(p0 && p1);
            FT r0 = calc_loc_on_line(bisector, *p0),
               r1 = calc_loc_on_line(bisector, *p1);
            int loc_id = -1;
            if (less(r0, r1)) {
                loc_id = locations.size();
                locations.emplace_back(*p0, this->cur_time);
                point->set_end(loc_id);
                if (prev->src_point->end_loc == -1)
                    cut_half_edge(next, prev->prev, 1, loc_id, point_pool, edge_pool);
                //point->branches[0] = prev->src_point->branches[0];
            }
            else if (less(r1, r0)) {
                loc_id = locations.size();
                locations.emplace_back(*p1, this->cur_time);
                point->set_end(loc_id);
                if (next->dest_point->end_loc == -1)
                    cut_half_edge(prev, next->next, 0, loc_id, point_pool, edge_pool);
                //point->branches[0] = next->dest_point->branches[0];
            }
            else {
                loc_id = locations.size();
                locations.emplace_back(*p1, this->cur_time);
                point->set_end(loc_id);
                if (next->dest_point->end_loc == -1)
                    cut_half_point(prev->prev, next->next, loc_id, point_pool, edge_pool);
                //point->branches[0] = next->dest_point->branches[0];
            }
            //point->branches[0]->in_branches.push_back(point);
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
}