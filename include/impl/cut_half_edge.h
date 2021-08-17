#pragma once
#include "CenterLineSolver.h"
namespace CenterLineSolver{
    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::cut_half_edge(
        EdgeData *base, EdgeData *edge, bool type, int location,
        std::vector<PointData *> &point_pool, std::vector<EdgeData *> &edge_pool)
    {
        //std::cout << "base->prev=" << base->prev->corr_line << std::endl;
        //std::cout << "base->next=" << base->next->corr_line << std::endl;
        //std::cout << "edge->prev=" << edge->prev->corr_line << std::endl;
        //std::cout << "edge->next=" << edge->next->corr_line << std::endl;
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
                inter_seg = new PointData(locations, location, base->to_vector(),
                                          edge->to_vector());
                point_pool.push_back(inter_seg);
                if (edge->src_point->end_loc == -1) {
                    edge->src_point->set_end(location);
                }
                //edge->src_point->branches[1].emplace_back(inter_seg, 0);
                //inter_seg->branches[0].emplace_back(edge->src_point, 1);

                //auto inter0 = CGAL::intersection(bisector, base->src_point->path()),
                //     inter1 = CGAL::intersection(bisector, edge->dest_point->path());
                auto inter0 = intersection(bisector, base->src_point->path(), line_loc[base->corr_line]),
                     inter1 = intersection(bisector, edge->dest_point->path(), line_loc[edge->corr_line]);
                if (!inter0 || !inter1) {
                    std::cout << "bisector = " << bisector << std::endl;
                    std::cout << "bisector_v = " << bisector.to_vector() << std::endl;
                    std::cout << "base = " << base->new_loc(this->line_loc, this->cur_time) << std::endl;
                    std::cout << "edge = " << edge->new_loc(this->line_loc, this->cur_time) << std::endl;
                    std::cout << "edge->dest = " << edge->dest_point->path() << std::endl;
                    std::cout << "inter=" << inter0 << std::endl
                              << inter1 << std::endl;
                    throw("cut_half_edge is_parallel,type=0: no intersection");
                }
                Point_2 *p0 = boost::get<Point_2>(&*inter0),
                        *p1 = boost::get<Point_2>(&*inter1); // ? �����forced_return
                if (!p0 || !p1) {
                    throw("cut_half_edge is_parallel,type=0: intersection not Point_2");
                }
                FT l0 = calc_loc_on_line(bisector, *p0),
                   l1 = calc_loc_on_line(bisector, *p1);
                int loc_id;
                //std::cout << "p0=" << *p0 << std::endl;
                //std::cout << "p1=" << *p1 << std::endl;
                if (less(l0, l1)) {
                    loc_id = locations.size();
                    locations.push_back(Location(*p1, this->cur_time));
                    inter_seg->set_end(loc_id);
                    if (edge->dest_point->end_loc == -1)
                        cut_half_edge(base, edge->next, 0, loc_id, point_pool, edge_pool);
                    // TODO: �Ƿ��п���ѭ������?
                    //inter_seg->branches[0] = edge->dest_point->branches[0];
                }
                else if (less(l1, l0)) {
                    loc_id = locations.size();
                    locations.push_back(Location(*p0, this->cur_time));
                    inter_seg->set_end(loc_id);
                    if (base->src_point->end_loc == -1)
                        cut_half_edge(edge, base->prev, 1, loc_id, point_pool, edge_pool);
                    //inter_seg->branches[0] = base->src_point->branches[0];
                }
                else {
                    loc_id = locations.size();
                    locations.push_back(Location(*p0, this->cur_time));
                    inter_seg->set_end(loc_id);
                    if (base->src_point->end_loc == -1)
                        cut_half_point(base->prev, edge->next, loc_id, point_pool, edge_pool);
                    //inter_seg->branches[0] = base->src_point->branches[0];
                }
                //inter_seg->branches[0]->in_branches.push_back(inter_seg);
            }
            else { // type == 1
                inter_seg = new PointData(locations, location, edge->to_vector(),
                                          base->to_vector());
                point_pool.push_back(inter_seg);
                if (edge->dest_point->end_loc == -1) {
                    edge->dest_point->set_end(location);
                }
                //edge->dest_point->branches[1].emplace_back(inter_seg, 0);
                //inter_seg->branches[0].emplace_back(edge->dest_point, 1);

                //auto inter0 = CGAL::intersection(bisector, base->dest_point->path()),
                //     inter1 = CGAL::intersection(bisector, edge->src_point->path());
                auto inter0 = intersection(bisector, base->dest_point->path(), line_loc[base->corr_line]),
                     inter1 = intersection(bisector, edge->src_point->path(), line_loc[edge->corr_line]);
                assert(inter0 && inter1);
                if(!inter0) throw("cut_half_edge.h, inter0 is null");
                if(!inter1) throw("cut_half_edge.h, inter1 is null");
                Point_2 *p0 = boost::get<Point_2>(&*inter0),
                        *p1 = boost::get<Point_2>(&*inter1);
                int loc_id = -1;
                assert(p0 && p1);
                FT r0 = calc_loc_on_line(bisector, *p0),
                   r1 = calc_loc_on_line(bisector, *p1);
                if (less(r0, r1)) {
                    loc_id = locations.size();
                    locations.push_back(Location(*p0, this->cur_time));
                    inter_seg->set_end(loc_id);
                    if (base->dest_point->end_loc == -1)
                        cut_half_edge(edge, base->next, 0, loc_id, point_pool, edge_pool);
                    //inter_seg->branches[0] = base->dest_point->branches[0];
                }
                else if (less(r1, r0)) {
                    loc_id = locations.size();
                    locations.push_back(Location(*p1, this->cur_time));
                    inter_seg->set_end(loc_id);
                    if (edge->src_point->end_loc == -1)
                        cut_half_edge(base, edge->prev, 1, loc_id, point_pool, edge_pool);
                    //inter_seg->branches[0] = edge->src_point->branches[0];
                }
                else {
                    loc_id = locations.size();
                    locations.push_back(Location(*p1, this->cur_time));
                    inter_seg->set_end(loc_id);
                    if (edge->src_point->end_loc == -1)
                        cut_half_point(edge->prev, base->next, loc_id, point_pool, edge_pool);
                    //inter_seg->branches[0] = edge->src_point->branches[0];
                }
            }
            return;
        }

        if (type == 0) {
            PointData *point = new PointData(locations, location, base->to_vector(),
                                             edge->to_vector());
            point_pool.push_back(point);
            if (edge->src_point->end_loc == -1) {
                edge->src_point->set_end(location); // ��ֹһ����
            }
            //edge->src_point->branches[1].emplace_back(point, 0);
            //point->branches[0].emplace_back(edge->src_point, 1);

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
            PointData *point = new PointData(locations, location, edge->to_vector(),
                                             base->to_vector());
            point_pool.push_back(point);

            if (edge->dest_point->end_loc == -1) {
                edge->dest_point->set_end(location); // ��ֹһ����
            }
            //edge->dest_point->branches[1].emplace_back(point, 0);
            //point->branches[0].emplace_back(edge->dest_point, 1);

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
}