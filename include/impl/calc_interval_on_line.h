#pragma once
#include "CenterLineSolver.h"
namespace CenterLineSolver{
    template <typename K, class Poly_with_holes, class Poly>
    inline bool CenterLineSolver<K, Poly_with_holes, Poly>::calc_interval_on_line(
        const Line_2 &bisector, const EdgeData &edge, const EdgeData &op_edge, Interval &interval,
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
        Ray_2 src_path = src->path(), dest_path = dest->path();
        if (is_debugging) {
            std::cout << "bisector_v = " << bisector.to_vector() << std::endl;
            std::cout << tmp_inner_product << " " << tmp_outer_product << std::endl;
        }
        // parallel
        if (parallel) {
            //auto a = CGAL::intersection(bisector, src_path);
            //auto b = CGAL::intersection(bisector, dest_path);
            auto a = intersection(bisector, src_path, line_loc[edge.corr_line]);
            auto b = intersection(bisector, dest_path, line_loc[edge.corr_line]);
            if (!a || !b)
                return false;
            if (tmp_inner_product < 0)
                std::swap(a, b);
            if (!(tmp_p = boost::get<Point_2>(&*a)) || !(tmp_q = boost::get<Point_2>(&*b)))
                return false;
            if (is_debugging)
                std::cout << "tmp_p, tmp_q" << *tmp_p << " " << *tmp_q << std::endl;
            FT loc_p = calc_loc_on_line(bisector, *tmp_p),
               loc_q = calc_loc_on_line(bisector, *tmp_q);
            if (is_debugging)
                std::cout << "loc_p, loc_q = " << loc_p << " " << loc_q << std::endl;
            // if(loc_p >= loc_q)
            if (!less(loc_p, loc_q))
                return false;
            interval.set(loc_p, loc_q, *tmp_p, *tmp_q, is_endpoint);
            return true;
        }
        else { // tmp_outer_product > 0
            //auto a = CGAL::intersection(bisector, src_path);
            //auto b = CGAL::intersection(bisector, dest_path);
            auto a = intersection(bisector, src_path, line_loc[edge.corr_line]);
            auto b = intersection(bisector, dest_path, line_loc[edge.corr_line]);
            Segment_2 new_seg = edge.new_loc(line_loc, this->cur_time);
            if (is_debugging)
                std::cout << "new_seg" << new_seg << std::endl;
            Point_2 new_src = new_seg.source(), new_dest = new_seg.target();
            Point_2 begin(0, 0), end(0, 0);
            auto o = CGAL::intersection(bisector, line_loc[edge.corr_line]);
            assert(bool(o));
            Point_2 origin = *boost::get<Point_2>(&*o);

            if (is_debugging) {
                std::cout << "bisector= " << bisector << std::endl;
                std::cout << new_src << " " << new_dest << " " << origin << std::endl;
            }
            FT src_loc = calc_loc_on_line(bisector, new_src);
            FT dest_loc = calc_loc_on_line(bisector, new_dest);
            FT origin_loc = calc_loc_on_line(bisector, origin);
            // ���߶��˶����򱳶�ֱ�ߣ������Ȼ���ཻ��˵���˶�����ƫת�������ߣ�˵�����ڱ�ƫת������λ�ߣ��������ڱ߲����Ƕ�λ�ߵ����ڱߣ���ô��Ȼ������ڱ߸��ܴ�����һ�����ײ�����Բ��ÿ��ǣ�����������������ͬ�࣬���������ߵ�һ��һ���ǿ�����һ�ߵ�һ���ȽӴ����ߣ�����Ӧ�����ڱ�һ�����ӽ��������ߵıߣ��ҽǶȸ��ӽ���
            if (is_debugging)
                std::cout << src_loc << " " << dest_loc << " " << origin_loc << std::endl;
            // if(src_loc <= origin_loc && dest_loc <= origin_loc)
            if (equal(src_loc, dest_loc))
                return false; // �Ѿ��˻�Ϊһ����
            if (!less(origin_loc, src_loc) && !less(origin_loc, dest_loc))
                return false;

            FT min_loc = CGAL::min(src_loc, dest_loc), max_loc = CGAL::max(src_loc, dest_loc);
            if (less(origin_loc, min_loc) || less(max_loc, origin_loc)) {
                if (tmp_inner_product < 0) {
                    std::swap(a, b);
                    std::swap(src, dest);
                    std::swap(src_path, dest_path);
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
                else if ((tmp_p = boost::get<Point_2>(&*a)) && in_box(*tmp_p)) {
                    begin = *tmp_p;
                }
                else { // ���������ཻ
                    //std::cout << "do not collide" << std::endl;
                    return false;
                }
                if (is_debugging)
                    std::cout << "begin = " << begin << std::endl;

                FT loc_p = calc_loc_on_line(bisector, begin), loc_q;

                if (!b)
                    interval.set(loc_p, begin, is_endpoint);
                else if ((tmp_q = boost::get<Point_2>(&*b)) && in_box(*tmp_q)) {
                    end = *tmp_q;
                    //std::cout << end << std::endl;
                    loc_q = calc_loc_on_line(bisector, end);
                    //std::cout << "loc_q = " << loc_q << std::endl;
                    interval.set(loc_p, loc_q, begin, end, is_endpoint);
                }
                else
                    interval.set(loc_p, begin, is_endpoint); // �߶κ�����û�н��㣬end��������������,��˴˴�˵�������㹻Զ
                return true;
            }
            // �߶�һ��ʼ���������ཻ
            // ��Ҫ�ų�һ��ʼ������Զ������
            else {
                tmp_p = &origin;
                // std::cout << tmp_p << std::endl;
                // if(src_loc == origin_loc)
                if (equal(src_loc, origin_loc))
                    is_endpoint = 1;
                // else if(dest_loc == origin_loc)
                else if (equal(dest_loc, origin_loc))
                    is_endpoint = 2;

                FT loc_p = origin_loc, loc_q;
                begin = *tmp_p;

                if (edge.prev == op_edge.next)
                    a = boost::none;
                if (edge.next == op_edge.prev)
                    b = boost::none; // ���һ����bisector�ཻ���������й����ڱߣ������ǵ���Ӧ�˵�һ���غ�
                if (!a && !b) {
                    interval.set(loc_p, begin, is_endpoint);
                }
                else if (a && b) {
                    // ��ǰis_endpointֻ�ж���origin_loc�Ƿ����src_loc��dest_loc
                    if (is_endpoint == 1) {
                        if ((tmp_r = boost::get<Ray_2>(&*a))) {
                            is_endpoint = -1;
                        }
                        // Զ��bisector
                        // else if(outer_product(src->speed(), bisector.to_vector()) > 0){
                        else if (less(0, outer_product(src->speed(), bisector.to_vector()),
                                      FT(eps) * calc_line_scale(bisector))) {
                            return false;
                        }
                        tmp_q = boost::get<Point_2>(&*b);
                        if (!tmp_q)
                            throw("calc_interval_on_line, tmp_p, a&&b,1: b no intersection");
                        end = *tmp_q;
                        if (in_box(end)) {
                            loc_q = calc_loc_on_line(bisector, end);
                            interval.set(loc_p, loc_q, begin, end, is_endpoint);
                        }
                        else
                            interval.set(loc_p, begin, is_endpoint);
                    }
                    else if (is_endpoint == 2) {
                        if ((tmp_r = boost::get<Ray_2>(&*b))) {
                            is_endpoint = -2;
                        }
                        // Զ��bisector
                        // else if(outer_product(dest->speed(), bisector.to_vector()) < 0){
                        else if (less(outer_product(dest->speed(), bisector.to_vector()), 0,
                                      FT(eps) * calc_line_scale(bisector))) {
                            return false;
                        }
                        tmp_q = boost::get<Point_2>(&*a);
                        if (!tmp_q)
                            throw("calc_interval_on_line, tmp_p, a&&b,2: b no intersection");
                        end = *tmp_q;
                        if (in_box(end)) {
                            loc_q = calc_loc_on_line(bisector, end);
                            interval.set(loc_p, loc_q, begin, end, is_endpoint);
                        }
                        else
                            interval.set(loc_p, begin, is_endpoint);
                    }
                    else {
                        // ��ʱ����һ�����ڲ�
                        Point_2 *tmp_src = boost::get<Point_2>(&*a);
                        Point_2 *tmp_dest = boost::get<Point_2>(&*b);
                        assert(tmp_src && tmp_dest); // �߶���bisector�������߶��ڲ�
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
                        if (in_box(end))
                            interval.set(loc_p, loc_q, begin, end, is_endpoint);
                        else
                            interval.set(loc_p, begin, is_endpoint);
                    }
                }
                // a��·���н����b��·���н���
                else {
                    if (a && (tmp_q = boost::get<Point_2>(&*a))) {
                        end = *tmp_q;
                    }
                    else if (b && (tmp_q = boost::get<Point_2>(&*b))) {
                        end = *tmp_q;
                    }
                    if (is_debugging) {
                        std::cout << "origin = " << origin << std::endl;
                        std::cout << "begin = " << begin << std::endl;
                        std::cout << "a = " << a << std::endl;
                        std::cout << "b = " << b << std::endl;
                    }
                    if (in_box(end)) {
                        loc_q = calc_loc_on_line(bisector, end);
                        interval.set(loc_p, loc_q, begin, end, is_endpoint);
                    }
                    else
                        interval.set(loc_p, begin, is_endpoint);
                }
            }
        }
        return true;
    }

}