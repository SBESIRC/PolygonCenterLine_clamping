#include "CenterLineSolver.h"
namespace CenterLineSolver{
    template <typename K, class Poly_with_holes, class Poly>
    inline void CenterLineSolver<K, Poly_with_holes, Poly>::process_event(const Event &event){
        FT last_time = 0;
        std::vector<PointData *> point_pool;
        std::vector<EdgeData *> edge_pool;
		this->cur_time = event.time;
		FT exact_time = CGAL::sqrt(event.time);
		for (size_t i = 0; i < line_data_pool.size(); ++i) {
			line_loc[i] = Line_2(line_data_pool[i].point(0) + line_speed[i] * exact_time,
									line_data_pool[i].to_vector());
		}
		std::cout << "time = sqrt(" << event.time << ") = " << exact_time
					<< std::endl;

		EdgeData *event_a = event.a, *event_b = event.b;
		/* 用于调试弧形墙
		if(event_b->corr_line == 322 || event_a->corr_line == 322){
			int l_a = event_a->corr_line, l_b = event_b->corr_line;
			Line_2 line_a = line_data_pool[event_a->corr_line], line_b = line_data_pool[event_b->corr_line];
			std::cout << "a = " << l_a << std::endl;
			start_debug();
			std::cout << "point="<<line_a.point(0) << " "<< std::endl;
			std::cout << "point="<<line_b.point(0) << std::endl;
			std::cout << "vect="<<line_a.to_vector() << " "<< std::endl;
			std::cout << "vect="<<line_b.to_vector() << " "<< std::endl;
			std::cout << "speed="<<line_speed[l_a] << std::endl;
			std::cout << "speed="<<line_speed[l_b] << std::endl;
			std::cout << "line_loc="<<line_loc[event_a->corr_line] << std::endl
					<< "line_loc="<<line_loc[event_b->corr_line] << std::endl;
			std::cout << "src = " << event_a->src_point->path() << std::endl;
			std::cout << "dest = " << event_a->dest_point->path() << std::endl;
			std::cout << "src = " << event_b->src_point->path() << std::endl;
			std::cout << "dest = " << event_b->dest_point->path() << std::endl;
			std::cout << event_a->prev->corr_line << " " << line_data_pool[event_a->prev->corr_line].to_vector() << std::endl;
			std::cout << line_data_pool[event_a->corr_line].to_vector() << std::endl;
			std::cout << event_a->next->corr_line << " " << line_data_pool[event_a->next->corr_line].to_vector() << std::endl;

			std::cout << event_b->prev->corr_line << " " << line_data_pool[event_b->prev->corr_line].to_vector() << std::endl;
			std::cout << event_b->next->corr_line << " " << line_data_pool[event_b->next->corr_line].to_vector() << std::endl;
		}
		*/
		Segment_2 seg0 = event_a->new_loc(line_loc, this->cur_time),
					seg1 = event_b->new_loc(line_loc, this->cur_time);
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
			if (it->src_point->end_loc == -1) {
				Line_2 line[3] = {line_loc[it->corr_line], line_loc[it->next->corr_line], line_loc[it->prev->corr_line]};
				auto inter0 = CGAL::intersection(line[1], line[2]);
				auto inter1 = CGAL::intersection(line[2], line[0]);
				auto inter2 = CGAL::intersection(line[0], line[1]);
				if (!inter0 || !inter1 || !inter2)
					throw("triangle: no intersection");
				Point_2 *pts[3] = {boost::get<Point_2>(&*inter0), boost::get<Point_2>(&*inter1), boost::get<Point_2>(&*inter2)};
				PointData *data[3] = {it->next->dest_point, it->src_point, it->dest_point};
				if (!pts[0] || !pts[1] || !pts[2])
					throw("triangle: intersection not Point_2");
				Point_2 barycenter = CGAL::ORIGIN + (((*pts[0] - CGAL::ORIGIN) + (*pts[1] - CGAL::ORIGIN) + (*pts[2] - CGAL::ORIGIN)) / 3);
				FT end_time = CGAL::squared_distance(barycenter, line[0]);
				int loc_id = locations.size();
				locations.emplace_back(barycenter, end_time);
				data[0]->set_end(loc_id);
				data[1]->set_end(loc_id);
				data[2]->set_end(loc_id);
			}
		}
		else if (event.type == EventType::PARALLEL) {
			/* 输出公共部分 TODO: 改成创建点 (branches) */
			Point_2 point_l, point_r;
			Vector_2 vector_l, vector_r;
			if (l0 < l1)
				point_l = seg1.target();
			else
				point_l = seg0.source();
			if (r0 < r1)
				point_r = seg0.target();
			else
				point_r = seg1.source();
			if (outer_product(event_a->to_vector(), event_b->to_vector()) > 0) {
				vector_l = event_a->to_vector();
				vector_r = event_b->to_vector();
			}
			else {
				vector_l = event_b->to_vector();
				vector_r = event_a->to_vector();
			}
			int start_loc = locations.size(), end_loc = start_loc + 1;
			locations.emplace_back(point_l, event.time);
			locations.emplace_back(point_r, event.time);
			PointData *inter_seg = new PointData(locations, start_loc, vector_l, vector_r);
			inter_seg->set_end(end_loc);
			point_pool.push_back(inter_seg);

			if (event_a->prev != event_b->next) {
				// if(l0 < l1)
				if (less(l0, l1)) {
					cut_half_edge(event_a, event_b->next, 0, start_loc, point_pool, edge_pool);
				}
				// else if(l0 > l1)
				else if (less(l1, l0)) {
					cut_half_edge(event_b, event_a->prev, 1, start_loc, point_pool, edge_pool);
				}
				else {
					cut_half_point(event_a->prev, event_b->next, start_loc, point_pool, edge_pool);
				}
				/* if(seg0.source()和seg1.target()相邻){ 直接跳过 }
else if(l0 < l1){ 用seg1.target()的位置建点，切割seg0，重连 }
else if(l0 > l1){ 用seg0.source()的位置建点，切割seg1，重连 }*/
			}
			else {
				event_a->src_point->set_end(start_loc);
				event_b->dest_point->set_end(start_loc);
				// branches
			}
			if (event_a->next != event_b->prev) {
				// if(r0 < r1)
				if (less(r0, r1))
					cut_half_edge(event_b, event_a->next, 0, end_loc, point_pool, edge_pool);
				// else if(r0 > r1)
				else if (less(r1, r0))
					cut_half_edge(event_a, event_b->prev, 1, end_loc, point_pool, edge_pool);
				else {
					cut_half_point(event_b->prev, event_a->next, end_loc, point_pool, edge_pool);
				}
				/* if(seg0.source()和seg1.target()相邻){ 直接跳过 }
if(r0 < r1){ 用seg0.target()的位置建点，切割seg1，重连 }
else if(r0 > r1){ 用seg1.source()的位置建点，切割seg0，重连 } */
			}
			else {
				event_a->dest_point->set_end(end_loc);
				event_b->src_point->set_end(end_loc);
				//branches;
			}
			//std::cout << "parallel: a=" << event_a->src_point->start_loc << " "
			//        << event_a->src_point->end_loc << std::endl;
			//std::cout << "parallel: a=" << event_a->dest_point->start_loc << " "
			//        << event_a->dest_point->end_loc << std::endl;
			//std::cout << "parallel: b=" << event_b->src_point->start_loc << " "
			//        << event_b->src_point->end_loc << std::endl;
			//std::cout << "parallel: b=" << event_b->dest_point->start_loc << " "
			//        << event_b->dest_point->end_loc << std::endl;
		}
		else {
			if (l0 > r0) {
				std::swap(seg0, seg1);
				std::swap(l0, r1);
				std::swap(l1, r0);
				std::swap(event_a, event_b);
			}
			int loc_id = locations.size();
			if (event_a->prev == event_b->next) {
				locations.emplace_back(seg0.source(), this->cur_time);
				cut_half_point(event_b, event_a, loc_id, point_pool, edge_pool);
			}
			// else if(l0 < l1){
			else if (less(l0, l1)) {
				locations.emplace_back(seg1.target(), this->cur_time);
				cut_edge(event_a, event_b, event_b->next, loc_id, point_pool, edge_pool);
			}
			// else if(l0 > l1){
			else if (less(l1, l0)) {
				locations.emplace_back(seg0.source(), this->cur_time);
				cut_edge(event_b, event_a->prev, event_a, loc_id, point_pool, edge_pool);
			}
			else {
				locations.emplace_back(CGAL::ORIGIN + (seg0.source() - seg1.target()) / 2, this->cur_time);
				cut_half_point(event_b, event_a, loc_id, point_pool, edge_pool);
				cut_half_point(event_a->prev, event_b->next, loc_id, point_pool, edge_pool);
			}
		}

		for (EdgeData *it : edge_pool) {
			append_collision_events(it);
		}
		edge_data_pool.insert(edge_data_pool.end(), edge_pool.begin(), edge_pool.end());
		point_data_pool.insert(point_data_pool.end(), point_pool.begin(), point_pool.end());
	}
}
