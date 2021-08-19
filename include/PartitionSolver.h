#include <iostream>
#include <vector>

#include <CGAL/approximated_offset_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/offset_polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Object.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/ch_melkman.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Bbox_2.h>

namespace CenterLine {
	template<typename K>
	struct PartitionSolver{
		using FT = typename K::FT;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Point_2 = CGAL::Point_2<K>;
        using Vector_2 = CGAL::Vector_2<K>;
        using Segment_2 = CGAL::Segment_2<K>;
        using Line_2 = CGAL::Line_2<K>;

		using Gps_traits_2 = CGAL::Gps_circle_segment_traits_2<K>;
		using Offset_polygon_2 = typename Gps_traits_2::Polygon_2;
		using Offset_polygon_with_holes_2 = typename Gps_traits_2::Polygon_with_holes_2;
		using Polygon_set_2 = CGAL::General_polygon_set_2<Gps_traits_2>;
		using X_monotone_curve_2 = typename Gps_traits_2::X_monotone_curve_2;
		using Intersect_2 = typename Gps_traits_2::Intersect_2;
		using Intersection_map = typename X_monotone_curve_2::Intersection_map;
		using One_root_point_2 = typename Gps_traits_2::Point_2;
		using CoordNT = typename Gps_traits_2::Point_2::CoordNT;
		// params
		Polygon_with_holes_2 space;
		std::vector<Segment_2> segments;
		std::vector<std::pair<FT, FT>> seg_dis;
		// results
		std::vector<Offset_polygon_with_holes_2> remainders;
		std::vector<Polygon_2> convex_remainders;
		std::vector<CGAL::Bbox_2> boxes;
		std::vector<std::vector<size_t>> remainders_in_box;
		std::vector<Polygon_with_holes_2> parts, cencerline_parts;

		const double err_bound = 1e-8;

		// https://doc.cgal.org/4.3/Boolean_set_operations_2/Boolean_set_operations_2_2circle_segment_8cpp-example.html#_a4
		static Offset_polygon_2 construct_poly(const Polygon_2 &origin){
			Offset_polygon_2 res;
			for(auto it = origin.edges_begin(); it != origin.edges_end();++it) {
				X_monotone_curve_2 seg(it->source(), it->target());
				res.push_back(seg);
			}
			return res;
		}
		static Offset_polygon_with_holes_2 construct_poly(const Polygon_with_holes_2 &origin){
			std::vector<Offset_polygon_2> holes;
			for(auto it = origin.holes_begin(); it != origin.holes_end();++it) holes.push_back(construct_poly(*it));
			return Offset_polygon_with_holes_2(construct_poly(origin.outer_boundary()), holes.begin(), holes.end());
		}
		template<typename NT0, typename NT1>
		static NT1 trans_nt(const NT0 &x){ return CGAL::to_double(x); }
		static Point_2 trans_point_traits(const One_root_point_2 &point){
			auto x = point.x(), y = point.y();
			return Point_2(trans_nt<CoordNT, FT>(x), trans_nt<CoordNT, FT>(y));
		}
		static void dump_points(const Offset_polygon_2 &origin, std::vector<Point_2> &pts){
			for(auto c_it = origin.curves_begin();c_it != origin.curves_end();++c_it) {
				pts.emplace_back(trans_point_traits(c_it->source()));
			}
		}
		static FT calc_loc_on_line(const Line_2 &line, const Point_2 &point) {
            return line.b() * point.x() - line.a() * point.y();
        }

		bool calc_remained_polygons(FT R){
			Offset_polygon_with_holes_2 origin = construct_poly(space);
			Polygon_set_2 U(origin), S;
			for(auto &seg : segments){
				Polygon_2 tmp;
				tmp.push_back(seg.source()); tmp.push_back(seg.target());
				Offset_polygon_with_holes_2 p = CGAL::approximated_offset_2(tmp, R, err_bound);
				S.join(p);
			}
			remainders.clear();
			U.difference(S);
			U.polygons_with_holes(std::back_inserter(remainders));
			convex_remainders.clear();
			for(auto &remainder : remainders){
				// get convex_hull of remainder
				std::vector<Point_2> p_set;
				dump_points(remainder.outer_boundary(), p_set);
				std::vector<Point_2> convex;
				CGAL::ch_melkman(p_set.begin(), p_set.end(), std::back_inserter(convex));
				convex_remainders.push_back(Polygon_2(convex.begin(), convex.end()));
			}
			return true;
		}
		size_t merge_ch(std::vector<size_t> &indices, size_t a, size_t b){
			if(convex_remainders[indices[b]].size() > convex_remainders[indices[a]].size()) std::swap(a, b);
			std::vector<Point_2> p_set(convex_remainders[indices[a]].vertices_begin(), convex_remainders[indices[a]].vertices_end());
			p_set.insert(p_set.end(), convex_remainders[indices[b]].vertices_begin(), convex_remainders[indices[b]].vertices_end());
			std::vector<Point_2> res;
			CGAL::convex_hull_2(p_set.begin(), p_set.end(), std::back_inserter(res));
			convex_remainders[indices[a]].clear();
			convex_remainders[indices[a]].insert(convex_remainders[indices[a]].vertices_end(), res.begin(), res.end());
			indices[b] = -1;
			return indices[a];
		}
		void merge_opposite_ch(){
			std::vector<size_t> ch_id;
			std::vector<std::pair<std::pair<FT, FT>, size_t>> range_info;
			for(size_t i = 0;i < convex_remainders.size();++i) ch_id.push_back(i);
			for(int i = 0;i < segments.size();++i){
				Segment_2 &seg = segments[i];
				Line_2 &line = seg.supporting_line();
				FT max_dis = CGAL::max(seg_dis[i].first, seg_dis[i].second);
				FT l = calc_loc_on_line(line, seg.source()), r = calc_loc_on_line(line, seg.target());
				if(l > r) std::swap(l, r);

				range_info.clear();
				//std::cout << "seg = " << seg << std::endl;
				//std::cout << max_dis << std::endl;
				for(size_t j = 0;j < ch_id.size();++j){
					Polygon_2 &poly = convex_remainders[ch_id[j]];
					FT mn = r, mx = l;
					for(auto e_it = poly.edges_begin();e_it != poly.edges_end();++e_it){
						FT dis0 = CGAL::squared_distance(e_it->source(), line);
						FT dis1 = CGAL::squared_distance(e_it->target(), line);
						FT loc0 = calc_loc_on_line(line, e_it->source());
						FT loc1 = calc_loc_on_line(line, e_it->target());
						if(dis0 > dis1){
							std::swap(dis0, dis1);
							std::swap(loc0, loc1);
						}
						//std::cout << "e_it=" << *e_it << std::endl;
						//std::cout << dis0 << " " << dis1 << " " << loc0 << " " << loc1 << std::endl;
						if(dis1 <= max_dis){
							if(loc0 > loc1) swap(loc0, loc1);
							mn = CGAL::min(mn, CGAL::max(l, loc0));
							mx = CGAL::max(mx, CGAL::min(r, loc1));
						}
						// sqrt(max_dis) = x * sqrt(dis0) + (1-x) * sqrt(dis1);
						else if(dis0 < max_dis){
							double mx_d = CGAL::sqrt(CGAL::to_double(max_dis));
							double d0 = CGAL::sqrt(CGAL::to_double(dis0)), d1 = CGAL::sqrt(CGAL::to_double(dis1));
							FT x = (d1 - mx_d) * (d0 + d1) / (dis1 - dis0);
							FT loc = x * loc0 + (1 - x) * loc1;
							if(loc0 > loc1) swap(loc0, loc1);
							loc1 = loc;
							mn = CGAL::min(mn, CGAL::max(l, loc0));
							mx = CGAL::max(mx, CGAL::min(r, loc1));
						}
					}
					if(mn < mx) range_info.push_back(std::make_pair(std::make_pair(mn, mx), j));
				}
				std::sort(range_info.begin(), range_info.end());
				std::cout << "range_info = " << std::endl;
				for(auto info : range_info){ std::cout << info.first.first << " "  << info.first.second << "\n" << convex_remainders[ch_id[info.second]]; }
				if(range_info.size() < 2) continue;
				FT cur_r = range_info[0].first.second;
				size_t cur_ch = 0;
				for(size_t i = 1;i < range_info.size();++i){
					if(range_info[i].first.first < cur_r){
						cur_ch = merge_ch(ch_id, cur_ch, i);
					}
					cur_r = CGAL::max(cur_r, range_info[i].first.second);
				}
				for(size_t i = 0;i < ch_id.size();++i) if(ch_id[i] == -1){
					std::swap(ch_id[i], ch_id[ch_id.size()-1]);
					ch_id.pop_back();
				}
			}
			std::vector<Polygon_2> ans;
			for(size_t i = 0;i < ch_id.size();++i) ans.push_back(convex_remainders[ch_id[i]]);
			convex_remainders.swap(ans);
		}

		void merge_boxes(){
			boxes.clear();
			for(size_t i = 0;i < convex_remainders.size();++i){
				boxes.push_back(convex_remainders[i].bbox());
				std::vector<size_t> tmp; tmp.push_back(i);
				remainders_in_box.push_back(tmp);
			}

			for(size_t i = boxes.size() - 1;i < boxes.size();--i){
				for(size_t j = i + 1;j < boxes.size();) {
					if(do_overlap(boxes[j], boxes[i])){
						boxes[i] += boxes[j];
						if(remainders_in_box[j].size() > remainders_in_box[i].size()) remainders_in_box[i].swap(remainders_in_box[j]);
						remainders_in_box[i].insert(remainders_in_box[i].begin(), remainders_in_box[j].begin(), remainders_in_box[j].end());
						std::swap(boxes[j], boxes[boxes.size() - 1]);
						remainders_in_box[j].swap(remainders_in_box[boxes.size() - 1]);
						boxes.pop_back();
						remainders_in_box.pop_back();
					}
					else ++j;
				}
			}
		}

		bool solve(FT R){
			calc_remained_polygons(R);
			std::copy(remainders.begin(), remainders.end(), std::ostream_iterator<Offset_polygon_with_holes_2>(std::cout, "\n"));
			merge_opposite_ch();
			merge_boxes();
			for(size_t i = 0;i < boxes.size();++i){
				auto &box = boxes[i];
				std::vector<size_t> &rems = remainders_in_box[i];
				Polygon_2 rect;
				rect.push_back(Point_2(box.xmin(), box.ymin()));
				rect.push_back(Point_2(box.xmax(), box.ymin()));
				rect.push_back(Point_2(box.xmax(), box.ymax()));
				rect.push_back(Point_2(box.xmin(), box.ymax()));
				std::vector<Polygon_with_holes_2> all_parts;
				CGAL::intersection(space, rect, std::back_inserter(all_parts));
				for(auto &part : all_parts){
					bool ok = false;
					for(int rem : rems){
						Point_2 p = *convex_remainders[rem].vertices_begin();
						if(CGAL::oriented_side(p, part) != CGAL::ON_NEGATIVE_SIDE){
							ok = true;
							break;
						}
					}
					if(ok) this->parts.push_back(part);
				}
			}
			CGAL::Polygon_set_2<K> U(space), S;
			for(auto &part : this->parts) S.join(part);
			U.difference(S);
			U.polygons_with_holes(std::back_inserter(this->cencerline_parts));
			return true;
		}
	}; // struct ParititonSolver
}
