#ifndef UCS_PARTITION_SOLVER_H
#define UCS_PARTITION_SOLVER_H
#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <cmath>

#include <CGAL/approximated_offset_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/offset_polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Object.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Polygon_set_2.h>

namespace CenterLine {
	const double PI = acos(-1.0);
	template<typename K>
	struct UcsPartitionSolver{
		using FT = typename K::FT;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Point_2 = CGAL::Point_2<K>;
        using Vector_2 = CGAL::Vector_2<K>;
        using Segment_2 = CGAL::Segment_2<K>;
        using Line_2 = CGAL::Line_2<K>;
		using Polygon_set_2 = CGAL::Polygon_set_2<K>;

        using Ss = CGAL::Straight_skeleton_2<K>;
        using Halfedge_iterator = typename Ss::Halfedge_iterator;
        using Halfedge_handle = typename Ss::Halfedge_handle;
        using Vertex_handle = typename Ss::Vertex_handle;
        using Vertex_iterator = typename Ss::Vertex_iterator;

		boost::shared_ptr<Ss> _skeleton;
		//results:
		std::vector<Polygon_with_holes_2> ucs_parts;
		std::vector<Vector_2> ucs_direction;

		UcsPartitionSolver(boost::shared_ptr<Ss> skeleton) : _skeleton(skeleton) {}

		static bool is_ans(Halfedge_handle edge){
			if(!edge->is_bisector()) return false;
			FT tm0 = edge->opposite()->vertex()->time(), tm1 = edge->vertex()->time();
			if(tm0 > tm1 || edge->opposite()->vertex()->id() > edge->vertex()->id()) return false;
			if(CGAL::squared_distance(edge->vertex()->point(), edge->opposite()->vertex()->point()) < 1e-2) return false;
			Halfedge_handle l = edge->defining_contour_edge(), r = edge->opposite()->defining_contour_edge();
			Vector_2 v_l = l->vertex()->point() - l->opposite()->vertex()->point();
			Vector_2 v_r = r->vertex()->point() - r->opposite()->vertex()->point();
			//std::cout << "l = " << l->opposite()->vertex()->point() << " " << l->vertex()->point() << std::endl;
			//std::cout << "r = " << r->opposite()->vertex()->point() << " " << r->vertex()->point() << std::endl;
			//std::cout << "v_l = " << v_l << "\nv_r = " << v_r << std::endl;
			FT inner = v_l.x() * v_r.x() + v_l.y() * v_r.y();
			FT outer = v_l.x() * v_r.y() - v_l.y() * v_r.x();
			FT absolute_scale = v_l.squared_length() * v_r.squared_length();
			return (inner < 0 && outer >= 0 && inner * inner * 4 >= absolute_scale);
		}

		static int ufs_find(std::unordered_map<int, int> &f, int x){
			if(!f.count(x)) return f[x] = x;
			else if(f[x] == x) return x;
			return f[x] = ufs_find(f, f[x]);
		}

		bool operator()(double eps){
// 对于每个face，找到所有edge到所有seg的最近距离，取最近的那个seg所在的ucs
			std::unordered_map<int, double> radian_of_halfedge;
			std::unordered_map<int, Segment_2> edge_of_halfedge;
			std::vector<std::pair<double, int>> radians_and_edge;
			std::unordered_map<int, int> ucs_of_edge;
			std::unordered_map<int, std::pair<double, int>> ucs_val;
			std::unordered_map<int, Polygon_set_2> polyset_of_ucs;
			std::vector<Halfedge_handle> ans_edges;

			const double eps_seg = 1 << 20; // 精确度1/(2^20)
			for(auto it = _skeleton->halfedges_begin();it != _skeleton->halfedges_end();++it){
				if(is_ans(it)){
					ans_edges.push_back(it);
					Vector_2 v = it->vertex()->point() - it->opposite()->vertex()->point();
					if(v.x() < 0) v = -v;
					double radian = atan2(CGAL::to_double(v.y()), CGAL::to_double(v.x()));
					if(radian < 0) radian += PI / 2;
					radian = (long long)(radian * eps_seg + 0.5) / eps_seg;
					radian_of_halfedge[it->id()] = radian;
					edge_of_halfedge[it->id()] = Segment_2(it->opposite()->vertex()->point(), it->vertex()->point());
					radians_and_edge.emplace_back(radian, it->id());

					//std::cout << "is_ans: " << it->id() << std::endl;
					//std::cout << edge_of_halfedge[it->id()] << ",\tradian=" << radian << std::endl;
				}
			}
			std::sort(radians_and_edge.begin(), radians_and_edge.end());
			// calc ucs_of_edge
			for(auto e_it : ans_edges) ucs_of_edge[e_it->id()] = e_it->id();
			for(size_t i = 1;i < radians_and_edge.size();++i) if(radians_and_edge[i].first < radians_and_edge[i-1].first + eps) {
				ucs_of_edge[ufs_find(ucs_of_edge, radians_and_edge[i].second)] = ufs_find(ucs_of_edge, radians_and_edge[i-1].second);
			}
			if(radians_and_edge.front().first + PI/2 < radians_and_edge.back().first + eps){
				ucs_of_edge[ufs_find(ucs_of_edge, radians_and_edge.front().second)] = ufs_find(ucs_of_edge, radians_and_edge.back().second);
			}
			for(auto it = ucs_of_edge.begin();it != ucs_of_edge.end();++it){
				ucs_of_edge[it->first] = ufs_find(ucs_of_edge, it->first);
				//std::cout << "ucs of id " << it->first << " = " << ucs_of_edge[it->first] << " , " << radian_of_halfedge[it->first] << std::endl;
				if(!ucs_val.count(it->first)) ucs_val[it->first] = std::make_pair(0, 0);
				ucs_val[it->first].first += radian_of_halfedge[it->first];
				++ucs_val[it->first].second;
			}
			// ucs_of_edge ok

			for(auto f_it = _skeleton->faces_begin();f_it != _skeleton->faces_end();++f_it){
				FT min_dis = -1, min_center_dis = -1;
				int ucs_id = -1;
				for(auto e_it : ans_edges){
					Segment_2 seg = edge_of_halfedge[e_it->id()];
					auto f_e_it = f_it->halfedge();
					bool ok = false;
					Vector_2 face_S(0, 0);
					int face_size = 0;
					FT face_min_dis = -1, center_dis;
					do{
						Segment_2 f_e_seg(f_e_it->opposite()->vertex()->point(), f_e_it->vertex()->point());
						FT dis = CGAL::squared_distance(seg, f_e_seg);
						if(!ok || dis < face_min_dis) face_min_dis = dis, ok = true;
						face_S += (f_e_it->vertex()->point() - CGAL::ORIGIN);
						++face_size;
						f_e_it = f_e_it->next();
					}while(f_e_it != f_it->halfedge());
					center_dis = CGAL::squared_distance(seg, CGAL::ORIGIN + face_S / face_size);
					//ucs_id = ucs_of_edge[e_it->id()]
					if(ucs_id == -1 || face_min_dis < min_dis){
						ucs_id = ucs_of_edge[e_it->id()];
						min_dis = face_min_dis;
						min_center_dis = center_dis;
					}
					else if(min_dis == 0 && center_dis < min_center_dis){
						ucs_id = ucs_of_edge[e_it->id()];
						min_center_dis = center_dis;
					}
				}
				Polygon_2 poly;
				auto f_e_it = f_it->halfedge();
				do{
					poly.push_back(f_e_it->vertex()->point());
					f_e_it = f_e_it->next();
				}while(f_e_it != f_it->halfedge());
				//std::cout << "poly=" << poly << std::endl << ucs_id << std::endl;
				polyset_of_ucs[ucs_id].join(poly);
			}

			for(auto pr : polyset_of_ucs){
				std::vector<Polygon_with_holes_2> polyset;
				pr.second.polygons_with_holes(std::back_inserter(polyset));
				double val = ucs_val[pr.first].first / ucs_val[pr.first].second;
				val = (long long)(val * eps_seg + .5) / eps_seg;
				Vector_2 vect(cos(val), sin(val));
				for(auto part : polyset){
					ucs_parts.push_back(part);
					ucs_direction.push_back(vect);
				}
			}

			return true;
		}
	};
}
#endif // UCS_PARTITION_SOLVER_H