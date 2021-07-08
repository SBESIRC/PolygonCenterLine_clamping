#ifndef POLYGON_CENTERLINE_TEST
#define POLYGON_CENTERLINE_TEST
#include <CenterLineViewer.h>

#include "PolygonCenterLine.h"
#include "parse.h"
using namespace CenterLine;
struct PolygonCenterLineTest{
	using K = CGAL::Exact_predicates_exact_constructions_kernel;
    using Point_2 = K::Point_2;
    using Polygon_2 = CGAL::Polygon_2<K>;
    using Segment_2 = CGAL::Segment_2<K>;
    using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

	struct Polygon_with_holes{
        Polygon_2 poly;
        std::vector<Polygon_2> holes;
        Polygon_with_holes() : poly(), holes() {}
        void Polygon_with_holes::init(const Parse::Block &block){
			std::vector<Point_2> pts;
			for(const Parse::Point & point : block.data[0]){
				pts.emplace_back(point.x, point.y);
			}
			poly = Polygon_2(pts.begin(), pts.end()-1);
			if(CGAL::orientation_2(poly.begin(), poly.end()) == CGAL::CLOCKWISE){
				poly.reverse_orientation();
			}
			holes.clear();
			for(size_t i = 1;i < block.data.size();++i){
				pts.clear();
				for(const Parse::Point&point : block.data[i]){
					pts.emplace_back(point.x, point.y);
				}
				holes.emplace_back(pts.begin(), pts.end() - 1);
			}
			for(Polygon_2 &hole : holes) if(CGAL::orientation_2(hole.begin(), hole.end()) == CGAL::COUNTERCLOCKWISE){
				hole.reverse_orientation();
			}
		}
    };

    std::vector<Polygon_with_holes_2> spaces;
    std::vector<Polygon_with_holes_2> walls;
    PolygonCenterLineTest(std::string filename){
		std::cout << filename << std::endl;
		std::vector<Parse::Block> blocks = Parse::parse_geojson(filename);
		Polygon_with_holes poly;
		for(Parse::Block & block : blocks){
			if(block.type == "Polygon" && block.cate == "Space"){
				poly.init(block);
				size_t spaces_id = spaces.size();
				spaces.emplace_back(poly.poly, poly.holes.begin(), poly.holes.end());
			}
			else if(block.type == "Polygon" && block.cate == "Wall"){
				//std::cout << "Wall cnt" << block.data.size() << std::endl;
				poly.init(block);
				size_t walls_id = walls.size();
				walls.emplace_back(poly.poly, poly.holes.begin(), poly.holes.end());
			}
			//else if(block.type == "LineString"){
			//    Point_2 last(block.data[0][0].x, block.data[0][0].y);
			//    for(size_t i = 1;i < block.data[0].size();++i){
			//        Point_2 t(block.data[0][i].x, block.data[0][i].y);
			//        segments.emplace_back(last, t);
			//        last = t;
			//    }
			//}
		}
		std::cout << "spaces cnt=" << spaces.size()<<std::endl;
		std::cout << "walls cnt=" << walls.size()<<std::endl;
		//std::cout << "segment cnt=" << segments.size()<<std::endl;
		// gmp.dll
	}
    /// calculate centerlines and show them
    void showCenterLine(){
		std::vector<Point_2> points;
        std::vector<std::pair<int, int>> segs, sub_segs;
        std::vector<Polygon_2> PolyParts;
        int argc = 1;
        const char* argv[2]={"t2_viewer","\0"};
        QApplication app(argc, const_cast<char**>(argv));

		PolygonCenterLine centerline;
        for(size_t space_id = 0;space_id < spaces.size();++space_id){
            Polygon_with_holes_2 &space = spaces[space_id];
			points.clear();
			segs.clear();
			sub_segs.clear();
			PolyParts.clear();

			// the polygon with holes is stored in PolyParts
			PolyParts.push_back(space.outer_boundary());
			PolyParts.insert(PolyParts.end(), space.holes_begin(), space.holes_end());

			if(!centerline.calcCenterLine(space)){
				//exit(1);
			}
			// TODO: 根据centerline中的答案得到points, segs...(在getCenterLine()中)

			auto result = centerline.centerline();
			auto sub_line = centerline.sub_centerline();

			for(auto &seg : result){
				points.push_back(seg.source());
				points.push_back(seg.target());
			}
			for(auto &seg : sub_line){
				points.push_back(seg.source());
				points.push_back(seg.target());
			}
			std::sort(points.begin(), points.end());
			auto unique_end = std::unique(points.begin(), points.end());
			points.erase(unique_end, points.end());

			for(auto &seg : result){
				int a = std::lower_bound(points.begin(), points.end(), seg.source()) - points.begin();
				int b = std::lower_bound(points.begin(), points.end(), seg.target()) - points.begin();
				segs.emplace_back(a, b);
			}
			for(auto &seg : sub_line){
				int a = std::lower_bound(points.begin(), points.end(), seg.source()) - points.begin();
				int b = std::lower_bound(points.begin(), points.end(), seg.target()) - points.begin();
				sub_segs.emplace_back(a, b);
			}

            std::string title = "output " + std::to_string(space_id);
			CGAL::CenterLineViewer<Polygon_2> mainwindow(app.activeWindow(), *PolyParts.begin(), title.c_str());
			mainwindow.drawPartitions(PolyParts);
			mainwindow.drawTree(points, segs);
			mainwindow.drawTree(points, sub_segs, CGAL::Color(0, 255, 127));
			mainwindow.show();
			app.exec();
			mainwindow.saveImage(title + ".png");
        }
	}
};
#endif