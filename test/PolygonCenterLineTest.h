#ifndef POLYGON_CENTERLINE_TEST
#define POLYGON_CENTERLINE_TEST
#include <string>

#include <CenterLineViewer.h>

#include "PolygonCenterLine.h"
using namespace CenterLine;
struct PolygonCenterLineTest{
	using K = CGAL::Exact_predicates_exact_constructions_kernel;
    using Point_2 = K::Point_2;
    using Polygon_2 = CGAL::Polygon_2<K>;
    using Segment_2 = CGAL::Segment_2<K>;
    using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

    //std::vector<Polygon_with_holes_2> spaces;
    //std::vector<Polygon_with_holes_2> walls;
	std::string geojson, filename;
    PolygonCenterLineTest(std::string fname){
		std::cout << fname << std::endl;
		std::ifstream in(fname);
		std::istreambuf_iterator<char>  beg(in), end;
		geojson = std::string(beg,  end);
		filename = fname;
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
        {
            Polygon_with_holes_2 &space = PolygonCenterLine::geojson_to_poly(geojson);
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

            std::string title = "output " + filename;
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