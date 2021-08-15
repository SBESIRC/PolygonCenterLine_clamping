#ifndef POLYGON_CENTERLINE_TEST
#define POLYGON_CENTERLINE_TEST
#include <string>

#include <CenterLineViewer.h>
#include <PartitionViewer.h>

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
	PolygonCenterLine centerline;

    PolygonCenterLineTest(std::string fname, int id){
		std::cout << fname << std::endl;
		std::ifstream in(fname);
		if(!in.is_open()) throw("input file not found");
		std::istreambuf_iterator<char>  beg(in), end;
		geojson = std::string(beg,  end);
		filename = std::to_string(id);
	}
    /// calculate centerlines and show them
    void showCenterLine(){
		std::vector<Point_2> points;
        std::vector<std::pair<int, int>> segs, sub_segs;
        std::vector<Polygon_2> PolyParts;
        int argc = 1;
        const char* argv[2]={"t2_viewer","\0"};
        QApplication app(argc, const_cast<char**>(argv));

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
			// TODO: ����centerline�еĴ𰸵õ�points, segs...(��getCenterLine()��)

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

            std::string title = "output-" + filename;
			CGAL::CenterLineViewer<Polygon_2> mainwindow(app.activeWindow(), *PolyParts.begin(), title.c_str());
			mainwindow.drawPartitions(PolyParts);
			mainwindow.drawTree(points, segs);
			mainwindow.drawTree(points, sub_segs, CGAL::Color(0, 255, 127));
			mainwindow.show();
			app.exec();
			mainwindow.saveImage("test_results/" + title + ".png");
        }
	}
	void showPartition(double R = 0){
        int argc = 1;
        const char* argv[2]={"t2_viewer","\0"};
		std::vector<double> dis_set;
		if(R == 0) for(auto dis : centerline.point_contour_distance){
			if(dis > 10) dis_set.push_back(CGAL::to_double(dis - 10));
			dis_set.push_back(CGAL::to_double(dis + 10));
		}
		else dis_set.push_back(R);
		std::sort(dis_set.begin(), dis_set.end());
		for(double dis : dis_set){
			QApplication app(argc, const_cast<char**>(argv));
			std::string title = "partition-" + filename + "-" + std::to_string(dis);
			centerline.calcPartition(dis);

			CGAL::PartitionViewer<K> mainwindow(app.activeWindow(), centerline.relative_poly, title.c_str());
			CGAL::Color c(75, 160, 255, 127);
			for(const auto &remainder : centerline.rel_remainders){
				std::cout << "remainder" << std::endl << remainder << std::endl;
				mainwindow.drawPoly(remainder, c);
			}
			mainwindow.drawPoly(centerline.relative_segments, c, CGAL::Color(255, 255, 0, 127));

			std::string new_title = "partition-" + filename + "-" + std::to_string(dis) + "result";
			CGAL::PartitionViewer<K> new_window(app.activeWindow(), centerline.relative_poly, new_title.c_str());
			c = CGAL::Color(160, 75, 160, 127);
			//for(const auto &remainder : centerline.rel_convex_remainders){
			for(const auto &remainder : centerline.rel_parts){
				new_window.drawPoly(remainder, c);
			}
			new_window.drawPoly(centerline.relative_segments, c, CGAL::Color(255, 255, 0, 127));

			mainwindow.show();
			new_window.show();
			app.exec();
			mainwindow.saveImage("partition_results/" + title + ".png");
			new_window.saveImage("partition_results/" + new_title + ".png");
		}
	}
};
#endif