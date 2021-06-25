#ifndef POLYGON_CENTERLINE_H
#define POLYGON_CENTERLINE_H
#include <iostream>
#include "parse.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CenterLineViewer.h>
#include <string>
#include <getCenterLine.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <convertKernel.h>
#include <vector>
#include <CGAL/Gmpfr.h>

namespace CenterLine_ns{

class PolygonCenterLine{
    //https://github.com/CGAL/cgal/issues/1873
    using K = CGAL::Exact_predicates_exact_constructions_kernel;
    using Point_2 = K::Point_2;
    using Polygon_2 = CGAL::Polygon_2<K>;
    using Segment_2 = CGAL::Segment_2<K>;
    using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
    /// internal general polygon structure.
    struct Polygon_with_holes{
        Polygon_2 poly;
        std::vector<Polygon_2> holes;
        Polygon_with_holes() : poly(), holes() {}
        void Polygon_with_holes::init(const Parse::Block &block);
    };

    std::vector<Polygon_with_holes_2> spaces;
    std::vector<Polygon_with_holes_2> walls;
    std::vector<Segment_2> segments;
    public:
    PolygonCenterLine(std::string filename);
    /// calculate centerlines and show them
    void showCenterLine();
};
} // namespace CenterLine_ns
#endif // POLYGON_CENTERLINE_H

#ifdef PolygonCenterLine_Implementation
namespace CenterLine_ns{
    void PolygonCenterLine::Polygon_with_holes::init(const Parse::Block &block){
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
    } // Polygon_with_holes::init(const Parse::Block);

    PolygonCenterLine::PolygonCenterLine(std::string filename) : spaces(), walls() {
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
    } // PolygonCenterLine::PolygonCenterLine(std::string filename)

    void PolygonCenterLine::showCenterLine(){
        std::vector<Point_2> points;
        std::vector<std::pair<int, int>> segs, aug_segs;
        int id = 0;
        std::vector<Polygon_2> PolyParts;
        int argc = 1;
        const char* argv[2]={"t2_viewer","\0"};
        QApplication app(argc, const_cast<char**>(argv));

        //using K1 = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
        //using K1 = CGAL::Simple_cartesian<CGAL::Lazy_exact_nt<CORE::Expr>>;
        //using K1 = CGAL::Simple_cartesian<double>;
        using K1 = Fixed_kernel;
        // Gmpfr的计算精度
        CGAL::Gmpfr::set_default_precision(256);
		CGAL::set_pretty_mode(std::cout);
        //KernelConverter<CGAL::Epeck, Epect_sqrt, Lazy_gmpq_to_Expr_converter> to_sqrt_kernel;
        //KernelConverter<CGAL::Epeck, typename K1, Lazy_gmpq_to_double_converter> to_double_kernel;
        KernelConverter<CGAL::Epeck, typename K1, Lazy_gmpq_to_gmpfr_converter> kernel_converter;
        //for(Polygon_2 & poly : polygons){
        for(size_t space_id = 0;space_id < spaces.size();++space_id){
            Polygon_with_holes_2 &space = spaces[space_id];
            points.clear();
            segs.clear();
            aug_segs.clear();
            PolyParts.clear();
            CGAL::Polygon_with_holes_2<K1> new_space = kernel_converter.convert(space);
            std::vector<K1::Point_2> new_points, aug_points;
            std::vector<CGAL::Polygon_2<K1>> new_PolyParts;
            
            try{
            getCenterLine<K1, CGAL::Polygon_with_holes_2<K1>, CGAL::Polygon_2<K1>>(new_space, new_points, segs, aug_points, aug_segs, new_PolyParts);
            std::string title = "output " + std::to_string(id);
            CGAL::CenterLineViewer<CGAL::Polygon_2<K1>> mainwindow(app.activeWindow(), *new_PolyParts.begin(), title.c_str());
            mainwindow.drawPartitions(new_PolyParts);
            mainwindow.drawTree(new_points, segs);
            mainwindow.drawTree(aug_points, aug_segs, CGAL::Color(0, 255, 127));
            mainwindow.show();
            app.exec();
            mainwindow.saveImage(title + ".png");
            } catch(std::exception &e){
                std::cerr << "error = " << e.what() << std::endl;
            }

            

            ++id;
        }
    } // void PolygonCenterLine::showCenterLine(double interval)
}
#endif // PolygonCenterLine_Implementation