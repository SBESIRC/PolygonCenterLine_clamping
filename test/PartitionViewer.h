#pragma once
#include <string>

#include <CGAL/Qt/Basic_viewer_qt.h>
#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/Object.h>
#include <CGAL/General_polygon_set_2.h>

namespace CGAL {
	template<typename K>
    struct PartitionViewer : public Basic_viewer_qt {
        using Base = Basic_viewer_qt;
		using FT = typename K::FT;
		using Point_2 = CGAL::Point_2<K>;
        using Segment_2 = CGAL::Segment_2<K>;
		using Polygon_2 = CGAL::Polygon_2<K>;
		using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

		using Gps_traits_2 = CGAL::Gps_circle_segment_traits_2<K>;
		using Offset_polygon_2 = typename Gps_traits_2::Polygon_2;
		using Offset_polygon_with_holes_2 = typename Gps_traits_2::Polygon_with_holes_2;
		using Polygon_set_2 = CGAL::General_polygon_set_2<Gps_traits_2>;
		using X_monotone_curve_2 = typename Gps_traits_2::X_monotone_curve_2;
		using Intersect_2 = typename Gps_traits_2::Intersect_2;
		using Intersection_map = typename X_monotone_curve_2::Intersection_map;
		using One_root_point_2 = typename Gps_traits_2::Point_2;
		using CoordNT = typename Gps_traits_2::Point_2::CoordNT;

		Polygon_with_holes_2 space;

		PartitionViewer(QWidget *parent, const Polygon_with_holes_2 &sp, const char *title = "CenterLine Viewer") : Base(parent, title, true, true, true, false, false), space(sp)
        {
            compute_elements();
        }
		void compute_elements(){
			clear();
            //CGAL::Color c(75, 160, 255);
            //drawPoly(space, c);
            traversePoly(space);
		}
        void traversePoly(const Polygon_2 &poly, CGAL::Color c = CGAL::Color(0, 0, 0)){
            for(auto e_it = poly.edges_begin();e_it != poly.edges_end();++e_it){
                add_point(e_it->source());         // Add vertex
                add_segment(e_it->source(), e_it->target(), c); // Add segment with previous point
                add_point_in_face(e_it->target()); // Add point in face
            }
            add_point_in_face(poly.edges_begin()->source());
        }
        void traversePoly(const Polygon_with_holes_2 &poly, CGAL::Color c = CGAL::Color(0, 0, 0)){
            traversePoly(poly.outer_boundary(), c);
            for(auto h_it = poly.holes_begin(); h_it != poly.holes_end();++h_it){
                traversePoly(*h_it, c);
            }
        }
        void traversePoly(const Offset_polygon_2 &poly, CGAL::Color c = CGAL::Color(0, 0, 0)){
            Polygon_2 linear_poly = straighten_poly(poly);
            std::cout << "linear_poly=" << linear_poly;
            traversePoly(linear_poly, c);
        }
        void traversePoly(const Offset_polygon_with_holes_2 &poly, CGAL::Color c = CGAL::Color(0, 0, 0)){
            traversePoly(poly.outer_boundary(), c);
            for(auto h_it = poly.holes_begin(); h_it != poly.holes_end();++h_it){
                traversePoly(*h_it, c);
            }
        }
        void traversePoly(const std::vector<Segment_2> &segs, CGAL::Color c = CGAL::Color(0, 0, 0)){
            for(auto &seg : segs){
                add_point(seg.source());
                add_segment(seg.source(), seg.target(), c);
            }
        }

        template<typename NT0, typename NT1>
		static NT1 trans_nt(const NT0 &x){ return CGAL::to_double(x); }
		static Point_2 trans_point_traits(const One_root_point_2 &point){
			auto x = point.x(), y = point.y();
			return Point_2(trans_nt<CoordNT, FT>(x), trans_nt<CoordNT, FT>(y));
		}
        static Polygon_2 straighten_poly(const Offset_polygon_2 &origin){
			Polygon_2 res;
			for(auto c_it = origin.curves_begin(); c_it != origin.curves_end();++c_it){
				res.push_back(trans_point_traits(c_it->source()));
				if(c_it->is_circular()){
					CoordNT L, R;
					CGAL::Bbox_2 box = c_it->bbox();
					L = c_it->source().x(), R = c_it->target().x();
					Intersection_map intersect_map;
					Intersect_2 intersect_solver(intersect_map);
					for(int i = 1;i < 8;++i){
						std::vector<CGAL::Object> pts;
						std::pair<One_root_point_2, typename Gps_traits_2::Multiplicity> p;
						FT x = trans_nt<CoordNT, FT>((R * i + L * (8 - i)) / 8);
						X_monotone_curve_2 vert(Point_2(x, box.ymin() - 10), Point_2(x, box.ymax() + 10));
						intersect_solver(*c_it, vert, std::back_inserter(pts));
						if(!pts.empty() && assign(p, pts[0])){
							res.push_back(trans_point_traits(p.first));
						}
					}
				}
			}
			return res;
		}

        template<typename Poly>
        void drawPoly(const Poly &poly, CGAL::Color c, CGAL::Color seg_c = CGAL::Color(0, 0, 0)){
            face_begin(c);
            traversePoly(poly, seg_c);
            face_end();
        }

		virtual void keyPressEvent(QKeyEvent *e)
        {
            // Test key pressed:
            //    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
            //    if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton)) { ... }

            // Call: * compute_elements() if the model changed, followed by
            //       * redraw() if some viewing parameters changed that implies some
            //                  modifications of the buffers
            //                  (eg. type of normal, color/mono)
            //       * update() just to update the drawing

            // Call the base method to process others/classicals key
            Base::keyPressEvent(e);
        }

		void saveImage(std::string filename)
        {
            qreal aspectRatio = width() / static_cast<qreal>(height());
            static ImageInterface *imageInterface = nullptr;
            static bool _expand_frustum;
            static qglviewer::SnapShotBackground _background;
            static double _oversampling;
            if (!imageInterface) {
                imageInterface = new ImageInterface(this, aspectRatio);
                imageInterface->imgWidth->setValue(width());
                imageInterface->imgHeight->setValue(height());

                if (imageInterface->exec() == QDialog::Rejected) {
                    return;
                }
                _expand_frustum = imageInterface->expandFrustum->isChecked();
                _background = qglviewer::SnapShotBackground(imageInterface->color_comboBox->currentIndex());
                _oversampling = imageInterface->oversampling->value();
            }
            QSize finalSize(width(), height());

            QImage *image = takeSnapshot(_background, finalSize, _oversampling, _expand_frustum);
            if (image) {
                image->save(QString::fromStdString(filename));
                delete image;
            }
        }
    };
} // End namespace CGAL
#endif