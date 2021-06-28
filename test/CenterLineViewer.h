#pragma once
#include <string>

#include <CGAL/Qt/Basic_viewer_qt.h>
#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/Polygon_2.h>
#include <CGAL/Random.h>

namespace CGAL {
    template <class P2>
    class CenterLineViewer : public Basic_viewer_qt {
        using Base = Basic_viewer_qt;
        using Point = typename P2::Point_2;

    public:
        /// Construct the viewer.
        /// @param ap2 the polygon to view
        /// @param title the title of the window
        // First draw: vertices; edges, faces; multi-color; no inverse normal
        CenterLineViewer(QWidget *parent, const P2 &ap2, const char *title = "CenterLine Viewer") : Base(parent, title, true, true, true, false, false), poly(ap2)
        {
            compute_elements();
        }

    protected:
        void compute_elements()
        {
            clear();

            if (poly.is_empty())
                return;

            Point prev = poly.vertex(poly.size() - 1);

            CGAL::Color c(75, 160, 255);
            face_begin(c);

            for (typename P2::Vertex_const_iterator i = poly.vertices_begin();
                 i != poly.vertices_end(); ++i) {
                add_point(*i);         // Add vertex
                add_segment(prev, *i); // Add segment with previous point
                add_point_in_face(*i); // Add point in face
                prev = *i;
            }

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

    protected:
        const P2 &poly;

    public:
        void drawPartitions(const std::vector<P2> &poly_list)
        {
            clear();
            for (const P2 &poly : poly_list) {
                Point prev = poly.vertex(poly.size() - 1);
                CGAL::Color c(75, 160, 255);
                face_begin(c);
                for (typename P2::Vertex_const_iterator i = poly.vertices_begin(); i != poly.vertices_end(); ++i) {
                    add_point(*i);         // Add vertex
                    add_segment(prev, *i); // Add segment with previous point
                    add_point_in_face(*i); // Add point in face
                    prev = *i;
                }
                face_end();
            }
        }

        void drawTree(const std::vector<Point> &point, const std::vector<std::pair<int, int>> &segment, CGAL::Color c = CGAL::Color(255, 255, 0))
        {
            //compute_elements();
            for (auto &v : point)
                add_point(v);
            for (auto &seg : segment)
                add_segment(point[seg.first], point[seg.second], c);
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

    template <class T, class C>
    static void drawPoly(const CGAL::Polygon_2<T, C> &ap2, const char *title = "Polygon_2 Basic Viewer")
    {
#if defined(CGAL_TEST_SUITE)
        bool cgal_test_suite = true;
#else
        bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

        if (!cgal_test_suite) {
            int argc = 1;
            const char *argv[2] = {"t2_viewer", "\0"};
            QApplication app(argc, const_cast<char **>(argv));
            SimplePolygon2ViewerQt<CGAL::Polygon_2<T, C>>
                mainwindow(app.activeWindow(), ap2, title);
            mainwindow.show();
            app.exec();
        }
    }

    // Specialization of draw function.
    template <typename T, typename C, typename P = typename T::Point_2>
    static void drawCenterLine(const std::vector<CGAL::Polygon_2<T, C>> &poly_list,
                               const std::vector<P> &points,
                               const std::vector<std::pair<int, int>> &segments,
                               const char *title = "CenterLine Viewer")
    {
#if defined(CGAL_TEST_SUITE)
        bool cgal_test_suite = true;
#else
        bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

        if (!cgal_test_suite) {
            int argc = 1;
            const char *argv[2] = {"t2_viewer", "\0"};
            QApplication app(argc, const_cast<char **>(argv));
            CenterLineViewer<CGAL::Polygon_2<T, C>>
                mainwindow(app.activeWindow(), *poly_list.begin(), title);
            mainwindow.drawPartitions(poly_list);
            mainwindow.drawTree(points, segments);
            mainwindow.show();
            app.exec();
            mainwindow.saveImage(string(title) + ".png");
        }
    }

    template <class T, class C>
    static void drawPartition(const std::vector<CGAL::Polygon_2<T, C>> &poly_list, const char *title = "Partition Viewer")
    {
#if defined(CGAL_TEST_SUITE)
        bool cgal_test_suite = true;
#else
        bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

        if (!cgal_test_suite) {
            int argc = 1;
            const char *argv[2] = {"t2_viewer", "\0"};
            QApplication app(argc, const_cast<char **>(argv));
            CenterLineViewer<CGAL::Polygon_2<T, C>>
                mainwindow(app.activeWindow(), *poly_list.begin(), title);
            mainwindow.drawPartitions(poly_list);
            mainwindow.show();
            app.exec();
        }
    }
} // End namespace CGAL
#endif
