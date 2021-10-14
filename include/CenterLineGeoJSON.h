#pragma once
#include "stdafx.h"

namespace CenterLine
{
    struct Point {
        double x, y;
        Point() : x(0), y(0) {}
        Point(double a, double b) : x(a), y(b) {}
    };

    struct Line {
        Point a, b;
        Line(Point p0, Point p1) : a(p0), b(p1) {}
    };

    struct Block {
        std::vector<std::vector<Point>> coords;
        std::string category;
        std::string sw;
        std::string privacy;
    };

    struct parseout {
        std::vector<Block> data;
        double minx = DBL_MAX, maxx = DBL_MIN, miny = DBL_MAX, maxy = DBL_MIN;
    };

    struct Fe {
        Point p;
        std::string type;
        int rspace;
    };

    template <typename K>
    struct valout {
        std::string type;
        std::vector<std::vector<typename K::Line_2>> oS;
    };

    void parse_geojson(const std::string& datastr, parseout& out);
    std::string out2str(std::vector<std::pair<Point, Point>> out);
    std::string out2str(const std::vector<Block> &ucs_blocks, const std::vector<Point> &corr_ucs);
    std::string out2str(const std::vector<Block> &rect_blocks, const std::vector<Block> &centerline_blocks);
}
