#pragma once
#include<vector>

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

//struct FireExt {
//    Point coords;
//    std::string category;
//};

struct parseout {
    std::vector<Block> data;
    //std::vector<FireExt> dataF;
    double minx = DBL_MAX, maxx = DBL_MIN, miny = DBL_MAX, maxy = DBL_MIN;
};

struct Fe {
    Point p;
    std::string type;
    int rspace;
};

template <typename K>
struct valout {
    //Point p;
    std::string type;
    std::vector<std::vector<typename K::Line_2>> oS;
    //General_Polygon_Set oS;
};
