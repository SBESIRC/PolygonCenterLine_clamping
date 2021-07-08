#pragma once
#include "CenterLineSolver.h"
namespace CenterLineSolver{
    template <typename K, class Poly_with_holes, class Poly>
    struct CenterLineSolver<K, Poly_with_holes, Poly>::Event {
        FT time;
        EdgeData *a, *b;
        EventType type;
        Event(FT t, EdgeData *e0, EdgeData *e1, EventType tp)
            : time(t), a(e0), b(e1), type(tp)
        {
            assert(e0);
            assert(e1);
        }
        bool is_active() const { return a->is_active() && b->is_active(); }
        bool operator<(const Event &b) const { return time > b.time; }
    };
}