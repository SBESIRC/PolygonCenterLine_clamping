#ifndef CAST_KERNEL_H
#define CAST_KERNEL_H
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Polygon_2.h>

namespace KernelConverter {
    using Epect_sqrt = CGAL::Simple_cartesian<CORE::Expr>;
    //using Fixed_nt = CGAL::Lazy_exact_nt<CGAL::Gmpfr>;
    using Fixed_nt = CGAL::Gmpfr;
    using Fixed_kernel = CGAL::Simple_cartesian<Fixed_nt>;
    using Simple_kernel = CGAL::Simple_cartesian<double>;

    // Precision=-1表示使用默认精度
    template <typename K, typename NewK, int Precision = -1>
    struct NumberConverter : public std::unary_function<K, NewK> {
        typename NewK operator()(const typename K &a) const
        {
            return NewK(CGAL::to_double(a));
        }
    };
    template <typename K>
    struct NumberConverter<K, double> : public std::unary_function<K, double> {
        double operator()(const K &a) const
        {
            return CGAL::to_double(a);
        }
    };
    template <>
    struct NumberConverter<CGAL::Lazy_exact_nt<CGAL::Gmpq>, Fixed_nt> : public std::unary_function<CGAL::Lazy_exact_nt<CGAL::Gmpq>, Fixed_nt> {
        Fixed_nt operator()(const CGAL::Lazy_exact_nt<CGAL::Gmpq> &a) const
        {
            return Fixed_nt(CGAL::to_double(a));
        }
    };
    template <int Precision>
    struct NumberConverter<CGAL::Lazy_exact_nt<CGAL::Gmpq>, Fixed_nt, Precision> : public std::unary_function<CGAL::Lazy_exact_nt<CGAL::Gmpq>, Fixed_nt> {
        Fixed_nt operator()(const CGAL::Lazy_exact_nt<CGAL::Gmpq> &a) const
        {
            return Fixed_nt(CGAL::to_double(a), Precision);
        }
    };
    // <CGAL::Epeck, Epect_sqrt, NumberConverter<A, B, 20> >
    template <typename From, typename To, typename Converter>
    struct KernelConverter : public CGAL::Cartesian_converter<From, To, Converter> {
        //using K = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
        using K = typename To;
        CGAL::Polygon_2<K> convert(const CGAL::Polygon_2<From> &poly)
        {
            std::vector<K::Point_2> tmp_pts;
            for (auto &point : poly) {
                tmp_pts.emplace_back(operator()(point));
            }
            return CGAL::Polygon_2<K>(tmp_pts.begin(), tmp_pts.end());
        }
        CGAL::Polygon_with_holes_2<K> convert(const CGAL::Polygon_with_holes_2<CGAL::Epeck> &poly)
        {
            CGAL::Polygon_2<K> bound = convert(poly.outer_boundary());

            std::vector<CGAL::Polygon_2<K>> holes;
            for (auto it = poly.holes_begin(); it != poly.holes_end(); ++it) {
                holes.emplace_back(convert(*it));
            }
            return CGAL::Polygon_with_holes_2<K>(bound, holes.begin(), holes.end());
        }
    };
}
#endif
