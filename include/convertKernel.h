#ifndef CAST_KERNEL_H
#define CAST_KERNEL_H
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Polygon_2.h>
using Lazy_expr = CGAL::Lazy_exact_nt<CORE::Expr>;
using Epect_sqrt = CGAL::Simple_cartesian<Lazy_expr>;
//using Fixed_nt = CGAL::Lazy_exact_nt<CGAL::Gmpfr>;
using Fixed_nt = CGAL::Gmpfr;
using Fixed_kernel = CGAL::Simple_cartesian<Fixed_nt>;
using Simple_kernel = CGAL::Simple_cartesian<double>;

struct Lazy_gmpq_to_Expr_converter : public std::unary_function< CGAL::Lazy_exact_nt<CGAL::Gmpq>, Lazy_expr>
{
    //CORE::Expr
    Lazy_expr
    operator()(const CGAL::Lazy_exact_nt<CGAL::Gmpq> &a) const
    {
      return Lazy_expr(::CORE::BigRat(exact(a).mpq()));
    }
};

struct Lazy_gmpq_to_double_converter : public std::unary_function<CGAL::Lazy_exact_nt<CGAL::Gmpq>, double>
{
    double operator () (const CGAL::Lazy_exact_nt<CGAL::Gmpq> &a) const{
        return CGAL::to_double(a);
    }
};

struct Lazy_gmpq_to_gmpfr_converter : public std::unary_function<CGAL::Lazy_exact_nt<CGAL::Gmpq>, Fixed_nt> {
    Fixed_nt operator () (const CGAL::Lazy_exact_nt<CGAL::Gmpq> &a) const{
        return Fixed_nt(CGAL::to_double(a), 20);
    }
};

// <CGAL::Epeck, Epect_sqrt, Lazy_gmpq_to_Expr_converter>
template<typename From, typename To, typename Converter>
struct KernelConverter : public CGAL::Cartesian_converter<From, To, Converter>
{
    //using K = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
    using K = typename To;
    CGAL::Polygon_2<K> convert(const CGAL::Polygon_2<From> &poly){
        std::vector<K::Point_2> tmp_pts;
        for(auto &point : poly) {
            tmp_pts.emplace_back(operator ()(point));
        }
        return CGAL::Polygon_2<K>(tmp_pts.begin(), tmp_pts.end());
    }
    CGAL::Polygon_with_holes_2<K> convert(const CGAL::Polygon_with_holes_2<CGAL::Epeck> &poly){
        CGAL::Polygon_2<K> bound = convert(poly.outer_boundary());

        std::vector<CGAL::Polygon_2<K>> holes;
        for(auto it = poly.holes_begin();it != poly.holes_end();++it){
            holes.emplace_back(convert(*it));
        }
        return CGAL::Polygon_with_holes_2<K>(bound, holes.begin(), holes.end());
    }
};

#endif
