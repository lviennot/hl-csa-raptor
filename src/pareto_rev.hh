#ifndef PARETO_REV_HH
#define PARETO_REV_HH

#include <assert.h>
#include <stdio.h>
#include <limits>
#include <vector>
#include <utility>

/**
 * 2D Pareto set in a sorted vector in reverse order, i.e. the points 
 * still form an increasing curve, but are read in decreasing order of x coord:
 * while scanning the vector, x decreases  while y increases.
 *
 */


template<typename W>

class pareto_rev {
    // 2D Pareto_Rev set

public:

    struct point {
        W x;
        W y;
        point() : x(1000000), y(1000000) {}
        point(W xx, W yy) : x(xx), y(yy) {}
        point add(W dx, W dy) {
            return point(x + dx, y + dy);
        }
        point(const point &p) : x(p.x), y(p.y) {}
    };
    
    std::vector<point> pts;
    typename std::vector<point>::const_reverse_iterator last_bellow;

    pareto_rev(size_t cap=128) { pts.reserve(cap); last_bellow = pts.rend(); }

    pareto_rev(const pareto_rev &p) : pts(p.pts) {}

    pareto_rev& operator=(const pareto_rev& p) {
        if (this != &p) {
            pts = p.pts;
            last_bellow = p.pts.rbegin();
        }
        return *this;
    }

    size_t size() { return pts.size(); }

    void clear() { pts.clear(); last_bellow = pts.rend(); }

    void check() {
        W x, y;
        bool first = true;
        for (auto p : pts) {
            if (first) { first = false; }
            else { assert(x >= p.x && y <= p.y && (x > p.x || y < p.y)); }
            x = p.x; y = p.y;
        }
    }

    W smallest_x_bellow(W x_max, W y_max) { // returns the minimal x coord. of a point bellow y_max, x_max if none exists
        auto pts_end = pts.crend();
        for (auto p = pts.crbegin(); p != pts_end; ++p) {
            if (p->y <= y_max) return p->x;
        }
        return x_max;
    }

    W min_x(W dft) const {
        if (pts.empty()) return dft;
        return pts.back().x;
    }

    bool add(W x, W y) { // returns true if point was added
        bool insert = false;
        auto p = pts.rbegin(), pts_end = pts.rend();
        auto pos = pts_end;
        point tmp(x, y);
        size_t n_dom = 0;
        for ( ; p != pts_end; ++p) {
            if (p->x <= x && p->y <= y) return false; // dominated
            if (p->x < x) continue;
            // x <= p->x
            if (y <= p->y) {// *p is dom (y < p->y if x == p->x as x,y not dom) 
                if ( ! insert) {
                    std::swap(*p, tmp);
                    insert = true;
                } else if (pos == pts_end) {
                    pos = p;
                }
            } else { // p->y <= y => x < p->x as x,y not dom
                break;
            }
        }
        if ( ! insert) {
            pts.insert(p.base(), tmp); // base of rev it is after
            return true;
        } else if (pos != pts_end) { // delete after p up to pos
            pts.erase(p.base(), pos.base());
        }
        return true;
    }

    bool del_dominated(W x, W y) { // returns true if x,y dominates some elt
        auto p = pts.rbegin(), pts_end = pts.rend();
        auto pos = pts_end;
        size_t n_dom = 0;
        bool first_dom = true;
        for ( ; p != pts_end ; ++p) {
            if (p->x < x && p->y <= y) return false; // dominated
            if (p->x < x) continue;
            // x <= p->x
            if (y <= p->y) { // *p is dom
                ++n_dom;
                if (first_dom) {
                    pos = p;
                    first_dom = false;
                }
            } else {
                break;
            }
        }
        if (first_dom) {
            return false;
        } else {
            pts.erase(p.base(), pos.base());
            return true;
        }
    }

    void print(std::ostream &cout = std::cout) {
        for (auto p : pts) {
            cout << p.x <<","<< p.y <<" ";
        }
        cout <<"\n";
    }

    bool dominates(W x, W y) {
        auto pts_end = pts.rend();
        for (auto p = pts.rbegin(); p != pts_end; ++p) {
            if (p->x <= x && p->y <= y) return true; // dominated
            if (p->x >= x) return false; // (y < p->y if x == p->x)
        }
        return false;
    }

    bool contains(W x, W y) {
        return std::any_of(pts.begin(), pts.end(),
                           [x,y](const point p) {
                               return p.x == x && p.y == y;
                           });
    }

};




namespace unit {

    void pareto_rev_test(size_t n, int rnd = 100) {
        std::cerr <<"pareto_rev_test: sizeof(int)="<< sizeof(int)
                  <<" ---------- \n";
        pareto_rev<int> ps(n);
        std::vector<pareto_rev<int>::point> dom;
        for (int i = 0; i < n; ++i) {
            int x = rand() % rnd, y = rand() % rnd;
            std::cout <<"add "<< x <<","<< y <<" ";
            bool dominated = ps.dominates(x, y);
            std::cout <<" dom="<< dominated <<"\n";
            if ( ! ps.add(x, y)) {
                ps.print();
                assert(dominated);
                dom.push_back(pareto_rev<int>::point(x,y));
            } else {
                ps.print();
                assert( ! dominated);
            }
            ps.check();
            for (auto p : dom) assert(ps.dominates(p.x, p.y));
            y = - rand() % rnd;
            x = ps.smallest_x_bellow(rnd, y);
            std::cout <<"smallest_x_bellow ("<< y <<") = "<< x <<"\n";
            for (auto p : ps.pts) assert(p.y > y || p.x >= x || x == rnd);
            if (i % (n/2) == n/8) {
                dom.clear();
                int x = rand() % (rnd/4+1), y = rand() % (rnd/4+1);
                pareto_rev<int> xy(1);
                xy.add(x,y);
                int n_dom = 0;
                for (auto p : ps.pts) if (xy.dominates(p.x, p.y)) ++n_dom;
                int ps_size = ps.size();
                bool dom = ps.del_dominated(x, y);
                std::cout <<"............del_dominated "<< x <<" "<< y <<" : "
                          << n_dom <<"\n";
                assert(dom == (n_dom > 0));
                assert(n_dom == ps_size - ps.size());
            }
        }
        ps.add(9,4);
        ps.print();
        ps.check();
        for (auto p : dom) assert(ps.dominates(p.x, p.y));
    }
    
}
    


#endif // PARETO_REV_HH
