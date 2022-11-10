#pragma once

#include <assert.h>
#include <vector>
#include <queue>
#include <set>

#include "timetable.hh"


typedef timetable::ST ST;
typedef timetable::S S;
typedef timetable::R R;
typedef timetable::T T;
    
typedef int TR; // trips

struct temp_edge {
    TR trip; // for debugging, cannot use it to vary minimum waiting time in double_scan
    ST from, to;
    T dep, arr;
    int index;
    temp_edge(TR tr, ST u, ST v, T d, T a, int i)
        : trip(tr), from(u), to(v), dep(d), arr(a), index(i) {}
};



template<class C> // cost type, must override < and + operators, and have a constructor from a temp_edge

class double_scan {

private:
    const timetable &ttbl;



    std::vector<temp_edge> e_dep; // temp_edges ordered by from station and dep time
    std::vector<std::pair<R, int> > trip_route;

    typedef int E; // index of a temp_edge in e_dep

    std::vector<E> e_arr; // temp_edges ordered by arrival time

    TR n_tr;
    E n_edg;
    
    std::vector<C> edge_cost;
    std::vector<C> best_cost;
    std::vector<E> parent;
    std::vector<C> min_cost;
    // TODO: std::vector<E> prev_in_trip; // next temp_edge in the trip

    struct e_dep_indexes {
        E first, last; // first/last edge from a station
        E left, right; // indexes of edges from a station currently reachable
    };  
    std::vector<e_dep_indexes> st_indexes; // Indexes of edges from a given station in e_dep

    struct interval {
        E left, right;
        C cost;
        E parent;
    };
    std::vector<std::vector<interval> > st_intervals;
    

    static const E not_a_temp_edge = -1;

public:
    double_scan(const timetable &tt)
        : ttbl(tt), n_tr(0), n_edg(0) {
        
        // create temp_edges and store them in e_dep:
        for (R r = 0; r < ttbl.n_r; ++r) {
            n_tr += tt.trips_of[r].size();
            for (int i = 0; i < tt.trips_of[r].size(); ++i) {
                n_edg += tt.trips_of[r][i].size() - 1;
            }
        }
        e_dep.reserve(n_edg);
        trip_route.reserve(n_tr);
        int i_tr = 0;
        for (R r = 0; r < ttbl.n_r; ++r) {
            const std::vector<S> &stops = tt.route_stops[r];
            for (int i = 0; i < tt.trips_of[r].size(); ++i) {
                for (int j = 1; j < tt.trips_of[r][i].size(); ++j) {
                    e_dep.emplace_back(i_tr,
                                       tt.stop_station[stops[j-1]],
                                       tt.stop_station[stops[j]],
                                       tt.trips_of[r][i][j-1].second,
                                       tt.trips_of[r][i][j].first,
                                       j); // 1 if first edge of a trip
                }
                trip_route.emplace_back(r, i);
                ++i_tr;
            }
        }
        std::sort(e_dep.begin(), e_dep.end(),
                  [](const temp_edge &c, const temp_edge &d) {
                      if (c.from != d.from) return c.from < d.from;
                      return c.dep < d.dep;
                  });

        // create e_arr:
        e_arr.reserve(n_edg);
        for (E i = 0; i < n_edg; ++i) { e_arr.push_back(i); }
        std::sort(e_arr.begin(), e_arr.end(),
                  [this](E i, E j) {
                      const temp_edge &c = e_dep[i];
                      const temp_edge &d = e_dep[j];
                      if (c.arr != d.arr) return c.arr < d.arr;
                      // c.arr = d.arr
                      // Be careful to 0 delay temp_edges:
                      // Heuristic for zero-acyclicity (not sufficient):
                      if (c.trip == d.trip) return c.index < d.index;
                      // 0 delay temp_edge after:
                      if (c.dep != c.arr) return true;
                      if (d.dep != d.arr) return false;
                      return c.trip < d.trip;
                  });
        
        // costs and parent associated to edges:
        edge_cost.reserve(n_edg);
        best_cost.reserve(n_edg);
        parent.reserve(n_edg);
        for (temp_edge edg : e_dep) {
            edge_cost.emplace_back(edg);
            best_cost.emplace_back();
            parent.push_back(not_a_temp_edge);
        }

    }


    T earliest_arrival_time(const ST src, const ST dst, const T t_dep,
                            const T min_waiting_time,
                            const T max_waiting_time) {
        assert(0 <= t_dep && t_dep <= 3600*48);
        assert(max_waiting_time > 0);

        // initialize
        for (int i = 0; i < n_edg; ++i) { parent[i] = not_a_temp_edge; }

        // scan temp_edges by non-decreasing arrival times:
        for (E i : e_arr) {
            const temp_edge &e = e_dep[i];
        }
    }
};


// Standard costs:
    
struct cost_eat {

    T arr;

    cost_eat() : arr(std::numeric_limits<T>::max()) {
    }
    
    cost_eat(const temp_edge &edg) : arr(edg.arr) {
    }
    
    friend bool operator<(const cost_eat &e, const cost_eat &f) {
        return e.arr < f.arr;
    }
    
    friend cost_eat operator+(cost_eat e, cost_eat f) {
        return f;
    }
};
