#pragma once

#include <assert.h>
#include <vector>
#include <queue>
#include <set>

#include "timetable.hh"

template<class C> // cost type, must override < and + operators, and have a constructor from a connection

class double_scan {

public:
    typedef timetable::ST ST;
    typedef timetable::S S;
    typedef timetable::R R;
    typedef timetable::T T;
    
    typedef int TR; // trips
    TR n_tr;
    
    struct connection {
        TR trip; // for debugging, cannot use it to vary minimum waiting time in double_scan
        ST from, to;
        T dep, arr;
        int index;
        connection(TR tr, S u, S v, T d, T a, int i)
            : trip(tr), from(u), to(v), dep(d), arr(a), index(i) {}
    };

private:
    const timetable &ttbl;



    std::vector<connection> e_dep; // connections ordered by from station and dep time
    std::vector<std::pair<R, int> > trip_route;

    typedef int CI; // index of a connection in e_dep

    std::vector<CI> e_arr; // connections ordered by arrival time

    std::vector<cost> edge_cost;
    std::vector<cost> best_cost;
    std::vector<CI> parent;

   

public:
    double_scan(const timetable &tt)
        : ttbl(tt), n_tr(0) {
        
        // create connections :
        int n_conn = 0;
        for (R r = 0; r < ttbl.n_r; ++r) {
            n_tr += tt.trips_of[r].size();
            for (int i = 0; i < tt.trips_of[r].size(); ++i) {
                n_conn += tt.trips_of[r][i].size() - 1;
            }
        }
        conn.reserve(n_conn);
        trip_route.reserve(n_tr);
        int i_tr = 0;
        for (R r = 0; r < ttbl.n_r; ++r) {
            const std::vector<S> &stops = tt.route_stops[r];
            for (int i = 0; i < tt.trips_of[r].size(); ++i) {
                for (int j = 1; j < tt.trips_of[r][i].size(); ++j) {
                    conn.emplace_back(i_tr,
                                      tt.stop_station[stops[j-1]],
                                      tt.stop_station[stops[j]],
                                      tt.trips_of[r][i][j-1].second,
                                      tt.trips_of[r][i][j].first,
                                      j);
                }
                trip_route.emplace_back(r, i);
                ++i_tr;
            }
        }
        
        // compute costs :
        edge_cost.reserve(e_dep.size());
        for (connection conn : e_dep) {
            edge_cost.emplace_back(conn);
        }
    }
}
