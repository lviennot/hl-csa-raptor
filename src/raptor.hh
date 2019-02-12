#ifndef RAPTOR_HH
#define RAPTOR_HH

#include <assert.h>
#include <vector>
#include <queue>
#include <set>

#include "timetable.hh"
#include "connection_scan.hh"
#include "pareto_rev.hh"


class raptor {
private:
    const timetable &ttbl;

    typedef timetable::ST ST;
    typedef timetable::S S;
    typedef timetable::R R;
    typedef timetable::T T;

    static const int ntrips_max = 48; // max number of trips in a journey

    const int not_stop_index = -1;

    typedef pareto_rev<T> pset;
    typedef pset::point point; // x = est.arr.time, y = walking.time
    
    std::vector<T> st_eat; // earliest arrival time at station, hub
    std::vector<T> stop_prev_dep; // dep time of previous trip at a stop
    std::vector<pset> all_pareto, incr_pareto, tmp_pareto, //eat vs walking time
        prev_dep_pareto, dst_pareto; //eat vs walking time at stop/destination
    //std::vector<T> eat; // earliest arrival time at stop
    int dst_pareto_kmax; // max nb of trips in dst pareto set
    
    struct parent_t {
        bool trip;
        union {
            S stop;
            ST station;
        };
        T eat;
        T dist;
        parent_t() : trip(false), station(-1),
                     eat(std::numeric_limits<int>::max() / 2),
                     dist(std::numeric_limits<int>::max() / 2) {}
        void clear() {
            trip = false; stop = -1;
            eat=std::numeric_limits<int>::max() / 2;
            dist=std::numeric_limits<int>::max() / 2;
        }
    };
    std::vector<std::vector<parent_t> > parent; // not_stop_index for src
    std::vector<ST> n_trips; // number of trips for current eat
    
    std::vector<ST> improved_stations, improved_hubs;
    std::vector<bool> station_has_improved, hub_has_improved;
    std::vector<bool> stop_has_improved;
    std::vector<R> improved_routes;
    std::vector<int> route_has_improved_from, route_has_improved_to;


public:
    raptor(const timetable &tt)
        : ttbl(tt),
          st_eat(tt.n_h+1), // a dummy dst for one to all queries
          n_trips(tt.n_h), stop_prev_dep(tt.n_s),
          all_pareto(tt.n_h), incr_pareto(tt.n_h),
          tmp_pareto(tt.n_s), dst_pareto(ntrips_max + 1),
          dst_pareto_kmax(0),
          prev_dep_pareto(tt.n_s),
          station_has_improved(tt.n_st, false),
          hub_has_improved(tt.n_h, false),
          stop_has_improved(tt.n_s, false),
          route_has_improved_from(tt.n_r),
          route_has_improved_to(tt.n_r)
    {
        for (R r = 0; r < ttbl.n_r; ++r) {
            route_has_improved_from[r] = not_stop_index;
        }
        improved_stations.reserve(tt.n_h);
        improved_routes.reserve(tt.n_r);
        improved_hubs.reserve(tt.n_h);

        parent.reserve(ntrips_max + 1);
        for (int i = 0; i <= ntrips_max; ++i) {
            parent.push_back(std::vector<parent_t>(tt.n_h));
        }
    }

    bool update_eat_trip(ST st, T t, T dt, S par, int k) {
        if (t < st_eat[st]) {
            st_eat[st] = t;
            n_trips[st] = k;
            parent[k][st].trip = true;
            parent[k][st].stop = par;
            parent[k][st].eat = t;
            parent[k][st].dist = dt;
            return true;
        }
        return false;
    }
    
    bool update_eat_walk(ST st, T t, T dt, ST par, int k) {
        if (t < st_eat[st]) {
            st_eat[st] = t;
            n_trips[st] = k;
            parent[k][st].trip = false;
            parent[k][st].station = par;
            parent[k][st].eat = t;
            parent[k][st].dist = dt;
            return true;
        }
        return false;
    }
    
    T earliest_arrival_time(const ST src, const ST dst, const T t_dep,
                            const bool use_hubs = true,
                            const bool use_transfers = false,
                            const T min_chg_bef = 60, // chg time before trip
                            const T min_chg_aft = 0, // chg time after trip
                            const int k_max = ntrips_max,
                            const bool at_least_one_trip = false,
                            connection_scan *earliest_only_csa = nullptr) {

        assert(k_max <= ntrips_max);
        
        // initialize
        for (int i = 0; i <= ttbl.n_h; ++i) { st_eat[i] = ttbl.t_max; }
        for (int i = 0; i < ttbl.n_h; ++i) { n_trips[i] = ntrips_max + 1000; }
        for (int i = 0; i < ttbl.n_s; ++i) { stop_prev_dep[i] = ttbl.t_max; }

        // update helper (first phase)
        auto reach_station_trip = [this, dst](ST st, T t, T dt,
                                              R r, S par, int k) {
            if (t < st_eat[dst] // target pruning
                // bad idea: && t + pll_lb.distance(st, dst) <= st_eat[dst]
                && update_eat_trip(st, t, dt, par, k)
                ) {
                if ( ! station_has_improved[st]) {
                    improved_stations.push_back(st);
                    station_has_improved[st] = true;
                }
            }
        };

        // update helper (second phase)
        auto reach_station_walk =
            [this, dst, min_chg_bef](ST st, T t, T dt,
                                 S par, int k, bool self_walk = false) {
            if (t < st_eat[dst] // target pruning
                // bad idea: && t + pll_lb.distance(st, dst) <= st_eat[dst]
                && (self_walk || update_eat_walk(st, t, dt, par, k))
                && (st != dst) // could happen with at_least_one_trip=true
                ) {
                if (st >= ttbl.n_st) return;
                for (S u : ttbl.station_stops[st]) {
                    if (t + min_chg_bef <= stop_prev_dep[u]) {
                        R r = ttbl.stop_route[u].first;
                        /* OPTIM. TO FIX: if (station_has_improved[st]
                            && parent[k][st].trip
                            && r == ttbl.stop_route[parent[k][st].stop].first) {
                            continue;
                            }*/
                        stop_has_improved[u] = true;
                        int i = ttbl.stop_route[u].second;
                        int i_prev = route_has_improved_from[r];
                        if (i_prev == not_stop_index) {
                            improved_routes.push_back(r);
                            route_has_improved_from[r] = i;
                            route_has_improved_to[r] = i;
                        } else {
                            if (i < i_prev) route_has_improved_from[r] = i;
                            if (i > route_has_improved_to[r])
                                route_has_improved_to[r] = i;
                        }
                    }
                }
            }
        };

        assert(improved_routes.empty());
        assert(improved_stations.empty());

        reach_station_walk(src, t_dep, 0, not_stop_index, 0);

        if (use_transfers) {
            for (auto f : ttbl.transfers[src]) {
                reach_station_walk(f.dst, t_dep + f.wgt, f.wgt, src, 0);
            }
        }

        if (use_hubs) {
            /* Only ok for pure arrival time (not pareto set with k): */
            if (earliest_only_csa != nullptr){
                assert(min_chg_aft == 0);
                st_eat[dst] = 1 + earliest_only_csa
                    ->earliest_arrival_time(src, dst, t_dep,
                                            false, true, min_chg_bef,k_max);
                if (st_eat[dst] > ttbl.t_max) st_eat[dst] = ttbl.t_max;
            }
            for (auto e : ttbl.outhubs[src]) {
                if (t_dep + e.wgt >= st_eat[dst]) break; // target prun
                reach_station_walk(e.dst, t_dep + e.wgt, e.wgt, src, 0);
            }
            if (dst < ttbl.n_st) for (auto e : ttbl.rev_inhubs_id[dst]) {
                reach_station_walk(dst, st_eat[e.dst] + e.wgt, e.wgt, e.dst, 0);
            }
            for (auto e : ttbl.outhubs[src]) {
                for (auto f : ttbl.inhubs[e.dst]) {
                    T t = st_eat[e.dst] + f.wgt;
                    if (t >= st_eat[dst]) break; // targetpr
                    reach_station_walk(f.dst, t, f.wgt, e.dst, 0);
                }
            }
        }

        if (at_least_one_trip) {
            for (int i = 0; i < ttbl.n_h; ++i) { st_eat[i] = ttbl.t_max; }
            //st_eat[dst] = ttbl.t_max;
            st_eat[src] = t_dep;
        }

        int k = 1;
        for (; k <= k_max; ++k) {

            // first phase
            if (improved_routes.empty()) { break; }
            for (const R r : improved_routes) {
                const std::vector<std::vector<std::pair<T,T> > > &trips
                    = ttbl.trips_of[r];
                const std::vector<S> &stops = ttbl.route_stops[r];
                int x_beg = route_has_improved_from[r];
                int x_end = stops.size();
                assert(x_beg >= 0);
                assert(x_beg < x_end);
                int x_last = route_has_improved_to[r];
                assert(x_last >= 0 && x_last < x_end);
                int y_end = trips.size();
                S par = not_stop_index;
                T par_eat = ttbl.t_max;
                for (int x = x_beg, y = y_end; x < x_end; ++x) {
                   //std::cerr <<" x="<< x <<" xend="<< x_end <<"\n";
                    S u = stops[x];
                    //std::cerr << "  u="<< u <<" n_s="<< ttbl.n_s
                    //          <<" "<< ttbl.stop_station.size() <<"\n";
                    ST st = ttbl.stop_station[u];
                    bool eat_improves = true;
                    if (stop_has_improved[u]) {
                        stop_has_improved[u] = false;
                        T eat_when_improved = parent[k-1][st].eat;
                        // Find first departing trip after eat_when_improved:
                        int y_prev = y;
                        if (x == x_beg) {// dicho for first
                            auto lower =
                              std::lower_bound(ttbl.stop_departures[u].begin(),
                                               ttbl.stop_departures[u].end(),
                                               eat_when_improved + min_chg_bef);
                            y = std::distance(ttbl.stop_departures[u].begin(),
                                              lower);
                        } else {
                            while (y-1 >= 0
                                   && ttbl.stop_departures[u][y-1]
                                   >= eat_when_improved + min_chg_bef) {
                                --y;
                            }
                        }
                        if (y < y_prev) {
                            par = u;
                            par_eat = eat_when_improved; // TODO: measure waiting time, this does not work: ttbl.stop_departures[u][y];
                            eat_improves = false;
                        } 
                    }
                    // stop_prev_dep:
                    if (y < y_end) {
                        if (y == 0) stop_prev_dep[u] = - ttbl.t_max;
                        else stop_prev_dep[u] = trips[y-1][x].second;
                    }
                    // target pruning:
                    if (x >= x_last && (y == y_end
                         || trips[y][x].first + min_chg_aft >= st_eat[dst])){
                        break;
                    }
                    // Improve eat:
                    if (y < y_end && eat_improves) {
                        T eat = st_eat[st];
                        T arr = trips[y][x].first + min_chg_aft;
                        if (arr < eat) {
                            assert(par != not_stop_index);
                            reach_station_trip(st, arr, arr-par_eat, r, par, k);
                        }
                    }
                }
                route_has_improved_from[r] = not_stop_index;
            }
            improved_routes.clear();

            // second phase
            //std::cerr <<"k="<< k <<" "<< improved_stations.size() <<" stations\n";
            if (improved_stations.empty()) { break; }
            if (use_transfers) {
                // Assumes transfer graph is transitively closed.
                for (ST st : improved_stations) {
                    for (auto e : ttbl.transfers[st]) { 
                        T t = st_eat[st] + e.wgt;
                        if (t >= st_eat[dst]) break; // targpr
                        reach_station_walk(e.dst, t, e.wgt, st, k,
                                           e.dst == st); // sf (self walk)
                    }
                }
            }
            if (use_hubs) {
                // Assumes OutHub x InHub is transitively closed.
                for (ST st : improved_stations) {
                    if ( ! hub_has_improved[st]) {
                        hub_has_improved[st] = true;
                        improved_hubs.push_back(st);
                    }
                    for (auto e : ttbl.outhubs[st]) {
                        T t = st_eat[st] + e.wgt;
                        if (t >= st_eat[dst]) break; // target pruning
                        if (update_eat_walk(e.dst, t, e.wgt, st, k)) {
                            if ( ! hub_has_improved[e.dst]) {
                                hub_has_improved[e.dst] = true;
                                improved_hubs.push_back(e.dst);
                            }
                        }
                    }
                }
                for (ST h : improved_hubs) {
                    for (auto f : ttbl.inhubs[h]) {
                        T t = st_eat[h] + f.wgt;
                        if (t >= st_eat[dst]) break; // targpr
                        if (f.dst < ttbl.n_st) {
                            reach_station_walk(f.dst, t, f.wgt, h, k,
                                               h < ttbl.n_st && f.dst == h);//sf
                        }
                    }
                    hub_has_improved[h] = false;
                }
                improved_hubs.clear();
            }
            for (ST st : improved_stations) {
                    station_has_improved[st] = false;
            }
            improved_stations.clear();

        }

        if (k > k_max) {
            for (const R r : improved_routes) {
                const std::vector<S> &stops = ttbl.route_stops[r];
                assert(route_has_improved_from[r] != not_stop_index);
                for (int x=route_has_improved_from[r];
                     x <= route_has_improved_to[r]; ++x) {
                    stop_has_improved[stops[x]] = false;
                }
                route_has_improved_from[r] = not_stop_index;
            }
            improved_routes.clear();
            for (ST st : improved_stations) {
                station_has_improved[st] = false;
            }
            improved_stations.clear();
        }

        if (at_least_one_trip && st_eat[dst] < ttbl.t_max) {
            // Fix st_eat for enabling print_journey(dst):
            ST st = dst;
            int k = n_trips[st];
            assert(k > 0);
            while (true) {
                const parent_t &par = parent[k][st];
                ST par_st = par.trip? ttbl.stop_station[par.stop] : par.station;
                if (par_st == not_stop_index) break;
                if (k >= 1) assert(st_eat[st] <= par.eat);
                else st_eat[st] = par.eat;
                st = par_st;
                k = par.trip ? k-1 : k;
                assert(k >= 0);
            }
        }

        //if ( st_eat[dst] < ttbl.t_max && ! use_hubs)
        //    std::cout <<"ntrips "<< n_trips[dst] <<" k "<< k <<"\n";

        return st_eat[dst];
    }

    int eat_one_to_all(const ST src, const T t_dep,
                            const bool use_hubs = true,
                            const bool use_transfers = false,
                            const T min_chg_bef = 60, // chg time before trip
                            const T min_chg_aft = 0, // chg time after trip
                        const int k_max = ntrips_max) {
        earliest_arrival_time(src, ttbl.n_h, t_dep, use_hubs, use_transfers,
                              min_chg_bef, min_chg_aft);
        size_t nvis = 0;
        for (ST i = 0; i < ttbl.n_st; ++i) {
            if (st_eat[i] < ttbl.t_max) ++nvis;
        }
        return nvis;
    }

    T eat_to(ST st) const { return st_eat[st]; }
    int nb_trips_to(ST st) const { return n_trips[st]; }
    
    void print_journey(ST dst,
                       std::ostream &cout = std::cout,
                       T min_chg_bef = 60, T min_chg_aft = 0,
                       int k = -1,
                       bool reverse = true) {
        //cout <<"chg_bef="<< min_chg_bef <<" chg_aft="<< min_chg_aft <<"\n";
        assert(st_eat[dst] < ttbl.t_max);
        if (k == -1) { // main call
            k = n_trips[dst];
            cout <<" to "<< dst <<"="<< ttbl.station_id[dst]
                 <<" : EAT="<< st_eat[dst]
                 <<" ("<< k <<" trips)"
                 <<" :\n";
        }
        const parent_t &par = parent[k][dst];
        ST par_st = par.trip ? ttbl.stop_station[par.stop] : par.station;
        if (par_st == not_stop_index) return;
        int par_k = par.trip ? k-1 : k;
        T par_eat = parent[par_k][par_st].eat;
        if ( ! reverse) print_journey(par_st, cout,
                                      min_chg_bef, min_chg_aft, par_k);
        if (par.trip) {
            assert(k > 0);
            par_st = ttbl.stop_station[par.stop];
            // Check trip exists:
            R r = ttbl.stop_route[par.stop].first;
            const std::vector<std::vector<std::pair<T,T> > > &trips
                    = ttbl.trips_of[r];
            int par_x = ttbl.stop_route[par.stop].second;
            S dst_stp = not_stop_index;
            for (auto s : ttbl.station_stops[dst]) {
                if (ttbl.stop_route[s].first == r
                    && ttbl.stop_route[s].second > par_x) {//looping trips exist
                    dst_stp = s; break;
                }
            }
            int x = ttbl.stop_route[dst_stp].second;
            // find last trip arriving at t:
            int y = -1;
            while (y+1 < ttbl.stop_arrivals[dst_stp].size()
                   && ttbl.stop_arrivals[dst_stp][y+1] +min_chg_aft <= par.eat){
                ++y;
            }
            //cout <<"x="<< x <<" par_x="<< par_x <<"\n";
            for (int i = x; i >= par_x; --i) {
                S s = ttbl.route_stops[r][i];
                ST st = ttbl.stop_station[s];
                cout <<"     "<< st <<"="<< ttbl.hub_id[st]
                     <<" (stop "<< s <<" idx="<< i <<") at "
                     << ttbl.trips_of[r][y][i].first <<">="<< st_eat[st] 
                     <<", dep. "<< ttbl.trips_of[r][y][i].second <<"\n";
            }
            // */
            cout <<"trip "<< k <<"="<< r <<"["<< y <<"/"<< trips.size() <<"]"
                 <<" from "<< par_st <<"="<< ttbl.hub_id[par_st]
                 <<" (stop "<< par.stop <<")"
                 <<" at "<< par_eat
                 <<" to "<< dst <<"="<< ttbl.hub_id[dst]
                 <<" at "<< par.eat <<">="<< st_eat[dst] 
                 <<" ("<< par.dist <<"s)\n";
            assert(ttbl.trips_of[r][y][x].first + min_chg_aft == par.eat);
            assert(par_eat + min_chg_bef <= ttbl.trips_of[r][y][par_x].second);
            assert(par_eat + par.dist == par.eat);
        } else {
            cout <<"walk "<< k <<" from "<< par_st <<"="<< ttbl.hub_id[par_st]
                 <<" at "<< par_eat
                 <<" to "<< dst <<"="<< ttbl.hub_id[dst]
                 <<" at "<< par.eat <<"<="<< st_eat[dst] 
                 <<" ("<< par.dist
                 <<"s>="<< ttbl.walking_time_str(par_st, dst) <<")\n";
            /* Check walking time, only with use_hubs:
            assert(par.dist == // chg_time_aft excet for walk from root
                   (parent[k][par_st].stop == not_stop_index ? min_chg_aft : 0)
                   + ttbl.walking_time(par_st, dst));
            */
            assert(par_eat + par.dist == par.eat);
            // */
        }
        //cout <<"par "<< par_st <<" eat "<< par_eat <<"\n";
        assert(par_eat + par.dist == par.eat);
        if (reverse) print_journey(par_st, cout,
                                   min_chg_bef, min_chg_aft, par_k);
    }

    void print_longest_transfer(ST src, T t_dep, ST dst, int k,
                                std::ostream &cout = std::cout) {
        std::tuple<T, ST, ST> lg = longest_transfer(dst, n_trips[dst]);
        cout <<"longest "<< std::get<0>(lg)
             <<" from "<< std::get<1>(lg)
             <<"="<< ttbl.hub_id[std::get<1>(lg)]
             <<" to "<< std::get<2>(lg)
             <<"="<< ttbl.hub_id[std::get<2>(lg)]
             <<" in "<< src <<"="<< ttbl.hub_id[src]
             <<" at "<< t_dep
             <<" -> "<< dst <<"="<< ttbl.hub_id[dst]
             <<"\n";
    }

    std::tuple<T, ST, ST> longest_transfer(ST dst, int k = -1) {
        if (k == -1) { k = n_trips[dst]; }
        std::tuple<T, ST, ST> m{0, dst, dst}, prev{0, dst, dst};
        bool last = true, prev_is_walk = false;
        // We follow the journey backward:
        while (k > 0) { // don't consider walk at begin
            const parent_t &par = parent[k][dst];
            ST par_st = par.trip ? ttbl.stop_station[par.stop] : par.station;
            if (par.trip) {
                last = false; // don't consider walk at end
                prev_is_walk = false;
            } else { // walk
                T w = par.dist;
                if (prev_is_walk) prev = std::make_tuple(std::get<0>(prev) + w,
                                                         par_st,
                                                         std::get<2>(prev));
                else prev = std::make_tuple(w, par_st, dst);
                if ((! last) && std::get<0>(prev) > std::get<0>(m)) m = prev;
                prev_is_walk = true;
            }
            dst = par_st;
            k = par.trip ? k-1 : k;
        }
        return m;
    }

    
    void print_missing_transfers(ST dst,
                                 std::ostream &cout = std::cout, int k = -1) {
        if (k == -1) { k = n_trips[dst]; }
        T prev_wlk = 0;
        ST prev_dst = dst;
        bool last = true, prev_is_walk = false;
        // We follow the journey backward:
        while (k > 0) { // don't consider walk at begin
            const parent_t &par = parent[k][dst];
            int par_k = par.trip ? k-1 : k;
            ST par_st = par.trip ? ttbl.stop_station[par.stop] : par.station;
            if (par.trip) {
                last = false; // don't consider walk at end
                prev_is_walk = false;
            } else { // walk
                T wlk = par.dist;
                ST s = par_st, d = prev_is_walk ? prev_dst : dst;
                if (( ! last) && parent[par_k][s].trip
                              && ! ttbl.rev_transfers_id.has_edge(d, s)) {
                    cout <<"missing_transfer "<< ttbl.station_id[s]
                         <<" "<< ttbl.station_id[d]
                         <<" "<< (prev_is_walk ? wlk + prev_wlk : wlk) <<"\n";
                }
                prev_is_walk = true;
                prev_wlk = wlk;
            }
            prev_dst = dst;
            dst = par_st;
            k = par_k;
        }
    }


    pset profile(raptor &rev_rpt,
                 const ST src, const ST dst,
                 const T t_beg, const T t_end, // departure in [t_beg,t_end)
                 const bool use_hubs = true,
                 const bool use_transfers = false,
                 const T min_chg = 60, // before
                 const int k_max = ntrips_max) {
        
        T _walk_time = ttbl.t_max;
        if (use_transfers) {
            for (auto e : ttbl.transfers[src]) {
                if (e.dst == dst) _walk_time = e.wgt; 
            }
        }
        if (use_hubs) {
            _walk_time = ttbl.walking_time(src, dst);
        }
        const T walk_time = _walk_time;
        
        bool walking_is_faster = false;
        int n_tr = 0, n_dom_walk = 0, n_walk = 0;
        T prev_arr = - ttbl.t_max;
        pset profile;
        for (T t = t_beg; t < t_end; ) {
            const T arr = earliest_arrival_time(src, dst, t,
                             use_hubs, use_transfers, min_chg, 0, k_max, true);
            const int k_arr = n_trips[dst];
            /*
            int pssize = earliest_walk_pareto(src, dst, t,
                                              use_hubs, use_transfers,
                                              min_chg, k_max);
            T arr = pareto_eat_walk_to(dst).min_x(ttbl.t_max);
            const int k_arr = dst_pareto_kmax;
            */

            if (arr >= ttbl.t_max) {
                if (walking_is_faster) {++n_dom_walk;} else {++n_walk; ++n_tr;}
                break;
            }
            // For trip interval time: if (arr >= t_end) break;
            
            const T dep = - rev_rpt.earliest_arrival_time(dst, src, - arr,
                              use_hubs, use_transfers, 0, min_chg, k_max, true);
            const int k_arr_rev = n_trips[dst];

            if (dep >= t_end) break;

            /*
            std::cerr <<"\n   "<< src <<" "<< dst <<" "<< t
                      <<", hubs="<< use_hubs
                      <<", walktime="<< walk_time <<" ---\n";
            std::cerr << t <<" "<< dep <<" "<< arr <<" arr_w="<< (t + walk_time)
                      <<" : "<< k_arr <<"\n";
            print_journey(dst, std::cerr, min_chg, 0);
            std::cerr <<" --------\n";
            rev_rpt.print_journey(src, std::cerr, 0, min_chg);
            std::cerr <<" ---------\n";
            // */            
            assert(dep >= t || arr >= t + walk_time);
            assert(arr > prev_arr);


            /*
            const T arr2 = earliest_arrival_time(src, dst, dep,
                                    use_hubs, use_transfers, min_chg, 0, k_max);
            print_journey(dst, std::cerr, min_chg, 0);
            std::cerr <<" --------\n";

            if (arr != arr2) {
                std::cerr <<"\n"<< arr <<" "<< arr2 <<"\n";
                std::cerr <<" for "<< src <<" -> "<< dst <<" at "<< t <<"\n";
                std::cerr << t <<" "<< dep <<" "<< arr <<" w="<< (t + walk_time)
                          <<" : "<< k_arr <<" "<< k_arr_rev <<"\n";
            }
            assert(arr2 <= arr);
            assert(arr == arr2 || arr - dep >= walk_time);
            */


            // Update profile
            ++n_tr;
            if (arr > dep + walk_time) {
                if (walking_is_faster) {
                    ++n_dom_walk;
                } else {
                    ++n_walk;
                    //std::cout <<"walk("<< walk_time <<") ";
                }
                walking_is_faster = true;
                t = arr - walk_time + 1;
            } else {
                //std::cout << dep <<","<< arr <<" ";
                assert(profile.add(arr, - dep)); // later departure is better
                walking_is_faster = false;
                t = dep + 1;
            }
            
            prev_arr = arr;
        }

        /*
        std::cout <<"  src="<< src <<" dst="<< dst
                  <<" walk="<< walk_time
                  <<" wlk="<< n_walk <<" dom="<< n_dom_walk <<"   ";
        // */
        return profile;
    }

    


    // Return the number of paths found
    int earliest_walk_pareto(const ST src, const ST dst, const T t_dep,
                             const bool use_hubs = true,
                             const bool use_transfers = false,
                             const T min_chg_time = 60, // before
                             const int k_max = ntrips_max) {

        assert(k_max <= ntrips_max);

        for (int i = 0; i < ttbl.n_h; ++i) { all_pareto[i].clear(); }
        for (int i = 0; i < ttbl.n_st; ++i) { incr_pareto[i].clear(); }
        for (int i = 0; i < ttbl.n_s; ++i) { tmp_pareto[i].clear(); }
        for (int i = 0; i < ttbl.n_s; ++i) { prev_dep_pareto[i].clear(); }
        for (int i = 0; i <= k_max; ++i) { dst_pareto[i].clear(); }
        dst_pareto_kmax = 0;
        
        auto reach_station_trip = [this, dst](ST st, T t, T w) {
            if ((! all_pareto[dst].dominates(t, w)) // target prun
                && all_pareto[st].add(t, w)) {
                assert(incr_pareto[st].add(t, w));
                if ( ! station_has_improved[st]) {
                    improved_stations.push_back(st);
                    station_has_improved[st] = true;
                }
            }
        };

        auto reach_station_walk =
            [this, dst, t_dep, min_chg_time](ST st, T t, T w,
                                             bool self_walk=false){
            if ((! all_pareto[dst].dominates(t, w)) // target prun
                && (self_walk || all_pareto[st].add(t, w))) {
                for (S u : ttbl.station_stops.at(st)) {
                    //std::cerr << st <<" "<< u <<" : "<< t <<","<< w <<" : "<<;
                    //tmp_pareto[u].print(std::cerr);
                    if (//( ! prev_dep_pareto[u].dominates(t+min_chg_time-1, w))
                         tmp_pareto[u].add(t, w)) {
                        R r = ttbl.stop_route[u].first;
                        //assert(tmp_pareto[u].add(t, w));
                        stop_has_improved[u] = true;
                        int i = ttbl.stop_route[u].second;
                        int i_prev = route_has_improved_from[r];
                        if (i_prev == not_stop_index) {
                            improved_routes.push_back(r);
                            route_has_improved_from[r] = i;
                            route_has_improved_to[r] = i;
                        } else {
                            if (i < i_prev) route_has_improved_from[r] = i;
                            if (i > route_has_improved_to[r])
                                route_has_improved_to[r] = i;
                        }
                    }
                }
            }
        };

        assert(improved_routes.empty());
        assert(improved_stations.empty());

        reach_station_walk(src, t_dep, 0);

        if (use_transfers) {
            for (auto f : ttbl.transfers[src]) {
                reach_station_walk(f.dst, t_dep + f.wgt, f.wgt);
            }
        }
        if (use_hubs) {
            for (auto e : ttbl.outhubs[src]) {
                incr_pareto[e.dst].add(t_dep + e.wgt, e.wgt);
            }
            for (auto f : ttbl.rev_inhubs_id[dst]) {
                assert(incr_pareto[f.dst].pts.size() <= 1);
                for (auto p : incr_pareto[f.dst].pts) {
                    all_pareto[dst].add(p.x+f.wgt, p.y+f.wgt);
                }
            }
            for (auto e : ttbl.outhubs[src]) {
                incr_pareto[e.dst].clear();
                for (auto f : ttbl.inhubs[e.dst]) {
                    if (all_pareto[dst].dominates(t_dep + e.wgt+f.wgt,
                                                e.wgt+f.wgt)) break; // targpr
                    reach_station_walk(f.dst, t_dep + e.wgt+f.wgt, e.wgt+f.wgt);
                }
            }
        }

        dst_pareto[0] = all_pareto[dst];

        int k = 1;
        for (; k <= k_max; ++k) {

            // first phase: follow trips from improved routes
            if (improved_routes.empty()) { break; }
            for (const R r : improved_routes) {
                const std::vector<std::vector<std::pair<T,T> > > &trips
                    = ttbl.trips_of[r];
                const std::vector<S> &stops = ttbl.route_stops[r];
                int x_end = stops.size();
                int y_end = trips.size();
                int y = y_end;
                while (route_has_improved_from[r] != not_stop_index) {
                    assert(route_has_improved_from[r] >= 0
                           && route_has_improved_from[r] < stops.size());
                    int x_beg = route_has_improved_from[r];
                    S fst = stops[x_beg];
                    route_has_improved_from[r] = not_stop_index;
                    const std::vector<point> &pts = tmp_pareto[fst].pts;
                    assert(pts.size() > 0);
                    if (y < y_end) {
                        T lst = pts[0].x;
                        if (ttbl.stop_departures[fst][y] < lst + min_chg_time) {
                            y = y_end;
                        } else {
                            ++y; // so that y < y_prev bellow
                        }
                    }
                    for (auto p : pts) {
                        T arr = p.x, wlk = p.y;
                        int y_prev = y;
                        while (y-1 >= 0
                               && ttbl.stop_departures[fst][y-1]
                                  >= arr + min_chg_time) {
                            --y;
                        }
                        if (y < y_prev) { // improve arrival times along trip
                            for (int x = x_beg + 1; x < x_end; ++x) {
                                S u = stops[x];
                                arr = trips[y][x].first;
                                ST st = ttbl.stop_station[u];
                                reach_station_trip(st, arr, wlk);
                                /*slower:if (reach_station_trip(st, arr, wlk)) {
                                    if (y == 0) {
                                        prev_dep_pareto[u].add(0, wlk);
                                    } else {
                                        prev_dep_pareto[u]
                                            .add(trips[y-1][x].second, wlk);
                                    }
                                }*/
                                if (stop_has_improved[u]) {
                                    tmp_pareto[u].del_dominated(arr, wlk);
                                    if (tmp_pareto[u].pts.size() == 0) {
                                        stop_has_improved[u] = false;
                                        if (route_has_improved_from[r] == x) {
                                            route_has_improved_from[r] =
                                                not_stop_index;
                                        }
                                    } else {
                                        if (route_has_improved_from[r]
                                              == not_stop_index) {
                                            route_has_improved_from[r] = x;
                                        }
                                        if (tmp_pareto[u]
                                            .dominates(arr - min_chg_time,
                                                       wlk)) {
                                            break;
                                        }
                                    }
                                }
                            }
                        } else if (y == y_end) {
                            for (int x = x_beg + 1;
                                 x <= route_has_improved_to[r]; ++x) {
                                S u = stops[x];
                                ST st = ttbl.stop_station[u];
                                if (stop_has_improved[u]) {
                                    route_has_improved_from[r] = x;
                                    break;
                                }
                            }
                        }
                    }
                    tmp_pareto[fst].clear();
                    stop_has_improved[fst] = false;
                }
                assert(route_has_improved_from[r] == not_stop_index);
            }
            improved_routes.clear();

            // second phase: walk from improved stations
            if (improved_stations.empty()) { break; }
            if (use_transfers) {
                for (ST st : improved_stations) {
                    for (auto p : incr_pareto[st].pts) {
                        reach_station_walk(st, p.x, p.y, true);
                        for (auto f : ttbl.transfers[st]) {
                            if (f.dst != st) {
                                reach_station_walk(f.dst, p.x+f.wgt, p.y+f.wgt);
                            }
                        }
                    }
                }
            }
            if (use_hubs) {
                for (ST st : improved_stations) {
                    if ( ! hub_has_improved[st]) {
                        hub_has_improved[st] = true;
                        improved_hubs.push_back(st);
                    }
                    for (auto p : incr_pareto[st].pts) {
                        for (auto e : ttbl.outhubs[st]) {
                            if (all_pareto[dst].dominates(p.x+e.wgt, p.y+e.wgt))
                                break; // target prun
                            if (e.dst != st
                                && all_pareto[e.dst].add(p.x+e.wgt, p.y+e.wgt)){
                                if ( ! hub_has_improved[e.dst]) {
                                    hub_has_improved[e.dst] = true;
                                    improved_hubs.push_back(e.dst);
                                }
                                incr_pareto[e.dst].add(p.x+e.wgt, p.y+e.wgt);
                            }
                        }
                    }
                    station_has_improved[st] = false;
                }
                improved_stations.clear();
                for (ST h : improved_hubs) {
                    for (auto p : incr_pareto[h].pts) {
                        if (h < ttbl.n_st) {
                            reach_station_walk(h, p.x, p.y, true);
                        }
                        for (auto f : ttbl.inhubs[h]) {
                            if (all_pareto[dst].dominates(p.x+f.wgt, p.y+f.wgt))
                                break; // target prun
                            if (f.dst != h && f.dst < ttbl.n_st) {
                                reach_station_walk(f.dst, p.x+f.wgt, p.y+f.wgt);
                            }
                        }
                    }
                    hub_has_improved[h] = false;
                    incr_pareto[h].clear();
                }
                improved_hubs.clear();
            } else {
                for (ST st : improved_stations) {
                    station_has_improved[st] = false;
                    incr_pareto[st].clear();
                }
                improved_stations.clear();
            }
            
            dst_pareto[k] = all_pareto[dst];
            if (dst_pareto[k].size() > 0) dst_pareto_kmax = k;
        }
        if (k <= k_max) {
            dst_pareto[k] = all_pareto[dst];
            if (dst_pareto[k].size() > 0) dst_pareto_kmax = k;
        }

        if (k > k_max) {
            for (const R r : improved_routes) {
                const std::vector<S> &stops = ttbl.route_stops[r];
                assert(route_has_improved_from[r] != not_stop_index);
                for (int x=route_has_improved_from[r];
                     x <= route_has_improved_to[r]; ++x) {
                    if (stop_has_improved[stops[x]]) {
                        stop_has_improved[stops[x]] = false;
                        tmp_pareto[stops[x]].clear();
                    }
                }
                route_has_improved_from[r] = not_stop_index;
            }
            improved_routes.clear();
            for (ST st : improved_stations) {
                station_has_improved[st] = false;
                incr_pareto[st].clear();
            }
            improved_stations.clear();
        }
 
        size_t ps_size = 0, m_size = 0;
        for (int i = 0; i <= std::min(k, k_max); ++i) {
            m_size = std::max(m_size, dst_pareto[i].pts.size());
            if (i == 0) { ps_size += dst_pareto[i].pts.size(); }
            else {
                for (auto p : dst_pareto[i].pts) {
                    if ( ! dst_pareto[i-1].dominates(p.x, p.y)) ++ps_size;
                }
            }
        }
        /*
        std::cout <<"k "<< k <<"\n";
        std::cout <<"ps_size "<< ps_size <<"\n";
        std::cout <<"m_size "<< m_size <<"\n";
        */

        return ps_size;
        /*
        if (all_pareto[dst].pts.size() > 0) {
            return all_pareto[dst].pts[0].x;
        } else {
            return ttbl.t_max;
        }
        */
    }

    pset pareto_eat_walk_to(ST st) {
        return all_pareto[st];
    }

    T test(int n_q, T t_beg, T t_end) {
        assert(t_beg < t_end);
        uint64_t sum = 0, n_reached = 0;
        for (int i = 0; i < n_q; ++i) {
            ST src = rand() % ttbl.n_st;
            ST dst = rand() % ttbl.n_st;
            T t = t_beg + rand() % (t_end - t_beg);
            T arr = earliest_arrival_time(src, dst, t);
            if (arr < ttbl.t_max) {
                ++n_reached;
                sum += arr - t;
            }
        }
        std::cerr << n_reached <<" reached\n";
        return (T) (sum / n_reached);
    }
};


#endif // RAPTOR_HH
