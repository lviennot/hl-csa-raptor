#ifndef TIMETABLE_HH
#define TIMETABLE_HH

#include <assert.h>
#include <limits>
#include <iostream>
#include <zlib.h>
#include <algorithm>
#include <vector>
#include <set>
#include <cstdarg>
#include <cmath>
#include <string>
#include <unordered_map>
#include <map>

#include "mgraph.hh"
#include "traversal.hh"
//#include "pruned_landmark_labeling.hh"
#include "file_util.hh"

/*******
 * Data structures for representing timetables and walking transfers of
 * a public transit network.
 *
 * What is usually called a stop in GTFS format is called station here.
 * A station can have multiple stops, each belonging to exactly one route.
 * We ensure that all trips in a route have same sequence of stops, and that
 * a trip never overpass another trip of the same route (this may require
 * to split some routes in several routes, but rarely in practice).
 */

class timetable {
public:
    typedef int ST; // stations
    typedef int S; // stops
    typedef int R; // routes
    typedef int T; // time

    const int t_min = 0, t_max = 1000000000;//std::numeric_limits<int>::max()/2;

    typedef std::string id;
    typedef mgraph<ST, T> graph;

    size_t n_st, n_tr, n_s, n_r, n_h; // number of stations, transfer nodes, stops, routes, hubs
    std::vector<std::vector<S> > station_stops; // stops of a station
    std::vector<id> station_id, hub_id; // id from the gtfs data
    std::vector<ST> stop_station;
    std::vector<std::pair<R, int> > stop_route; // route at a stop, and index of stop in the route sequence
    std::vector<std::vector<T> > stop_departures; // departure time
    std::vector<std::vector<T> > stop_arrivals; // arrival time
    std::vector<std::vector<S> > route_stops; // stop sequence of a route
    std::vector<std::vector<std::vector<std::pair<T,T> > > > trips_of; // trips of a route : arrival/departure times at each stop
    graph transfers, inhubs, rev_inhubs, outhubs, lowerboundgraph; // weight sorted
    graph rev_transfers_id, rev_inhubs_id, outhubs_id; // auxiliary graphs
    
    std::map<id, ST> id_to_station, id_to_hub;

public:
    timetable(std::string day, std::string date,
              std::string calendar, std::string calendar_dates,
              std::string tripsfile, std::string stop_times,
              std::string transfersfile, bool symmetrize = false)
        : n_st(0), n_tr(0), n_s(0), n_r(0), n_h(0)
    {
        auto services = services_at(day, date, calendar, calendar_dates);
        std::cerr << services.size() <<" services\n";
        auto trips = trips_served(services, tripsfile);
        std::cerr << trips.size() <<" trips\n";
        auto trip_seq = trips_sequences(stop_times, trips);
        int events = 0;
        for (auto seq : trip_seq) events += seq.second.size();
        std::cerr << events <<" events\n";
        auto transf = station_transfers(transfersfile);
        build(trip_seq, transf, symmetrize);
        build_auxiliary_graphs();
    }

    timetable(std::string stop_times,
              std::string transfersfile, bool symmetrize = false)
        : n_st(0), n_tr(0), n_s(0), n_r(0), n_h(0)
    {
        auto trip_seq = trips_sequences_oneday(stop_times);
        int events = 0;
        for (auto seq : trip_seq) events += seq.second.size();
        std::cerr << events <<" events\n";
        auto transf = station_transfers_2(transfersfile);
        build(trip_seq, transf, symmetrize);
        build_auxiliary_graphs();
    }

    timetable(std::string stop_times,
              std::string inhubsfile, std::string outhubsfile,
              std::string transfersfile, bool symmetrize = false,
              std::string walkingfile="", T t_from=0, T t_to=0)
        : n_st(0), n_tr(0), n_s(0), n_r(0), n_h(0)
    {
        auto trip_seq = trips_sequences_oneday(stop_times);
        int events = 0;
        for (auto seq : trip_seq) events += seq.second.size();
        std::cerr << events <<" events\n";
        auto transf = station_transfers_2(transfersfile);
        build(trip_seq, transf, symmetrize);
        // hubs :
        assert(n_st + n_tr == station_id.size());
        hub_id.reserve(n_st + n_tr);
        for (int st = 0; st < n_st + n_tr; ++st) {
            hub_id.push_back(station_id[st]);
            id_to_hub[station_id[st]] = st;
        }
        n_h = hub_id.size();
        auto get_hub = [this](std::string h) -> ST {
            auto search = id_to_hub.find(h);
            if (search == id_to_hub.end()) {
                assert(n_h == hub_id.size());
                hub_id.push_back(h);
                id_to_hub[h] = n_h;
                return n_h++;
            }
            return search->second;
        };
        // in-hubs :
        std::vector<graph::edge> edg;
        std::vector<bool> seen(n_st, false);
        auto rows = read_tuples(inhubsfile, 4);
        edg.reserve(rows.size());
        for (auto r : rows) {
            auto st = id_to_station.find(r[2]);
            if (st != id_to_station.end()) {
                ST h = get_hub(r[1]);
                ST s = st->second;
                seen[s] = true;
                T delay = std::stoi(r[3]);
                edg.push_back(graph::edge(h, s, delay));
            }
        }
        for (int st = 0; st < n_st; ++st)
            if ( ! seen[st]) edg.push_back(graph::edge(st, st, 0));
        std::sort(edg.begin(), edg.end(), [](const graph::edge &e,
                                             const graph::edge &f) {
                      return e.wgt < f.wgt;
                  });
        inhubs.set_edges(edg, n_h);
        std::cerr << inhubs.n() <<" in-hubs, avg in-degree "
                  << (inhubs.m() / n_st)
                  <<", avg out-degree "<< inhubs.m() / n_h
                  <<", max out-degree "<< inhubs.max_degree()
                  <<"\n";
        // out-hubs :
        edg.clear();
        for (int st = 0; st < n_st; ++st) seen[st] = false;
        rows = read_tuples(outhubsfile, 4);
        edg.reserve(rows.size());
        for (auto r : rows) {
            auto st = id_to_station.find(r[1]);
            if (st != id_to_station.end()) {
                ST s = st->second;
                auto h = id_to_hub.find(r[2]);
                if (h != id_to_hub.end()) {
                    T delay = std::stoi(r[3]);
                    edg.push_back(graph::edge(s, h->second, delay));
                }
            }
        }
        for (int st = 0; st < n_st; ++st)
            if ( ! seen[st]) edg.push_back(graph::edge(st, st, 0));
        std::sort(edg.begin(), edg.end(), [](const graph::edge &e,
                                             const graph::edge &f) {
                      return e.wgt < f.wgt;
                  });
        outhubs.set_edges(edg, n_h);
        std::cerr << outhubs.n() <<" out-hubs, avg out-degree "
                  << (outhubs.m() / n_st)
                  <<", max out-degree "<< outhubs.max_degree()
                  <<", avg in-degree "<< outhubs.m() / n_h
                  <<"\n";
        std::cerr << n_h <<" hubs\n";

        // Check weight sort:
        for (ST u : outhubs) {
            T dt = 0;
            for (auto e : outhubs[u]) { assert(dt <= e.wgt); dt = e.wgt; }
        }
        for (ST u : inhubs) {
            T dt = 0;
            for (auto e : inhubs[u]) { assert(dt <= e.wgt); dt = e.wgt; }
        }
        
        // lower-bound graph
        if (walkingfile != "") {
            size_t n_lb = n_h;
            std::map<id, ST> id;
            auto get_node = [this, &n_lb, &id](std::string u) -> ST {
                auto search = id_to_hub.find(u);
                if (search != id_to_hub.end()) {
                    return search->second;
                } // else
                search = id.find(u);
                if (search != id.end()) {
                    return search->second;
                } // else
                id[u] = n_lb;
                return n_lb++;
            };

            std::vector<graph::edge> edg;
            auto rows = read_tuples(walkingfile, 3);
            edg.reserve(rows.size() + 100000);
            for (auto r : rows) {
                ST u = get_node(r[0]), v = get_node(r[1]);
                T delay = std::stoi(r[2]);
                edg.push_back(graph::edge(u, v, delay));
            }

            // fastest trip for each connection
            if (t_to == 0) t_to = t_max;
            for (R r = 0; r < n_r; ++r) {
                const std::vector<std::vector<std::pair<T,T> > > &trips
                    = trips_of[r];
                const std::vector<S> &stops = route_stops[r];
                for (int x = 1; x < stops.size(); ++x) {
                    ST u = stop_station[stops[x-1]];
                    ST v = stop_station[stops[x]];
                    T t = t_max;
                    for (int y = 0; y < trips.size(); ++y) {
                        if (trips[y][x-1].second >= t_from
                            && trips[y][x].first <= t_to) {
                            T dt = trips[y][x].first - trips[y][x-1].second;
                            t = std::min(t, dt);
                        }
                    }
                    if (t < t_max) {
                        edg.push_back(graph::edge(u, v, t));
                    }
                }
            }

            lowerboundgraph.set_edges(edg);
            std::cerr <<"lower-bound graph "<< lowerboundgraph.n()
                      <<" nodes, "<< lowerboundgraph.m() <<" edges\n";
        }

        build_auxiliary_graphs();
    }

    // Hack for last departure time requests through EAT query:
    void reverse_time() {
        for (R r = 0; r < n_r; ++r) {
            // stop sequence :
            rev_vector(route_stops[r]);
            for (int i = 0; i < route_stops[r].size(); ++i) {
                S s = route_stops[r][i];
                stop_route[s] = std::make_pair(r, i);
            }
            //trips :
            rev_vector(trips_of[r]);
            for (int a = 0; a < trips_of[r].size(); ++a) {
                std::vector<std::pair<T,T> > &trip = trips_of[r][a];
                rev_vector(trip);
                for (int i = 0; i < trip.size(); ++i) {
                    T rev_dep = - trip[i].first, rev_arr = - trip[i].second;
                    trip[i] = std::make_pair(rev_arr, rev_dep);
                    S s = route_stops[r][i];
                    stop_departures[s][a] = rev_dep;
                    stop_arrivals[s][a] = rev_arr;
                }
            }
        }
        // graphs
        transfers = transfers.reverse();
        transfers.sort_neighbors_by_weight();
        graph tmpin = inhubs;
        inhubs = outhubs.reverse();
        inhubs.sort_neighbors_by_weight();
        outhubs = tmpin.reverse();
        outhubs.sort_neighbors_by_weight();
        check();
        build_auxiliary_graphs();
    }

private:
    void build(std::unordered_map<id, std::vector<std::tuple<T, T, id, int> > > &trip_seq,
               std::vector<std::tuple<id, id, T> > &transf,
               bool symmetrize) {
        auto create_station = [this](id st) {
            if (id_to_station.find(st) == id_to_station.end()) {
                id_to_station[st] = n_st++;
                station_id.push_back(st);
                station_stops.push_back(std::vector<S>{});
            }
        };

        std::vector<std::vector<std::tuple<T, T, id, int> > > inter;
        inter.reserve(trip_seq.size());
        for (auto seq : trip_seq) {
            auto s = seq.second;
            // sort trip according to stop_sequence:
            std::sort(s.begin(), s.end(),
                      [](const std::tuple<T, T, id, int> &a,
                         const std::tuple<T, T, id, int> &b) {
                          return std::get<3>(a) < std::get<3>(b);
                      });
            // check time increases along trip:
            for (int i = 0; i < s.size(); ++i) {
                const std::tuple<T, T, id, int> &a = s[i];
                assert(std::get<0>(a) <= std::get<1>(a));
                if (i+1 < s.size()) {
                    const std::tuple<T, T, id, int> &b = s[i+1];
                    if (std::get<1>(a) > std::get<0>(b)) {
                        std::cerr << "decr time in "<< seq.first <<" : "
                                  << std::get<2>(a) <<","<< std::get<3>(a) <<"\n";
                    }
                    assert(std::get<1>(a) <= std::get<0>(b));
                }
            }
            inter.push_back(s);
        }
        // sort trips
        std::sort(inter.begin(), inter.end(),
                  [](const std::vector<std::tuple<T, T, id, int> > &a,
                     const std::vector<std::tuple<T, T, id, int> > &b) {
                      return std::get<1>(a.at(0)) < std::get<1>(b.at(0));
                  });
        // construct timetables
        std::map<std::vector<id>, R> stations_to_route, stations_to_route2;
        int n_overpass = 0;
        for (auto s : inter) {
            // partition trips into routes with same stop sequence
            std::vector<id> stations;
            stations.reserve(s.size());
            for (int i = 0; i < s.size(); ++i) {
                const std::tuple<T, T, id, int> &a = s[i];
                stations.push_back(std::get<2>(a));
            }
            R rte = -1;
            bool overpass = false, new_route = true;
            auto search = stations_to_route.find(stations);
            if (search != stations_to_route.end()) {
                new_route = false;
                rte = search->second;
                assert(route_stops[rte].size() == s.size());
                const auto &prev = trips_of[rte].back();
                assert(prev[0].second <= std::get<1>(s[0]));
                for (int i = 0; i < s.size(); ++i) {
                    T arr = std::get<0>(s[i]);
                    T dep = std::get<1>(s[i]);
                    if (arr < prev[i].first || dep < prev[i].second) {
                        overpass = true;
                        break;
                    }
                }
            }
            if (overpass) {
                search = stations_to_route2.find(stations);
                if (search != stations_to_route2.end()) {
                    new_route = false;
                    rte = search->second;
                    assert(route_stops[rte].size() == s.size());
                } else {
                    new_route = true;
                    ++n_overpass;
                }
            }
            if (new_route) {
                rte = n_r++;
                if (overpass) { stations_to_route2[stations] = rte; }
                else { stations_to_route[stations] = rte; }
                // create stations:
                for (auto st : stations) {
                    create_station(st);
                }
                // create a stop for each station in the new route:
                std::vector<S> stops(stations.size());
                for (int i = 0; i < stations.size(); ++i) {
                    ST st = id_to_station[stations[i]];
                    S s = n_s++;
                    station_stops[st].push_back(s);
                    stop_station.push_back(st);
                    stops[i] = s;
                    stop_route.push_back(std::make_pair(rte, i));
                }
                route_stops.push_back(stops);
                trips_of.push_back(std::vector<std::vector<std::pair<T,T> > >{});
            }
            // create trip in route table:
            std::vector<std::pair<T,T> > trp{s.size()};
            for (int i = 0; i < s.size(); ++i) {
                T arr = std::get<0>(s[i]);
                T dep = std::get<1>(s[i]);
                trp[i] = std::make_pair(arr, dep);
            }
            trips_of[rte].push_back(trp);
        }
        std::cerr << n_st <<" stations, "<< n_s <<" stops, "
                  << n_r <<" routes ("<< n_overpass <<" for overpasses)\n";

        //for (ST u = 0; u < n_st; ++u) std::cout << station_id[u] <<"\n";
        
        n_overpass = 0;
        for(R rte = 0; rte < n_r; ++rte) {
            std::sort(trips_of[rte].begin(), trips_of[rte].end(),
                      [](const std::vector<std::pair<T,T> > &p,
                         const std::vector<std::pair<T,T> > &q) {
                          for (int i = 0; i < p.size(); ++i) {
                              if (p[i].second < q[i].second) // departs before
                                  return true;
                              if (p[i].second > q[i].second) // departs after
                                  return false;
                          }
                          return false; // equality
                      });
            for(int i = 1; i < trips_of[rte].size() ; ++i) {
                // check trip does not overpass:
                bool overpass = false;
                const auto &prev = trips_of[rte][i-1];
                auto &trp = trips_of[rte][i];
                for (int i = 0; i < trp.size(); ++i) {
                    if (trp[i].second < prev[i].second // departs before
                        || trp[i].first < prev[i].first) { // arrives before
                            overpass = true;
                            // quick fix :
                            trp[i].second = std::max(trp[i].second,
                                                     prev[i].second);
                            trp[i].first = std::max(trp[i].first,
                                                    prev[i].first);
                    }
                }
                if(overpass) {
                    ++n_overpass;
                    std::cerr << "overpass in " << rte <<" : ";
                    for (S u : route_stops[rte]) {
                        std::cerr << station_id[stop_station[u]] <<" ";
                    }
                    std::cerr <<"\n";
                }
            }
        }
        if (n_overpass > 0) {
            std::cerr <<"WARNING: timetable modified to fix "
                      << n_overpass <<" remaining overpasses\n";
        }

        stop_departures.reserve(n_s);
        stop_arrivals.reserve(n_s);
        for (S s = 0; s < n_s; ++s) {
            R rte = stop_route[s].first;
            int is = stop_route[s].second;
            assert(is < route_stops[rte].size());
            int ntrips = trips_of[rte].size();
            std::vector<T> deps, arrs;
            deps.reserve(ntrips);
            arrs.reserve(ntrips);
            for (int i = 0; i < ntrips; ++i) {
                assert(trips_of[rte][i].size() == route_stops[rte].size());
                arrs.push_back(trips_of[rte][i][is].first);
                deps.push_back(trips_of[rte][i][is].second);
            }
            stop_arrivals.push_back(arrs);
            stop_departures.push_back(deps);
        }

        // transfers :
        std::vector<graph::edge> st_transf;
        st_transf.reserve(transf.size());
        assert(n_st == station_id.size());
        for (auto tr : transf) {
            id from = std::get<0>(tr);
            id to = std::get<1>(tr);
            ST st_from, st_to;
            if (id_to_station.find(from) == id_to_station.end()) {
                st_from = n_st + n_tr++;
                station_id.push_back(from);
                id_to_station[from] = st_from;
            } else { st_from = id_to_station[from]; }
            if (id_to_station.find(to) == id_to_station.end()) {
                st_to = n_st + n_tr++;                
                station_id.push_back(to);
                id_to_station[to] = st_to;
            } else { st_to = id_to_station[to]; }
            T delay = std::get<2>(tr);
            st_transf.push_back(graph::edge(st_from, st_to, delay));
        }
        transfers.set_edges(st_transf, n_st + n_tr);
        if (symmetrize) { // add missing reverse links
            graph r = transfers.reverse();
            for (ST u : r) {
                for (auto e : r[u]) {
                    if ( ! r.has_edge(e.dst, u)) {
                        st_transf.push_back(graph::edge(u, e.dst, e.wgt));
                    }
                }
            }
            transfers.set_edges(st_transf, n_st + n_tr);
        }
        //transfers = transfers.simple();
        std::cerr << transfers.m() << " transfers with "
                  << n_st <<" + "<< n_tr <<" nodes ";
        size_t n_asym = transfers.asymmetry(false);
        if (n_asym == 0) {
            n_asym = transfers.asymmetry(true);
            std::cerr <<"(symmetric: "
                      << n_asym <<" reverse weights differ)\n";
        } else {
            std::cerr <<"(asymmetric: "
                      << n_asym <<" reverse links miss)\n";
        }
        transfers.sort_neighbors_by_weight();

        check();

    }

    
    void build_auxiliary_graphs() {
        size_t asym = transfers.asymmetry(false); 
        std::cerr << transfers.m() <<" transfers: "
                  << asym << " reverse links miss, "
                  << transfers.asymmetry(true) <<" reverse weights differ\n";
        std::cerr <<"  avg degree "<< (1.0*transfers.m()/transfers.n())
                  <<", max degree "<< transfers.max_degree() <<"\n"; 

        // transitive closure of transfer graph:
        std::vector<graph::edge> transf;
        traversal<graph> trav(transfers.n());
        //std::cout << "from_stop_id,to_stop_id,min_transfer_time\n";
        for (ST st = 0; st < n_st; ++st) {
            trav.clear();
            trav.dijkstra(transfers, st);
            for (int i = 0; i < trav.nvis(); ++i) { // include self loop
                ST ot = trav.visit(i);
                T t = trav.dist(ot);
                if (ot < n_st) {
                    transf.push_back(graph::edge(st, ot, t));
                    //std::cout << station_id[st] <<","<< station_id[ot] <<","<< t <<"\n";
                }
            }
        }
        transfers.set_edges(transf, n_st);
        asym = transfers.asymmetry(false); 
        std::cerr << transfers.m() <<" transitive transfers: "
                  << asym << " reverse links miss, "
                  << transfers.asymmetry(true) <<" reverse weights differ\n";
        std::cerr <<"  avg degree "<< (1.0*transfers.m()/transfers.n())
                  <<", max degree "<< transfers.max_degree() <<"\n"; 
        assert(asym <= 100); // strange network otherwise

        rev_transfers_id = transfers.reverse(); // sort by ID

        if (outhubs.n() > 0 || inhubs.n() > 0) {
            rev_inhubs_id = inhubs.reverse(); // sorted by ID, not weight
            std::cerr <<"rev_inhubs avg degree "
                      <<(1.0*rev_inhubs_id.m()/rev_inhubs_id.n())
                      <<", max degree "<< rev_inhubs_id.max_degree() <<"\n";

            rev_inhubs = rev_inhubs_id;
            rev_inhubs.sort_neighbors_by_weight();
            
            // Check hub distances vs transfers
            outhubs_id = outhubs.reverse().reverse(); // ID sorted
            
            for (auto u : transfers) {
                for (auto e : transfers[u]) {
                    T t = walking_time(u, e.dst);
                    //std::cerr << ttbl.hub_id[u] <<" -> "<< ttbl.hub_id[e.dst]
                    //          <<": "<< t <<","<< e.wgt <<"\n";
                    assert(t <= e.wgt);
                }
            }
        }


        // ---------------- lower bound graph ---------------------        
        //pll_lb(tt.lowerboundgraph)
        //pll_lb.print_stats(std::cerr);
    }

public:
    
    T walking_time(ST u, ST v) const {
        assert(u < n_st && v < n_st);
        auto uh = outhubs_id[u].begin();
        auto ue = outhubs_id[u].end();
        auto hv = rev_inhubs_id[v].begin();
        auto ve = rev_inhubs_id[v].end();
        T t = t_max;
        while (uh != ue && hv != ve) {
            if (uh->dst < hv->dst) { ++uh; }
            else if (uh->dst > hv->dst) { ++hv; }
            else { // uh->dst == hv-->dst
                T t_uhv = uh->wgt + hv->wgt;
                if (t_uhv <t) t = t_uhv;
                ++uh; ++hv;
            }
        }
        return t;
    }
    
    std::string walking_time_str(ST u, ST v) const {
        T t = t_max;
        ST best_hub = -1;
        if (u >= n_st) {
            t = rev_inhubs_id.edge_weight(v, u);
            best_hub = u;
        } else if (v >= n_st) {
            t = outhubs_id.edge_weight(u, v);
            best_hub = v;
        } else {
            auto uh = outhubs_id[u].begin();
            auto ue = outhubs_id[u].end();
            auto hv = rev_inhubs_id[v].begin();
            auto ve = rev_inhubs_id[v].end();
            while (uh != ue && hv != ve) {
                if (uh->dst < hv->dst) { ++uh; }
                else if (uh->dst > hv->dst) { ++hv; }
                else { // uh->dst == hv-->dst
                    T t_uhv = uh->wgt + hv->wgt;
                    if (t_uhv <t) { t = t_uhv; best_hub = uh->dst; }
                    ++uh; ++hv;
                }
            }
        }
        if (t == t_max) return std::string("infinity");
        return std::to_string(t) + " (hub=" + std::to_string(best_hub) +")";
    }
    

    void check() {
        assert(station_stops.size() == n_st);
        assert(stop_station.size() == n_s);
        for (ST st = 0; st < n_st; ++st) {
            for (S u : station_stops[st]) { assert(stop_station[u] == st); }
        }
        assert(trips_of.size() == n_r);
        assert(route_stops.size() == n_r);
        assert(stop_route.size() == n_s);
        assert(stop_departures.size() == n_s);
        assert(stop_arrivals.size() == n_s);
        for (S u = 0; u < n_s; ++u) {
            R r = stop_route[u].first;
            int i = stop_route[u].second;
            assert(stop_departures[u].size() == trips_of[r].size());
            assert(stop_arrivals[u].size() == trips_of[r].size());
            assert(i < route_stops[r].size());
            assert(route_stops[r][i] == u);
            for (int p = 0; p < trips_of[r].size(); ++p) {
                assert(trips_of[r][p][i].first == stop_arrivals[u][p]);
                assert(trips_of[r][p][i].second == stop_departures[u][p]);
            }
        }
        for (R r = 0; r < n_r; ++r) {
            for (S u : route_stops[r]) {
                assert(u >= 0 && u < n_s);
                assert(stop_route[u].first == r);
            }
        }
    }

    
    static std::set<id> services_at(std::string day, std::string date,
                                    std::string calendar,
                                    std::string calendar_dates) {
        std::set<id> services;
        // read calendar.txt file
        for (auto r : read_csv(calendar, 4,
                               "service_id", day.c_str(),
                               "start_date", "end_date")) {
            if (r[1] == "1" && r[2] <= date && date <= r[3]) {
                services.insert(r[0]);
            }
        }
        // read calendar_dates.txt file
        if (std::string{""} != calendar_dates) {
            for (auto r : read_csv(calendar_dates, 3,
                                   "service_id", "date", "exception_type")) {
                if (r[1] == date) {
                    if (r[2] == "1") services.insert(r[0]);
                    else if (r[2] == "2") services.erase(r[0]);
                    else assert(false);
                }
            }
        }
        return services;
    }

    static std::set<id> trips_served(const std::set<id> &services,
                                     std::string tripsfile) {
        std::set<id> trips;
        // read trips.txt file
        for (auto r : read_csv(tripsfile, 2, "service_id", "trip_id")) {
            if (services.find(r[0]) != services.end()) trips.insert(r[1]);
        }
        return trips;
    }

    static std::vector<std::tuple<id, id, T> >
    station_transfers(std::string transfersfile) {
        std::vector<std::tuple<id, id, T> > transf;
        for (auto r : read_csv(transfersfile, 4,
                               "from_stop_id", "to_stop_id", "transfer_type",
                               "min_transfer_time")) {
            assert(r[2] == "2");
            T t = std::stoi(r[3]);
            transf.push_back(std::make_tuple(r[0], r[1], t));
        }
        return transf;
    }

    static std::vector<std::tuple<id, id, T> >
    station_transfers_2(std::string transfersfile) {
        std::vector<std::tuple<id, id, T> > transf;
        for (auto r : read_csv(transfersfile, 3,
                               "from_stop_id", "to_stop_id",
                               "min_transfer_time")) {
            T t = std::stoi(r[2]);
            transf.push_back(std::make_tuple(r[0], r[1], t));
        }
        return transf;
    }

    static std::unordered_map<id, std::vector<std::tuple<T, T, id, int> > >
    trips_sequences(std::string stop_times,
                    const std::set<id> &trips = std::set<id>()) {
        std::unordered_map<id, std::vector<std::tuple<T, T, id, int> > >
            trip_seq;
        for (auto r : read_csv(stop_times, 5,
                                "trip_id", "arrival_time", "departure_time",
                                "stop_id", "stop_sequence")) {
            std::string trp = r[0];
            if (trips.empty() || trips.find(trp) != trips.end()) {
                T arr = time_of_string(r[1]);
                T dep = time_of_string(r[2]);
                id stp = r[3];
                int seq = std::stoi(r[4]);
                trip_seq[trp].push_back(std::make_tuple(arr, dep, stp, seq));
            }
        }
        return trip_seq;
    }

    static std::unordered_map<id, std::vector<std::tuple<T, T, id, int> > >
    trips_sequences_oneday(std::string stop_times) {
        std::unordered_map<id, std::vector<std::tuple<T, T, id, int> > >
            trip_seq;
        for (auto r : read_csv(stop_times, 5,
                                "trip_id", "arrival_time", "departure_time",
                                "stop_id", "stop_sequence")) {
            std::string trp = r[0];
            {
                T arr = std::stoi(r[1]);
                T dep = std::stoi(r[2]);
                id stp = r[3];
                int seq = std::stoi(r[4]);
                trip_seq[trp].push_back(std::make_tuple(arr, dep, stp, seq));
            }
        }
        return trip_seq;
    }

    static T time_of_string(const std::string &str) {
        auto v = split(str, ':');
        if (v.size() != 3) std::cerr << "buggy time : "<< str <<"\n";
        assert(v.size() == 3);
        int h = std::stoi(v[0]), m = std::stoi(v[1]), s = std::stoi(v[2]);
        return h * 3600 + m * 60 + s;
    }
    
public:

    static T decimeters_to_seconds(const int d) {
        return (T)(std::lround(3600.0 * d / (4 * 10000))); // 4 Km/h
    }
};


#endif // TIMETABLE_HH
