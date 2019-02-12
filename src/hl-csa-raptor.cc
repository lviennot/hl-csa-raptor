#include <iostream>

#include "timetable.hh"
#include "raptor.hh"
#include "connection_scan.hh"
#include "logging.hh"
#include "file_util.hh"

void usage_exit (char **argv) {
    auto paragraph = [](std::string s, int width=80) -> std::string {
        std::string acc;
        while (s.size() > 0) {
            int pos = s.size();
            if (pos > width) pos = s.rfind(' ', width);
            std::string line = s.substr(0, pos);
            acc += line + "\n";
            s = s.substr(pos);
        }
        return acc;
    };
    
    std::cerr <<"Usage: "<< argv[0]
              <<" [options] gtfs_dir\n"
              <<"\n";
    exit(1);
}

std::vector<std::string> get_args(int argc, char **argv) {
    std::vector<std::string> a;
    for (int i = 0; i < argc; ++i) {
        std::string s(argv[i]);
        if (s[0] == '-') continue; // option
        a.push_back(s);
    }
    return a;
}

bool has_opt(int argc, char **argv, std::string opt) {
    assert(opt[0] == '-');
    for (int i = 0; i < argc; ++i) {
        std::string s(argv[i]);
        if (s == opt) return true;
    }
    return false;
}

std::string get_opt(int argc, char **argv,
                    std::string opt_prefix, std::string dft) {
    int len = opt_prefix.size();
    assert(opt_prefix[0] == '-'
           && (opt_prefix[len-1] == '=' || opt_prefix[len-1] == ':'));
    for (int i = 0; i < argc; ++i) {
        std::string s(argv[i]);
        if (s.size() >= len && s.substr(0, len) == opt_prefix) {
            return s.substr(len);
        }
    }
    return dft;
}

int get_int_opt(int argc, char **argv,
                std::string opt_prefix, int dft) {
    return std::stoi(get_opt(argc, argv, opt_prefix, std::to_string(dft)));
}

typedef pareto_rev<int> pset;

int main (int argc, char **argv) {
    logging main_log("--");

    // ------------------------ usage -------------------------
    std::vector<std::string> args = get_args(argc, argv);
    if (args.size() != 2) {
        usage_exit(argv);
    }
    std::string dir{args[1]};
    std::cerr <<"--------------- "<< dir <<" ---------------------\n";
    dir += "/";

    // ------------------------ time -------------------------
    main_log.cerr() << "start\n";
    double t = main_log.lap();

    //std::cerr <<"first rands: "<< rand() <<" "<< rand() <<" "<< rand() <<"\n";


    // ------------------------- load timetable ----------------------
    /*
    timetable ttbl{args[1], args[2],
            dir+"calendar.txt", dir+"calendar_dates.txt",
            dir+"trips.txt", dir+"stop_times.txt", dir+"transfers.txt", true};
    */
    timetable ttbl{dir+"stop_times.csv.gz",
            dir+"walk_and_transfer_inhubs.gr.gz",
            dir+"walk_and_transfer_outhubs.gr.gz",
            dir+"transfers.csv.gz", true};
    //dir+"walking_and_transfers.gr", t_from, t_to};
    std::cerr << ttbl.n_r <<" routes, "<< ttbl.n_st <<" sations, "
              << ttbl.n_s <<" stops\n";
    main_log.cerr(t) << "timetable\n";
    t = main_log.lap();
    //exit(0);

    // --------------- earliest arrival time through Raptor and CSA ---------
    raptor rpt(ttbl);
    main_log.cerr(t) << "raptor initialized\n";
    t = main_log.lap();

    timetable rev_ttbl(ttbl);
    rev_ttbl.check();
    rev_ttbl.reverse_time();
    raptor rev_rpt(rev_ttbl);
    main_log.cerr(t) << "rev raptor initialized\n";
    t = main_log.lap();

    connection_scan csa(ttbl);
    main_log.cerr(t) << "csa initialized\n";
    t = main_log.lap();

    const bool hub=true, trf=false;
    const int chg=std::stoi(get_opt(argc, argv, "-min-change-time=", "60")),
        km=48;

    uint64_t sum = 0, n_ok = 0;


    // ------------------ read queries -------------------
    int n_q = 0;
    std::vector<std::tuple<int, int, int> > queries;
    if (get_opt(argc, argv, "-query-file=", "queries-unif.csv") != "") {
        auto rows = read_csv
            (dir + get_opt(argc, argv, "-query-file=", "queries-unif.csv"), 3,
             "source", "destination", "departure_time");
             //"source", "target", "time");
        int n_q_max = std::stoi(get_opt(argc, argv, "-nq=", "10000"));
        for (auto r : rows) {
            if (n_q >= n_q_max) break;
            if (has_opt(argc, argv, "-1") && n_q >= 1) break;
            if (has_opt(argc, argv, "-10") && n_q >= 10) break;
            if (has_opt(argc, argv, "-100") && n_q >= 100) break;
            int src = ttbl.id_to_station[r[0]];
            int dst = ttbl.id_to_station[r[1]];
            int t = std::stoi(r[2]);
            queries.push_back(std::make_tuple(src, dst, t));
            ++n_q;
        }
        main_log.cerr(t) << n_q << " queries\n";
        t = main_log.lap();
    }
    // */


    // make andom successful queries
    if (get_opt(argc, argv, "-random-queries=", "") != "") {
        n_q = std::stoi(get_opt(argc, argv, "-random-queries=", "1000"));
        int max_delay = std::stoi(get_opt(argc, argv, "-max-delay=",
                                          std::to_string(ttbl.t_max)));
        int t_beg = 0*3600, t_end = 24*3600;
        int n_try = 0, n_err = 0;
        while (queries.size() < n_q) {
            ++n_try;
            int src = rand() % ttbl.n_st;
            int dst = rand() % ttbl.n_st;
            int t = t_beg + rand() % (t_end - t_beg);
            int arr = csa.earliest_arrival_time(src, dst, t, false, true,
                                                 chg);
            if (arr <= t + max_delay) {
                ++n_ok;
                queries.push_back(std::make_tuple(src, dst, t));
            }
        }
        main_log.cerr(t) << n_q <<" Random queries, success rate : "
                         << (n_q*100/n_try) <<"%\n";
        
        t = main_log.lap();
    }
    // */


    //* Print a journey
    if (has_opt(argc, argv, "-journey")) {
        int src = std::stoi(get_opt(argc, argv, "-src=", "-1"));
        if (src == -1)
            src = ttbl.id_to_station[get_opt(argc, argv, "-src-id=", "")];
        int dst = std::stoi(get_opt(argc, argv, "-dst=", "-1"));
        if (dst == -1)
            dst = ttbl.id_to_station[get_opt(argc, argv, "-dst-id=", "")];
        int t = std::stoi(get_opt(argc, argv, "-t=", "0"));
        bool hubs = has_opt(argc, argv, "-hubs");

        int arr = rpt.earliest_arrival_time(src, dst, t, hubs, ! hubs, chg);
        std::cout <<" -------- "<< (hubs ? "HL_" : "") <<"Raptor "
                  <<"from "<< src << "=" << ttbl.station_id[src] <<" at "<< t;
        rpt.print_journey(dst, std::cout, chg);

        arr = rpt.earliest_arrival_time(src, dst, t, hubs, ! hubs,
                                        chg, 0, km, true);
        std::cout <<" -------- "<< (hubs ? "HL_" : "") <<"Raptor_trips "
                  <<"from "<< src << "=" << ttbl.station_id[src] <<" at "<< t;
        rpt.print_journey(dst, std::cout, chg);

        int dep = - rev_rpt.earliest_arrival_time(dst, src, - arr,
                                                hubs, ! hubs, 0, chg, km, true);
        std::cout <<" -------- "<< (hubs ? "HL_" : "") <<"RaptorREV_trips "
                  <<"from "<< dst << "=" << ttbl.station_id[dst] <<" at "<< t;
        rev_rpt.print_journey(src, std::cout, 0, chg);
        
        arr = csa.earliest_arrival_time(src, dst, t, hubs, ! hubs, chg);
        std::cout <<" -------- "<< (hubs ? "HL_" : "") <<"CSA "
                  <<"from "<< src << "=" << ttbl.station_id[src] <<" at "<< t;
        csa.print_journey(dst, hubs, ! hubs, std::cout, chg);
    }

    //* Arrival times
    if (has_opt(argc, argv, "-arrival-times")) {
        std::cout <<"src,dst,tdep,eat,eat_Unrestricted_Walking,eat_Walk_Only\n";
        n_ok = 0;
        for (auto q : queries) {
            int src = std::get<0>(q);
            int dst = std::get<1>(q);
            int t = std::get<2>(q);
            int arr1 = rpt.earliest_arrival_time(src, dst, t, false, true, chg);
            int arr2 = csa.earliest_arrival_time_opt(src, dst, t, false, true, chg);
            assert(arr1 == arr2); // can fail if chg == 0
            int arrHL1 = rpt.earliest_arrival_time(src, dst, t, hub, trf, chg);
            int arrHL2 = csa.earliest_arrival_time_opt(src, dst, t, hub, trf, chg);
            assert(arrHL1 == arrHL2); // can fail if chg == 0
            std::cout << ttbl.station_id[src] <<","<< ttbl.station_id[dst]
                      <<","<< t <<","<< arr1 <<","<< arrHL1
                      <<","<< (t + ttbl.walking_time(src, dst)) <<"\n";
            std::cout.flush();
            ++n_ok;
        }
        main_log.cerr(t) << n_ok << " arrival times\n";
        t = main_log.lap();
    }
    // */

    //* Missing transfers
    if (has_opt(argc, argv, "-missing-transfers")) {
        bool hubs = true;
        for (auto q : queries) {
            int src = std::get<0>(q);
            int t = std::get<2>(q);
            int nvis = rpt.eat_one_to_all(src, t, hubs, ! hubs, chg);
            for (int i = 0; i < ttbl.n_st; ++i) {
                if (rpt.eat_to(i) < ttbl.t_max) {
                    rpt.print_missing_transfers(i, std::cout);
                }
            }
        }
    }

    //* One to all arrival times
    if (has_opt(argc, argv, "-one-to-all")) {
        bool hubs = has_opt(argc, argv, "-hubs");
        std::srand(std::time(0));
        n_ok = 0;
        std::cout <<"src,tdep,nb_reached,rank,dst,travel_time,walk_time\n";
        for (auto q : queries) {
            int src = std::get<0>(q);
            int t = std::get<2>(q);
            int nvis = rpt.eat_one_to_all(src, t, hubs, ! hubs, chg);
            std::vector<std::tuple<int, int, int> > eat(ttbl.n_st);
            for (int i = 0; i < ttbl.n_st; ++i) {
                eat[i] = std::make_tuple(i, rpt.eat_to(i),
                                         ttbl.walking_time(src, i)) ;
            }
            std::sort(eat.begin(), eat.end(),
                      [](const std::tuple<int, int, int>&a,
                         const std::tuple<int, int, int>&b){
                          if (std::get<1>(a) != std::get<1>(b))
                              return std::get<1>(a) < std::get<1>(b);
                          return std::get<2>(a) < std::get<2>(b);
                      });
            for (int r = 0; r < ttbl.n_st; ++r) {
                int st = std::get<0>(eat[r]), arr = std::get<1>(eat[r]),
                    walk = std::get<2>(eat[r]);
                std::cout << ttbl.station_id[src]<<","<< t <<","<< nvis
                          <<","<< r <<","<< ttbl.station_id[st]
                          <<","<< (arr - t) <<","<< walk <<"\n";
            }
            std::cout.flush();
            ++n_ok;
        }
        main_log.cerr(t) << n_ok << " one to all arrival times\n";
        t = main_log.lap();
    }
    // */

    //* Generate queries by walking time rank
    if (has_opt(argc, argv, "-rank")) {
        bool hubs = true, trf = true;
        bool walk_rank = ! has_opt(argc, argv, "-travel-rank");
        bool all_ranks = has_opt(argc, argv, "-all-ranks");
        std::srand(std::time(nullptr));
        n_ok = 0;
        std::cout <<"source,destination,departure_time"
                  <<",log2_of_station_rank,station_rank,walk_time\n";
        std::vector<bool> seen(ttbl.n_st, false);
        for (auto q : queries) {
            if (n_ok >= std::stoi(get_opt(argc, argv, "-nq=", "1000"))) break;
            int src = std::get<0>(q);
            int dst = std::get<1>(q);
            int t = std::get<2>(q);
            if (seen[src]) continue;
            seen[src] = true;
            int nvis = rpt.eat_one_to_all(src, t, ! trf, trf, chg);
            std::vector<std::tuple<int, int, int> > eat(ttbl.n_st);
            for (int i = 0; i < ttbl.n_st; ++i) {
                eat[i] = std::make_tuple(i, ttbl.walking_time(src, i),
                                         rpt.eat_to(i));
            }
            std::sort(eat.begin(), eat.end(),
                      [src, walk_rank](const std::tuple<int, int, int>&a,
                                       const std::tuple<int, int, int>&b){
                          if (walk_rank) {
                              if (std::get<1>(a) != std::get<1>(b))
                                  return std::get<1>(a) < std::get<1>(b);
                              if (std::get<0>(a) == src) return true;
                              if (std::get<0>(b) == src) return false;
                              return std::get<0>(a) < std::get<0>(b);
                          } else {
                              if (std::get<2>(a) != std::get<2>(b))
                                  return std::get<2>(a) < std::get<2>(b);
                              if (std::get<1>(a) != std::get<1>(b))
                                  return std::get<1>(a) < std::get<1>(b);
                              if (std::get<0>(a) == src) return true;
                              if (std::get<0>(b) == src) return false;
                              return std::get<0>(a) < std::get<0>(b);
                          }
                      });
            int nvis_unrestr = rpt.eat_one_to_all(src, t, hubs, ! hubs, chg);
            int r_dst = -1;
            for (int r = 0; r < ttbl.n_st; ++r) {
                if (std::get<0>(eat[r]) == dst) r_dst = r;
            }
            int log_n_st = 0, rl = ttbl.n_st - 1;
            while (rl > 1) { rl /= 2; ++log_n_st; }
            int log_rnd = std::rand() % (log_n_st + 1);
            int rbase=2, log=1;
            // at least 2 min:
            while (std::get<1>(eat[rbase]) < 120) { rbase *= 2; ++log; }
            // reachable:
            int last_reach = (int)ttbl.n_st - 1;
            while (std::get<1>(eat[last_reach]) >= ttbl.t_max) { --last_reach; }
            // for all ranks:
            for ( ; rbase < ttbl.n_st; rbase *= 2, ++log){
                if (all_ranks || log == log_rnd) {
                    int interv_len = std::min(rbase, last_reach - rbase);
                    int r = rbase + (std::rand() % interv_len);
                    if (rbase <= r_dst && r_dst < 2*rbase) r = r_dst;
                    int st = std::get<0>(eat[r]);
                    int arr = std::get<2>(eat[r]), walk = std::get<1>(eat[r]);
                    int unrestr = rpt.eat_to(st);
                    std::cout
                        << ttbl.station_id[src] <<","<< ttbl.station_id[st]
                        <<","<< t <<","<< log <<","<< r <<","<< walk
                        //<<","<< (arr == ttbl.t_max ? ttbl.t_max : arr - t)
                     //<<","<<(unrestr == ttbl.t_max ? ttbl.t_max : unrestr - t)
                        <<"\n";
                }
            }
            std::cout.flush();
            ++n_ok;
        }
        main_log.cerr(t) << n_ok << " eat times rank\n";
        t = main_log.lap();
    }
    // */

    //* Generate queries with times and walking-time station-rank
    if (has_opt(argc, argv, "-rerank")) {
        bool hubs = false;
        std::srand(std::time(nullptr));
        n_ok = 0;
        std::cout <<"source,destination,departure_time"
                  <<",log2_of_station_rank,station_rank,walk_time\n";
        int src_prev = -1, t = 0, nvis = 0;
        std::vector<std::tuple<int, int> > eat(ttbl.n_st);
        for (auto q : queries) {
            int src = std::get<0>(q);
            int dst = std::get<1>(q);
            int t_req = std::get<2>(q);
            if (src != src_prev) {
                bool preserve_time = ! has_opt(argc, argv, "-gen-times");
                t = preserve_time ? t_req : std::rand() % (24*3600);
                nvis = rpt.eat_one_to_all(src, t, hubs, ! hubs, chg);
                for (int i = 0; i < ttbl.n_st; ++i) {
                    eat[i] = std::make_tuple(i, ttbl.walking_time(src, i));
                }
                std::sort(eat.begin(), eat.end(),
                          [src](const std::tuple<int, int>&a,
                             const std::tuple<int, int>&b){
                              if (std::get<1>(a) != std::get<1>(b))
                                  return std::get<1>(a) < std::get<1>(b);
                              if (std::get<0>(a) == src) return true;
                              if (std::get<0>(b) == src) return false;
                              return std::get<0>(a) < std::get<0>(b);
                          });
                src_prev = src;
            }
            assert(std::get<0>(eat[0]) == src);
            assert(rpt.eat_to(src) == t);
            int r = -1, walk = 0;
            for (int i = 0; i < ttbl.n_st; ++i) {
                if (std::get<0>(eat[i]) == dst) {
                    r = i;
                    walk = std::get<1>(eat[i]);
                    break;
                }
            }
            assert(r > 0);
            int log = 0, rl = r;
            while (rl > 1) { rl /= 2; ++log; }
            int arr = rpt.eat_to(dst);
            std::cout << ttbl.station_id[src] <<","<< ttbl.station_id[dst]
                      <<","<< t <<","<< log <<","<< r <<","<< walk
                      //<<","<< (arr == ttbl.t_max ? ttbl.t_max : arr - t)
                      <<"\n";
            std::cout.flush();
            ++n_ok;
        }
        main_log.cerr(t) << n_ok << " rerank\n";
        t = main_log.lap();
    }
    // */
    



    if (has_opt(argc, argv, "-exit")) exit(0);


    auto cout_avg_time = [n_q,&t,&main_log](std::string s) {
        std::cout << s <<"_avg_time "
                  << (main_log.lap() - t) * 1000.0 / n_q <<"\n";
    };
    t = main_log.lap();
    
    
    // go Raptor restricted walk
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int arr = rpt.earliest_arrival_time(src, dst, t, false, true, chg);
        //std::cout << arr <<"\n";
        //assert(arr < ttbl.t_max);
        if (arr < ttbl.t_max) {
            sum += arr - t;
            ++n_ok;
        }
    }
    cout_avg_time("Raptor");
    main_log.cerr(t) << n_q << " Raptor queries done, avg_time = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();


    // go HLRaptor
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int arr = rpt.earliest_arrival_time(src, dst, t, hub, trf, chg);
        //if (arr < ttbl.t_max) rpt.print_journey(dst);
        //std::cout << src <<","<< dst <<","<< t <<" : "<< (arr - dep) <<"\n";
        //std::cout << arr <<"\n";
        //assert(arr < ttbl.t_max);
        if (arr < ttbl.t_max) {
            sum += arr - t;
            ++n_ok;
        }
    }
    cout_avg_time("HLRaptor");
    main_log.cerr(t) << n_q << " HLRaptor queries done, avg_time = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();
    
    // go CSA
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int arr = csa.earliest_arrival_time_opt(src, dst, t, false, true, chg);
        //std::cout << arr <<"\n";
        //assert(arr < ttbl.t_max);
        if (arr < ttbl.t_max) {
            sum += arr - t;
            ++n_ok;
        }
    }
    cout_avg_time("CSA");
    main_log.cerr(t) << n_q << " CSA queries done, avg_time = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();


    // go HLCSA
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int arr = csa.earliest_arrival_time_opt(src, dst, t, hub, trf, chg);
        //assert(arr < ttbl.t_max);
        if (arr < ttbl.t_max) {
            sum += arr - t;
            ++n_ok;
        }
    }
    cout_avg_time("HLCSA");
    main_log.cerr(t) << n_q << " HLCSA queries done, avg_time = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();

    
    if (has_opt(argc, argv, "-exit-after-csa")) exit(0);

    
    //* go Pareto
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int npath = rpt.earliest_walk_pareto(src, dst, t, false, true, chg);
        if (npath > 0) {
            sum += npath;
            ++n_ok;
        }
    }
    cout_avg_time("mcRaptor");
    main_log.cerr(t) << n_q << " Pareto Raptor queries done, avg_npaths = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();
    // */
    
    //* go HLPareto
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        //int arr = rpt.earliest_arrival_time(src, dst, t, hub, trf, chg);
        int npath = rpt.earliest_walk_pareto(src, dst, t, hub, trf, chg);
        if (npath > 0) {
            sum += npath;
            ++n_ok;
        }
    }
    cout_avg_time("HLmcRaptor");
    main_log.cerr(t) << n_q << " Pareto HLRaptor queries done, avg_npaths = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();
    // */



    // go last departure Raptor
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int arr = rpt.earliest_arrival_time(src, dst, t, false, true, 0);
        int dep = - rev_rpt.earliest_arrival_time(dst, src, - arr, false, true, 0);
        int arr2 = rpt.earliest_arrival_time(src, dst, dep, false, true, 0);
        //std::cout << t <<" "<< dep <<" "<< arr <<" "<< arr2 <<"\n";
        //assert(arr < ttbl.t_max);
        if (arr < ttbl.t_max) {
            assert(arr == arr2);
            sum += dep - t;
            ++n_ok;
        }
    }
    //cout_avg_time("lastdepRaptor");
    main_log.cerr(t) << n_q << " last dep Raptor queries done, avg_dep_time = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();
    // */


    // int range_2h = has_opt(argc, argv, "-2h-range");
    for (int range_2h = 1; range_2h != -1; --range_2h) {
    
    // go profile Raptor
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int t_beg = range_2h ? t : 0;
        int t_end = range_2h ? t + 7200 : 24 * 3600;
        std::cout << src <<" "<< dst
                  <<" ["<< t_beg <<","<< t_end <<") : ";
        pset prof = rpt.profile(rev_rpt, src, dst, t_beg, t_end,
                                 false, true, chg);
        int ntrips = prof.size();
        std::cout << ntrips <<"\n";
        sum += ntrips;
        ++n_ok;
    }
    cout_avg_time("prRaptor");
    main_log.cerr(t) << n_q << " Profile Raptor queries done, avg_ntrips = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();

    
    // go profile HLRaptor
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int t_beg = range_2h ? t : 0;
        int t_end = range_2h ? t + 7200 : 24 * 3600;
        std::cout << src <<" "<< dst
                  <<" ["<< t_beg <<","<< t_end <<") : ";
        pset prof = rpt.profile(rev_rpt, src, dst, t_beg, t_end,
                                 hub, trf, chg);
        int ntrips = prof.size();
        std::cout << ntrips <<"\n"; std::cout.flush();
        if (src == std::stoi(get_opt(argc, argv, "-src=", "-1"))
            && dst == std::stoi(get_opt(argc, argv, "-dst=", "-1")))
            prof.print();
        sum += ntrips;
        ++n_ok;
    }
    cout_avg_time("HLprRaptor");
    main_log.cerr(t) << n_q << " Profile HLRaptor queries done, avg_ntrips = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();
    // */    


    // go profile CSA
    sum = 0, n_ok = 0;
    bool prescan = has_opt(argc, argv, "-csa-profile-prescan");
    int t_get = std::stoi(get_opt(argc, argv, "-t-beg=", "0"));
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int t_beg = range_2h ? t : 0;
        int t_end = range_2h ? t + 7200 : 24 * 3600;
        std::cout << src <<"(="<< ttbl.station_id[src] <<") "
                  << dst <<"(="<< ttbl.station_id[dst] <<") "
                  <<" ["<< t_beg <<","<< t_end<<") : ";
        pset prof = csa.profile(src, dst, t_beg, t_end, false, true, chg,
                                0, km, prescan);
        int ntrips = prof.size();
        std::cout << ntrips <<"\n";
        sum += ntrips;
        ++n_ok;
    }
    cout_avg_time("prCSA");
    main_log.cerr(t) << n_q << " Profile CSA queries done, avg_ntrips = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();
    // */    


    // go profile HLCSA
    sum = 0, n_ok = 0;
    for (auto q : queries) {
        int src = std::get<0>(q);
        int dst = std::get<1>(q);
        int t = std::get<2>(q);
        int t_beg = range_2h ? t : 0;
        int t_end = range_2h ? t + 7200 : 24 * 3600;
        std::cout << src <<" "<< dst
                  <<" ["<< t_beg <<","<< t_end <<") : ";
        pset prof = csa.profile(src, dst, t_beg, t_end, hub, trf, chg,
                                0, km, prescan);
        int ntrips = prof.size();
        std::cout << ntrips <<"\n";
        if (src == std::stoi(get_opt(argc, argv, "-src=", "-1"))
            && dst == std::stoi(get_opt(argc, argv, "-dst=", "-1")))
            prof.print();
        sum += ntrips;
        ++n_ok;
    }
    cout_avg_time("HLprCSA");
    main_log.cerr(t) << n_q << " Profile HLCSA queries done, avg_ntrips = "
                     << (sum / n_ok)
                     << "  "<< n_ok <<"/"<< queries.size() <<" ok\n";
    t = main_log.lap();
    // */

    }

    // ------------------------ end -------------------------
    main_log.cerr() << "end "<< dir <<"\n";

}
