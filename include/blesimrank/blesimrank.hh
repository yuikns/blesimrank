#ifndef BLESIMRANK_BLESIMRANK_HH
#define BLESIMRANK_BLESIMRANK_HH

#include <cstring>
#include <cstdlib>
#include <cmath>

#include "graph.hh"
#include <string>
#include <vector>
#include <set>
#include <map>
#include <sys/stat.h>

#include "argcv/random/random.hh"
#include "argcv/type/heapq.hh"

using namespace argcv::random;

using argcv::type::heapq;

class Blesimrank {
public:
    Blesimrank(const std::string &dataset, int _t = 11, int _p = 11, int _q = 5, int _r1 = 100,
               int _r2 = 10000, int _k = 20, double _theta = 0.01)
        : dataset(dataset),
          network_filename("data/" + dataset + ".graph"),
          output_filename("result/" + dataset + ".ssr"),
          _t(_t),
          _p(_p),
          _q(_q),
          _r1(_r1),
          _r2(_r2),
          _k(_k),
          _theta(_theta),
          _c(0.6) {
        printf(
            "[Blesimrank] dataset: %s, T:%d, P:%d, Q:%d, R1:%d, R2:%d, K:%d, "
            "theta:%f\n",
            dataset.c_str(), T(), P(), Q(), R1(), R2(), K(), theta());
        init();
    }

    ~Blesimrank() {}

    void init() {
        FILE *fp = fopen(network_filename.c_str(), "r");
        char line[100];
        fgets(line, 1000, fp);
        read_ds_size(line, _n, _m);
        printf("[Blesimrank] init ... N: %d, M: %d, R1=%d, R2=%d\n", N(), M(), R1(), R2());
        fflush(NULL);
        _g1 = new BGrapgh(N(), R1());
        printf("[Blesimrank] G1 inited\n");
        fflush(NULL);

        size_t count = 0;
        while (fgets(line, 100, fp)) {
            int a, b;
            double c;
            sscanf(line, "%d\t%d\t%lf", &a, &b, &c);
            // printf("[Blesimrank] get node : a:%d b:%d c:%f\n", a, b, c);
            _g1->get_node(a)->weight += c;
            _g1->get_node(b)->weight += c;
            _g1->get_node(a)->edges.push_back(GEdge(b, c, c));
            _g1->get_node(b)->edges.push_back(GEdge(a, c, c));
            //_g1->get_node(a)->edges.insert(std::make_pair(b, GEdge(c, c)));
            //_g1->get_node(b)->edges.insert(std::make_pair(a, GEdge(c, c)));

            count++;
            if (count % 1000 == 0) {
                printf("[Blesimrank] load data %zu\n", count);
            }
        }
        fclose(fp);
        fprintf(stdout, "[Blesimrank] load data finished, node size : %zu , edge size: %zu \n", _g1->nsize(),
                _g1->esize() / 2);
    }

    void read_ds_size(char *line, int &n, int &m) {
        char *k = strtok(line, "\t");
        if (k == NULL) {
            fprintf(stderr, "[read_ds_size] error unknown val : %s\n", line);
            return;
        }
        n = atoi(k);
        k = strtok(NULL, "\t");
        if (k == NULL) {
            fprintf(stderr, "[read_ds_size] erroe unknown val : %s\n", line);
            return;
        }
        m = atoi(k);
    }

    void trans_pros() {
        for (int i = 0; i < _g1->nsize(); i++) {
            GNode *_node = _g1->get_node(i);
            double sum_weight = _node->weight;
            if (_node->edges.size() > 0) {
                _node->edges[0].accum_weight = _node->edges[0].accum_weight / sum_weight;
                for (int j = 1; j < _node->edges.size(); j++) {
                    _node->edges[i].accum_weight = _node->edges[i - 1].accum_weight
                                                   + _node->edges[i].accum_weight / sum_weight;
                }
            }
        }
    }

    std::vector<size_t> randomWalk(BGrapgh *_g, int v) {
        std::vector<size_t> path;
        int c_node = v;
        int ct = 0;
        GNode *_node = _g->get_node(c_node);
        for (ct = 0; ct < _t; ct++) {
            std::vector<GEdge> &e = _node->edges;
            if (_node->edges.size() > 0) {
                double r = random_double();
                int j;
                for (j = 0; j < _node->edges.size(); j++) {
                    if (e[j].accum_weight > r) break;
                }
                c_node = e[j].id;
                path.push_back(c_node);
                _node = _g->get_node(c_node);
            } else {
                break;
            }
        }
        return path;
    }

    void gen_bipartitle_graph() {
        _g2 = new BGrapgh(_g1);  // fork a new node&edges instance
        for (int v = 0; v < _g1->nsize(); v++) {
            for (int i = 0; i < _p; i++) {
                std::vector<std::vector<size_t>> paths;
                for (int q = 0; q < _q; q++) {
                    paths.push_back(randomWalk(_g1, v));
                }
                for (int t = 0; t < _t; t++) {
                    std::map<size_t, bool> hitter;
                    for (int q = 0; q < _q; q++) {
                        if (paths[q].size() > t) {
                            hitter[paths[q][t]] = hitter.find(paths[q][t]) == hitter.end() ? false : true;
                        }
                    }
                    for (std::map<size_t, bool>::const_iterator it = hitter.begin(); it != hitter.end();
                         it++) {
                        if (it->second) {
                            _g2->get_node(v)->edges.push_back(GEdge(it->first, 0, 0));
                            _g2->get_node(it->first)->edges.push_back(GEdge(v, 0, 0));  // is it necessary ?
                        }
                    }
                }
            }
            if (v % 1000 == 0) {
                printf("[gen_bipartitle_graph]::ing node %d\n", v);
            }
        }
        printf("[gen_bipartitle_graph]::done\n");
    }

    void cal_gamma_bound() {
        for (int v = 0; v < _g1->nsize(); v++) {
            GNode *_node = _g1->get_node(v);
            std::vector<GEdge> &e = _node->edges;
            _node->gamma = 0;
            if (e.size() > 0) {
                std::vector<std::vector<size_t>> paths;
                for (int r = 0; r < _r1; r++) {
                    paths.push_back(randomWalk(_g1, v));
                }
                for (int t = 0; t < _t; t++) {
                    std::map<size_t, size_t> hitter;
                    for (int r = 0; r < _r1; r++) {
                        if (paths[r].size() > t)
                            hitter[paths[r][t]]
                                = hitter.find(paths[r][t]) == hitter.end() ? 1 : hitter[paths[r][t]] + 1;
                    }
                    double mu = 0;
                    for (std::map<size_t, size_t>::const_iterator it = hitter.begin(); it != hitter.end();
                         it++) {
                        size_t path_number = it->second;
                        mu += (1 - _c) * path_number * path_number / (_r1 * _r1);
                    }

                    double cgamma = sqrt(mu);
                    if (cgamma > _node->gamma) {
                        _node->gamma = cgamma;
                    }
                }
            }
        }
    }

    void cal_beta_bound() { fprintf(stderr, "NOT USED\n"); }

    void preproc() {
        trans_pros();
        gen_bipartitle_graph();
        cal_gamma_bound();
    }

    double simrank(int v, int u, int _r) {
        std::vector<std::vector<size_t>> paths_v;
        std::vector<std::vector<size_t>> paths_u;
        for (int r = 0; r < _r; r++) {
            paths_v.push_back(randomWalk(_g1, v));
            paths_u.push_back(randomWalk(_g1, u));
        }
        double eta = 0;
        for (int t = 0; t < _t; t++) {
            std::map<size_t, size_t> hitter_v;
            std::map<size_t, size_t> hitter_u;
            for (int r = 0; r < _r1; r++) {
                if (paths_v[r].size() > t)
                    hitter_v[paths_v[r][t]]
                        = hitter_v.find(paths_v[r][t]) == hitter_v.end() ? 1 : hitter_v[paths_v[r][t]] + 1;
                if (paths_u[r].size() > t)
                    hitter_u[paths_u[r][t]]
                        = hitter_u.find(paths_u[r][t]) == hitter_u.end() ? 1 : hitter_u[paths_u[r][t]] + 1;
            }
            for (std::map<size_t, size_t>::const_iterator it = hitter_v.begin(); it != hitter_v.end(); it++) {
                size_t alpha = it->second;
                std::map<size_t, size_t>::const_iterator it_u = hitter_v.find(it->first);
                if (it_u != hitter_v.end()) {
                    size_t beta = it_u->second;
                    eta += pow(_c, t) * (1 - _c) * alpha * beta / (_r * _r);
                }
            }
        }
        return eta;
    }

    static int pair_compare_by_value(std::pair<size_t, double> a, std::pair<size_t, double> b) {
        return a.second > b.second ? 1 :
                                   //(a.second == b.second ? (int_compare(a.first,b.first)) : -1);
                   (a.second == b.second ? 0 : -1);
    }

    template <typename T>
    void vector_reverse(std::vector<T> &v) {
        size_t len = v.size();
        for (size_t i = 0; i < len / 2; i++) {
            T t = v[i];
            v[i] = v[len - 1 - i];
            v[len - 1 - i] = t;
        }
    }

    std::vector<std::pair<size_t, double>> query(int v) {  // cands and scores
        GNode *_node = _g2->get_node(v);
        std::vector<GEdge> &e = _node->edges;
        std::vector<std::pair<size_t, double>> cands;
        if (e.size() == 0) {
            return cands;
        } else {
            std::set<size_t> cset;
            for (int i = 0; i < e.size(); i++) {
                std::vector<GEdge> &ce = _g2->get_node(e[i].id)->edges;
                for (int j = 0; j < ce.size(); j++) {
                    if (ce[j].id != v) {
                        cset.insert(ce[j].id);
                    }
                }
            }

            heapq<std::pair<size_t, double>> hq(_k, pair_compare_by_value);
            for (std::set<size_t>::const_iterator it = cset.begin(); it != cset.end(); it++) {
                if (_g1->get_node(*it)->gamma >= _theta) {
                    int id = *it;
                    double score = simrank(v, id, 100);
                    hq.push(std::make_pair(id, score));
                }
            }
            std::pair<size_t, double> tval;
            while (hq.pop(tval)) {
                cands.push_back(tval);
            }
            vector_reverse(cands);
        }
        return cands;
    }

    void save() {
        mkdir("result", 0755);
        FILE *fp = fopen(output_filename.c_str(), "w");
        for (int i = 0; i < _g1->nsize(); i++) {
            std::vector<std::pair<size_t, double>> r = query(i);
            for (std::vector<std::pair<size_t, double>>::const_iterator it = r.begin(); it != r.end(); it++) {
                fprintf(fp,"%zu:%f\t",it->first,it->second);
            }
            fprintf(fp,"\n");
            if(i % 100 == 0 ) {
                printf("saving : %d \n",i);
            }
        }
        fclose(fp);
        printf("[save] all done\n");
    }

    const int T() { return _t; }
    const int P() { return _p; }
    const int Q() { return _q; }
    const int R1() { return _r1; }
    const int R2() { return _r2; }
    const int K() { return _k; }
    const double theta() { return _theta; }
    const int N() { return _n; }
    const int M() { return _m; }

private:
    const std::string dataset;
    const std::string network_filename;
    const std::string output_filename;
    int _t;
    int _p;
    int _q;
    int _r1;
    int _r2;
    int _k;
    double _theta;
    double _c;
    int _n;
    int _m;
    BGrapgh *_g1;
    BGrapgh *_g2;
};

#endif  //  BLESIMRANK_BLESIMRANK_HH
