#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <array>
#include <algorithm>

using namespace std;

static const double NEG_INF = -1e300;

static inline double logp(double p) { return (p <= 0.0) ? NEG_INF : log(p); }

// stable log-sum-exp for 3 values
static inline double logsumexp3(double a, double b, double c) {
    double m = max(a, max(b, c));
    if (m <= NEG_INF/2) return NEG_INF;
    return m + log(exp(a - m) + exp(b - m) + exp(c - m));
}

static inline double logsumexp_vec(const vector<double>& v) {
    double m = NEG_INF;
    for (double x : v) m = max(m, x);
    if (m <= NEG_INF/2) return NEG_INF;
    double s = 0.0;
    for (double x : v) s += exp(x - m);
    return m + log(s);
}

static string readCString(ifstream &input) {
    string result;
    char ch;
    while (input.get(ch)) {
        if (ch == '\0') break;
        result.push_back(ch);
    }
    return result;
}

static inline void trim_cr(string &s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

static inline int obs_state(char c1, char c2) {
    bool g1 = (c1=='-');
    bool g2 = (c2=='-');
    if (!g1 && !g2) return 1; // Match
    if ( g1 && !g2) return 2; // Insertion
    if (!g1 &&  g2) return 3; // Deletion
    return -1; // "--"
}

static inline double emit_log(int state, const string& o,
                              const vector<unordered_map<string,double>>& e,
                              double floor_prob) {
    auto it = e[state].find(o);
    double p = (it == e[state].end()) ? floor_prob : it->second;
    return logp(p);
}

static void load_initial_values(const string &path,
                                vector<vector<double>> &A,
                                vector<unordered_map<string,double>> &E) {
    const int S = 5;
    A.assign(S, vector<double>(S, 0.0));
    E.assign(S, unordered_map<string,double>());

    ifstream file(path);
    if (!file.is_open()) {
        cerr << "Error opening file: " << path << "\n";
        exit(1);
    }

    for (int i=0;i<S;++i) {
        for (int j=0;j<S;++j) {
            string tok;
            if (!(file >> tok)) {
                cerr << "Unexpected EOF while reading A\n";
                exit(1);
            }
            trim_cr(tok);
            A[i][j] = stod(tok);
        }
    }

    string key;
    while (file >> key) {
        trim_cr(key);
        string colon;
        double p;
        if (!(file >> colon >> p)) break;
        trim_cr(colon);
        trim_cr(key);

        if (key.size()>=2 && key[0]=='-') E[2][key] = p;
        else if (key.size()>=2 && key[1]=='-') E[3][key] = p;
        else E[1][key] = p;
    }
}

// save A/E in same format (A 5x5 then emissions)
static void save_values(const string& path,
                        const vector<vector<double>>& A,
                        const vector<unordered_map<string,double>>& E) {
    ofstream out(path);
    if (!out.is_open()) {
        cerr << "Error opening output: " << path << "\n";
        exit(1);
    }
    out.setf(std::ios::fixed); out.precision(10);

    const int S = 5;
    for (int i=0;i<S;++i) {
        for (int j=0;j<S;++j) {
            out << A[i][j];
            if (!(i==S-1 && j==S-1)) out << "\t";
        }
        out << "\n";
    }

    // Emissions: dump in stable order for reproducibility
    auto dump_map = [&](const unordered_map<string,double>& mp) {
        vector<pair<string,double>> v(mp.begin(), mp.end());
        sort(v.begin(), v.end(), [](auto& a, auto& b){ return a.first < b.first; });
        for (auto& kv : v) out << kv.first << " : " << kv.second << "\n";
    };

    dump_map(E[1]);
    dump_map(E[2]);
    dump_map(E[3]);
    out.close();
}

// Build observations and colType from aligned strings.
static void build_obs(const string& s1, const string& s2,
                      vector<string>& O, vector<int>& colType) {
    size_t L0 = min(s1.size(), s2.size());
    O.clear(); colType.clear();
    O.reserve(L0); colType.reserve(L0);

    for (size_t j=0;j<L0;++j) {
        char c1=s1[j], c2=s2[j];
        if (c1=='-' && c2=='-') continue;
        int t = obs_state(c1,c2);
        if (t<0) continue;
        string obs; obs.push_back(c1); obs.push_back(c2);
        O.push_back(obs);
        colType.push_back(t);
    }
}

// Forward-Backward for one sequence of observations.
// Returns log-likelihood and accumulates expected counts into global accumulators.
static double forward_backward_one(const vector<string>& O,
                                  const vector<int>& colType,
                                  const vector<vector<double>>& A,
                                  const vector<unordered_map<string,double>>& E,
                                  // accumulators:
                                  array<double,3>& init_acc,                 // Begin->(M/I/D)
                                  array<double,3>& end_acc,                  // (M/I/D)->End
                                  array<array<double,3>,3>& trans_acc,       // (M/I/D)->(M/I/D)
                                  unordered_map<string,double>& emitM_acc,
                                  unordered_map<string,double>& emitI_acc,
                                  unordered_map<string,double>& emitD_acc,
                                  double emit_floor=1e-12,
                                  double trans_floor=1e-15)
{
    const int BEGIN=0, M=1, I=2, D=3, END=4;
    const int K=3; // M,I,D as 0,1,2 indices in DP arrays

    int T = (int)O.size();
    if (T==0) return 0.0;

    auto trans_log = [&](int from, int to){
        // Optionally enforce structure:
        if (from == END) return NEG_INF;    // nothing after End
        if (to == BEGIN) return NEG_INF;    // never go to Begin
        double p = A[from][to];
        if (p<=0.0) p = trans_floor;
        return logp(p);
    };

    auto emit_masked = [&](int state, int t)->double{
        // state is 1/2/3 (M/I/D)
        if (colType[t] != state) return NEG_INF;
        return emit_log(state, O[t], E, emit_floor);
    };

    // alpha[t][k], beta[t][k] for k in {M,I,D} mapped to {0,1,2}
    vector<array<double,3>> alpha(T), beta(T);

    // --- Forward init ---
    alpha[0][0] = trans_log(BEGIN,M) + emit_masked(M,0);
    alpha[0][1] = trans_log(BEGIN,I) + emit_masked(I,0);
    alpha[0][2] = trans_log(BEGIN,D) + emit_masked(D,0);

    // --- Forward recursion ---
    for (int t=1;t<T;++t) {
        // to M
        alpha[t][0] = logsumexp3(
            alpha[t-1][0] + trans_log(M,M),
            alpha[t-1][1] + trans_log(I,M),
            alpha[t-1][2] + trans_log(D,M)
        ) + emit_masked(M,t);

        // to I
        alpha[t][1] = logsumexp3(
            alpha[t-1][0] + trans_log(M,I),
            alpha[t-1][1] + trans_log(I,I),
            alpha[t-1][2] + trans_log(D,I)
        ) + emit_masked(I,t);

        // to D
        alpha[t][2] = logsumexp3(
            alpha[t-1][0] + trans_log(M,D),
            alpha[t-1][1] + trans_log(I,D),
            alpha[t-1][2] + trans_log(D,D)
        ) + emit_masked(D,t);
    }

    // --- Termination ---
    double loglik = logsumexp3(
        alpha[T-1][0] + trans_log(M,END),
        alpha[T-1][1] + trans_log(I,END),
        alpha[T-1][2] + trans_log(D,END)
    );

    // --- Backward init ---
    beta[T-1][0] = trans_log(M,END);
    beta[T-1][1] = trans_log(I,END);
    beta[T-1][2] = trans_log(D,END);

    // --- Backward recursion ---
    for (int t=T-2;t>=0;--t) {
        // from M at time t -> next state at t+1
        beta[t][0] = logsumexp3(
            trans_log(M,M) + emit_masked(M,t+1) + beta[t+1][0],
            trans_log(M,I) + emit_masked(I,t+1) + beta[t+1][1],
            trans_log(M,D) + emit_masked(D,t+1) + beta[t+1][2]
        );
        // from I
        beta[t][1] = logsumexp3(
            trans_log(I,M) + emit_masked(M,t+1) + beta[t+1][0],
            trans_log(I,I) + emit_masked(I,t+1) + beta[t+1][1],
            trans_log(I,D) + emit_masked(D,t+1) + beta[t+1][2]
        );
        // from D
        beta[t][2] = logsumexp3(
            trans_log(D,M) + emit_masked(M,t+1) + beta[t+1][0],
            trans_log(D,I) + emit_masked(I,t+1) + beta[t+1][1],
            trans_log(D,D) + emit_masked(D,t+1) + beta[t+1][2]
        );
    }

    // --- Accumulate expected counts ---
    // gamma[t][k] = P(state=k at t | O)
    // in log: loggamma = alpha+beta-loglik
    // Begin->state at t=0:
    {
        double lgM = alpha[0][0] + beta[0][0] - loglik;
        double lgI = alpha[0][1] + beta[0][1] - loglik;
        double lgD = alpha[0][2] + beta[0][2] - loglik;
        init_acc[0] += exp(lgM);
        init_acc[1] += exp(lgI);
        init_acc[2] += exp(lgD);
    }

    // emissions expected counts
    for (int t=0;t<T;++t) {
        double lgM = alpha[t][0] + beta[t][0] - loglik;
        double lgI = alpha[t][1] + beta[t][1] - loglik;
        double lgD = alpha[t][2] + beta[t][2] - loglik;

        // Only one of these will be non-NEG_INF because of masking, but we keep generic.
        double gM = (lgM<=NEG_INF/2) ? 0.0 : exp(lgM);
        double gI = (lgI<=NEG_INF/2) ? 0.0 : exp(lgI);
        double gD = (lgD<=NEG_INF/2) ? 0.0 : exp(lgD);

        if (gM>0) emitM_acc[O[t]] += gM;
        if (gI>0) emitI_acc[O[t]] += gI;
        if (gD>0) emitD_acc[O[t]] += gD;
    }

    // transitions expected counts: xi[t][k->u] for t=0..T-2
    for (int t=0;t<T-1;++t) {
        // helper to compute log xi for a particular (from,to)
        auto lxi = [&](int fromState, int toState, int fromIdx, int toIdx)->double{
            // fromState/toState are 1/2/3 ; fromIdx/toIdx are 0/1/2
            return alpha[t][fromIdx]
                 + trans_log(fromState,toState)
                 + emit_masked(toState,t+1)
                 + beta[t+1][toIdx]
                 - loglik;
        };

        // accumulate all 9 transitions among M/I/D
        double vMM = lxi(M,M,0,0);
        double vMI = lxi(M,I,0,1);
        double vMD = lxi(M,D,0,2);

        double vIM = lxi(I,M,1,0);
        double vII = lxi(I,I,1,1);
        double vID = lxi(I,D,1,2);

        double vDM = lxi(D,M,2,0);
        double vDI = lxi(D,I,2,1);
        double vDD = lxi(D,D,2,2);

        // add exp(log xi)
        trans_acc[0][0] += (vMM<=NEG_INF/2)?0.0:exp(vMM);
        trans_acc[0][1] += (vMI<=NEG_INF/2)?0.0:exp(vMI);
        trans_acc[0][2] += (vMD<=NEG_INF/2)?0.0:exp(vMD);

        trans_acc[1][0] += (vIM<=NEG_INF/2)?0.0:exp(vIM);
        trans_acc[1][1] += (vII<=NEG_INF/2)?0.0:exp(vII);
        trans_acc[1][2] += (vID<=NEG_INF/2)?0.0:exp(vID);

        trans_acc[2][0] += (vDM<=NEG_INF/2)?0.0:exp(vDM);
        trans_acc[2][1] += (vDI<=NEG_INF/2)?0.0:exp(vDI);
        trans_acc[2][2] += (vDD<=NEG_INF/2)?0.0:exp(vDD);
    }

    // end transitions (state at T-1 -> End)
    {
        double lgM = alpha[T-1][0] + trans_log(M,END) - loglik;
        double lgI = alpha[T-1][1] + trans_log(I,END) - loglik;
        double lgD = alpha[T-1][2] + trans_log(D,END) - loglik;
        end_acc[0] += exp(lgM);
        end_acc[1] += exp(lgI);
        end_acc[2] += exp(lgD);
    }

    return loglik;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: baum_welch <input.bin> <initial_values.txt> <iters> [max_pairs]\n";
        return 1;
    }
    string binfile = argv[1];
    string initfile = argv[2];
    int iters = stoi(argv[3]);
    long long max_pairs = (argc>=5) ? stoll(argv[4]) : -1; // -1 means all

    vector<vector<double>> A;
    vector<unordered_map<string,double>> E;
    load_initial_values(initfile, A, E);

    for (int iter=1; iter<=iters; ++iter) {
        // accumulators reset each iter
        array<double,3> init_acc{0,0,0};
        array<double,3> end_acc{0,0,0};
        array<array<double,3>,3> trans_acc{};
        for (int i=0;i<3;++i) for (int j=0;j<3;++j) trans_acc[i][j]=0.0;

        unordered_map<string,double> emitM_acc, emitI_acc, emitD_acc;

        ifstream in(binfile, ios::binary);
        if (!in.is_open()) {
            cerr << "Error opening " << binfile << "\n";
            return 1;
        }

        vector<string> O;
        vector<int> colType;

        long long pairs=0;
        double total_loglik = 0.0;

        while (true) {
            int cleanedCount;
            if (!in.read(reinterpret_cast<char*>(&cleanedCount), sizeof(int))) break;
            if (cleanedCount <= 0) continue;

            for (int i=0;i<cleanedCount;++i) {
                string n1 = readCString(in);
                if (!in) { cerr << "EOF reading n1\n"; return 1; }
                string n2 = readCString(in);
                if (!in) { cerr << "EOF reading n2\n"; return 1; }
                string s1 = readCString(in);
                if (!in) { cerr << "EOF reading s1\n"; return 1; }
                string s2 = readCString(in);
                if (!in) { cerr << "EOF reading s2\n"; return 1; }

                build_obs(s1,s2,O,colType);
                double ll = forward_backward_one(O,colType,A,E,
                                                init_acc,end_acc,trans_acc,
                                                emitM_acc,emitI_acc,emitD_acc);
                total_loglik += ll;

                pairs++;
                if (max_pairs > 0 && pairs >= max_pairs) break;
            }
            if (max_pairs > 0 && pairs >= max_pairs) break;
        }
        in.close();

        // --- M-step: update A and E from expected counts ---
        const int BEGIN=0, M=1, I=2, D=3, END=4;

        vector<vector<double>> Anew = A;
        vector<unordered_map<string,double>> Enew = E;

        // Update Begin row: Begin -> {M,I,D}
        {
            double s = init_acc[0] + init_acc[1] + init_acc[2];
            if (s > 0) {
                Anew[BEGIN][M] = init_acc[0]/s;
                Anew[BEGIN][I] = init_acc[1]/s;
                Anew[BEGIN][D] = init_acc[2]/s;
            }
            // forbid Begin->Begin and Begin->End explicitly
            Anew[BEGIN][BEGIN] = 0.0;
            Anew[BEGIN][END]   = 0.0;
        }

        // Update transitions among M/I/D and to End
        auto norm_row = [&](int fromState, int fromIdx) {
            double s = trans_acc[fromIdx][0] + trans_acc[fromIdx][1] + trans_acc[fromIdx][2] + end_acc[fromIdx];
            if (s <= 0) return;
            Anew[fromState][M]   = trans_acc[fromIdx][0]/s;
            Anew[fromState][I]   = trans_acc[fromIdx][1]/s;
            Anew[fromState][D]   = trans_acc[fromIdx][2]/s;
            Anew[fromState][END] = end_acc[fromIdx]/s;
            Anew[fromState][BEGIN] = 0.0; // forbid to Begin
        };

        norm_row(M,0);
        norm_row(I,1);
        norm_row(D,2);

        // End row stays zero
        for (int j=0;j<5;++j) Anew[END][j]=0.0;

        // Update emissions: normalize expected emission counts for each state
        auto norm_emit = [&](int state, unordered_map<string,double>& acc) {
            double s = 0.0;
            for (auto& kv : acc) s += kv.second;
            if (s <= 0) return;
            Enew[state].clear();
            for (auto& kv : acc) Enew[state][kv.first] = kv.second / s;
        };

        norm_emit(M, emitM_acc);
        norm_emit(I, emitI_acc);
        norm_emit(D, emitD_acc);

        A.swap(Anew);
        E.swap(Enew);

        cerr << "Iter " << iter << ": pairs=" << pairs
             << " total_loglik=" << total_loglik << "\n";

        // write intermediate result each iter
        string outname = "updated_values_iter" + to_string(iter) + ".txt";
        save_values(outname, A, E);
    }

    cout << "Done. Wrote updated_values_iter*.txt\n";
    return 0;
}
