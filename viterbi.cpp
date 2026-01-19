#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <array>
#include <algorithm>

using namespace std;

// NEG_INF koristimo kao "log(0)" (praktično -beskonačno).
// Kad ne želimo dopustiti neku opciju u Viterbiju (npr. krivo stanje za opažanje),
// vratimo NEG_INF da taj put nikad ne pobijedi.
static const double NEG_INF = -1e300;

// Pretvorba vjerojatnosti u log-vjerojatnost.
// U Viterbiju zbrajamo logove umjesto da množimo vjerojatnosti (stabilnije numerički).
static inline double logp(double p) { return (p <= 0.0) ? NEG_INF : log(p); }

// -------------------- BINARNI FORMAT: čitanje null-terminated stringova --------------------
// Writer je zapisivao stringove kao bytes + '\0' terminator.
// Ova funkcija čita bajtove dok ne naiđe na '\0'.
static string readCString(ifstream &input)
{
    string result;
    char ch;
    while (input.get(ch))
    {
        if (ch == '\0')
            break;
        result.push_back(ch);
    }
    return result;
}

// -------------------- MAPIRANJE KOLONE PORAVNANJA U TIP STANJA (M/I/D) --------------------
// Za pair-HMM, svaki stupac poravnanja ima prirodan tip:
// - oba nisu '-'  -> Match (1)
// - prvi je '-'   -> Insertion (2)   (opažanje "-X")
// - drugi je '-'  -> Deletion  (3)   (opažanje "X-")
// - oba '-'       -> ignoriramo ("--")
static inline int obs_state(char c1, char c2)
{
    bool g1 = (c1 == '-');
    bool g2 = (c2 == '-');
    if (!g1 && !g2)
        return 1; // Match
    if (g1 && !g2)
        return 2; // Ins
    if (!g1 && g2)
        return 3; // Del
    return -1; // "--"
}

// -------------------- EMISIJE: log P(o | state) --------------------
// E je vektor mapa: e[state][obs] = vjerojatnost emisije
// - state=1: Match emisije poput "AA","CT",...
// - state=2: Insertion emisije poput "-A","-C",...
// - state=3: Deletion emisije poput "A-","C-",...
//
// floor_prob služi kao "mala vjerojatnost" ako ključ ne postoji u mapi (da ne pukne model).
static inline double emit_log(int state, const string &o,
                              const vector<unordered_map<string, double>> &e,
                              double floor_prob)
{
    auto it = e[state].find(o);
    double p = (it == e[state].end()) ? floor_prob : it->second;
    return logp(p);
}

// Ako je initial_values.txt nastao na Windowsu, može imati '\r' na kraju tokena.
// Ovo uklanja taj '\r' da ključevi budu točni (npr. "AA" umjesto "AA\r").
static inline void trim_cr(string &s)
{
    if (!s.empty() && s.back() == '\r')
        s.pop_back();
}

// -------------------- UČITAVANJE PARAMETARA HMM-a (A i E) --------------------
// Učitavamo iz initial_values.txt:
// 1) prvih 25 brojeva (5x5) -> tranzicijska matrica A
// 2) zatim emisije: "KEY : PROB"
//
// A dimenzije 5x5 za stanja: Begin, Match, Ins, Del, End.
static void load_initial_values(const string &path,
                                vector<vector<double>> &a,
                                vector<unordered_map<string, double>> &e)
{
    const int S = 5; // Begin, Match, Ins, Del, End
    a.assign(S, vector<double>(S, 0.0));
    e.assign(S, unordered_map<string, double>());

    ifstream file(path);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << path << "\n";
        exit(1);
    }

    // (1) Čitanje matrice A: 25 brojeva (row-major)
    for (int i = 0; i < S; ++i)
    {
        for (int j = 0; j < S; ++j)
        {
            string tok;
            if (!(file >> tok))
            {
                cerr << "Unexpected EOF while reading transition matrix A\n";
                exit(1);
            }
            trim_cr(tok);
            a[i][j] = stod(tok);
        }
    }

    // (2) Čitanje emisija: <key> : <prob>
    // Ključevi su npr. "AA", "CT", "-A", "A-"
    string key;
    while (file >> key)
    {
        trim_cr(key);

        string colon;
        double p;
        if (!(file >> colon >> p))
            break;
        trim_cr(colon);
        trim_cr(key);

        // Ovdje klasificiramo emisiju u jedno od 3 emitirajuća stanja:
        // "-A" -> state 2 (Insertion)
        // "A-" -> state 3 (Deletion)
        // "AA" -> state 1 (Match)
        if (key.size() >= 2 && key[0] == '-')
            e[2][key] = p; // Insertion
        else if (key.size() >= 2 && key[1] == '-')
            e[3][key] = p; // Deletion
        else
            e[1][key] = p; // Match
    }
}

// -------------------- STRUKTURA REZULTATA ZA JEDAN PAR --------------------
struct Summary
{
    string n1, n2;
    int L = 0;            // broj kolona opažanja (bez "--")
    int nM = 0, nI = 0, nD = 0;       // koliko je koraka Viterbi put bio u M/I/D
    int blocksI = 0, blocksD = 0;     // koliko "blokova" insertion/deletion (broj ulazaka u I/D segmente)
    double logScore = NEG_INF;        // log-score najboljeg puta
};

// -------------------- VITERBI ZA JEDAN PAR PORAVNANIH SEKVENCI --------------------
// Ulaz: poravnate sekvence s1 i s2 (s '-' gapovima)
// Iz njih gradimo niz opažanja O[t] (string dužine 2) i colType[t] (1/2/3)
// Zatim radimo Viterbi DP nad stanjima (M/I/D).
static Summary viterbi_one(const string &n1, const string &n2,
                           const string &s1, const string &s2,
                           const vector<vector<double>> &A,
                           const vector<unordered_map<string, double>> &E,
                           double emit_floor = 1e-12,
                           double trans_floor = 1e-15)
{
    // Indeksi stanja u A:
    const int BEGIN = 0, M = 1, I = 2, D = 3, END = 4;

    size_t L0 = min(s1.size(), s2.size());

    // (A) Pretvori poravnanje u opažanja O i tip kolone colType.
    // O[t] je npr. "AA", "-A", "C-".
    // colType[t] je 1/2/3 (M/I/D).
    vector<string> O;
    vector<int> colType; // 1=Match, 2=Ins, 3=Del
    O.reserve(L0);
    colType.reserve(L0);

    for (size_t j = 0; j < L0; ++j)
    {
        char c1 = s1[j], c2 = s2[j];

        // "--" kolone ignoriramo (nema informacije)
        if (c1 == '-' && c2 == '-')
            continue;

        int t = obs_state(c1, c2);
        if (t < 0)
            continue;

        string obs;
        obs.push_back(c1);
        obs.push_back(c2);

        O.push_back(obs);
        colType.push_back(t);
    }

    Summary out;
    out.n1 = n1;
    out.n2 = n2;
    out.L = (int)O.size();
    if (out.L == 0)
        return out;

    // (B) Tranzicijska log-prob funkcija: log A[from][to],
    // uz floor za nule (da ne dobijemo -inf i da model ostane robustan).
    auto trans = [&](int from, int to){
        if (from == END) return NEG_INF;      // ništa nakon End
        if (to == BEGIN) return NEG_INF;      // nitko ne ide u Begin
        double p = A[from][to];
        if (p <= 0.0) p = trans_floor;
        return logp(p);
    };

    // (C) Ključna stvar: MASKIRANE EMISIJE
    // Dozvoljavamo da se opažanje emitira samo iz stanja koje ima smisla:
    // - ako je colType[t]==D, onda samo state D smije emitirati "X-"
    // - ostali dobiju -inf
    auto emit_masked = [&](int state, int t) -> double {
        if (colType[t] != state)
            return NEG_INF;
        return emit_log(state, O[t], E, emit_floor);
    };

    // (D) Viterbi DP:
    // dp[t][k] = najbolja log-vjerojatnost za prvi t+1 opažanja,
    //           završavajući u stanju k (0=M, 1=I, 2=D).
    vector<array<double, 3>> dp(out.L);
    // back[t][k] = iz kojeg je prethodnog stanja došlo najbolje rješenje
    //             (0=M, 1=I, 2=D).
    vector<array<int, 3>> back(out.L);

    // (E) INICIJALIZACIJA za t=0 (prvo opažanje):
    dp[0][0] = trans(BEGIN, M) + emit_masked(M, 0);
    dp[0][1] = trans(BEGIN, I) + emit_masked(I, 0);
    dp[0][2] = trans(BEGIN, D) + emit_masked(D, 0);
    back[0] = {-1, -1, -1};

    // (F) REKURZIJA: za svaki t od 1 do L-1 računamo najbolji put do M/I/D
    for (int t = 1; t < out.L; ++t)
    {
        // ---- prijelaz u M (dp[t][0]) ----
        {
            double best = dp[t - 1][0] + trans(M, M);
            int arg = 0;

            double cand = dp[t - 1][1] + trans(I, M);
            if (cand > best) { best = cand; arg = 1; }

            cand = dp[t - 1][2] + trans(D, M);
            if (cand > best) { best = cand; arg = 2; }

            dp[t][0] = best + emit_masked(M, t);
            back[t][0] = arg;
        }

        // ---- prijelaz u I (dp[t][1]) ----
        {
            double best = dp[t - 1][0] + trans(M, I);
            int arg = 0;

            double cand = dp[t - 1][1] + trans(I, I);
            if (cand > best) { best = cand; arg = 1; }

            cand = dp[t - 1][2] + trans(D, I);
            if (cand > best) { best = cand; arg = 2; }

            dp[t][1] = best + emit_masked(I, t);
            back[t][1] = arg;
        }

        // ---- prijelaz u D (dp[t][2]) ----
        {
            double best = dp[t - 1][0] + trans(M, D);
            int arg = 0;

            double cand = dp[t - 1][1] + trans(I, D);
            if (cand > best) { best = cand; arg = 1; }

            cand = dp[t - 1][2] + trans(D, D);
            if (cand > best) { best = cand; arg = 2; }

            dp[t][2] = best + emit_masked(D, t);
            back[t][2] = arg;
        }
    }

    // (G) TERMINACIJA: dodaj prijelaz u End i izaberi najbolji završni state
    double best = dp[out.L - 1][0] + trans(M, END);
    int last = 0;

    double cand = dp[out.L - 1][1] + trans(I, END);
    if (cand > best) { best = cand; last = 1; }

    cand = dp[out.L - 1][2] + trans(D, END);
    if (cand > best) { best = cand; last = 2; }

    out.logScore = best;

    // (H) BACKTRACKING: rekonstruiraj najbolji put stanja unatrag
    vector<int> path(out.L);
    path[out.L - 1] = last;

    for (int t = out.L - 1; t > 0; --t)
    {
        path[t - 1] = back[t][path[t]];
        if (path[t - 1] < 0)
            path[t - 1] = 0;
    }

    // (I) SUMARIZACIJA: prebroji koliko puta je put bio u M/I/D i koliko blokova I/D
    int prev = -1;
    for (int t = 0; t < out.L; ++t)
    {
        int st = path[t]; // 0=M,1=I,2=D (ovo su indeksi dp/back, ne A indeksi)

        if (st == 0) out.nM++;
        else if (st == 1) out.nI++;
        else out.nD++;

        // "blok" znači ulazak u I (ili D) iz ne-I (ne-D)
        if (st == 1 && prev != 1) out.blocksI++;
        if (st == 2 && prev != 2) out.blocksD++;
        prev = st;
    }

    return out;
}

int main(int argc, char *argv[])
{
    // Očekuje: input.bin, initial_values.txt, opcionalno output.csv
    if (argc < 3)
    {
        cerr << "Usage: viterbi <input.bin> <initial_values.txt> [output.csv]\n";
        return 1;
    }

    string binfile = argv[1];
    string initfile = argv[2];
    string csvfile = (argc >= 4) ? argv[3] : string("viterbi_summary.csv");

    // (1) Učitaj HMM parametre A i E
    vector<vector<double>> A;
    vector<unordered_map<string, double>> E;
    load_initial_values(initfile, A, E);

    // (2) Otvori binarni dataset (iz writer.cpp)
    ifstream in(binfile, ios::binary);
    if (!in.is_open())
    {
        cerr << "Error opening " << binfile << "\n";
        return 1;
    }

    // (3) Otvori CSV izlaz
    ofstream out(csvfile);
    if (!out.is_open())
    {
        cerr << "Error opening " << csvfile << "\n";
        return 1;
    }

    // Header za CSV: dodali smo gap_rate = fracI + fracD
    out << "name1,name2,L,nM,nI,nD,fracM,fracI,fracD,gap_rate,blocksI,blocksD,logScore\n";

    long long total = 0;
    int printed = 0;

    // (4) Čitaj binarni format:
    // svaka "grupa" počinje s cleanedCount (int), nakon toga ide cleanedCount zapisa:
    // name1\0 name2\0 seq1\0 seq2\0
    while (true)
    {
        int cleanedCount;
        if (!in.read(reinterpret_cast<char *>(&cleanedCount), sizeof(int)))
            break;
        if (cleanedCount <= 0)
            continue;

        for (int i = 0; i < cleanedCount; ++i)
        {
            // Čitanje jednog zapisa (4 stringa)
            string n1 = readCString(in);
            if (!in) return 1;
            string n2 = readCString(in);
            if (!in) return 1;
            string s1 = readCString(in);
            if (!in) return 1;
            string s2 = readCString(in);
            if (!in) return 1;

            total++;

            // (5) Viterbi dekodiranje za taj par
            Summary s = viterbi_one(n1, n2, s1, s2, A, E);

            // (6) Izračun frakcija i gap_rate
            double fracM = s.L ? (double)s.nM / s.L : 0.0;
            double fracI = s.L ? (double)s.nI / s.L : 0.0;
            double fracD = s.L ? (double)s.nD / s.L : 0.0;
            double gap_rate = fracI + fracD;

            // (7) Upis u CSV
            out << s.n1 << "," << s.n2 << ","
                << s.L << "," << s.nM << "," << s.nI << "," << s.nD << ","
                << fracM << "," << fracI << "," << fracD << ","
                << gap_rate << ","
                << s.blocksI << "," << s.blocksD << ","
                << s.logScore << "\n";

            // (8) Ispis prvih 20 u konzolu (sanity check)
            if (printed < 20)
            {
                printed++;
                cout << "[" << printed << "] "
                     << s.n1 << " vs " << s.n2
                     << " | L=" << s.L
                     << " | frac(M/I/D)=" << fracM << "/" << fracI << "/" << fracD
                     << " | blocks(I/D)=" << s.blocksI << "/" << s.blocksD
                     << " | logScore=" << s.logScore << "\n";
            }
        }
    }

    in.close();
    out.close();

    // (9) Kratki sanity print da emisije postoje (opcijski)
    cerr << "Check emissions: "
         << "AA=" << (E[1].count("AA") ? E[1].at("AA") : -1)
         << " -A=" << (E[2].count("-A") ? E[2].at("-A") : -1)
         << " A-=" << (E[3].count("A-") ? E[3].at("A-") : -1)
         << "\n";

    cout << "Done. Decoded alignments: " << total << "\n";
    cout << "CSV written: " << csvfile << "\n";
    return 0;
}
