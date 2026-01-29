#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>

using namespace std;

static const double NEG_INF = -1e300;
static inline double logp(double p) { return (p <= 0.0) ? NEG_INF : log(p); }

// -------------------- UČITAVANJE HMM PARAMETARA --------------------
static inline void trim_cr(string &s)
{
    if (!s.empty() && s.back() == '\r')
        s.pop_back();
}

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

    // Učitaj tranzicijsku matricu A (5x5)
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

    // Učitaj emisije: <key> : <prob>
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

        // Klasifikacija:
        // "-A" -> Insertion (state 2)
        // "A-" -> Deletion (state 3)
        // "AA" -> Match (state 1)
        if (key.size() >= 2 && key[0] == '-')
            e[2][key] = p;
        else if (key.size() >= 2 && key[1] == '-')
            e[3][key] = p;
        else
            e[1][key] = p;
    }
}

// -------------------- EMISIJSKA VJEROJATNOST --------------------
// Za Match: kombinira oba znaka (npr. "AA", "CT")
// Za Ix: kombinira samo x znak (npr. "-A", "-G")
// Za Iy: kombinira samo y znak (npr. "A-", "G-")
static inline double emit_log(int state, char x, char y,
                              const vector<unordered_map<string, double>> &E,
                              double floor_prob = 1e-12)
{
    string obs;
    if (state == 1)
    { // Match
        obs.push_back(x);
        obs.push_back(y);
    }
    else if (state == 2)
    { // Insertion in x (gap u y)
        obs.push_back('-');
        obs.push_back(y);
    }
    else if (state == 3)
    { // Deletion from x (gap u x)
        obs.push_back(x);
        obs.push_back('-');
    }
    else
    {
        return NEG_INF;
    }

    auto it = E[state].find(obs);
    double p = (it == E[state].end()) ? floor_prob : it->second;
    return logp(p);
}

// -------------------- PAIR-HMM VITERBI ALIGNMENT --------------------
// Stvara poravnanje od neporavnatih sekvenci x i y
struct AlignmentResult
{
    string aligned_x;
    string aligned_y;
    double logScore;
    int nM, nIx, nIy;
};

static AlignmentResult viterbi_align(const string &x, const string &y,
                                     const vector<vector<double>> &A,
                                     const vector<unordered_map<string, double>> &E,
                                     double trans_floor = 1e-15)
{
    const int BEGIN = 0, M = 1, Ix = 2, Iy = 3, END = 4;

    int lenX = (int)x.size();
    int lenY = (int)y.size();

    // DP tablice: dp[i][j][state] = najbolja log-vjerojatnost do pozicije (i,j) u stanju state
    // state: 0=M, 1=Ix, 2=Iy
    vector<vector<vector<double>>> dp(lenX + 1,
                                      vector<vector<double>>(lenY + 1, vector<double>(3, NEG_INF)));
    vector<vector<vector<int>>> back(lenX + 1,
                                     vector<vector<int>>(lenY + 1, vector<int>(3, -1)));

    // Tranzicijska funkcija
    auto trans = [&](int from, int to) -> double
    {
        if (from == END)
            return NEG_INF;
        if (to == BEGIN)
            return NEG_INF;
        double p = A[from][to];
        if (p <= 0.0)
            p = trans_floor;
        return logp(p);
    };

    // INICIJALIZACIJA: BEGIN -> M/Ix/Iy
    dp[0][0][0] = trans(BEGIN, M);
    dp[0][0][1] = trans(BEGIN, Ix);
    dp[0][0][2] = trans(BEGIN, Iy);

    // REKURZIJA
    for (int i = 0; i <= lenX; ++i)
    {
        for (int j = 0; j <= lenY; ++j)
        {
            // M(i,j) - konzumira oba: x[i-1] i y[j-1]
            if (i > 0 && j > 0)
            {
                double emit = emit_log(M, x[i - 1], y[j - 1], E);

                // Dolazi iz M(i-1,j-1)
                double scoreM = dp[i - 1][j - 1][0] + trans(M, M) + emit;
                // Dolazi iz Ix(i-1,j-1)
                double scoreIx = dp[i - 1][j - 1][1] + trans(Ix, M) + emit;
                // Dolazi iz Iy(i-1,j-1)
                double scoreIy = dp[i - 1][j - 1][2] + trans(Iy, M) + emit;

                if (scoreM >= scoreIx && scoreM >= scoreIy)
                {
                    dp[i][j][0] = scoreM;
                    back[i][j][0] = 0; // iz M
                }
                else if (scoreIx >= scoreIy)
                {
                    dp[i][j][0] = scoreIx;
                    back[i][j][0] = 1; // iz Ix
                }
                else
                {
                    dp[i][j][0] = scoreIy;
                    back[i][j][0] = 2; // iz Iy
                }
            }

            // Ix(i,j) - konzumira samo x[i-1], gap u y
            if (i > 0)
            {
                double emit = emit_log(Ix, x[i - 1], 'X', E); // 'X' placeholder, koristi se '-'

                // Dolazi iz M(i-1,j) - otvara gap
                double scoreM = dp[i - 1][j][0] + trans(M, Ix) + emit;
                // Dolazi iz Ix(i-1,j) - produži gap
                double scoreIx = dp[i - 1][j][1] + trans(Ix, Ix) + emit;

                if (scoreM >= scoreIx)
                {
                    dp[i][j][1] = scoreM;
                    back[i][j][1] = 0; // iz M
                }
                else
                {
                    dp[i][j][1] = scoreIx;
                    back[i][j][1] = 1; // iz Ix
                }
            }

            // Iy(i,j) - konzumira samo y[j-1], gap u x
            if (j > 0)
            {
                double emit = emit_log(Iy, 'X', y[j - 1], E);

                // Dolazi iz M(i,j-1) - otvara gap
                double scoreM = dp[i][j - 1][0] + trans(M, Iy) + emit;
                // Dolazi iz Iy(i,j-1) - produži gap
                double scoreIy = dp[i][j - 1][2] + trans(Iy, Iy) + emit;

                if (scoreM >= scoreIy)
                {
                    dp[i][j][2] = scoreM;
                    back[i][j][2] = 0; // iz M
                }
                else
                {
                    dp[i][j][2] = scoreIy;
                    back[i][j][2] = 2; // iz Iy
                }
            }
        }
    }

    // TERMINACIJA: pronađi najbolje završno stanje
    double bestScore = dp[lenX][lenY][0] + trans(M, END);
    int lastState = 0;

    double scoreIx = dp[lenX][lenY][1] + trans(Ix, END);
    if (scoreIx > bestScore)
    {
        bestScore = scoreIx;
        lastState = 1;
    }

    double scoreIy = dp[lenX][lenY][2] + trans(Iy, END);
    if (scoreIy > bestScore)
    {
        bestScore = scoreIy;
        lastState = 2;
    }

    // BACKTRACKING: rekonstruiraj poravnanje
    string aligned_x, aligned_y;
    int i = lenX, j = lenY;
    int state = lastState;

    int countM = 0, countIx = 0, countIy = 0;

    while (i > 0 || j > 0)
    {
        if (state == 0)
        { // Match
            if (i > 0 && j > 0)
            {
                aligned_x.push_back(x[i - 1]);
                aligned_y.push_back(y[j - 1]);
                countM++;
                int prevState = back[i][j][0];
                i--;
                j--;
                state = prevState;
            }
            else
            {
                break;
            }
        }
        else if (state == 1)
        { // Ix - gap u y
            if (i > 0)
            {
                aligned_x.push_back(x[i - 1]);
                aligned_y.push_back('-');
                countIx++;
                int prevState = back[i][j][1];
                i--;
                state = prevState;
            }
            else
            {
                break;
            }
        }
        else if (state == 2)
        { // Iy - gap u x
            if (j > 0)
            {
                aligned_x.push_back('-');
                aligned_y.push_back(y[j - 1]);
                countIy++;
                int prevState = back[i][j][2];
                j--;
                state = prevState;
            }
            else
            {
                break;
            }
        }
        else
        {
            break;
        }
    }

    // Poravnanje je naopako (backtracking ide unatrag)
    reverse(aligned_x.begin(), aligned_x.end());
    reverse(aligned_y.begin(), aligned_y.end());

    AlignmentResult result;
    result.aligned_x = aligned_x;
    result.aligned_y = aligned_y;
    result.logScore = bestScore;
    result.nM = countM;
    result.nIx = countIx;
    result.nIy = countIy;

    return result;
}

// -------------------- MAIN --------------------
int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        cerr << "Usage: viterbi2 <seq1> <seq2> <initial_values.txt> [output.txt]\n";
        cerr << "Example: viterbi2 ACGT AGCT initial_values.txt alignment.txt\n";
        return 1;
    }

    string seq1 = argv[1];
    string seq2 = argv[2];
    string initfile = argv[3];
    string outfile = (argc >= 5) ? argv[4] : "alignment_out.txt";

    // Učitaj HMM parametre
    vector<vector<double>> A;
    vector<unordered_map<string, double>> E;
    load_initial_values(initfile, A, E);

    cout << "Aligning sequences:\n";
    cout << "X: " << seq1 << " (len=" << seq1.size() << ")\n";
    cout << "Y: " << seq2 << " (len=" << seq2.size() << ")\n\n";

    // Poravnaj sekvence
    AlignmentResult result = viterbi_align(seq1, seq2, A, E);

    // Ispis rezultata
    cout << "Alignment result:\n";
    cout << "X: " << result.aligned_x << "\n";
    cout << "Y: " << result.aligned_y << "\n";
    cout << "Log-score: " << result.logScore << "\n";
    cout << "Match states: " << result.nM << "\n";
    cout << "Insertion (gap in Y): " << result.nIx << "\n";
    cout << "Deletion (gap in X): " << result.nIy << "\n";

    // Spremi u file
    ofstream out(outfile);
    if (out.is_open())
    {
        out << "# Pairwise alignment using HMM Viterbi\n";
        out << "# X: " << seq1 << "\n";
        out << "# Y: " << seq2 << "\n";
        out << "# Log-score: " << result.logScore << "\n";
        out << "#\n";
        out << result.aligned_x << "\n";
        out << result.aligned_y << "\n";
        out.close();
        cout << "\nAlignment saved to: " << outfile << "\n";
    }

    return 0;
}
