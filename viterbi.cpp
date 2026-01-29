#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <omp.h>

using namespace std;

static const double NEG_INF = -1e300;
static inline double logp(double p) { return (p <= 0.0) ? NEG_INF : log(p); }

static inline void trim_cr(string &s)
{
    if (!s.empty() && s.back() == '\r')
        s.pop_back();
}

static void load_initial_values(const string &path,
                                vector<vector<double>> &a,
                                vector<unordered_map<string, double>> &e)
{
    a.assign(5, vector<double>(5, 0.0));
    e.assign(5, unordered_map<string, double>());

    ifstream file(path);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << path << "\n";
        exit(1);
    }

    // Matrica A
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
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

    // Matrica E
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

        // -X -> Insertion
        // X- -> Deletion
        // XX -> Match
        if (key[0] == '-')
            e[2][key] = p;
        else if (key[1] == '-')
            e[3][key] = p;
        else
            e[1][key] = p;
    }
}

static inline double emit_log(int state, char x, char y,
                              const vector<unordered_map<string, double>> &E,
                              double floor_prob = 1e-12)
{
    string obs;
    if (state == 1) // Match
    {
        obs.push_back(x);
        obs.push_back(y);
    }
    else if (state == 2) // Insertion
    {
        obs.push_back('-');
        obs.push_back(y);
    }
    else if (state == 3) // Deletion
    {
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

    // dp[i][j][state] = najbolja log-vjerojatnost do pozicije (i,j) u stanju state
    vector<vector<vector<double>>> dp(lenX + 1,
                                      vector<vector<double>>(lenY + 1, vector<double>(3, NEG_INF)));
    vector<vector<vector<int>>> back(lenX + 1,
                                     vector<vector<int>>(lenY + 1, vector<int>(3, -1)));

    // lambda funkcija koja daje log-vjerojatnost prijelaza
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

    // Inicijalizacija
    dp[0][0][0] = trans(BEGIN, M);
    dp[0][0][1] = trans(BEGIN, Ix);
    dp[0][0][2] = trans(BEGIN, Iy);

    for (int i = 0; i <= lenX; ++i)
    {
        for (int j = 0; j <= lenY; ++j)
        {
            // Match
            if (i > 0 && j > 0)
            {
                double emit = emit_log(M, x[i - 1], y[j - 1], E);

                double scoreM = dp[i - 1][j - 1][0] + trans(M, M) + emit;
                double scoreIx = dp[i - 1][j - 1][1] + trans(Ix, M) + emit;
                double scoreIy = dp[i - 1][j - 1][2] + trans(Iy, M) + emit;

                if (scoreM >= scoreIx && scoreM >= scoreIy)
                {
                    dp[i][j][0] = scoreM;
                    back[i][j][0] = M;
                }
                else if (scoreIx >= scoreIy)
                {
                    dp[i][j][0] = scoreIx;
                    back[i][j][0] = Ix;
                }
                else
                {
                    dp[i][j][0] = scoreIy;
                    back[i][j][0] = Iy;
                }
            }

            // Insertion
            if (i > 0)
            {
                double emit = emit_log(Ix, x[i - 1], '-', E);

                double scoreM = dp[i - 1][j][0] + trans(M, Ix) + emit;
                double scoreIx = dp[i - 1][j][1] + trans(Ix, Ix) + emit;

                if (scoreM >= scoreIx)
                {
                    dp[i][j][1] = scoreM;
                    back[i][j][1] = M;
                }
                else
                {
                    dp[i][j][1] = scoreIx;
                    back[i][j][1] = Ix;
                }
            }

            // Deletion
            if (j > 0)
            {
                double emit = emit_log(Iy, '-', y[j - 1], E);

                double scoreM = dp[i][j - 1][0] + trans(M, Iy) + emit;
                double scoreIy = dp[i][j - 1][2] + trans(Iy, Iy) + emit;

                if (scoreM >= scoreIy)
                {
                    dp[i][j][2] = scoreM;
                    back[i][j][2] = M;
                }
                else
                {
                    dp[i][j][2] = scoreIy;
                    back[i][j][2] = Iy;
                }
            }
        }
    }

    double bestScore = dp[lenX][lenY][0] + trans(M, END);
    int lastState = M;

    double scoreIx = dp[lenX][lenY][1] + trans(Ix, END);
    if (scoreIx > bestScore)
    {
        bestScore = scoreIx;
        lastState = Ix;
    }

    double scoreIy = dp[lenX][lenY][2] + trans(Iy, END);
    if (scoreIy > bestScore)
    {
        bestScore = scoreIy;
        lastState = Iy;
    }

    // Rekonstrukcija poravnanja
    string aligned_x, aligned_y;
    int i = lenX, j = lenY;
    int state = lastState;

    int countM = 0, countIx = 0, countIy = 0;

    while (i > 0 || j > 0)
    {
        if (state == M)
        { // Match
            if (i > 0 && j > 0)
            {
                aligned_x.push_back(x[i - 1]);
                aligned_y.push_back(y[j - 1]);
                countM++;
                state = back[i][j][M];
                i--;
                j--;
            }
            else
            {
                break;
            }
        }
        else if (state == Ix)
        { // Insert
            if (i > 0)
            {
                aligned_x.push_back(x[i - 1]);
                aligned_y.push_back('-');
                countIx++;
                state = back[i][j][Ix];
                i--;
            }
            else
            {
                break;
            }
        }
        else if (state == Iy)
        { // Delete
            if (j > 0)
            {
                aligned_x.push_back('-');
                aligned_y.push_back(y[j - 1]);
                countIy++;
                state = back[i][j][Iy];
                j--;
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

    // Poravnanje je naopako
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

string readCString(ifstream &input)
{
    string result;
    char ch;
    // citaj do NULL ili EOF
    while (input.get(ch))
    {
        if (ch == '\0')
            break;
        result.push_back(ch);
    }
    return result;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Usage: viterbi <input.bin> <initial_values.txt> [output.txt]\n";
        cerr << "Example: viterbi input.bin initial_values.txt alignment.txt\n";
        return 1;
    }

    string inputFile = argv[1];
    string initfile = argv[2];
    string outfile = (argc >= 4) ? argv[3] : "alignment_out.txt";

    vector<vector<double>> A;
    vector<unordered_map<string, double>> E;
    load_initial_values(initfile, A, E);

    ifstream file(inputFile, ios::binary);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << inputFile << endl;
        return 1;
    }

    vector<tuple<string, string, string, string>> sequencePairs;
    int cleanedCount;
    while (true)
    {
        if (!file.read(reinterpret_cast<char *>(&cleanedCount), sizeof(int)))
            break;
        if (cleanedCount <= 0)
            continue;
        for (int i = 0; i < cleanedCount; ++i)
        {
            string name1 = readCString(file);
            if (!file)
            {
                cerr << "Unexpected EOF reading name1\n";
                exit(1);
            }
            string name2 = readCString(file);
            if (!file)
            {
                cerr << "Unexpected EOF reading name2\n";
                exit(1);
            }
            string seq1 = readCString(file);
            if (!file)
            {
                cerr << "Unexpected EOF reading seq1\n";
                exit(1);
            }
            string seq2 = readCString(file);
            if (!file)
            {
                cerr << "Unexpected EOF reading seq2\n";
                exit(1);
            }
            sequencePairs.push_back(make_tuple(name1, name2, seq1, seq2));
        }
    }
    file.close();

    vector<AlignmentResult> results(sequencePairs.size());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)sequencePairs.size(); ++i)
    {
        string name1 = get<0>(sequencePairs[i]);
        string name2 = get<1>(sequencePairs[i]);
        string seq1 = get<2>(sequencePairs[i]);
        string seq2 = get<3>(sequencePairs[i]);

        // Ukloni gapove
        seq1.erase(remove(seq1.begin(), seq1.end(), '-'), seq1.end());
        seq2.erase(remove(seq2.begin(), seq2.end(), '-'), seq2.end());

        AlignmentResult result = viterbi_align(seq1, seq2, A, E);
        results[i] = result;
    }

    ofstream out(outfile);
    for (size_t i = 0; i < sequencePairs.size(); ++i)
    {
        string seq1 = get<2>(sequencePairs[i]);
        string seq2 = get<3>(sequencePairs[i]);
        seq1.erase(remove(seq1.begin(), seq1.end(), '-'), seq1.end());
        seq2.erase(remove(seq2.begin(), seq2.end(), '-'), seq2.end());

        if (out.is_open())
        {
            out << "X: " << seq1 << "\n";
            out << "Y: " << seq2 << "\n";
            out << "Log-score: " << results[i].logScore << "\n";
            out << "\n";
            out << results[i].aligned_x << "\n";
            out << results[i].aligned_y << "\n\n";
        }
    }
    out.close();
    cout << "\nAlignment saved to: " << outfile << "\n";

    return 0;
}
