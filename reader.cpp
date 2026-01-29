// reader.cpp
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <map>
#include <ctime>

using namespace std;

namespace std
{
    template <>
    struct hash<pair<string, string>>
    {
        size_t operator()(const pair<string, string> &p) const noexcept
        {
            size_t h1 = hash<string>{}(p.first);
            size_t h2 = hash<string>{}(p.second);
            return h1 ^ (h2 << 1);
        }
    };
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

double safe_div(double num, double den)
{
    if (den <= 0.0)
        return 0.0;
    return num / den;
}

void read_and_predict_hmm(const string &input,
                          unordered_map<pair<string, string>, pair<string, string>, hash<pair<string, string>>> &alignments,
                          vector<vector<double>> &a,
                          vector<unordered_map<string, double>> &e)
{
    ifstream file(input, ios::binary);
    if (!file.is_open())
    {
        cerr << "Greška pri otvaranju datoteke: " << input << endl;
        exit(1);
    }

    long long beg_match = 0, beg_ins = 0, beg_del = 0;
    long long match_match = 0, match_ins = 0, match_del = 0, match_end = 0;
    long long ins_match = 0, ins_ins = 0, ins_del = 0, ins_end = 0;
    long long del_match = 0, del_del = 0, del_ins = 0, del_end = 0;
    long long beg = 0, ins = 0, del = 0, match = 0;

    if (a.size() < 4)
        a.assign(4, vector<double>(5, 0.0));
    else
    {
        for (auto &row : a)
            if (row.size() < 5)
                row.assign(5, 0.0);
    }
    if (e.size() < 4)
        e.assign(4, unordered_map<string, double>());

    while (true)
    {
        int cleanedCount;
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

            alignments[{name1, name2}] = {seq1, seq2};

            size_t L = min(seq1.length(), seq2.length());
            for (size_t j = 0; j < L; ++j)
            {
                bool g1 = (seq1[j] == '-');
                bool g2 = (seq2[j] == '-');

                if (j == 0)
                {
                    if (!g1 && !g2)
                    {
                        beg_match++;
                        match++;
                        e[1][seq1.substr(j, 1) + seq2.substr(j, 1)] += 1.0;
                    }
                    else if (g1 && !g2)
                    {
                        beg_ins++;
                        ins++;
                        e[2]["-" + seq2.substr(j, 1)] += 1.0;
                    }
                    else if (!g1 && g2)
                    {
                        beg_del++;
                        del++;
                        e[3][seq1.substr(j, 1) + "-"] += 1.0;
                    }
                }
                else
                {
                    if (!g1 && !g2)
                    {
                        match++;
                        bool pg1 = (seq1[j - 1] == '-');
                        bool pg2 = (seq2[j - 1] == '-');
                        if (!pg1 && !pg2)
                            match_match++;
                        else if (pg1 && !pg2)
                            ins_match++;
                        else if (!pg1 && pg2)
                            del_match++;
                        if (j == L - 1)
                            match_end++;
                        e[1][seq1.substr(j, 1) + seq2.substr(j, 1)] += 1.0;
                    }
                    else if (g1 && !g2)
                    {
                        ins++;
                        bool pg1 = (seq1[j - 1] == '-');
                        bool pg2 = (seq2[j - 1] == '-');
                        if (!pg1 && !pg2)
                            match_ins++;
                        else if (pg1 && !pg2)
                            ins_ins++;
                        else if (!pg1 && pg2)
                            del_ins++;
                        if (j == L - 1)
                            ins_end++;
                        e[2]["-" + seq2.substr(j, 1)] += 1.0;
                    }
                    else if (!g1 && g2)
                    {
                        del++;
                        bool pg1 = (seq1[j - 1] == '-');
                        bool pg2 = (seq2[j - 1] == '-');
                        if (!pg1 && !pg2)
                            match_del++;
                        else if (pg1 && !pg2)
                            ins_del++;
                        else if (!pg1 && pg2)
                            del_del++;
                        if (j == L - 1)
                            del_end++;
                        e[3][seq1.substr(j, 1) + "-"] += 1.0;
                    }
                }
            }
        }
        beg += cleanedCount;
    }

    a[0][0] = 0.0;
    a[0][1] = safe_div((double)beg_match, (double)beg);
    a[0][2] = safe_div((double)beg_ins, (double)beg);
    a[0][3] = safe_div((double)beg_del, (double)beg);
    a[0][4] = 0.0;

    a[1][0] = 0.0;
    a[1][1] = safe_div((double)match_match, (double)(match));
    a[1][2] = safe_div((double)match_ins, (double)(match));
    a[1][3] = safe_div((double)match_del, (double)(match));
    a[1][4] = safe_div((double)match_end, (double)(match));

    a[2][0] = 0.0;
    a[2][1] = safe_div((double)ins_match, (double)(ins));
    a[2][2] = safe_div((double)ins_ins, (double)(ins));
    a[2][3] = safe_div((double)ins_del, (double)(ins));
    a[2][4] = safe_div((double)ins_end, (double)(ins));

    a[3][0] = 0.0;
    a[3][1] = safe_div((double)del_match, (double)(del));
    a[3][2] = safe_div((double)del_ins, (double)(del));
    a[3][3] = safe_div((double)del_del, (double)(del));
    a[3][4] = safe_div((double)del_end, (double)(del));

    // normalize emissions
    if (match > 0)
        for (auto &kv : e[1])
            kv.second /= (double)match;
    if (ins > 0)
        for (auto &kv : e[2])
            kv.second /= (double)ins;
    if (del > 0)
        for (auto &kv : e[3])
            kv.second /= (double)del;

    file.close();
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Usage: reader <input.bin>\n";
        return 1;
    }
    string input = argv[1];

    unordered_map<pair<string, string>, pair<string, string>> alignments;
    vector<vector<double>> a(5, vector<double>(5, 0.0));
    vector<unordered_map<string, double>> e(4);

    time_t start = time(nullptr);

    read_and_predict_hmm(input, alignments, a, e);

    cout << "Transition matrix A (rows Begin/Match/Ins/Del/End -> columns Begin/Match/Ins/Del/End):\n";
    for (auto &row : a)
    {
        for (auto &v : row)
            cout << v << "\t";
        cout << "\n";
    }

    vector<string> states = {"Begin", "Match", "Insertion", "Deletion", "End"};
    for (int i = 1; i <= 3; ++i)
    {
        cout << "\nEmissions for state " << states[i] << ":\n";
        for (auto &kv : e[i])
            cout << kv.first << " : " << kv.second << "\n";
    }

    cout << "\nTotal alignments read: " << alignments.size() << "\n";

    time_t end = time(nullptr);
    cout << "Time taken: " << difftime(end, start) << " seconds\n;" << endl;
    return 0;
}
