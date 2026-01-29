/*
 Petar 'PetarV' Velickovic
 Algorithm: Needleman-Wunsch
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>
#include <fstream>
#include <unistd.h>

#define MAX_N 3501

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

int n, m;
int match_score, mismatch_score, gap_score;
string A, B;
int dp[MAX_N][MAX_N];

/*
 Needleman-Wunsch algorithm for determining the optimal alignment between two strings
 assuming a given score for hits, gaps and mismatches.
 Complexity: O(n * m) time, O(n * m) memory
*/

inline int needleman_wunsch()
{
    for (int i = 0; i <= n; i++)
        dp[i][0] = dp[0][i] = i * gap_score;
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= m; j++)
        {
            int S = (A[i - 1] == B[j - 1]) ? match_score : mismatch_score;
            dp[i][j] = max(dp[i - 1][j - 1] + S, max(dp[i - 1][j] + gap_score, dp[i][j - 1] + gap_score));
        }
    }
    return dp[n][m];
}

inline pair<string, string> get_optimal_alignment()
{
    string retA, retB;
    stack<char> SA, SB;
    int ii = n, jj = m;
    while (ii != 0 || jj != 0)
    {
        if (ii == 0)
        {
            SA.push('-');
            SB.push(B[jj - 1]);
            jj--;
        }
        else if (jj == 0)
        {
            SA.push(A[ii - 1]);
            SB.push('-');
            ii--;
        }
        else
        {
            int S = (A[ii - 1] == B[jj - 1]) ? match_score : mismatch_score;
            if (dp[ii][jj] == dp[ii - 1][jj - 1] + S)
            {
                SA.push(A[ii - 1]);
                SB.push(B[jj - 1]);
                ii--;
                jj--;
            }
            else if (dp[ii - 1][jj] > dp[ii][jj - 1])
            {
                SA.push(A[ii - 1]);
                SB.push('-');
                ii--;
            }
            else
            {
                SA.push('-');
                SB.push(B[jj - 1]);
                jj--;
            }
        }
    }
    while (!SA.empty())
    {
        retA += SA.top();
        retB += SB.top();
        SA.pop();
        SB.pop();
    }
    return make_pair(retA, retB);
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
    if (argc < 5)
    {
        printf("Usage: %s <match_score> <mismatch_score> <gap_score> <input.bin>\n", argv[0]);
        return 1;
    }
    match_score = atoi(argv[1]);
    mismatch_score = atoi(argv[2]);
    gap_score = atoi(argv[3]);
    string input = argv[4];
    ifstream file(input, ios::binary);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << input << endl;
        exit(1);
    }
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
            cout << "Alignment between " << name1 << " and " << name2 << ":\n";
            A = seq1;
            B = seq2;
            A.erase(remove(A.begin(), A.end(), '-'), A.end());
            B.erase(remove(B.begin(), B.end(), '-'), B.end());
            // cout << "A: " << A << "\nB: " << B << "\n";
            n = A.length();
            m = B.length();
            // printf("n = %d, m = %d\n", n, m);
            printf("%d\n", needleman_wunsch());
            // pair<string, string> alignment = get_optimal_alignment();
            // printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());
            // sleep(1);
        }
    }
    return 0;
}