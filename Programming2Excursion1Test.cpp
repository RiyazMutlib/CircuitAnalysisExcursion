/*
EEL4837 - Excursion 1: Circuit Analysis Tool
Single-file solution, visually organized into 3 parts:
  Jhon - Input & Data Preparation
  Riyaz - Matrix Construction (MNA setup)
  Christian - Solver + Output + Integration

Matches assignment I/O and constraints (No extra credit yet).
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
using namespace std;


// Jhon - INPUT & DATA PREPARATION

struct Branch
{
    string label;
    int src;
    int dst;
    double value;
    bool isVoltage;
    int vsrcIndex; // Added because it is needed later for MNA: index of this voltage source among all voltage sources (0..M-1)
};
vector<Branch> branches;
map<int, int> node_to_index;
set<int> unique_nodes;
struct SparseMatrix
 {
    int rows, cols;
    vector<vector<pair<int, double>>> data; // row -> list of (col, value)

    SparseMatrix(int r, int c) : rows(r), cols(c), data(r) {}

    void set(int row, int col, double val)
    {
        if (abs(val) < 1e-12)
            return; // treat as zero, do not store
        for (auto& entry : data[row])
        {
            if (entry.first == col)
            {
                entry.second = val;
                return;
            }
        }
        data[row].emplace_back(col, val); // add new entry
    }
    double get(int row, int col)const
    {
        for (const auto& entry : data[row])
        {
            if (entry.first == col)
                return entry.second;
        }
        return 0.0; // default zero
    }
    // convert from dense matrix
    void fromDense(const vector<vector<double>>& A)
    {
		// determine size
        int rows = (int)A.size();
		int cols = (int)A[0].size();
		// clear existing data
        data.clear();
		data.resize(rows);
        for( int i = 0; i< rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
				if (abs(A[i][j]) > 1e-12)   
                    data[i].emplace_back(j, A[i][j]);
            }
        }
    }
};

void readNetlist()
{
    // read file
    ifstream fin("netlist.txt");
    if (!fin.is_open())
    {
        cout << "Error! Could not open netlist.txt\n";
        return;
    }
    Branch b;
    while (fin >> b.label >> b.src >> b.dst >> b.value)
    {
        // classify components: voltage and resistor
        if (b.label[0] == 'V')
            b.isVoltage = true;
        else
            b.isVoltage = false;

        branches.push_back(b);
        // store nodes
        unique_nodes.insert(b.src);
        unique_nodes.insert(b.dst);
    }
    fin.close();
}
// remove ground node (node 0) if it exists, as it will be the reference node in MNA
void removeGround()
{
    if (unique_nodes.count(0))
        unique_nodes.erase(0);
}
// assign an index to each unique node for matrix construction
void NodeIndex()
{
    int index = 0;
    for (auto node : unique_nodes)
    {
        node_to_index[node] = index;
        index++;
    }
}
// count voltage sources and resistors for MNA matrix size
int VoltageSources()
{
    int count = 0;
    for (auto& b : branches)
        if (b.isVoltage)
            count++;
    return count;
}
// count resistors for MNA matrix size
int Resistors()
{
    int count = 0;
    for (auto& b : branches)
        if (!b.isVoltage)
            count++;
    return count;
}
//  debug functions that print the branches, node mapping, and matrices for verification
void printBranches()
{
    for (auto& b : branches)
    {
        cout << b.label << " "
            << b.src << " "
            << b.dst << " "
            << b.value << " "
            << b.isVoltage << endl;
    }
}
// print node to index mapping for verification
void printNodeMap()
{
    for (auto& p : node_to_index)
    {
        cout << p.first << " -> " << p.second << endl;
    }
}
// print a matrix for verification
void printMatrix(const vector<vector<double>>& A)
{
    for (auto& row : A)
    {
        for (double val : row)
            cout << val << " ";
        cout << endl;
    }
}


// Riyaz - MATRIX CONSTRUCTION (MNA setup)

void buildMNASystem(vector<vector<double>>& A, vector<double>& rhs, int& N, int& M, int& B)
{
    // Determine system size:
    // N = number of non-ground nodes (unknown node potentials)
    // M = number of voltage sources (unknown source currents)
    // B = number of branches (for later output)
    N = (int)node_to_index.size();
    M = VoltageSources();
    B = (int)branches.size();

    // Assign an index to each voltage source so we know where its current lives in the solution
    int vsCounter = 0;
    for (auto& br : branches)
    {
        if (br.isVoltage)
            br.vsrcIndex = vsCounter++;
        else
            br.vsrcIndex = -1;
    }

    int dim = N + M; // unknowns: [e1..eN, iV1..iVM]

    // Create MNA matrix A and RHS rhs
    A.assign(dim, vector<double>(dim, 0.0));
    rhs.assign(dim, 0.0);

    // Stamp resistors:
    // For resistor between nodes a and b with conductance g = 1/R
    // add g to diagonals, -g to off-diagonals (if both not ground)
    for (const auto& br : branches)
    {
        if (br.isVoltage)
            continue;

        double R = br.value;
        double g = 1.0 / R;

        int i = (br.src == 0) ? -1 : node_to_index[br.src];
        int j = (br.dst == 0) ? -1 : node_to_index[br.dst];

        if (i != -1)
            A[i][i] += g;
        if (j != -1)
            A[j][j] += g;
        if (i != -1 && j != -1)
        {
            A[i][j] -= g;
            A[j][i] -= g;
        }
    }

    // Stamp voltage sources:
    // For voltage source k between src->dst with value Vs:
    // KCL coupling: node rows connect to column (N+k)
    // Constraint row: row (N+k) enforces V(src) - V(dst) = Vs
    for (const auto& br : branches)
    {
        if (!br.isVoltage)
            continue;

        int k = br.vsrcIndex; // 0..M-1
        int rowcol = N + k;   // row/col index for this voltage source current unknown

        int i = (br.src == 0) ? -1 : node_to_index[br.src];
        int j = (br.dst == 0) ? -1 : node_to_index[br.dst];

        // KCL coupling (B block)
        if (i != -1)
            A[i][rowcol] += 1.0;
        if (j != -1)
            A[j][rowcol] -= 1.0;

        // Constraint equation (C block)
        if (i != -1)
            A[rowcol][i] += 1.0;
        if (j != -1)
            A[rowcol][j] -= 1.0;

        rhs[rowcol] = br.value;
    }
}


// ============================================================
// SOLVER (Gaussian Elimination)
// ============================================================

vector<double> gaussianSolve(vector<vector<double>> A, vector<double> b)
{
    int n = (int)A.size();

    for (int col = 0; col < n; col++)
    {
        // partial pivot
        int pivot = col;
        double best = fabs(A[col][col]);
        for (int r = col + 1; r < n; r++)
        {
            double v = fabs(A[r][col]);
            if (v > best)
            {
                best = v;
                pivot = r;
            }
        }

        if (best < 1e-15)
        {
            // inputs should be legal per spec, but guard against singular matrix
            return vector<double>(n, 0.0);
        }

        if (pivot != col)
        {
            swap(A[pivot], A[col]);
            swap(b[pivot], b[col]);
        }

        // eliminate
        for (int r = col + 1; r < n; r++)
        {
            double factor = A[r][col] / A[col][col];
            if (fabs(factor) < 1e-18)
                continue;

            for (int c = col; c < n; c++)
                A[r][c] -= factor * A[col][c];

            b[r] -= factor * b[col];
        }
    }

    // back substitution
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; j++)
            sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }
    return x;
}
// ===========================================================
// SOLVER (Sparse Gaussian Elimination)
// ===========================================================
vector<double> solveSparseSystem(SparseMatrix A, vector<double> b)
{
    int n =(int)A.rows;

    for (int i = 0; i < n; i++)
    {
        double pivot = A.get(i, i);

        for (int k = i + 1; k < n; k++)
        {
            double factor = A.get(k, i) / pivot;

            if (abs(factor) < 1e-12)
                continue;

            for (auto& col : A.data[i])
            {
                int j = col.first;
                double val = col.second;

                double newVal = A.get(k, j) - factor * val;
                A.set(k, j, newVal);
            }

            b[k] -= factor * b[i];
        }
    }

    vector<double> x(n);

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = b[i];

        for (auto& col : A.data[i])
        {
            int j = col.first;
            if (j > i)
                sum -= col.second * x[j];
        }

        x[i] = sum / A.get(i, i);
    }

    return x;
}
// helper to get node potential including ground
double nodeVoltage(int nodeLabel, const vector<double>& e)
{
    if (nodeLabel == 0)
        return 0.0;

    auto it = node_to_index.find(nodeLabel);
    if (it == node_to_index.end())
        return 0.0;

    return e[it->second];
}


// ============================================================
// OUTPUT & INTEGRATION
// ============================================================

int main()
{
    // ----  read + prep ----
    readNetlist();
    removeGround();
    NodeIndex();

    // Optional debugging
    // printBranches();
    // cout << endl;
    // printNodeMap();
    // cout << endl;

    // ----  build matrix system ----
    vector<vector<double>> A;
    vector<double> rhs;
    int N = 0, M = 0, B = 0;
    buildMNASystem(A, rhs, N, M, B);

    // ---- solve ----
        vector<double> x;
    if ((int)A.size() >= 4) {
        cout << "Using Sparse Matrix solver \n";
            SparseMatrix S((int)A.size(), (int)A[0].size());
            S.fromDense(A);
        x = solveSparseSystem(S, rhs);
    }
    else {
        cout << "Using Dense Matrix solver \n";
            x = gaussianSolve(A, rhs);
    }

    // Extract node potentials e1..eN
    vector<double> e(N, 0.0);
    for (int i = 0; i < N; i++)
        e[i] = x[i];

    // Extract voltage source currents iV1..iVM
    vector<double> iV(M, 0.0);
    for (int k = 0; k < M; k++)
        iV[k] = x[N + k];

    // Compute branch voltages and currents in netlist order
    vector<double> vBranch(B, 0.0);
    vector<double> iBranch(B, 0.0);

    for (int idx = 0; idx < B; idx++)
    {
        const auto& br = branches[idx];

        // branch voltage drop: v = V(src) - V(dst)
        double v = nodeVoltage(br.src, e) - nodeVoltage(br.dst, e);
        vBranch[idx] = v;

        // branch current:
        // - resistor: i = v / R
        // - voltage source: i is the solved unknown for that source
        if (br.isVoltage)
            iBranch[idx] = iV[br.vsrcIndex];
        else
            iBranch[idx] = v / br.value;
    }

    // Write output.txt exactly as required: e's then v's then i's, single row, 3 decimals :contentReference[oaicite:1]{index=1}
    ofstream fout("output.txt");
    if (!fout.is_open())
    {
        cerr << "Error: could not create output.txt\n";
        return 1;
    }

    fout << fixed << setprecision(3);

    bool first = true;
    auto writeVal = [&](double val)
        {
            if (!first)
                fout << ' ';
            first = false;

            // avoid -0.000
            if (fabs(val) < 0.0005)
                val = 0.0;

            fout << val;
        };

    // e1 e2 ... eN (ordered by node_to_index creation, which follows set ordering)
    for (int i = 0; i < N; i++)
        writeVal(e[i]);

    // v1 v2 ... vB (netlist order)
    for (int i = 0; i < B; i++)
        writeVal(vBranch[i]);

    // i1 i2 ... iB (netlist order)
    for (int i = 0; i < B; i++)
        writeVal(iBranch[i]);

    fout << "\n";
    fout.close();

    return 0;
}
