/*
EEL4837 – Excursion 1: Circuit Analysis Tool
Single-file solution, visually organized into 3 parts:
  Jhon – Input & Data Preparation
  Riyaz – Matrix Construction (MNA setup)
  Christian – Solver + Output + Integration

Matches assignment I/O and constraints (Extra Credit Include: Sparse Matrix Implementation
 and sparse matrix solver).
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
#include <cstdio>
using namespace std;


// Jhon – INPUT & DATA PREPARATION

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


// Riyaz – MATRIX CONSTRUCTION (MNA setup)

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


// Christian – SOLVER + OUTPUT + INTEGRATION

// Gaussian elimination with partial pivoting to solve A x = rhs
vector<double> solveSystem(vector<vector<double>> A, vector<double> b)
{
    int n = A.size();
    // Augment matrix with b
    for (int i = 0; i < n; i++) {
        A[i].push_back(b[i]);
    }

    // Forward elimination with partial pivoting
    for (int i = 0; i < n; i++) {
        // Find pivot
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);

        if (abs(A[i][i]) < 1e-12) {
            cout << "Warning: Matrix may be singular!" << endl;
        }

        for (int k = i + 1; k < n; k++) {
            double c = -A[k][i] / A[i][i];
            for (int j = i; j <= n; j++) {
                if (i == j)
                    A[k][j] = 0.0;
                else
                    A[k][j] += c * A[i][j];
            }
        }
    }

    // Back substitution
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = A[i][n];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}
// Sparse matrix gaussinan solver
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
// This matches the exact sample output shown in the assignment PDF
string formatNumber(double val)
{
    char buf[64];
    sprintf(buf, "%.3f", val);
    string s = buf;
    size_t dot = s.find('.');
    if (dot != string::npos)
    {
        s.erase(s.find_last_not_of('0') + 1, string::npos);
        if (!s.empty() && s.back() == '.')
            s.pop_back();
    }
    return s;
}

// Write the final output to output.txt in the exact required format
void writeOutput(const vector<double>& nodePotentials, 
                 const vector<double>& branchVoltages, 
                 const vector<double>& branchCurrents)
{
    ofstream fout("output.txt");
    if (!fout.is_open()) {
        cout << "Error! Could not create output.txt\n";
        return;
    }

    // e1 e2 ... en
    for (double e : nodePotentials) {
        fout << formatNumber(e) << " ";
    }
    // v1 v2 ... vn
    for (double v : branchVoltages) {
        fout << formatNumber(v) << " ";
    }
    // i1 i2 ... in
    for (double i : branchCurrents) {
        fout << formatNumber(i) << " ";
    }
    fout << endl;
    fout.close();
}

int main()
{
    readNetlist();
    removeGround();
    NodeIndex();

    vector<vector<double>> A;
    vector<double> rhs;
    int N, M, B;

    buildMNASystem(A, rhs, N, M, B);

    // Solve the system
        vector<double> x;
    if ((int)A.size() >= 15) {
        cout << "Using Sparse Matrix solver \n";
            SparseMatrix S((int)A.size(), (int)A[0].size());
            S.fromDense(A);
        x = solveSparseSystem(S, rhs);
    }
    else {
        cout << "Using Dense Matrix solver \n";
            x = solveSystem(A, rhs);
    }

    // Build potentials map (ground is 0.0)
    map<int, double> potentials;
    potentials[0] = 0.0;
    for (auto& p : node_to_index) {
        potentials[p.first] = x[p.second];
    }

    // Compute branch voltages and currents
    vector<double> branchVoltages;
    vector<double> branchCurrents;

    for (const auto& br : branches) {
        double v = potentials[br.src] - potentials[br.dst];
        double i_val;
        if (br.isVoltage) {
            i_val = x[N + br.vsrcIndex];   // uses MNA sign convention (matches sample)
        } else {
            i_val = (abs(br.value) > 1e-12) ? v / br.value : 0.0;
        }
        branchVoltages.push_back(v);
        branchCurrents.push_back(i_val);
    }

    // Collect node potentials in matrix order
    vector<double> nodePotentials;
    for (auto& p : node_to_index) {
        nodePotentials.push_back(potentials[p.first]);
    }

    // Write output file
    writeOutput(nodePotentials, branchVoltages, branchCurrents);

    cout << "Program completed successfully. Check output.txt" << endl;
    // printBranches();     // uncomment for debugging
    // printNodeMap();

    return 0;
}
