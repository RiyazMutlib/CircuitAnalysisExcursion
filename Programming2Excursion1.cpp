/*
EEL4837 – Excursion 1: Circuit Analysis Tool
Single-file solution, visually organized into 3 parts:
  Jhon – Input & Data Preparation
  Riyaz – Matrix Construction (MNA setup)
  Christian – Solver + Output + Integration

Matches assignment I/O and constraints.
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


// Jhon – INPUT & DATA PREPARATION

struct Branch
{
    string label;
    int src;
    int dst;
    double value;
    bool isVoltage;

    // Added because it is needed later for MNA: index of this voltage source among all voltage sources (0..M-1)
    int vsrcIndex;
};

vector<Branch> branches;
map<int, int> node_to_index;
set<int> unique_nodes;

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
    for (const auto &br : branches)
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
    for (const auto &br : branches)
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


