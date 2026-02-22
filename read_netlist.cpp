#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
using namespace std;

struct Branch
{
    string label;
    int src;
    int dst;
    double value;
    bool isVoltage;
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
    for (auto &b : branches)
        if (b.isVoltage)
            count++;
    return count;
}
// count resistors for MNA matrix size
int Resistors()
{
    int count = 0;
    for (auto &b : branches)
        if (!b.isVoltage)
            count++;
    return count;
}
//  debug functions that print the branches, node mapping, and matrices for verification
void printBranches()
{
    for (auto &b : branches)
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
    for (auto &p : node_to_index)
    {
        cout << p.first << " -> " << p.second << endl;
    }
}
// print a matrix for verification
void printMatrix(const vector<vector<double>> &A)
{
    for (auto &row : A)
    {
        for (double val : row)
            cout << val << " ";
        cout << endl;
    }
}

int main()
{

    readNetlist();
    removeGround();
    NodeIndex();

    printBranches();
    cout << endl;
    printNodeMap();
    cout << endl;

    std::cout << "Total branches: " << branches.size() << std::endl;
    std::cout << "Total nodes: " << node_to_index.size() << std::endl;
    std::cout << "Voltage sources: " << VoltageSources() << std::endl;

    return 0;
}