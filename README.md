EEL4837 Programming for Electrical Engineers II – Spring 2026
Excursion 1: Circuit Analysis Tool

Team:
- Jhon (Input & Data Preparation)
- Riyaz (Matrix Construction / MNA setup)
- Christian Vincente (Solver + Output + Integration)

How to Compile:
    g++ *.cpp -o ex1

How to Execute / Run:
    ./ex1
    (netlist.txt must be in the same directory as the executable)

Input (netlist.txt):
    Each line follows the exact format shown in class and the assignment:
        branch_label source_node destination_node value
    Example:
        V1 1 0 5
        R1 1 2 10
        R2 2 0 20
    • Voltage sources start with V
    • Resistors start with R
    • Node 0 is always ground (reference)
    • No empty lines, no extra spaces

Output (output.txt):
    A single line containing:
    e1 e2 … en   v1 v2 … vn   i1 i2 … in
    • e  = node potentials (voltages at non-ground nodes)
    • v  = voltage drop across each branch (source_node - destination_node)
    • i  = current through each branch
    Numbers are shown with up to 3 decimal places and no trailing zeros/decimal points (e.g., 5 not 5.000).
    The order of e, v, and i follows the order of nodes and branches in the netlist.

Bonus points:
    We are going for bonus points(sparse matrix data structure and the corresponding sparse matrix solver).

Additional Notes:
- Uses only standard C++ libraries (no external matrix libraries).
- Implements Modified Nodal Analysis (MNA) + Gaussian elimination with partial pivoting.
- Output formatting exactly matches the sample in the assignment PDF.
- Ground node is always the 0 V reference.
- Code is cleanly split by contributor for readability and maintainability.
