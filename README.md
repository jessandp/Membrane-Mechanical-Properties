# Membrane-Mechanical-Properties
This repository contains all files pertaining to identifying the mechanical properties of a Molecular Dynamics simulated membrane.
    Membrane Builder - CHARMM-GUI
    MD Simulation - gromacs
    Visualization - VMD
    Analysis - R-Studio / Pycharm

    
# CHOLESTEROL'S POSITION
----- Objectives -----
1. Import data: Phosphorus (upper/lower leaflet) and Cholesterol (head/tail)
   Open a file and read contents line by line
2. Identify information needed from the file: System size, atom type, etc. 
   Search for strings in file
3. Convert data type: _____IDENTIFY DATA TYPE_____
4. Create a new file: "py_pdb.py"
   Readable / writable by python 
5. Identify the number of snapshots and separate them 
   Use the function: _____IDENTIFY FUNCTION_____
6. Identify the average of each coordinate of the snapshots
   - Phosphorus average
   - Cholesterol head average
   - Cholesterol tail average
7. Subtract the Cholesterol from the Phosphorus average 
   This value will give the average position of cholesterol relative to the position of phosphorus. Phosphorus is representing the outermost upper and lower leaflets.
----- Key Points -----
   + When working with file paths use: os.path
   + Use the function: readlines()
   + Use if statements to find strings within a file
   + To separate elements in a string use the function: split()
   + Choose data type best suited for running analysis
   

# BENDING MODULUS
----- Objectives -----
1. Import data: Phosphorus (upper/lower leaflet) and Cholesterol (head/tail)
   Open a file and read contents line by line
2. Identify information needed from the file: System size, atom type, etc. 
   Search for strings in file
3. Convert data type: _____IDENTIFY DATA TYPE_____
4. Create a new file: "py_pdb.py"
   Readable / writable by python 
5. Identify the number of snapshots and separate them 
   Use the function: _____IDENTIFY FUNCTION_____
6. _____CONTINUE OBJECTIVES_____

----- Key Points -----
   + When working with file paths use: os.path
   + Use the function: readlines()
   + Use if statements to find strings within a file
   + To separate elements in a string use the function: split()
   + Choose data type best suited for running analysis
    
References 
1. "File Parsing"
   http://education.molssi.org/python-scripting-biochemistry/chapters/File_Parsing.html
2. ""