#ifndef GRID_H
#define GRID_H

#include "molecule.h"
#include <vector>
using namespace std;

enum BoundaryType { Box, Periodic };

class Grid 
{
public:
	//the constructor get the system size vector and a range parameter to indicate the number of grid cells
	//to include in every direction, around the central cell of the molecule, during the pair potential calculation
	Grid(vector<double> & sys_sizes, BoundaryType bc, int range);
	~Grid();
	void RegisterMol(Molecule* mol_ptr);
	void RemoveMol(Molecule* mol_ptr);

	struct Nbr
	{
		vector<Molecule*> nbr_vec;
		vector<vector<int>> shift;
	};
	Nbr getNbr(vector<double> location, bool shift); //if shift==1 the function calculate the shift field

	static int mod(int a, int b); //TODO :move outside to an external class
	int grid_mol_num;
private:
	int grid_range;
	BoundaryType grid_bc;
	vector<int> grid_size;
	vector<vector<vector<Molecule*>>>    grid_2D; //2D: array containing vectors of pointers, to the molecules in each grid cell  
	vector<vector<vector<vector<Molecule*>>>>    grid_3D; //3D: array containing vectors of pointers, to the molecules in each grid cell  
	
	vector<int> getGridPoint(vector<double> loc);
};

#endif // GRID_H