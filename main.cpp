#include "./include/mol_sys.h"
#include "./include/molecule.h"
#include "./include/defined.h"
#include "./include/init.h"

using namespace std;

int main(int argc, char* argv[])
{
	//initializing all the vectors:
	vector<double> sys_sizes = Init::get_sys_sizes();
	vector<int> molecules_in_each_directions = Init::get_mols_each_dir();
	vector<vector<int> > colloid_molecules = Init::get_col_indices();
	vector<double> temperature_range = Init::get_temp_range();
	vector<double> initial_spin = Init::get_init_spin();
	vector<Molecule> molecules = Init::get_molecules(molecules_in_each_directions, initial_spin);

	//add randomization to the initial location and orientation of the vectors:
	Init::add_randomization(molecules, colloid_molecules, molecules_in_each_directions, sys_sizes);
	cout << "finished initialization, start colloing:" << endl;	
//here we can call to some function that will modify the system if user want to.
	Mol_Sys * lc_system = new Mol_Sys(sys_sizes, molecules, temperature_range);
	lc_system->start_cooling();
	delete lc_system;

	cout << "finished successfully" << endl;
}

