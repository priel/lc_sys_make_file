#include "grid.h"

#ifndef INIT_H
#define INIT_H

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <random>

#include "./defined.h"
#include "molecule.h"
#include "grid.h"

#define SIZE_OF_ARRAY(A) (sizeof(A) / sizeof(A[0]))

using namespace std;
/// this is a static class (constructor in private) contating all the function to initialize the system.

class Init
{
    public:
		static vector<int> get_mols_each_dir();

		static vector<double> get_sys_sizes();

		static vector< vector<int> > get_col_indices();

		static vector<double> get_temp_range();

		static vector<double> get_init_spin();

		static vector<Molecule> get_molecules(const vector<int> & molecules_in_each_directions, const vector<double> & initial_spin);

		static BoundaryType get_boundary_condition();

		static int get_range();

		static void add_randomization(vector<Molecule> & molecules, const vector< vector<int> > & colloid_molecules,
			const vector<int> & molecules_in_each_directions, const vector<double> & sys_sizes);

    protected:

    private:
		//static class, private constuctor.
		Init();
		virtual ~Init();
};

#endif // INIT_H
