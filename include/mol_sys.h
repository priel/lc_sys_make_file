#ifndef MOL_SYS_H
#define MOL_SYS_H

#include <ctime>
#include <cstdlib>
#include <random>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <string>


#include "./mol_sys.h"
#include "./defined.h"
#include "./molecule.h"
#include "./model.h"
#include "./file_writer.h"

using namespace std;

class Mol_Sys
{
    public:

        /// this is some kind of custom constructor where all the parameters are pre-defined.
        Mol_Sys(vector<double> & sys_sizes, vector<Molecule> & mols, vector<double> temperature_range);

        ///nothing has made with new, nothing to delete.
        ~Mol_Sys();

        vector<double> m_sys_sizes; ///array in the length of dimensions which determine the x,y,(z) of the system.

		vector<Molecule> m_molecules; /// ///pointer to the first molecule array.

		vector<double> m_temperature_range; ///hold the range of temperature we want to check (Starting from 0 to max_temp -1;
		unsigned int m_current_index_temp;

		Model* m_model; ///hold model calculated parameters

		File_Writer* m_file_writer; ///class in charge of files outputs

		/// have all the pairs of potential for example pair_potential[0][1] has the potential between molecule 0 and 1.
		///the last column has the sum potential of this
		vector< vector<double> > m_pair_potentials;

        /** example of potentials as a matrix for m_molecules.size()=4:
            0   1   2   3
        0  N/A N/A N/A N/A
        1  2.0 N/A N/A N/A
        2  1.0 2.1 N/A N/A
        3  2.2 0.8 0.5 N/A   //this is not really needed and just here for simplicity.
        */

        // We could have saved more memory by defining only the triangle, however it left that way for simplicity of the code.
		// if memory will become an issue we can change this.

        /** the public functions: */


        /// update the system potential based on the pair potential
        double get_sys_potential();

		/// in future will use some module how to cool the system.
		/// currently will just perform x monte carlos for each temperature from the array.
        void start_cooling();
        
		///get the potential of all the pairs with the index
        double get_all_pair_potential_of_index(unsigned int index);
		
		/// in charge of updating the system.
        void update_sys(Molecule &mol_chosen, unsigned int index, double* potential, double tot_pot_update);
        
		///doing NUMBER_OF_STEPS times monte carlo steps
        void monte_carlo();
        


    protected:

    private:
		Mol_Sys(); /// no default constructr should be in use, arguments must be provided.
};

#endif // MOL_SYS_H
