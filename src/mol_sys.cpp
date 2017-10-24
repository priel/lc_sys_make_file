#include "./../include/mol_sys.h"

Mol_Sys::Mol_Sys(vector<double> & sys_sizes, vector<Molecule> & mols, vector<double> temperature_range, BoundaryType bc, int range)
	:m_sys_sizes(sys_sizes), m_molecules(mols), m_temperature_range(temperature_range),m_bc(bc),m_range(range)
{
	m_model = new Model();
	m_grid = new Grid(sys_sizes,m_bc,m_range);
	m_file_writer = new File_Writer();
	for (unsigned int i = 0; i < m_molecules.size(); i++)
	{
		vector<double> pot(i); //generate a vector in size of i (0 for the first and so on).
		m_pair_potentials.push_back(pot);
	}

	for (unsigned int i = 0; i < m_molecules.size(); i++) {
		for (unsigned int j = 0; j < i; j++)
			m_pair_potentials[i][j] = m_molecules[i].potential(&m_molecules[j], m_model);
	}

	for (unsigned int i = 0; i < m_molecules.size(); i++) {
		m_molecules[i].ID = i;
		m_grid->RegisterMol(&m_molecules[i]);
	}

#ifdef DEBUG
	vector<Molecule*> nbr_ptr;
	for (unsigned int i = 0; i < m_molecules.size(); i++) {
		nbr_ptr = m_grid->getNbr(m_molecules[i].m_location, 0).nbr_vec;
		cout << "I am:" << m_molecules[i].ID << endl;
		cout << "My nbr ar:" << endl;
		for (vector<Molecule *>::iterator it = nbr_ptr.begin(); it != nbr_ptr.end(); it++) {
			if ((*it) == NULL) {
				cout << "NULL  " ;
			}
			else {
				cout << (*it)->ID << "  ";
			}
		}
		cout << endl;
	}
#endif //DEBUG
}


Mol_Sys::~Mol_Sys()
{
	delete m_file_writer;
	delete m_model;
	delete m_grid;
}



double Mol_Sys::get_sys_potential()
{
	///update the system potential based on the pair potential of all the molecules.

	/// this function is less performance sensitive so we dont mind doing things twice:
	double potential = 0.0;
	vector<Molecule*> nbr_vec;
	for (unsigned int i = 0; i < m_molecules.size(); i++)
	{
		nbr_vec = m_grid->getNbr(m_molecules[i].m_location, 0).nbr_vec;
		potential += get_all_pair_potential_of_index(i, nbr_vec);
	}
	///since we calculated the pair potential i,j for each pair twice, once for i and once for j we need to divide by two
	return potential / 2;
}

void Mol_Sys::start_cooling()
{
#ifdef SHOW_TEMP_TIMMING
	clock_t prev, curr;
	double duration;
	curr = clock();
#endif // SHOW_TEMP_TIMMING


	double sys_potential = get_sys_potential();
	m_file_writer->make_model_directory();
	m_file_writer->write_state2xyz(m_molecules, m_temperature_range[0], sys_potential);

	/// in future will use some module how to cool the system.
	/// currently will just perform x monte carlos for each temperature from the array.
	for (m_current_index_temp = 0; m_current_index_temp < m_temperature_range.size(); m_current_index_temp++)
	{

		/// need to add print of the system here to xyz.
		/// first implement just a simple print
		monte_carlo();
		sys_potential = get_sys_potential();
		m_file_writer->write_state2xyz(m_molecules, m_temperature_range[m_current_index_temp], sys_potential);

#ifdef SHOW_TEMP_TIMMING
		prev = curr;
		curr = clock();
		duration = (curr - prev) / (double)CLOCKS_PER_SEC;
		cout << "temperature index " << m_current_index_temp << " took " << duration << " secs. potential = " << sys_potential << endl;
#endif // SHOW_TEMP_TIMMING
	}
	m_file_writer->write_list_file();
}

double Mol_Sys::get_all_pair_potential_of_index(unsigned int index, vector<Molecule*> nbr_vec)
{
	///sum the pair potentials:
	double potential = 0.0;
	int row;
	int col;
	int nbr_ID;
	for (vector<Molecule*>::iterator it = nbr_vec.begin(); it != nbr_vec.end(); it++) {
		nbr_ID = (*it)->ID;
		if (nbr_ID == index)
			continue;
		row = (index > nbr_ID) ? index : nbr_ID;
		col = (index < nbr_ID) ? index : nbr_ID;
		potential += m_pair_potentials[row][col];
	}

	return potential;
}

void Mol_Sys::update_sys(Molecule &mol_chosen, vector<Molecule*> nbr_vec, vector<double> potential)
{
	///inputs:
	/// Molecule mol_chosen: molecule to update.
	/// int index: the index of the molecule to update.
	/// double* potential: an array of all the new pair potential (size of m_molecules_size)
	///         potential[index] is undefined and should not use!!
	///
	/// the function update the system with the new vars.

	//remove "old" molecule from grid
	int mol_ID = mol_chosen.ID;
	m_grid->RemoveMol(&m_molecules[mol_ID]);

	///update the molecule
	m_molecules[mol_ID] = mol_chosen;

	//register "new" molecule to grid
	m_grid->RegisterMol(&m_molecules[mol_ID]);

	///update the pair potential:
	int nbr_ID;
	int row;
	int col;
	int nbr_cnt = 0;
	for (vector<Molecule*>::iterator it = nbr_vec.begin(); it != nbr_vec.end(); it++) {
		if (*it == NULL) {
			continue;
		}
		else {
			nbr_ID = (*it)->ID;
			if (nbr_ID == mol_ID) {
				continue;
			}
			row = (mol_ID > nbr_ID) ? mol_ID : nbr_ID;
			col = (mol_ID < nbr_ID) ? mol_ID : nbr_ID;
			m_pair_potentials[row][col] = potential[nbr_cnt];
			nbr_cnt++;
		}

	}

}

void Mol_Sys::monte_carlo()
{
	/// for each step:
	///	choose molecule randomly
	///	change location and spin with gauss distribution
	///	calculate dE
	///	if dE<0 take the step
	///	if not take the step with probability of e^-dE/Kb*T */

	int num_mol_chosen;
	vector<Molecule*> curr_nbr;
	vector<Molecule*> temp_nbr;
	vector<vector<int>> temp_nbr_shift;
	int grid_cell_count;
	double prob, dE, current_total_pot, suggested_location, suggested_spin, spin_norm, temp_total_pot;
	vector<double> potential;
	Molecule mol_chosen;

	///initiate the random generators:
	srand((unsigned int)time(0));
	std::default_random_engine loc_gen((unsigned int)time(0));
	std::normal_distribution<double> loc_dist(0.0, STD_LOCATION);

	std::default_random_engine spin_gen((unsigned int)time(0));
	std::normal_distribution<double> spin_dist(0.0, STD_SPIN);

	std::random_device  rand_dev;
	std::mt19937 generator(rand_dev());

	std::uniform_real_distribution<> distr_double(0, 1);
	std::uniform_int_distribution<int> distr_int(0, m_molecules.size()-1);

	/*std::random_device;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0, 1); */

	for (int i = 0; i < NUMBER_OF_STEPS; i++)
	{
		spin_norm = 0;

		///choose molecule:
		num_mol_chosen = distr_int(generator);

		///change location around gauss dist:
		mol_chosen = m_molecules[num_mol_chosen];
#ifdef DONT_MOVE_COLS
		if (mol_chosen.m_mol_type == col)
			continue;
#endif // DONT_MOVE_COLS



		for (unsigned int j = 0; j < mol_chosen.m_location.size(); j++)
		{
#ifdef DEBUG
			unsigned int counter = 0;
#endif //DEBUG
			///it's actually multivariate normal distribution where E=loc, std=std given, and no correlation between the axis.
			if (m_bc == Box) {
				do
				{
#ifdef DEBUG
					counter++;
					if (counter > 500)
					{
						cout << "you probably gave location for molecule such that it excceed dimension " << j << "\n";
						exit(EXIT_FAILURE);
					}

#endif //DEBUG
					suggested_location = mol_chosen.m_location[j] + loc_dist(loc_gen);
				} while ((suggested_location > m_sys_sizes[j]) || (suggested_location < 0));
				mol_chosen.m_location[j] = suggested_location;
			}
			else if (m_bc == Periodic) {
				suggested_location = mol_chosen.m_location[j] + loc_dist(loc_gen);
				mol_chosen.m_location[j] = Grid::mod(suggested_location, m_sys_sizes[j]);
			}
		}

		for (unsigned int j = 0; j < mol_chosen.m_spin.size(); j++)
		{
			//suggest the spin and normalize it.
			suggested_spin = mol_chosen.m_spin[j] + spin_dist(spin_gen);
			mol_chosen.m_spin[j] = suggested_spin;
			spin_norm += (suggested_spin*suggested_spin);
		}
		spin_norm = sqrt(spin_norm);
		for (unsigned int j = 0; j < mol_chosen.m_spin.size(); j++)
			mol_chosen.m_spin[j] /= spin_norm;


		/// we now have the location vector and the spin vector suggested, now we have to calculate dE for them
		///since all changed is this 1 molecule we will:
		/// calculate row of for the potential done by this molecule
		temp_nbr = m_grid->getNbr(mol_chosen.m_location, 1).nbr_vec;
		temp_nbr_shift = m_grid->getNbr(mol_chosen.m_location, 1).shift;
		temp_total_pot = 0;
		grid_cell_count = 0;
		potential.resize(0);
		for (vector<Molecule*>::iterator it = temp_nbr.begin(); it != temp_nbr.end(); it++)
		{
			//NULL value separate between nbr's in the same grid cell
			//untill we reach NULL all nbr's are in the same grid cell and the same shift is relavant to all.
			if ((*it) == NULL) { 
				grid_cell_count++;
				continue;
			}
			if ((*it)->ID == num_mol_chosen) continue;
			potential.push_back(mol_chosen.potential((*it), m_model, temp_nbr_shift[grid_cell_count]));
			temp_total_pot += potential.back();
		}
		curr_nbr = m_grid->getNbr(m_molecules[num_mol_chosen].m_location, 0).nbr_vec;
		current_total_pot = get_all_pair_potential_of_index(num_mol_chosen, curr_nbr);
		dE = current_total_pot - temp_total_pot;
		if (temp_total_pot <= current_total_pot)
		{
			update_sys(mol_chosen, temp_nbr , potential);
		}
		else
		{
			prob = distr_double(generator);

			if (prob < exp(dE / (m_temperature_range[m_current_index_temp] * K_B)))
			{
				update_sys(mol_chosen, temp_nbr , potential);
			}
		}

		//no need to delete mol_chosen since it wan't created by new it limited to this scope.
	}
}

