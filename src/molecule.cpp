#include "./../include/molecule.h"
//oshri's commit from VS

Molecule::Molecule(std::vector<double> loc, std::vector<double> spin, Mol_Type mol_type)
                    : m_location(loc), m_spin(spin), m_mol_type(mol_type){}

Molecule::Molecule()
{
	m_mol_type = lc; //default molecule is lc
	m_location.resize(DIMENSIONS);
	m_spin.resize(DIMENSIONS);
}  

Molecule::~Molecule(){}

//TODO :use internal function in overloading
double Molecule::potential(const Molecule * mol, const Model * model, const vector<int>& shift)
{
  double potential;

    //this is the most performance sensitive. try to be as minimal as possible
    //try to do as much as possible in compile time.
    double norm_r, sigma0, chi, alpha, first, second, sigma, R;
    double epsilon, epsilon_ni, epsilon_miu, epsilon0, alpha_tag, chi_tag;
    double dot_spin1_nr, dot_spin2_nr, dot_spin1_spin2;
	double r0, r1;
#if DIMENSIONS == 3
	double r2;
#endif // DIMENSIONS

	if (shift.size() == 0) {
		r0 = mol->m_location[0] - m_location[0];
		r1 = mol->m_location[1] - m_location[1];
#if DIMENSIONS == 3
		r2 = mol->m_location[2] - m_location[2];
#endif // DIMENSIONS%
	}
	else {
		r0 = mol->m_location[0] - m_location[0] - shift[0];
		r1 = mol->m_location[1] - m_location[1] - shift[1];
#if DIMENSIONS == 3
		r2 = mol->m_location[2] - m_location[2] - shift[2];
#endif // DIMENSIONS%
	}

    #if DIMENSIONS == 2
        norm_r = sqrt(r0*r0 + r1*r1);
        dot_spin1_nr = (m_spin[0]*r0 + m_spin[1]*r1) / norm_r;
        dot_spin2_nr = (mol->m_spin[0]*r0 + mol->m_spin[1]*r1) / norm_r;
        dot_spin1_spin2 = (mol->m_spin[0]*m_spin[0] + mol->m_spin[1]*m_spin[1]);
    #elif DIMENSIONS ==3
		norm_r = sqrt(r0*r0 + r1*r1 + r2*r2);
        dot_spin1_nr = (m_spin[0]*r0 + m_spin[1]*r1 + m_spin[2]*r2) / norm_r;
        dot_spin2_nr = (mol->m_spin[0]*r0 + mol->m_spin[1]*r1 + mol->m_spin[2]*r2) / norm_r;
        dot_spin1_spin2 = (mol->m_spin[0]*m_spin[0] + mol->m_spin[1]*m_spin[1] + mol->m_spin[2]*m_spin[2]);
    #endif // DIMENSIONS

    if ((m_mol_type == 0)) // both are LC
    {
        if (mol->m_mol_type == 0)
        {
            sigma0 = model->m_sigma_00;
            chi = model->m_chi_00;
            alpha = model->m_alpha_00;

            epsilon0 = EPSILON0_0;
            alpha_tag = model->m_alpha_tag_00;
            chi_tag = model->m_chi_tag_00;
        }
        else
        {
            sigma0 = model->m_sigma_01;
            chi = model->m_chi_01;
            alpha = model->m_alpha_01; //in the original

            epsilon0 = sqrt(EPSILON0_0 * EPSILON0_1);
            alpha_tag = model->m_alpha_tag_01;
            chi_tag = model->m_chi_tag_01;
        }

    }
    else
    {
        if (mol->m_mol_type == 0)
        {
            sigma0 = model->m_sigma_01;
            chi = model->m_chi_01;
            alpha = 1 / model->m_alpha_01; //in the original

            epsilon0 = sqrt(EPSILON0_0 * EPSILON0_1);
            alpha_tag = model->m_alpha_tag_01;
            chi_tag = model->m_chi_tag_01;
        }
        else
        {
            sigma0 = model->m_sigma_11;
            chi = model->m_chi_11;
            alpha = model->m_alpha_11;

            epsilon0 = EPSILON0_1;
            alpha_tag = model->m_alpha_tag_11;
            chi_tag = model->m_chi_tag_11;
        }
    }

    first = (dot_spin1_nr * alpha + dot_spin2_nr / alpha);
    first *=first;
    first /= (1.0 + chi * dot_spin1_spin2);

    second = (dot_spin1_nr * alpha - dot_spin2_nr / alpha);
    second *= second;
    second /= (1.0 - chi * dot_spin1_spin2);
#ifdef DEBUG
	if ((1.0 - chi * (first + second) / 2.0) < 0)
	{
		cout << "chi has been negative again" << endl;
		exit(EXIT_FAILURE);
	}
#endif
    sigma = sigma0 / sqrt(1.0 - chi * (first + second) / 2.0);
    R = sigma0 / (norm_r - sigma + sigma0);

#ifdef DEBUG
	if ((1.0 - (chi * chi) * (dot_spin1_spin2 * dot_spin1_spin2)) < 0)
	{
		cout << "chi square with doted has been negative" << endl;
		exit(EXIT_FAILURE);
	}
#endif
    epsilon_ni = 1.0 / sqrt(1.0 - (chi * chi) * (dot_spin1_spin2 * dot_spin1_spin2));
    // these first and second have nothing to do with the above first and second.
    first = alpha_tag * dot_spin1_nr + (1/alpha_tag)*dot_spin2_nr; // I assume that there was a problem in the paper and it should be alpha and not alpha inpower of -1
    first *= first;
    first /= (1.0 + chi_tag * dot_spin1_spin2);

    second = alpha_tag * dot_spin1_nr - (1/alpha_tag)*dot_spin2_nr;
    second *= second;
    second /= (1.0 - chi_tag * dot_spin1_spin2);

    epsilon_miu = 1.0 - chi_tag * (first + second) / 2.0;

    epsilon = epsilon0 * pow(epsilon_ni, NI) * pow(epsilon_miu, MIU);

    potential = (4 * epsilon * (pow(R,12) - pow(R,6)));
    return potential;

}
