#ifndef MOLECULE_H
#define MOLECULE_H
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "./defined.h"
#include "./model.h"

using namespace std;

class Molecule
{
    public:

        Molecule();

        Molecule(std::vector<double> loc, std::vector<double> spin, Mol_Type mol_type);


        /** Default destructor
         the destructor will always delete the vectors */
        virtual ~Molecule();

        double potential(const Molecule* mol, const Model* model);

        std::vector<double> m_location;
        std::vector<double> m_spin;

        Mol_Type m_mol_type;


    protected:

    private:
};

#endif // MOLECULE_H
