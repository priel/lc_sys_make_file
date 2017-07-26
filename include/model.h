#ifndef MODEL_H
#define MODEL_H

#include <cmath>
#include "./defined.h"

/// this is a container to hold all the variables of the system.
/// this should be hold by the system in order to calculate the potential.
/// the length of the molecules also hold here for performance.
/// another approach would be that every molecule with hold their structure.

class Model
{

    public:
        Model();
        virtual ~Model();

        double m_sigma_00;
        double m_sigma_11;
        double m_sigma_01;

        double m_chi_00;
        double m_chi_11;
        double m_chi_01;

        double m_alpha_00;
        double m_alpha_11;
        double m_alpha_01;


        double m_chi_tag_00;
        double m_chi_tag_11;
        double m_chi_tag_01;

        double m_alpha_tag_00;
        double m_alpha_tag_11;
        double m_alpha_tag_01;


    protected:

    private:
};

#endif // MODEL_H
