#include "./../include/model.h"



Model::Model()
{
    m_sigma_00 = D_0;
    m_sigma_11 = D_1;
    m_sigma_01 = sqrt(pow(D_0, 2) + pow(D_1, 2) ) / 2.0;

    m_chi_00 = (pow(L_0/D_0, 2) - 1) / (pow(L_0/D_0, 2) + 1);
    m_chi_11 = (pow(L_1/D_1, 2) - 1) / (pow(L_1/D_1, 2) + 1);
    m_chi_01 = sqrt( (pow(L_0, 2) - pow(D_0, 2)) * (pow(L_1, 2) - pow(D_1, 2) )
               /(pow(L_1, 2) + pow(D_0, 2) ) / (pow(L_0, 2) + pow(D_1, 2)) );

    m_alpha_00 = 1.0;
    m_alpha_11 = 1.0;
    m_alpha_01 = pow(((pow(L_0, 2) - pow(D_0, 2) ) * (pow(L_1, 2) + pow(D_0, 2) )
               /  (pow(L_1, 2) - pow(D_1, 2) ) / (pow(L_0, 2) + pow(D_1, 2) ) ) ,0.25); // i=0,j=1

    m_chi_tag_00 = (1 - pow(EE_DIV_ES_0 ,( 1 / MIU )) / (1 + pow(EE_DIV_ES_0,(1/ MIU)) ) );
    m_chi_tag_11 = (1 - pow(EE_DIV_ES_1 ,( 1 / MIU )) / (1 + pow(EE_DIV_ES_1,(1/ MIU)) ) );
    m_chi_tag_01 = (1 - pow(EE_DIV_ES_0 * EE_DIV_ES_1, (1 / (2 * MIU)) ) / (1 + pow(EE_DIV_ES_0 * EE_DIV_ES_1, 1 / (2 * MIU)) ) );

    m_alpha_tag_00 = sqrt(1 / (1 + pow(EE_DIV_ES_0, (1/MIU)) ) );
    m_alpha_tag_11 = sqrt(1 / (1 + pow(EE_DIV_ES_1, (1/MIU)) ) );
    m_alpha_tag_01 = sqrt(1 / (1 + pow(EE_DIV_ES_0 * EE_DIV_ES_1, 1/(2 * MIU)) ) );

}

Model::~Model()
{
    //dtor
}
