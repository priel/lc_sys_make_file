#include "./../include/file_writer.h"



File_Writer::File_Writer()
{
	m_xyz_file_counter = 0;
}

File_Writer::~File_Writer()
{
    //dtor
}

void File_Writer::make_model_directory()
/// create folder in runs dir with the model name.
{
	string model_name = MODEL_NAME;
	string original_dir = ".//runs//" + model_name;
	string dir_to_create = original_dir;
	int nError = 0, i = 0, nErrorLinux = 0;

	do
	{
#if defined _MSC_VER
		nError = _mkdir(dir_to_create.c_str());
#elif defined __GNUC__
		nError = mkdir(dir_to_create.c_str(), 0777);
#endif
#if defined DEBUG
		if (i>300)
		{
			cout << "couldnt craete file after 300 tries, this could be because there is no runs directory. exiting.." << endl;
			cout << "GNUC defined, error is " << nError << "runs erro r" << nErrorLinux << endl;
			exit(EXIT_FAILURE);
		}
#endif
		i++;
		dir_to_create = original_dir + "_" + to_string(i);
	} while (nError != 0);

	i--;
	dir_to_create = original_dir + "_" + to_string(i);

	//copy the model file to the the model directory
	std::ifstream  src(".//include//defined.h", std::ios::binary);
	string defined_file = dir_to_create + "//defined_file.h";
	std::ofstream  dst(defined_file, std::ios::binary);

	dst << src.rdbuf();

	//make output directory for all the xyz files.
	m_output_dir = dir_to_create + "//output";
#if defined _MSC_VER
	nError = _mkdir(m_output_dir.c_str());
#elif defined __GNUC__
	nError = mkdir(m_output_dir.c_str(), 0777);
#endif
	if (nError != 0)
	{
		cout << "couldnt create output directory, exiting..." << endl;
		exit(EXIT_FAILURE);
	}
}

void File_Writer::write_state2xyz(const vector<Molecule> & molecules, double temperature, double potential)
{
	//craeting file in format: lqs_sys_0002
	ofstream xyz_file;
	stringstream fn_suffix;
	string file_name = "lqs_sys_";
	fn_suffix << setfill('0') << setw(4) << m_xyz_file_counter;
	file_name = m_output_dir + "//" + file_name + fn_suffix.str() + ".xyz";
	xyz_file.open(file_name);

	int num_col_mol = 0;
	int num_of_molecules = molecules.size();

	for (int i = 0; i < num_of_molecules; i++)
		if (molecules[i].m_mol_type == col)
			num_col_mol++;

	//turns out that AVIZ can parse any format, so leaving default format as is.
	/*
	//setting format of double for writing:
	xyz_file.precision(6);
	xyz_file << fixed;
	*/

	int line_counter = num_of_molecules;
#if DIMENSIONS == 2
	line_counter += (DEGREES_IN_CIRCLE / DIF_ANGLES_COL_REPRESENTATION) * num_col_mol - 1;
#elif DIMENSIONS == 3
	line_counter += (DEGREES_IN_CIRCLE / DIF_ANGLES_COL_REPRESENTATION) * (DEGREES_IN_CIRCLE / DIF_ANGLES_COL_REPRESENTATION / 2) * num_col_mol - 1;
#endif
	xyz_file << line_counter << endl;
	xyz_file << "Liquid Crystals with Colloide in temperature=" << temperature << ", and potential=" << potential << endl;
	for (int i = 0; i < num_of_molecules; i++)
	{
		if (molecules[i].m_mol_type == col)
		{
			int lines_to_draw_x = 360 / DIF_ANGLES_COL_REPRESENTATION;
#if	DIMENSIONS == 2
			double phi;
			for (int j = 0; j < (DEGREES_IN_CIRCLE / DIF_ANGLES_COL_REPRESENTATION); j++)
			{
				phi = (2 * PI * DIF_ANGLES_COL_REPRESENTATION / DEGREES_IN_CIRCLE) * j;
				xyz_file << "d 0.0 " << molecules[i].m_location[0] << " " << molecules[i].m_location[1] << " 0.0 "
					<< cos(phi) << " " << sin(phi) << endl;
			}
#elif DIMENSIONS == 3
			double phi, theta;
			for (int j = 0; j < (DEGREES_IN_CIRCLE / DIF_ANGLES_COL_REPRESENTATION); j++)
			{
				phi = (2 * PI * DIF_ANGLES_COL_REPRESENTATION / DEGREES_IN_CIRCLE) * j;
				for (int k = 0; k < (DEGREES_IN_CIRCLE / DIF_ANGLES_COL_REPRESENTATION) / 2; k++)
				{
					theta = (2 * PI * DIF_ANGLES_COL_REPRESENTATION / DEGREES_IN_CIRCLE) * k;
					xyz_file << "d " << molecules[i].m_location[0] << " " << molecules[i].m_location[1] << " " << molecules[i].m_location[2]
						<< " " << sin(theta) * cos(phi) << " " << sin(theta) * sin(phi) << " " << cos(theta) << endl;
				}
			}
#endif
		}
		else
		{
#if DIMENSIONS == 2 //ommit x:
			xyz_file << "Sp 0.0 " << molecules[i].m_location[0] << " " << molecules[i].m_location[1] << " 0.0 "
				<< molecules[i].m_spin[0] << " " << molecules[i].m_spin[1] << endl;
#elif DIMENSIONS == 3
			xyz_file << "Sp " << molecules[i].m_location[0] << " " << molecules[i].m_location[1] << " " << molecules[i].m_location[2]
				<< " " << molecules[i].m_spin[0] << " " << molecules[i].m_spin[1] << " " << molecules[i].m_spin[2] << endl;
#endif
		}
	}
	xyz_file.close();
	m_xyz_file_counter++;
}

void File_Writer::write_list_file()
{
	ofstream list_file;
	string file_name = m_output_dir + "//xyz.list";
	list_file.open(file_name);
	stringstream xyz_suffix;
	string xyz_name = "lqs_sys_";
	string line;

	for (int i = 0; i < m_xyz_file_counter; i++)
	{
		xyz_suffix.clear();
		xyz_suffix.str(string());
		xyz_suffix << setfill('0') << setw(4) << i;
		line = xyz_name + xyz_suffix.str() + ".xyz";
		list_file << line << endl;
	}
	list_file.close();
}
