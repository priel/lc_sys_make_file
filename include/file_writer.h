#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>


#if defined _MSC_VER
#include <direct.h>
#elif defined __GNUC__
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "./defined.h"
#include "molecule.h"

#define DEGREES_IN_CIRCLE 360

///this class in charge of all the files.
///
using namespace std;

class File_Writer
{
    public:
        File_Writer();
        virtual ~File_Writer();

		string m_output_dir;

		int m_xyz_file_counter;

		/**public functions**/
		///creting folders for the new model and copy the defined to it.
		void make_model_directory();

		///writing xyz file of the state of molecules (location and orientation).
		void write_state2xyz(const vector<Molecule> & molecules, double temperature, double potential);

		///writing a list of all files for Aviz
		void write_list_file();

    protected:

    private:
};

#endif // FILE_WRITER_H
