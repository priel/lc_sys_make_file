#include "./../include/grid.h"

Grid::Grid(vector<double> & sys_sizes, BoundaryType bc, int range)
	:_grid_bc(bc), _grid_range(range)
{
	m_grid_mol_num = 0;
	vector<int> grid_size(sys_sizes.begin(), sys_sizes.end());
	if (grid_size.size() == 2) {
		_grid_2D.resize(grid_size[0]);
		for (int i = 0; i < grid_size[0]; i++) {
			_grid_2D[i].resize(grid_size[1]);
		}
	}
	else if (grid_size.size() == 3) {
		_grid_3D.resize(grid_size[0]);
		for (int i = 0; i < grid_size[0]; i++) {
			_grid_3D[i].resize(grid_size[1]);
			for (int j = 0; j < grid_size[1]; j++) {
				_grid_3D[i][j].resize(grid_size[2]);
			}
		}
	}
	else {
		cout << "cannot initiallize grid. dimension must be 2 or 3" << endl;
		exit(EXIT_FAILURE);
	}

	Grid::_grid_size = grid_size;
}

Grid::~Grid() {
}


void Grid::RegisterMol(Molecule* mol_ptr) {
	vector<int> grid_pt = getGridPoint(mol_ptr->m_location);
	if (grid_pt.size() == 2) {
		_grid_2D[grid_pt[0]][grid_pt[1]].push_back(mol_ptr);
		m_grid_mol_num++;
	}
	else if (grid_pt.size() == 3) {
		_grid_3D[grid_pt[0]][grid_pt[1]][grid_pt[2]].push_back(mol_ptr);
		m_grid_mol_num++;
	}
	else {
		cout << "cannot initiallize grid. dimension must be 2 or 3" << endl;
		exit(EXIT_FAILURE);
	}
}

void Grid::RemoveMol(Molecule* mol_ptr) {
	int size1;
	int size2;
	vector<int> grid_pt = getGridPoint(mol_ptr->m_location);
	if (grid_pt.size() == 2) {
		for (vector<Molecule*>::iterator it = _grid_2D[grid_pt[0]][grid_pt[1]].begin();
			it != _grid_2D[grid_pt[0]][grid_pt[1]].end(); it++) {
			if ((*it)->ID == mol_ptr->ID) {
				size1 = _grid_2D[grid_pt[0]][grid_pt[1]].size();
				_grid_2D[grid_pt[0]][grid_pt[1]].erase(it);
				size2 = _grid_2D[grid_pt[0]][grid_pt[1]].size();
				if (size2 != size1 - 1) {
					cout << "ERROR";
				}
				m_grid_mol_num--;
				break;
			}
		}
	}
	else if (grid_pt.size() == 3) {
		for (vector<Molecule*>::iterator it = _grid_3D[grid_pt[0]][grid_pt[1]][grid_pt[2]].begin();
			it != _grid_3D[grid_pt[0]][grid_pt[1]][grid_pt[2]].end(); it++) {
			if ((*it)->ID == mol_ptr->ID) {
				_grid_3D[grid_pt[0]][grid_pt[1]][grid_pt[2]].erase(it);
				m_grid_mol_num--;
				break;
			}
		}
	}
	else {
		cout << "cannot initiallize grid. dimension must be 2 or 3" << endl;
		exit(EXIT_FAILURE);
	}
}

Grid::Nbr Grid::getNbr(vector<double> location, bool shift) {
	Nbr res;
	vector<int> grid_pt = getGridPoint(location);
	vector<int> low_idx;
	vector<int> high_idx;
	vector<int> tmp_v;
	int low;
	int high;

	for (int i = 0; i < grid_pt.size(); i++) {
		low = (grid_pt[i] - _grid_range);
		high = (grid_pt[i] + _grid_range);
		if (_grid_bc == Box) {
			low = (low < 0) ? 0 : low;
			high = (high > _grid_size[i] - 1) ? _grid_size[i] - 1 : high;
		}
		low_idx.push_back(low);
		high_idx.push_back(high);
	}
	
	if (grid_pt.size() == 2) {
		for (int i = low_idx[0]; i <= high_idx[0]; i++) {
			for (int j = low_idx[1]; j <= high_idx[1]; j++) {
				// concatenation of vectors in neighbor grid cells 
				res.nbr_vec.insert(res.nbr_vec.end(), 
							   _grid_2D[mod(i,_grid_size[0])][mod(j,_grid_size[1])].begin(),
							   _grid_2D[mod(i,_grid_size[0])][mod(j,_grid_size[1])].end());
				if (shift) {
					res.nbr_vec.push_back(NULL); // NULL ptr separate between sequences of nbr with same shift
					tmp_v.push_back(i - mod(i, _grid_size[0]));
					tmp_v.push_back(j - mod(j, _grid_size[1]));
					res.shift.push_back(tmp_v);
				}
			}
		}
	}
	else if (grid_pt.size() == 3) {
		for (int i = low_idx[0]; i <= high_idx[0]; i++) {
			for (int j = low_idx[1]; j <= high_idx[1]; j++) {
				for (int k = low_idx[2]; k <= high_idx[2]; k++) {
					// concatenation of vectors in neighbor grid cells 
					res.nbr_vec.insert(res.nbr_vec.end(),
								   _grid_3D[mod(i, _grid_size[0])][mod(j, _grid_size[1])][mod(k, _grid_size[2])].begin(),
								   _grid_3D[mod(i, _grid_size[0])][mod(j, _grid_size[1])][mod(k, _grid_size[2])].end());
					if (shift) {
						res.nbr_vec.push_back(NULL); // NULL ptr separate between sequences of nbr with same shift
						tmp_v.push_back(i - mod(i, _grid_size[0]));
						tmp_v.push_back(j - mod(j, _grid_size[1]));
						tmp_v.push_back(k - mod(k, _grid_size[2]));
						res.shift.push_back(tmp_v);
					}
				}
			}
		}
	}
	return res;
}

vector<int> Grid::getGridPoint(vector<double> loc) {
	for (vector<double>::iterator it = loc.begin(); it != loc.end(); it++) {
		*it = round(*it);
	}
	vector<int> grid_pt(loc.begin(), loc.end());
	return grid_pt;
}

//implemnting modulu operation = euclidean reminder
int Grid::mod(int a, int b)
{
	int r = a % b;
	return r >= 0 ? r : r + std::abs(b);
}