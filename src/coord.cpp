# include "coord.hpp"
# include "types.hpp"

int divmod (int a, int b) {return  (a%b+b) % b;} //modular division

coord::coord(const int& x, const int& y, const int& z, const int& l){
	this->x = x;
	this->y = y;
	this->z = z;
	this->l = l;
};
	
coord::coord(const int&c, const coord& S){
	int L = S.x;
	int M = S.y;
	int N = S.z;
	*this = {c % (M*L*N) % (M*L) % L, c % (M*L*N) % (M*L) / L, c % (M*L*N) / (M*L), c / (M*L*N)};
}
int coord::hash(const coord& S){
	int L = S.x;
	int M = S.y;
	int N = S.z;
	return l*L*M*N + z*L*M + y*L + x;
}
void coord::operator=(const coord& c) {
	x = c.x;
	y = c.y;
	z = c.z;
	l = c.l;
}

//dual lattice (cubes)
//copied from 211, verified on 4/28/22
int coord::getFaceQubits(const coord& S, const int& face, const int& direction, const int& k){
	vector<int> helper {0,1,2,3};
	if (direction == 1) {
		iter_swap(helper.begin() + 0, helper.begin() + 3);
		iter_swap(helper.begin() + 1, helper.begin() + 2);
	}
	if (face == 0) { // x-y plane
		if (z % 2 == 0) {
			iter_swap(helper.begin() + 0, helper.begin() + 1);
			iter_swap(helper.begin() + 2, helper.begin() + 3);
		}
		if (k == helper[0]) {
			return coord(x, divmod(y + 1, S.y), z, 0).hash(S);//9
		} else if (k == helper[1]) {
			return coord(x, y, z, 0).hash(S);//8
		} else if (k == helper[2]) {
			return coord(divmod(x + 1, S.x), y, z, 1).hash(S);//10
		} else{
			return coord(x, y, z, 1).hash(S);//11
		}
	} else if (face == 1) { //z-x plane
		if (y % 2 == 0) {
			iter_swap(helper.begin() + 0, helper.begin() + 1);
			iter_swap(helper.begin() + 2, helper.begin() + 3);
		}
		if (k == helper[0]) {
			return coord(divmod(x + 1, S.x), y, z, 2).hash(S);//4
		} else if (k == helper[1]) {
			return coord(x, y, z, 2).hash(S);//6
		} else if (k == helper[2]) {
			return coord(x, y, z, 0).hash(S);//8
		} else{
			return coord(x, y, divmod(z+1, S.z), 0).hash(S);//0
		}
	} else { //z-y  plane
		if (x % 2 == 0) {
			iter_swap(helper.begin() + 0, helper.begin() + 1);
			iter_swap(helper.begin() + 2, helper.begin() + 3);
		}
		if (k == helper[0]) {
			return coord(x, y, z, 1).hash(S);//11
		} else if (k == helper[1]) {
			return coord(x, y, divmod(z+1, S.z), 1).hash(S);//3
		} else if (k == helper[2]) {
			return coord(x, y, z, 2).hash(S);//6
		} else{
			return coord(x, divmod(y + 1, S.x), z, 2).hash(S);//7
		}
	}
}

//dual lattice (cubes)
int coord::getFaceQubits(const coord& S, const int& k){ //new version
	switch(k) {
		case 10:
			return coord(divmod(x + 1, S.x), y, z, 1).hash(S); //10
		case 8:
			return coord(x, y, z, 0).hash(S);//8
		case 9:
			return coord(x, divmod(y + 1, S.y), z, 0).hash(S); //9
		case 0:
			return coord(x, y, divmod(z+1, S.z), 0).hash(S);//0
		case 6:
			return coord(x, y, z, 2).hash(S);//6		
		case 4:
			return coord(divmod(x + 1, S.x), y, z, 2).hash(S);//4
		case 7:
			return coord(x, divmod(y + 1, S.x), z, 2).hash(S);//７
		case 3:
			return coord(x, y, divmod(z+1, S.z), 1).hash(S);//３
		case 11:
			return coord(x, y, z, 1).hash(S);//11
	}
	return 0;
}

//print coord
ostream& operator<<(ostream& os, const coord& c) {
	os << "(" << c.x << "," << c.y << "," << c.z << "," << c.l << ")";
	return os;
}

subcoord::subcoord(const int& x, const int& y, const int& l){
	this->x = x;
	this->y = y;
	this->l = l;
}

subcoord::subcoord(const int& c, const subcoord& S){
	int L = S.x;
	int M = S.y;
	*this = {c % (M*L) % L, c % (M*L) / L, c % (M*L)};
}

int subcoord::hash(const subcoord& S){
	int L = S.x;
	int M = S.y;
	return l*L*M + y*L + x;
}

void subcoord::operator=(const subcoord& c) {
	x = c.x;
	y = c.y;
	l = c.l;
}