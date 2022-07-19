#include "vertex.hpp"

//################################## vertex ##############################

vertex::vertex(const int& c, const coord& S){
	this->c = c;
	coord C(c,S);
	partial ={
		coord(divmod(C.x-1, S.x), C.y, C.z, 0).hash(S), coord(C.x, C.y, C.z, 0).hash(S),
		coord(C.x, divmod(C.y-1, S.y), C.z, 1).hash(S), coord(C.x, C.y, C.z, 1).hash(S),
		coord(C.x, C.y, divmod(C.z-1, S.z), 2).hash(S), coord(C.x, C.y, C.z, 2).hash(S),
	};
	//each vertex contains 6 physical_qubits, left (Mx-), right(Mx+), up (Ly-), down (Ly+), less (Nz-), more (Nz+)
};

subvertex::subvertex(const int& c, const subcoord& S){
	int L = S.x;
	int M = S.y;
	x = c % L; //x subcoordinate
	y = c / L; //y subcoordinate
	partial = {y*L+divmod(x-1,M), c, M*L+divmod(y-1,L)*L+x, M*L+c
	};
	//each subvertex contains 4 physical_qubits, left (Mx-), right (Mx+), up, down
	//each vertex contains 6 physical_qubits, left (Mx-), right(Mx+), up (Ly-), down (Ly+), less (Nz-), more (Nz+)
};
