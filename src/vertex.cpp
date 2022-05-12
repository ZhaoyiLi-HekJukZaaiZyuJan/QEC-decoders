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

