// implementation of the vector printing functions

template <class T>
ostream& operator<<(ostream& os, const vector<T> vec) {
	for (int i = 0; i < vec.size(); i++) {
		os << vec[i] << "|";
	}
	return os;
}

template <class T, class S>
ostream& operator<<(ostream& os, const pair<T,S> aPair) {
	os << "(" << aPair.first << "," << aPair.second << ")";
	return os;
}

template <class T>
void printMatrix(const vector<T> vec) {
	for (int i=0; i<vec.size(); i++) {
		cout << vec[i] << endl;
	}
}