#include <vector>
#include <iostream>

template<typename T>
void print_vector(const std::vector<T> &__vec) {
	if (!__vec.size())
		std::cerr<<"print_vector() : Empty vector\n";
	std::cout<<"Printing vector :\n";
	for (int i=0; i<__vec.size(); i++) {
		std::cout<<"i : "<<i<<" :\n"<<__vec[i]<<"\n";
	}
}


