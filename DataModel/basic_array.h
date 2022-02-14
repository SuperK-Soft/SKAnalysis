/* vim:set noexpandtab tabstop=8 wrap */
#ifndef BasicArray_H
#define BasicArray_H
#include "type_name_as_string.h"

#include <typeinfo>
#include <iostream>
#include <vector>
#include <iterator>

/*
class generic_array {
	public:
	generic_array(){};
	~generic_array(){};
	protected:
	void SetParent(generic_array* parent_ptr){
		parent=parent_ptr;
	};
	generic_array* GetParent(){
		return parent;
	}
	generic_array* parent=NULL;
};
*/

template <class T, bool = (std::is_array<T>::value|
			    ((std::is_pointer<T>::value && 
			      std::is_array<typename std::remove_pointer<T>::type>::value) |
			     (std::is_pointer<T>::value && 
			      std::is_pointer<typename std::remove_pointer<T>::type>::value))
			  ) >
class basic_array /*: public generic_array*/ {
	typedef typename std::remove_pointer<T>::type W;
	
	public:
	basic_array(){
		//std::cout<<"default constructing a basic_array of primitives"<<std::endl;
	};
	basic_array(intptr_t addr_in, size_t sizein/*, generic_array* parent=NULL*/){
		//std::cout<<"constructing a basic array (of primitives) of size "<<sizein<<" and type "
		//	 << type_name<T>() <<std::endl;
		a_size=sizein;
		
		bool isptr = std::is_pointer<T>::value;
		W** app = reinterpret_cast<W**>(addr_in);
		//std::cout<<"isptr is "<<isptr<<", app is "<<app<<std::endl;
		addr =   (isptr) ? reinterpret_cast<W*>(app) :
				   reinterpret_cast<W*>(addr_in);
		//std::cout<<"addr is of type "<<type_name<W*>()<<" and is at "<<addr<<std::endl;
	}
	basic_array(intptr_t addr_in, std::vector<size_t> sizes/*, generic_array* parent=NULL*/){
		a_size = 0;
		if(sizes.size()){
			a_size=sizes.front();
		}
		//std::cout<<"constructing a basic array (of primitives) of size "<<a_size<<" and type "
		//	 << type_name<T>() <<std::endl;
		
		bool isptr = std::is_pointer<T>::value;
		W** app = reinterpret_cast<W**>(addr_in);
		//std::cout<<"isptr is "<<isptr<<", app is "<<app<<std::endl;
		addr =   (isptr) ? reinterpret_cast<W*>(app) :
				   reinterpret_cast<W*>(addr_in);
		//std::cout<<"addr is of type "<<type_name<W*>()<<" and is at "<<addr<<std::endl;
	}
	
	// pointer to array
	template<class X>
	basic_array(X ptr_in, typename std::enable_if<std::is_same<typename std::remove_extent<typename std::remove_pointer<X>::type>::type, T>::value, bool>::type potato= true){
		//std::cout<<"pointer constructor with arg type "<<type_name<X>()<<std::endl;
		a_size = sizeof(typename std::remove_pointer<X>::type)/sizeof(T);
		addr = &(*ptr_in)[0];
	}
	
	// reference to array
	template<typename X, size_t N>
	basic_array(X (&ref_in)[N]){
		//std::cout<<"reference constructor with arg type "<<type_name<decltype(ref_in)>()<<std::endl;
		//std::cout<<"X is of type "<<type_name<X>()<<std::endl;
		a_size = N; //sizeof(typename std::remove_pointer<X>::type)/sizeof(T);
		addr = &(const_cast<W*>(ref_in))[0];
	}
	
	W operator[](int i) const {
		return addr[i];
	}
	
	~basic_array(){};
	
	const W& at(int i) const {
		if(i<0||i>a_size){
			throw std::out_of_range ("out of range exception requesting element "+std::to_string(i)+" in "+__FILE__+"::"+std::to_string(__LINE__));
		}
		return addr[i];
	}
	const W& back() const {
		return addr[a_size-1];
	}
	const W* begin() const {
		return cbegin();
	}
	const W* cbegin() const {
		const W* theit = &addr[0];
		return theit;
	}
	const W* cend() const {
		const W* theit = &addr[a_size];
		return theit;
	}
	const W* crbegin() const {
		return std::reverse_iterator<W*>(&addr[a_size]);
	}
	const W* crend() const {
		return std::reverse_iterator<W*>(&addr[0]-1);
	}
	const W* data() const {
		return &addr[0];
	}
	bool empty() const {
		return (a_size==0);
	}
	const W* end() const {
		return cend();
	}
	const W& front() const {
		return addr[0];
	}
	W* rbegin() const {
		return crbegin();
	}
	W* rend() const {
		return crend();
	}
	int size() const {
		return a_size;
	}
	int max_size() const {
		return a_size;
	}
	int dimensions() const {
		return 1;
	}
	
//	fill - not implemented, this class is const!
//	swap - not implemented, this class is const!
	
	protected:
	W* addr;
	size_t a_size;
};

// implementation for higher-dimensional arrays
template < class T>
class basic_array<T,true> /*: public generic_array*/ {
	typedef typename std::remove_pointer<T>::type W;
	typedef typename std::remove_extent<W>::type U;  // raw element type
	typedef basic_array<U> V; // wrapped element type
	public:
	basic_array(){
		//std::cout<<"default constructing a basic_array of arrays"<<std::endl;
	};
	basic_array(intptr_t addr_in, size_t sizein_1, size_t sizein_2=0/*, generic_array* parent=NULL*/){
		//std::cout<<"constructing a basic array (of arrays) of size "<<sizein_1<<" and type "
		//	 << type_name<T>() <<std::endl;
		a_size = sizein_1;
		Init(addr_in, std::vector<size_t>{sizein_2});
		//std::cout<<"array of arrays is a wrapper around an array of "<<type_name<W>()
		//	 <<" objects. subarray is of type "<<type_name<decltype(subarray)>()
		//	 <<" and holds objects of type "<<type_name<decltype(subarray[0])>()
		//	 <<" or alternatively "<<type_name<decltype(subarray.at(0))>()
		//	 <<" while return type is "<<type_name<V>()<<std::endl;
	}
	basic_array(intptr_t addr_in, std::vector<size_t> sizes/*, generic_array* parent=NULL*/){
		if(sizes.size()==0){
			a_size = 0;
		} else {
			a_size = sizes.front();
			sizes.erase(sizes.begin());
		}
		//std::cout<<"constructing a basic array (of arrays) of size "<<a_size<<" and type "
		//	 << type_name<T>() <<std::endl;
		Init(addr_in, sizes);
		//std::cout<<"array of arrays is a wrapper around an array of "<<type_name<W>()
		//	 <<" objects. subarray is of type "<<type_name<decltype(subarray)>()
		//	 <<" and holds objects of type "<<type_name<decltype(subarray[0])>()
		//	 <<" or alternatively "<<type_name<decltype(subarray.at(0))>()
		//	 <<" while return type is "<<type_name<V>()<<std::endl;
	}
	
	template<class X>
	basic_array(X ref_in, typename std::enable_if<std::is_same<typename std::remove_extent<typename std::remove_pointer<X>::type>::type, T>::value, bool>::type potato= true){
		//std::cout<<"reference constructor with arg type "<<type_name<X>()<<std::endl;
		a_size = sizeof(typename std::remove_pointer<X>::type)/sizeof(T);
		Init(reinterpret_cast<intptr_t>(&ref_in[0]), std::vector<size_t>{});
	}
	
	bool Init(intptr_t addr_in, std::vector<size_t> sizes){
		bool isptr = std::is_pointer<decltype(addr_in)>::value;
		
		T** app = reinterpret_cast<T**>(addr_in);
		addr =   (isptr) ? reinterpret_cast<W*>(app[0]) :
				   reinterpret_cast<W*>(addr_in);
		//std::cout<<"type of array element 0 is: "
		//	 <<type_name<decltype(addr[0])>() <<std::endl;
		if(sizes.size()==0 || sizes.front()==0){
			sizes = std::vector<size_t>{sizeof(W)/sizeof(U)};
		}
		a_dimensions = sizes.size();
		subarray.clear();
		//std::cout<<"expect subarray to hold items of type "
		//	 <<type_name<W>()<<" while addr[0] is of type "
		//	 <<type_name<decltype(addr[0])>()<<std::endl;
		
		bool is_nested_ptr = (std::is_pointer<T>::value && 
			      std::is_pointer<typename std::remove_pointer<T>::type>::value);
		for(int i=0; i<a_size; ++i){
			W* subaddr = (is_nested_ptr) ? (W*)addr[i] :
						       (W*)&addr[i]; // : *(addr[i]);
			void* subaddr_v = (void*)subaddr;
			intptr_t subaddr_i = reinterpret_cast<intptr_t>(subaddr_v);
			subarray.emplace_back(subaddr_i, sizes/*,(generic_array*)this*/);
		}
		return true;
	}
	
	V operator[](int i) const {
		return subarray[i];
	}
	
	~basic_array(){};
	
	const V& at(int i) const {
		if(i<0||i>a_size){
			throw std::out_of_range ("out of range exception requesting element "+std::to_string(i)+" in "+__FILE__+"::"+std::to_string(__LINE__));
		}
		return subarray.at(i);
	}
	const V back() const {
		return subarray[a_size-1];
	}
	const V* begin() const {
		return cbegin();
	}
	const V* cbegin() const {
		return &subarray[0];
	}
	const V* cend() const {
		return &subarray[a_size];
	}
	const V* crbegin() const {
		return std::reverse_iterator<V*>(&subarray[a_size]);
	}
	const V* crend() const {
		return std::reverse_iterator<V*>(&subarray[0]-1);
	}
	const V* data() const {
		return addr;
	}
	bool empty() const {
		return (a_size==0);
	}
	const V* end() const {
		return cend();
	}
	const V front() const {
		return subarray[0];
	}
	V* rbegin() const {
		return crbegin();
	}
	V* rend() const {
		return crend();
	}
	int size() const {
		return a_size;
	}
	int max_size() const {
		return a_size;
	}
	int dimensions() const {
		//return (std::rank<W>::value)+1; // doesn't work properly with dynamic arrays
		return a_dimensions;
	}
//	fill - not implemented, this class is const!
//	swap - not implemented, this class is const!
	
	protected:
	W* addr;
	size_t a_size;
	size_t a_dimensions;
	std::vector<basic_array<U>> subarray;
};

#endif // define BasicArray_H

/* Usage:
	// STATICALLY SIZED ARRAYS
	// =======================
	
	// 1. DECLARE AN ARRAY
	// --------------------
	int anarray[4] = {5,6,7,8};
//	int anarray[] = {5,6,7,8};
//	int anarray[][3] = {{5,6,7},{55,66,77},{555,666,777},{5555,6666,7777}};
//	int anarray[4][3] = {{5,6,7},{55,66,77},{555,666,777},{5555,6666,7777}};
//	double anarray[][3][2] = {{{5.1,5.2},{6.1,6.2},{7.1,7.2}},{{50.1,50.2},{60.1,60.2},{70.1,70.2}},{{500.1,500.2},{600.1,600.2},{700.1,700.2}},{{5000.1,5000.2},{6000.1,6000.2},{7000.1,7000.2}}};
	
	// 2. BUILD THE WRAPPER
	// ---------------------
	// note that &anarray returns a pointer to an array of the correct type
	// (e.g. float(*)[3], which is a pointer to an array of 3 floats, which may be declared as:
	// float (*myarraypointer)[3] = &myarray;   // and used as
	// (*myarraypointer)[1]=20;
	// note this is NOT just a pointer to the first element (float* thing = &anarray[0])
	basic_array<std::remove_extent<decltype(anarray)>::type> myarray(&anarray);
	// float (&myarrayref)[3] = myarray;        // is a reference to an array of 3 floats
	// alternatively:
	// typedef int array_type[3];
	// array_type& ra = a;
	
	// DYNAMICALLY SIZED ARRAYS
	// ========================
	
	// 1. DECLARE AN ARRAY
	// --------------------
	// 1D
//	int* array = new int[4];
//	size_t dim = 4;
//	array[0]=5; array[1]=6; array[2]=7; array[3]=8;
	
	// 2D
//	int (*array)[3] = new int[4][3];
//	std::vector<size_t> dims{4,3};
//	for(int i=0; i<4; ++i){
//		for(int j=0; j<3; ++j){
//			array[i][j]=(i+1)*(j+1);
//		}
//	}
//	//delete [] array;    // cleanup
	// --
	// OR
	// --
//	size_t rowSize = 4;
//	//size_t colSize = 4;  // must be constexpr
//	std::vector<size_t> dims{rowSize,colSize};
//	int(*array)[colSize] = (int(*)[colSize]) new int[rowSize*colSize];
//	// (sizeof(array)/sizeof(array[0])) does not work here
//	for(int i=0; i<rowSize; ++i){
//		for(int j=0; j<colSize; ++j){
//			array[i][j]=(i+1)*(j+1);
//		}
//	}
//	// delete [] array   // cleanup
	// --
	// OR
	// --
//	size_t rowSize2 = 4;
//	size_t colSize2 = 3;
//	std::vector<size_t> dims{rowSize2,colSize2};
//	int **array = new int*[rowSize2];
//	for(int i = 0; i < rowSize2; ++i) {
//		array[i] = new int[colSize2];
//		std::cout<<"row "<<i<<" starts at "<<array[i]<<std::endl;
//	}
//	for(int i=0; i<4; ++i){
//		for(int j=0; j<3; ++j){
//			array[i][j]=(i+1)*(j+1);
//		}
//	}
//	//for(int i = 0; i < colSize2; ++i) { delete [] ary[i]; } delete [] ary;   // cleanup
	
	// 3D
//	size_t rowSize2 = 4;
//	size_t colSize2 = 3;
//	size_t aisleSize2 = 2;
//	std::vector<size_t> dims{rowSize2, colSize2, aisleSize2};
//	int ***array = new int**[rowSize2];
//	for(int i = 0; i < rowSize2; ++i) {
//		array[i] = new int*[colSize2];
//		for(int j=0; j < colSize2; ++j){
//			array[i][j] = new int[aisleSize2];
//			for(int k=0; k<aisleSize2; ++k){
//				array[i][j][k]=(i+1)*(j+1)*(k+1);
//				std::cout<<"array["<<i<<"]["<<j<<"]["<<k<<"]="<<array[i][j][k]<<", ";
//			}
//			std::cout<<std::endl;
//		}
//	}
//	// cleanup
//	//for (int i = 0; i < rowSize2; ++i){
//	//	for (int j = 0; j < colSize2; ++j){
//	//		delete[] A[i][j];
//	//	}
//	//	delete[] A[i];
//	//}
//	//delete[] A;
	
	// 2. BUILD THE WRAPPER
	// ---------------------
//	basic_array<std::remove_extent<decltype(anarray)>::type> myarray(ip, 4);    // 1D: specify size directly
//	basic_array<std::remove_extent<decltype(anarray)>::type> myarray(ip, 4, 3); // 2D: specify size directly
//	basic_array<std::remove_extent<decltype(anarray)>::type> myarray(ip, dims); // ANYD: pass vector of sizes
//	// rather than std::remove_extent<decltype(anarray)>::type, it may be easier in some cases to specify
//	// type directly. Use 'int*' for a 1D array, and continue adding '*'s for each further dimension:
//	// e.g. basic_array<int*> myarray(ip, 4), basic_array<int**>... for 2D, basic_array<int***> for 3D etc
	
	
	// USING THE WRAPPER
	// =================
	// use myarray as if it were a std::array<T, size>
	// or std::array<std::array<T, size1>, size2>  for 2D, or further nested for higher dimensions
	
	std::cout<<"myarray has size "<<myarray.size();
	std::cout<<" and dimensionality "<<myarray.dimensions()<<std::endl;
	std::cout<<", beginning "<<myarray.begin();
	std::cout<<", and end "<<myarray.end()<<std::endl;
	
	// 1D
//	std::cout<<", first element "<<myarray.front();
//	std::cout<<", last element "<<myarray.back();
//	for(auto it = myarray.begin(); it!=myarray.end(); ++it){
//		std::cout<<"myarray element "<<std::distance(myarray.begin(),it)<<" is "<<(*it)<<std::endl;
//	}
	// 2D
//	std::cout<<", first element has size "<<myarray.front().size()<<std::endl;
//	for(auto it = myarray.begin(); it!=myarray.end(); ++it){
//		std::cout<<"myarray element "<<std::distance(myarray.begin(),it)<<"[0] is "<<((*it)[0])<<std::endl;
//		int* mem = reinterpret_cast<int*>((*it)[0]);
//		std::cout<<"manual: "<<mem<<std::endl;
//	}
	// 3D
//	std::cout<<"first element has size "<<myarray.front().front().size();
//	std::cout<<std::endl<<std::endl;
//	std::cout<<"array dimensions are: "<<myarray.dimensions();
//	std::cout<<"->";
//	std::cout<<myarray.at(0).dimensions();
//	std::cout<<"->";
//	std::cout<<myarray.at(0).at(0).dimensions();
//	std::cout<<std::endl;
//	std::cout<<"getting top element 0"<<std::endl;
//	auto athing1  = myarray.at(0);
//	std::cout<<" it has type "<<type_name<decltype(athing1)>()<<std::endl;
//	std::cout<<"getting middle element 0"<<std::endl;
//	auto athing2 = myarray.at(0).at(0);
//	std::cout<<" it has type "<<type_name<decltype(athing2)>()<<std::endl;
//	std::cout<<"getting bottom element 0"<<std::endl;
//	auto athing3 = myarray.at(0).at(0).at(0);
//	std::cout<<" it has type "<<type_name<decltype(athing3)>()<<std::endl;
//	
//	std::cout<<"again with operator [], top:"<<std::endl;
//	auto athing11  = myarray[0];
//	std::cout<<" top has type "<<type_name<decltype(athing11)>()<<std::endl;
//	auto athing22  = myarray[0][0];
//	std::cout<<" middle has type "<<type_name<decltype(athing22)>()<<std::endl;
//	auto athing33  = myarray[0][0][0];
//	std::cout<<" bottom has type "<<type_name<decltype(athing33)>()<<std::endl;
//	for(auto it = myarray.begin(); it!=myarray.end(); ++it){
//		std::cout<<"myarray element "<<std::distance(myarray.begin(),it)<<"[0][0] is "
//			 <<((*it)[0][0])<<std::endl;
//	}

*/
