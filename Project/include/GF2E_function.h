#ifndef H_GF2EF
#define H_GF2EF

#include <NTL/mat_GF2E.h>
#include <NTL/mat_GF2.h>
#include <iostream>
#include <iomanip>
#include <set>
#include <cstdint>

NTL::mat_GF2 RandomMatGF2(unsigned long n, unsigned long m);
//Return a random n x m matrix in GF2

NTL::mat_GF2 random_invertible_mat_GF2(unsigned long n);
//Return a random invertible nxn matrix in GF2

NTL::vec_GF2 vecGF2FromUint8(uint8_t x);
uint8_t uint8FromMatGF2(NTL::mat_GF2 const & v);
NTL::mat_GF2 matGF2FromUint8(uint8_t x);
NTL::mat_GF2 matGF2FromUint8_reverse(uint8_t x);
uint16_t uint16FromMatGF2(NTL::mat_GF2 const & v);
NTL::mat_GF2 matGF2FromUint16(uint16_t x);
uint32_t uint32FromMatGF2(NTL::mat_GF2 const & v);
NTL::mat_GF2 matGF2FromUint32(uint32_t x);

std::set<unsigned int> reducedEchelonizeWithPivot(NTL::mat_GF2 & M);
//Compute a row-reduced echelon form of M
// --> remove all zeros rows after echelonization <-- (rank(M) = M.NumRows() at the end
//Return a set containing the index of the column used as pivot

std::set<unsigned int> echelonizeCustom(NTL::mat_GF2 & M, unsigned int nbPivotMax);

void compactPrint(NTL::mat_GF2 const & M, std::ostream& os);
void superCompactPrint(NTL::mat_GF2 const & M, std::ostream& os);

template<typename T> T mat_augment(T const & M1, T const & M2){
	//Return a matrix constructed as M = M1|M2
	//Number of rows of M is the min between number of rows of M1 and M2
	T M;
	M.SetDims( std::min(M1.NumRows(),M2.NumRows()),  M1.NumCols() + M2.NumCols());
	
	for(long i = 0; i < M.NumRows(); i++){
		auto & rowM = M[i];
		const auto & rowM1 = M1[i];
		const auto & rowM2 = M2[i];
		for(long j = 0; j <  M1.NumCols(); j++){
			rowM[j] = rowM1[j];
		}
		for(long j = M1.NumCols(); j <  M1.NumCols() + M2.NumCols(); j++){
			rowM[j] = rowM2[j-M1.NumCols()];
		}
	}
	
	return M;
}

template<typename T> void ref_mat_augment(T & M1, T const & M2){
	//Set M1 to M1|M2
	//Number of rows of M is the min between number of rows of M1 and M2
	//SetDims whith a bigger number of columns would completly reset M1, so we need a temp matrix M
	T M;
	M.SetDims( std::min(M1.NumRows(),M2.NumRows()),  M1.NumCols() + M2.NumCols());
	unsigned int M1_numcols = M1.NumCols();
	unsigned int M1M2_numcols = M1.NumCols() + M2.NumCols();
	
	for(long i = 0; i < M.NumRows(); i++){
		auto & rowM = M[i];
		auto & rowM1 = M1[i];
		const auto & rowM2 = M2[i];
		for(long j = 0; j <  M1_numcols; j++){
			rowM[j] = rowM1[j];
		}
		for(long j = M1_numcols; j < M1M2_numcols; j++){
			rowM[j] = rowM2[j-M1_numcols];
		}
	}
	
	M1 = std::move(M);
}


template<typename T> T mat_stack(T const & M1, T const & M2){
//Return a matrix constructed as M = M1
//									 M2
//Number of columns of M is the min between number of columns of M1 and M2

	T M;
	M.SetDims(M1.NumRows() + M2.NumRows(), std::min(M1.NumCols(),M2.NumCols()));
	
	for(long i = 0; i < M1.NumRows(); i++){
		M[i] = M1[i];
	}
	
	for(long i = M1.NumRows(); i < M1.NumRows() + M2.NumRows(); i++){
		M[i] = M2[i-M1.NumRows()];
	}
	
	return M;
}

template<typename T> void ref_mat_stack(T & M1, T const & M2){
//set M1 = M1
//		   M2
//Number of columns of M is the min between number of columns of M1 and M2
	long oldM1_numrows = M1.NumRows();
	M1.SetDims(M1.NumRows() + M2.NumRows(), std::min(M1.NumCols(),M2.NumCols()));
	
	// for(long i = oldM1_numrows; i < oldM1_numrows + M2.NumRows(); i++){
	// 	auto & rowM1 = M1[i];
	// 	const auto & rowM2 = M2[i-oldM1_numrows];
	// 	for(long j = 0; j < M1.NumCols(); j++){
	// 		rowM1[j] = rowM2[j];
	// 	}
	// }

	for(long i = oldM1_numrows; i < oldM1_numrows + M2.NumRows(); i++){
		M1[i] = M2[i-oldM1_numrows];
	}
}

NTL::GF2E GF2EFromLong(const long x);
//Return a GF2E element from binary representation of x

long LongFromGF2X(NTL::GF2X const & P);

template<typename T> unsigned int reducedEchelonize(T & M, unsigned int c){
	//reduced echelon form of M on the first c columns
	//return the index of the lastPivot
	unsigned int indPivot = 0;
	for(unsigned int j = 0; j < c; j++){
		//for each column < c

		//get the first row with a non-zero coefficient on column j
		bool foundPivot = false;
		for(unsigned int i = indPivot; i < M.NumRows(); i++){
			//row < j are already pivot
			if(M[i][j] != 0){
				if(!foundPivot){
					//if we didn't find the pivot before, here it is
					swap(M[i], M[indPivot]);
					foundPivot = true;
				}
				else{
					//else we already have our pivot on row j, so echelonize
					M[i] -= M[indPivot];
				}
			}
		}

		//and now reduce
		if(foundPivot){
			auto eqPivot = M[indPivot];
			for(int i = indPivot-1; i >= 0; i--){
				auto & rowM = M[i];
				if(rowM[j] != 0){
					rowM -= eqPivot;
				}
			}
			++indPivot;
		}
	}

	return indPivot-1;
}
#endif