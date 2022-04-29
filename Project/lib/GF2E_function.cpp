#include "GF2E_function.h"


using namespace std;
using namespace NTL;



mat_GF2 RandomMatGF2(unsigned long n, unsigned long m){
	//Return a random non zero n x m matrix in GF2
	while(1){
		mat_GF2 M;
		M.SetDims(n,m);
		for(unsigned int i = 0; i < n; i++){
			M[i] = random_vec_GF2(m);
		}
		if(!IsZero(M)){
			return M;
		}
	}
}


mat_GF2 random_invertible_mat_GF2(unsigned long n){
	//Return a random invertible matrix nxn in GF2
	while(1){
		mat_GF2 M = random_mat_GF2(n,n);
		if(determinant(M) != 0){
			return M;
		}
	}
}

uint8_t uint8FromMatGF2(NTL::mat_GF2 const & v){
	//return x the 8bit uint such that the binary form of x is the first 8 bit of the first column of v (v is supposed to be a 8x1 matrix)
	//返回x这个8位uint，即x的二进制形式是v的第一列的第一个8位（v应该是一个8x1的矩阵）。
	//uint8_t 就是 unsigned char，一个字节是char类型
	uint8_t x = 0;
	for(unsigned int i = 0; i < 8; i++){
		if(v[i][0] == 1){
			x ^= (1 << i);
		}
	}
	return x;
}

mat_GF2 matGF2FromUint8(uint8_t x){
	//return a 8x1 GF2 matrix built from the binary repr of x
	mat_GF2 v;
	v.SetDims(8,1);
	for(unsigned int i = 0; i < 8; i++){
		if(x & (1 << i)){
			v[i][0] = 1;
		}
	}
	return v;
}

mat_GF2 matGF2FromUint8_reverse(uint8_t x){
	//return a 1x8 GF2 matrix built from the binary repr of x
	mat_GF2 v,temp;
	vec_GF2 a;
	v.SetDims(8,1);
	for(unsigned int i = 0; i < 8; i++){
		if(x & (1 << i)){
			v[i][0] = 1;
		}
	}
	temp = transpose(v);
	a = reverse(temp[0]);
	VectorCopy(temp[0],a,8);//1x8
	
	return temp;
}



vec_GF2 vecGF2FromUint8(uint8_t x){
	//返回一个8bit的向量
	vec_GF2 v;
	v.SetLength(8);
	for(unsigned int i = 0; i < 8; i++){
		if(x & (1 << i)){
			v[i] = 1;
		}
	}
	v =  reverse(v);
	return v;
}

uint16_t uint16FromMatGF2(NTL::mat_GF2 const & v){
	//16*1 的矩阵
	uint16_t x = 0;
	for(unsigned int i = 0; i < 16; i++){
		if(v[i][0] == 1){
			x ^= (1 << i);
		}
	}
	return x;
}
mat_GF2 matGF2FromUint16(uint16_t x){
	mat_GF2 v;
	v.SetDims(16,1);
	for(unsigned int i = 0; i < 16; i++){
		if(x & (1 << i)){
			v[i][0] = 1;
		}
	}
	return v;
}

uint32_t uint32FromMatGF2(mat_GF2 const & v){
	uint32_t x = 0;
	for(unsigned int i = 0; i < 32; i++){
		if(v[i][0] == 1){
			x ^= (1 << i);
		}
	}
	return x;
}

mat_GF2 matGF2FromUint32(uint32_t x){
	mat_GF2 v;
	v.SetDims(32,1);
	for(unsigned int i = 0; i < 32; i++){
		if(x & (1 << i)){
			v[i][0] = 1;
		}
	}
	return v;
}




void compactPrint(NTL::mat_GF2 const & M, ostream& os){
//compact binary print of a NTL GF2 matrix
//打印NTL GF2矩阵
	for(unsigned int i = 0; i < M.NumRows(); i++){
		for(unsigned int j = 0; j < M.NumCols(); j++){
			if(M[i][j] == 0){
				os << 0;
			}
			else{
				os << 1;
			}
			if((j+1)%8==0){
				os << " ";
			}
		}
		os << endl;
	}
}

void superCompactPrint(NTL::mat_GF2 const & M, ostream& os){
//compact binary print of a NTL GF2 matrix, using hex notation 使用16进制打印矩阵
	//for example, a row 0001 1111 wil be printed as 0x1f
	for(unsigned int i = 0; i < M.NumRows(); i++){
		uint8_t x = 0; 
		for(unsigned int j = 0; j < M.NumCols(); j++){
			x <<= 1;
			if(M[i][j] == 1){
				x ^= 1;
			}

			if((j+1)%4 == 0){
				// printf("%01x",x);
				os << hex << setw(1) << setfill('0') << int(x)<< dec;
				x = 0;
			}
		}
		os << endl;
	}
}

std::set<unsigned int> echelonizeCustom(mat_GF2 & M, unsigned int nbPivotMax){
	//echelonize M on at most nbPivotMax pivots
	//return a set of column index of the pivots

	std::set<unsigned int> setPivots;
	unsigned int sysPos = 0;
	auto M_numrows = M.NumRows();
	auto M_numcols = M.NumCols();

	for(unsigned int j = 0; j < M_numcols; j++){
		bool elim = false;
		unsigned int pos = sysPos;
		while(pos != M_numrows){
			if(M[pos][j] != 0){
				//if coefficient j is non zero
				if(!elim){
					//if non zero coeff on column j was not found until now, use the current row as the pivot
					swap(M[pos],M[sysPos]);
					sysPos++;
					elim = true;
					setPivots.emplace(j);
				}
				else{
					//we had our pivot for the column j, so echelonize
					//the pivot is M[sysPos-1]
					M[pos] -= M[sysPos-1];
				}
			}
			++pos;
		}
		if(setPivots.size() == nbPivotMax){
			break;
		}
	}

	return setPivots;
}

GF2E GF2EFromLong(const long x){
	GF2X p;
	long i = 0;
	long loc_x = x;
	while(loc_x != 0){
		SetCoeff(p,i,(loc_x&1));
		loc_x >>= 1;
		i++;
	}
	
	return to_GF2E(p);
}

long LongFromGF2X(GF2X const & P){
	long x = 0;
	for(unsigned int i = 0; i <= deg(P); i++){
		if(coeff(P,i) == 1){
			x ^= (1 << i);
		}
	}

	return x;
}

std::set<unsigned int> reducedEchelonizeWithPivot(NTL::mat_GF2 & M){
//Compute a row-reduced echelon form of M
// --> remove all zeros rows after echelonization <-- (rank(M) = M.NumRows() at the end
//Return a set containing the index of the column used as pivot

	unsigned int indPivot = 0;
	auto const M_numrows = M.NumRows();
	auto const M_numcols = M.NumCols();
	std::set<unsigned int> setPivots;
	for(unsigned int j = 0; j < M_numcols; j++){
		//get the first row with a non-zero coefficient on column j
		bool foundPivot = false;
		for(unsigned int i = indPivot; i < M_numrows; i++){
			//row < indPivot are already pivot
			if(M[i][j] != 0){
				if(!foundPivot){
					//if we didn't find the pivot before, here it is (column j)
					swap(M[i], M[indPivot]);
					setPivots.emplace(j);
					foundPivot = true;
				}
				else{
					//else we already have our pivot on row indPivot, so echelonize
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

	//remove the zero lines (if any)
	M.SetDims(indPivot,M_numcols);

	return setPivots;
}