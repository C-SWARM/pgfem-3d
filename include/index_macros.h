/* HEADER */
#ifndef INDEX_MACROS_H
#define INDEX_MACROS_H

/** This file defines indexing macros for matrices/tensors */

#define idx_2(row,col) ((unsigned int) ((row)*3+(col)))

#define idx_2_gen(row,col,nrows,ncols) ((unsigned int) ((row)*(ncols) + (col)))

#define idx_3(i,j,k) ((unsigned int) ((i)*9 + (j)*3 + (k)))

#define idx_3_gen(i,j,k,I,J,K) ((unsigned int) ((i)*(J)*(K) + (j)*(K) + (k)))

#define idx_4(i,j,k,l) ((unsigned int) ((i)*27 + (j)*9 + (k)*3 + (l)))

#define idx_4_gen(i,j,k,l,I,J,K,L) ((unsigned int) ((i)*(J)*(K)*(L) \
						    + (j)*(K)*(L) \
						    + (k)*(L) + (l)))

#define idx_6(i,j,k,l,m,n) ((unsigned int) ((i)*243 + (j)*81 \
					    + (k)*27 + (l)*9 \
					    + (m)*3 + (n) ))

#define idx_K(a,b,w,g,nne,ndofn) ((unsigned int) ((((a)*(ndofn)+(b))*(nne)*(ndofn)\
						   + (w)*(ndofn) + (g))))

#define idx_K_gen(a,b,w,g,n_row,d_row,n_col,d_col) ((unsigned int) (((a)*(d_row)\
								     +(b))\
								    *(n_col)*(d_col)\
								    +(w)*(d_col)+(g)))

#endif /* #ifndef INDEX_MACROS_H */
