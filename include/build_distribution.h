#ifndef BUILD_DISTRIBUTION_H
#define BUILD_DISTRIBUTION_H


/** This function builds the distribution array for a distributed
    sparse matrix in CSR format. DomDof is a vector of size nproc
    containing the number of columns on the domain. Dist is a vector
    of size nproc+1 which after this function will contain the the
    index of the first column on the domain. The communicator is
    only included in this function for generality of determining the
    number of domains. This function should be run on all processors
    and will not require communication. */
void build_distribution(long *DomDof, int *Dist, int nproc);

#endif /* #ifndef BUILD_DISTRIBUTION_H */
