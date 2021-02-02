/**
 * Let $T$ be a threshold and let $[i..j]$ be an interval of the MS array such that  
 * $MS[i]>=T, MS[j]>=T$, and $MS[k]<T$ for all $k \in [i+1..j-1]$. Let $i',j'$ be the one-
 * -bits of the MS bitvector that correspond to positions $i,j$, and let $p$ be the first
 * zero after $i'$ and $q$ be the first one before $j'$. The program builds a DAG that 
 * encodes how to permute $[p..q]$ into a sequence of run-pairs $(nZeros_x,nOnes_x)$, 
 * where every $nZeros_x$ and $nOnes_x$ is positive, that minimizes the size of the 
 * encoding, where each $nZeros_x$ is delta-coded, and each $nOnes_x$ is encoded with 
 * respect to $nZeros_x$. The minimal encoding is equivalent to a shortest path in the 
 * DAG: the program computes all shortest paths in order to answer multiple queries later.
 */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>


/** 
 * Representation of the DAG
 */
uint8_t ALLOW_NEGATIVE_MS;
uint32_t THRESHOLD;
uint64_t N_ZEROS, N_ONES, N_NODES, N_ARCS;
uint64_t **inNeighbors;
uint64_t *rowLengths;
int64_t *lastInNeighbor;
uint64_t *cost;
int64_t *neighbor;


void allocateMemory() {
	const uint64_t CAPACITY = 10;  // Arbitrary
	uint64_t i;
	
	inNeighbors=(uint64_t **)malloc(N_NODES*sizeof(uint64_t *));
	for (i=0; i<N_NODES; i++) inNeighbors[i]=(uint64_t *)malloc(CAPACITY*sizeof(uint64_t *));
	rowLengths=(uint64_t *)malloc(N_NODES*sizeof(uint64_t));
	for (i=0; i<N_NODES; i++) rowLengths[i]=CAPACITY;
	lastInNeighbor=(int64_t *)malloc(N_NODES*sizeof(uint64_t));
	cost=(uint64_t *)malloc(N_NODES*sizeof(uint64_t));
	neighbor=(int64_t *)malloc(N_NODES*sizeof(int64_t));
	for (i=0; i<N_NODES; i++) lastInNeighbor[i]=-1;
}


void deallocateMemory() {
	uint64_t i;
	
	for (i=0; i<N_NODES; i++) free(inNeighbors[i]);
	free(rowLengths);
	free(lastInNeighbor);
	free(cost);
	free(neighbor);
}


void clean() {
	uint64_t i;
	
	for (i=0; i<N_NODES; i++) lastInNeighbor[i]=-1;
	for (i=0; i<N_NODES; i++) cost[i]=ULLONG_MAX;
	for (i=0; i<N_NODES; i++) neighbor[i]=-1;
}


/**
 * @param nZeros,nOnes >=1.
 */
uint64_t coordinates2id(const uint64_t nZeros, const uint64_t nOnes, const uint64_t totalNZeros, const uint64_t totalNOnes) {
	return (nZeros-1)*totalNOnes+(nOnes-1);
}


void id2coordinates(const uint64_t nodeID, const uint64_t totalNZeros, const uint64_t totalNOnes, uint64_t out[]) {
	out[0]=1+nodeID/totalNOnes;
	out[1]=1+nodeID%totalNOnes;
}


void addArc(const uint64_t nZerosFrom, const uint64_t nOnesFrom, const uint64_t nZerosTo, const uint64_t nOnesTo, const uint64_t totalNZeros, const uint64_t totalNOnes) {
	const uint64_t fromNode = coordinates2id(nZerosFrom,nOnesFrom,totalNZeros,totalNOnes);
	const uint64_t toNode = coordinates2id(nZerosTo,nOnesTo,totalNZeros,totalNOnes);
	
	lastInNeighbor[toNode]++;
	if (lastInNeighbor[toNode]==rowLengths[toNode]) {
		inNeighbors[toNode]=(uint64_t *)realloc(inNeighbors[toNode],(rowLengths[toNode]<<1)*sizeof(uint64_t));
		rowLengths[toNode]<<=1;
	}
	inNeighbors[toNode][lastInNeighbor[toNode]]=fromNode;
	N_ARCS++;
}


uint64_t length(uint64_t n) {
  uint64_t b = 0;
  while (n!=0) { b++; n>>=1; }
  return b;
}


uint64_t deltaCodeLength(const uint64_t value) {
	uint64_t len = length(value);
	uint64_t llen = length(len);
	return len+llen+llen-2;
}


/**
 * Encodes $nOnes$ as a difference with respect to $nZeros$, if $nZeros>threshold$.
 * If $nZeros<=threshold$, the procedure just returns $nOnes$.
 *
 * @param nOnes > 0;
 * @param base first integer > 0 to be used as the destination of the projection.
 */
uint64_t project(const uint64_t nZeros, const uint64_t nOnes, const uint64_t base, const uint64_t threshold) {
	uint64_t delta;

	if (nZeros<=threshold) return nOnes;
	if (nOnes==nZeros) return base;
	if (nOnes<nZeros) {
        delta=nZeros-nOnes;
        return base-1+(delta<<1);
	}
	delta=nOnes-nZeros;
	if (delta<=nZeros-1) return base+(delta<<1);
	return base+((nZeros-1)<<1)+(delta-nZeros+1);
}


/**
 * @return number of bits for encoding the run pair $(nZeros,nOnes)$.
 */
uint64_t encodingSize(const uint64_t nZeros, const uint64_t nOnes) {
	return deltaCodeLength(nZeros)+deltaCodeLength(project(nZeros,nOnes,1,0/*Forces to use always diff. encoding*/));
}


void buildDAG(const uint32_t initialMSValue) {
	const uint32_t maxMSValue = THRESHOLD-1;  // Max MS value in a permuted range
	uint8_t canBeFirst;
	int64_t i, j, p, q;
	uint64_t nZeros, nOnes;
	uint64_t fromNode, minCost, newCost;
	int64_t msValue, minNeighbor, fromI, toI, fromP, toP;
	uint64_t out[5];
	
	// Building the DAG
	printf("Building DAG for initialMSValue=%u, N_NODES=%lu, N_ZEROS=%lu, N_ONES=%lu... ",initialMSValue,N_NODES,N_ZEROS,N_ONES);
	for (j=1; j<=N_ONES; j++) {
		msValue=initialMSValue-j;
		if (ALLOW_NEGATIVE_MS) fromI=1;
		else {
			fromI=j-initialMSValue;
			if (fromI<1) fromI=1;
			else if (fromI>N_ZEROS) break;
		}
		toI=maxMSValue-msValue;
		if (toI>N_ZEROS) toI=N_ZEROS;
		toP=toI+1;  // This is because we want:
		//
		// (1) initialMSValue-q+p <= maxMSValue   // last one of pair (p,q)
		// (2) MS(i,j)+(p-i)-1 <= maxMSValue        // first one of pair (p,q)
		//     
		// (1) p <= maxMSValue-initialMSValue+q
		// (2) p <= maxMSValue-(initialMSValue-j+i)+i+1
		//     p <= maxMSValue-initialMSValue+j+1
		//
		if (toP>N_ZEROS) toP=N_ZEROS;
		for (i=fromI; i<=toI; i++) {
			for (q=j+1; q<=N_ONES; q++) {
				if (ALLOW_NEGATIVE_MS) fromP=i+1;
				else {
					fromP=q-initialMSValue;
					if (fromP<i+1) fromP=i+1;
					else if (fromP>N_ZEROS) break;
				}
				for (p=fromP; p<=toP; p++) addArc(i,j,p,q,N_ZEROS,N_ONES);
			}
		}
	}
	printf("done.  N_ARCS=%lu \n",N_ARCS);
	
	// Computing all single-source shortest paths in the DAG
	printf("Computing shortest paths... ");
	for (i=0; i<N_NODES; i++) {
		id2coordinates(i,N_ZEROS,N_ONES,out);
		nZeros=out[0]; nOnes=out[1];
		canBeFirst=1;
		if (initialMSValue+nZeros>maxMSValue+1) canBeFirst=0;
		else {
			if (!ALLOW_NEGATIVE_MS && initialMSValue+nZeros<nOnes) canBeFirst=0;
		}
		if (canBeFirst) minCost=encodingSize(nZeros,nOnes);
		else minCost=ULLONG_MAX;
		minNeighbor=-1;
		for (j=0; j<=lastInNeighbor[i]; j++) {
			fromNode=inNeighbors[i][j];
			if (cost[fromNode]==ULLONG_MAX) continue;
			id2coordinates(fromNode,N_ZEROS,N_ONES,out);
			newCost=cost[fromNode]+encodingSize(nZeros-out[0],nOnes-out[1]);
			if (newCost<minCost) {
				minCost=newCost;
				minNeighbor=fromNode;
			}				
		}
		cost[i]=minCost; neighbor[i]=minNeighbor;
	}
	printf("done \n");
}


/**
 * @param argv in the following order: 
 * threshold
 * initialMSValue_first (< threshold)
 * initialMSValue_last (< threshold)
 * maxNZeros 
 * maxNOnes 
 * allowNegativeMS 
 * saveAll 0=saves just the $neighbor$ array;
 * destinationDir
 */
int main(int argc, char **argv){
	THRESHOLD=(uint32_t)atoi(argv[1]);
	const uint32_t INITIAL_MS_VALUE_FROM = (uint64_t)atol(argv[2]);
	const uint32_t INITIAL_MS_VALUE_TO = (uint64_t)atol(argv[3]);
	N_ZEROS=(uint64_t)atol(argv[4]);
	N_ONES=(uint64_t)atol(argv[5]);
	ALLOW_NEGATIVE_MS=(uint8_t)atoi(argv[6]);
	const uint8_t SAVE_ALL = (uint8_t)atoi(argv[7]);
	const char *DIRECTORY = argv[8];
	N_NODES=N_ZEROS*N_ONES;

	uint32_t i, j;
	FILE *file;
	char *buffer = (char *)malloc(100*sizeof(char));
	
	allocateMemory();
	for (i=INITIAL_MS_VALUE_FROM; i<=INITIAL_MS_VALUE_TO; i++) {
		N_ARCS=0;
		clean();
		buildDAG(i);
		sprintf(buffer,"%s/table-diff-%d-%d-%lu-%lu-%d",DIRECTORY,THRESHOLD,i,N_ZEROS,N_ONES,ALLOW_NEGATIVE_MS);
		file=fopen(buffer,"w");
		fwrite(&N_NODES,sizeof(uint64_t),1,file);
		fwrite(neighbor,sizeof(int64_t),N_NODES,file);
		if (SAVE_ALL) {
			fwrite(cost,sizeof(uint64_t),N_NODES,file);
			fwrite(lastInNeighbor,sizeof(int64_t),N_NODES,file);
			for (j=0; j<N_NODES; j++) fwrite(inNeighbors[j],sizeof(uint64_t),lastInNeighbor[j]+1,file);
		}
		fclose(file);
	}
	free(buffer);
	deallocateMemory();
}