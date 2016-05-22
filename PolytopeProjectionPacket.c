#include <stdio.h>
#include "math.h"
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <float.h>

//prototipi
void projection2(double *uPointer, double *resPointer, int length);
double absoluteNorm(double *pointer, int length);
void projectionIntoBinarySubspace(double *pointer, double *outputPointer, int length);
double projectionIntoBinary(double value);
double absolute(double value);
void selection_sort(double v[],int nmax);
void multiplyMatrix(double *pointerA, double *pointerB, double *pointerC, int ma, int na, int mb, int nb);

//main
void mexFunction(
        int nlhs,               // number of outputs
        mxArray *plhs[],        // outputs vector
        int nrhs,               // number of inputs
        const mxArray *prhs[]   // inputs vector
        ) {
        
    const mxArray *matrixIn = prhs[0];
    const mxArray *d = prhs[1];
    int m = (int)mxGetM(matrixIn);
    int n = (int)mxGetN(matrixIn);
    int md = (int)mxGetM(d);
    int nd = (int)mxGetN(d);


    if(m > 1 && n > 1) 		mexErrMsgTxt("Function only takes vectors");
    if(m == 1 && n > 1)		mexErrMsgTxt("Function takes only vectors in row");
    if(m == 0 || n == 0) 	mexErrMsgTxt("Function does not take empty matrix"); 
    if(m == 1 && n == 1) 	mexErrMsgTxt("Function does not take scalars");
    //if(md != m && nd != n)	mexErrMsgTxt("Incompatible vectors length");
    
    
    if(m > 1 && n == 1) {
        plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
	    double *in;
        double *out;
        double *dp;
        in = mxGetPr(prhs[0]);
        dp = mxGetPr(prhs[1]);
        out = mxGetPr(plhs[0]);
        
        //projection2(in, out, m);
        
        int prev, curr;
        prev = 0;
        curr = (int)*dp;
        projection2(in, out, (int)*dp);
        for (int i = 1; i < md; i++) {
        	prev = curr;
        	curr += (int)*(dp+i);
        	projection2(in+prev, out+prev, (int)*(dp+i));
        }
        
    }
    
    /*
    
    z(1:d(1),1) = PolytopeProjection(v(1:d(1)),1);
    currPos = d(1);
    for j=2:m
        prevPos = currPos;
        currPos = currPos + d(j);
        z(prevPos+1:currPos,1) = PolytopeProjection(v(prevPos+1:currPos),1);
    end
    
    */
    
    
    return;
    
}

//function
void projection2(double *uPointer, double *resPointer, int length) {
	double z[length];
	double *zPointer = &z[0];
	projectionIntoBinarySubspace(uPointer, zPointer, length);
	
	double v[length];
	for (int i = 0; i < length; i++) {
		if (z[i] > 0.5) {
			v[i] = 1;
		} else {
			v[i] = 0;
		}
	}
	
	double norm = absoluteNorm(&v[0], length);
	if (!((int)(norm)%2)) {
		
		int j=0;
		double min;
		min = absolute(z[0]-0.5);
		
		for (int i = 1; i < length; i++) {
			if (absolute(z[i]-0.5)<min) {
				j=i;
				min = absolute(z[i]-0.5);
			}
		}
		v[j] = 1-v[j];
	}
	
	double zv[length];
	double *zvPointer = &zv[0];
	double tau[length];
	double *tauPointer = &tau[0];
	int tLength = 0;
	for (int i = 0; i < length; i++) {
		zv[i]  = z[i]-v[i];
	}
	if (absoluteNorm(zvPointer, length) >= 1) {
		for (int y=0; y< length; y++) *(resPointer+y) = z[y];
		return;
	}
	
	for (int i = 0; i< length; i++) {
		tau[i] = *(uPointer+i) * (1-v[i]) + (1-*(uPointer+i)) * v[i];
		if (tau[i]<0) tLength++;
	}	
	
	
	double t[tLength];
	int k=0;
	for (int i = 0; i < length; i++) {
		if (tau[i]<0) t[k++]=-tau[i];
	}
	selection_sort(t, tLength);
	double tauBin[length];
	double *tauBinPointer = &tauBin[0];
	projectionIntoBinarySubspace(tauPointer, tauBinPointer, length);
	double delta = 1 - absoluteNorm(tauBinPointer, length);
	
	
	double N = length-tLength;
	
	for (int n = 0; n < tLength; n++) {
		if (delta/N <= t[n]) {
			break;
		} else {
			delta = delta + t[n];
			N = N+1;
		}
	}	
	
	double lambda = delta/N;
	
	double uu[length];
	double *uuPointer = &uu[0];
	for (int i = 0; i < length; i++) {
		uu[i] = *(uPointer+i) - lambda*(2*v[i]-1);
	}
	
	projectionIntoBinarySubspace(uuPointer, resPointer, length);
	
}

//utility
double absoluteNorm(double *pointer, int length) {
	double sum = 0;
	for (int i = 0; i < length; i++) {
		double *currentPointer = pointer+i;
		sum += (*currentPointer > 0 ? *currentPointer : -(*currentPointer));
	}
	return sum;
}

void projectionIntoBinarySubspace(double *pointer, double *outputPointer, int length) {
	double *currentPointer;
	for (int i = 0; i < length; i=i+1) {
		currentPointer = pointer+i;
		//printf("%f	->", *currentPointer); 
		*(outputPointer+i)= projectionIntoBinary(*currentPointer);
		//printf("	%f\n", *(outputPointer+i));
	}
}

double projectionIntoBinary(double value) {
	if (value > 1) value = 1;
	if (value < 0) value = 0;
	return value; 
}

double absolute(double value) {
	return value<0 ? -value : value;
}

void selection_sort(double v[],int nmax){
    int i,j,min;
    double t;
    
    for (i=0; i<nmax-1; i++) {
		
        min = i;

        for (j=i+1; j<nmax; j++) {
            if (v[min]> v[j]) {
                min = j;
            }
        }
        
        if (min != i) {
			t = v[i];
			v[i] = v[min];
			v[min] = t;        
		}
		
	}
}

void multiplyMatrix(double *pointerA, double *pointerB, double *pointerC, int ma, int na, int mb, int nb) {
	
	if (nb==1) {
		for (int i = 0; i < ma; i++) {
		for (int j = 0; j < nb; j++) {
			double current=0;
			for (int k = 0; k < na; k++) {
				current = current + (*(pointerA + na*i + k))*(*(pointerB + j + nb*k));
			}
			*(pointerC + i + j*ma) = current;	
		}
	}
	} else {
		for (int i = 0; i < ma; i++) {
		for (int j = 0; j < nb; j++) {
			double current=0;
			for (int k = 0; k < na; k++) {
				current = current + (*(pointerA + na*i + k))*(*(pointerB + j + nb*k));
			}
			*(pointerC + j + i*ma) = current;	
		}
	}
	}
}
