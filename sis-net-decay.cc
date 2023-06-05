/* Decay rate for SIS on networks via direct numerical simulations.
Undirected and unweighted networks assumed.

When using this code, please cite the following paper:

Naoki Masuda, Victor M. Preciado, Masaki Ogura.
Analysis of the susceptible-infected-susceptible epidemic dynamics in networks via the non-backtracking matrix.
IMA Journal of Applied Mathematics, 85, 214-230 (2020).

The input file, the edge list, should be in the following format:
    The first row of the input file should be "# number-of-nodes number-of-edges"
    Starting from the second row, the first two columns have two nodes to form an edge, and the third column has the edge weight.
    However, we ignore the edge weight (i.e., third column, starting from the second row).
    The node index starts from 1, not from 0. 
    The node index should be consecutive.
    If users do not like these rules on the input file, it should be easy to rewrite the read-file part of the following code.

 The output is cout.
    1st column: infection rate (i.e., lambda)
    2nd column: decay rate
    3rd column: Pearson correlation coefficient in the linear regression between the time and the fraction of infectious nodes
    4th column: number of time points used for the linear regression
*/

#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cstring>
#include <cmath> // sqrt

#include "mt19937ar.c"
/* Mersenne Twister to generate random variates on [0,1].
   mt19937ar.c is available at 
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c

   If one prefers to use the in-built random number generator,
   (genrand_int32()+0.5)/4294967296.0 should be replaced by (double)rand()/RAND_MAX
   and
   init_genrand(time(NULL)) should be replaced by, e.g., srand(time(NULL)) */

int compare_sir_net (const void *a, const void *b) {
  return (*(int*)a-*(int*)b);
}

int main (int argc, char **argv) {

    if (argc != 2 && argc != 3) {
        cerr << "Usage: sis-net-decay.out infile_name (or data_type; either has to be specified) (infection-rate)" << endl;
        cerr << "infection-rate is optional" << endl; 
    exit(8);
    }

    init_genrand(time(NULL)); // seed the random number generator

    int i, j, vs, ve; // counters
    int nV; // # nodes
    int nE; // # bidirectional edges
    char dummy[8];

    ifstream fin; // input file
    double lam=0.1; // infection rate
    if (argc == 3)
        lam = atof(argv[2]);

    int lam_samples=30; // # lambda values
    double lam_min, lam_max;
    int trials=10000; // not multiplied by # nodes

    int data = atoi(argv[1]); // data_type (or input file name)
    if (data >= 1 && data <= 8) {
        if (data==1) {
            cerr << "regular random graph" << endl;
            fin.open("rrg-n100k6-decayrate.mat");
            lam_min = 0.0;
            lam_max = 0.28;
        } else if (data==2) {
            cerr << "BA model" << endl;
            fin.open("ba-n100m3-decayrate.mat");
            lam_min=0.0;
            lam_max=0.2;
        } else if (data==3) { 
            cerr << "cycle graph" << endl;
            fin.open("cycle-n100.mat");
            lam_min=0.0;
            lam_max=1.6;
        } else if (data==4) {
            cerr << "LFR benchmark network" << endl;
            fin.open("lfr-n100k6-decayrate.mat");
            lam_min=0.0;
            lam_max=0.22;
        } else if (data==5) {
            cerr << "dolphin network" << endl;
            fin.open("dolphins.mat");
            lam_min=0.0;
            lam_max=0.35;
        } else if (data==6) {
            cerr << "LCC of the network of network scientists" << endl;
            fin.open("netscience-lcc.mat");
            lam_min=0.0;
            lam_max=0.24;
        } else if (data==7) {
            cerr << "email network" << endl;
            fin.open("email-arenas.mat");
            lam_min=0.0;
            lam_max=0.066; // With 0.07, the decay rate is negative and corr = NaN.
        } else if (data==8) {
            cerr << "LCC of the Hamsterster netework" << endl;
            fin.open("hamsterster-lcc.mat");
            lam_min=0.0;
            lam_max=0.03;
        }
    } else {
        fin.open(argv[1]); // input edge list
    }
    if (!fin) {
        cerr << "sis-net-decay.out: cannot open input edge list file" << endl;
        exit(8);
    }

    // read input file
    fin >> dummy >> nV >> nE;
    cerr << nV << " nodes, " << nE << " edges" << endl; // undirected net assumed
    int *accum_k; // cumulative degree
    accum_k = (int *)malloc(nV*sizeof(int)); // allocate memory. The length of array accum_k is nV, the number of nodes
    for (i=0 ; i<nV ; i++) accum_k[i] = 0; // initialization
    int *E; // edge list
    E = (int *)malloc(4*nE*sizeof(int)); // allocate memory
    for (i=0 ; i<nE ; i++) {
        fin >> vs >> ve >> dummy; // read an edge from the input file. Ignore the third column, which is the edge weight.
        vs--; ve--; // node index of the input file starts from 1. So, with this, the node index starts from 0 to nV-1
        E[4*i+3]=E[4*i]=vs;
        E[4*i+2]=E[4*i+1]=ve; // a directed edge from vs to ve and a directed edge from ve to vs. If the original network is a directed network, users should modify these lines
        accum_k[vs]++; // degree of node vs. Undirected graph assumed
        accum_k[ve]++; // degree of node ve. Undirected graph assumed
    }
    fin.close();
    // loading of network data done

    // preprocessing of the edge list
    for (i=1 ; i<nV ; i++) accum_k[i] += accum_k[i-1]; // accum_k[i] is the sum of degrees from node 0 to node i
    qsort(E, 2*nE, 2*sizeof(int), compare_sir_net); // Sort the edge list according to the node index in the first column.
    // nE is the number of undirected edges. So, 2*nE is the number of directed edges.
    for (i=0 ; i<2*nE ; i++) 
        E[i] = E[2*i+1]; // reuse E[]. Because accum_k has the information about the degree of all nodes, we do not need to keep the original E[0], ..., E[2*nE-1]
    // preprocessing of the edge list done


    if (trials >= 2)
        cout << "# lam decay-rate corr number-of-data-points" << endl;
    double t, t_prev; // time
    double t_max = 50000;
    int st[nV]; // state of the nodes. 0: susceptible, 1: infected
    int Inow[nV]; // index of infected nodes
    int accum_kInow[nV]; // accumated deg of infected nodes
    int kInow;
    int upper, lower; // working variable
    int ind, tr; // counters
    // If one wants to do x times for each index patient, set trials = x*n;
    int tmp_max_ind; // working variable
    double ra1; // random variable
    double rate; // total state-transition rate 
    double skip = 1.0; // time spacing for measuring the fraction of infectious nodes
    int nI; // number of infected nodes
    double *p_I; // fraction of infectious nodes
    p_I = (double *)malloc((int)(t_max/skip)*sizeof(double)); // allocate memory
    for (i=1 ; i < (int)(t_max/skip) ; i++) // will set p_I[0] below
        p_I[i] = 0.0; // p_I[i] = fraction of infected nodes at t = i * skip

    // for calculating the decay rate by linear regression after running all the SIS simulations
    int nX; // # samples used for linear regression to calculate the decay rate
    double tmp[5];
    double corr, slope, intercept;

    for (ind=0 ; ind<lam_samples ; ind++) { // for each infection rate value

        lam =(lam_samples>=2)?
            lam_min + (lam_max-lam_min)*(double)ind/(lam_samples-1) : lam; // set the infection rate value
        cerr << "lam = " << lam;

        for (tr=0 ; tr<trials ; tr++) { // for each simulation for a fixed infection rate value

            // initialization 
            nI=nV; // # infected nodes
            for (i=0 ; i<nV ; i++) {
                st[i] = 1; // initially all nodes are infectious
                Inow[i] = i; // list of the index of the infectious nodes
                accum_kInow[i] = accum_k[i]; // sum of the degree of infectious nodes
            }
            t=0; // initialize the time

            // dynamics
            while (nI>0 && t < t_max) { // while at least one infectious node exists or max time is not reached

            	rate = lam*accum_kInow[nI-1] + nI; // total event rate
                t_prev = t;
            	t += -1.0/rate * log ((genrand_int32()+0.5)/4294967296.0); // increment time according to the Gillespie's direct algorithm
                if ((int)(t_prev/skip) < (int)(t/skip)) { // check the status of the dynamics when t is a multiple of "skip"
                    if ((int)(t/skip) >= (int)(t_max/skip)) // t > t_max reached
                        tmp_max_ind = (int)(t_max/skip)-1;
                    else
                        tmp_max_ind = (int)(t/skip);   
                    for (i = (int)(t_prev/skip) + 1 ; i <= tmp_max_ind ; i++)
                        p_I[i] += nI; // p_I[i] = fraction of infectious nodes at time t = i * skip
                }

                // determine the node to be updated
                ra1 = (genrand_int32()+0.5)/4294967296.0 * rate; // ra1 is uniformly distributed on (0, rate)
                if (ra1<nI) { // recovery
                    ve = (int)ra1; // 0 <= i < nI
                    if (st[Inow[ve]] != 1) {
                	    cerr << "sis-net-decay: st[ve]=1 violated" << endl; // error caught
                	    exit(8);
                	  }
                    st[Inow[ve]] = 0; // I -> S
                    kInow = (ve>=1)? accum_kInow[ve]-accum_kInow[ve-1] : accum_kInow[0];
                    for (j=ve ; j<nI-1 ; j++) {
                        accum_kInow[j] = accum_kInow[j+1]-kInow; 
                        Inow[j]=Inow[j+1];
                    }
                    nI--;
            	} else { // infection
                    ra1 -= nI;
                    // determine a tentative starting I node
                    lower=0;
                    upper=nI-1;
                    vs=upper;
                    while (lower != vs || upper != vs) { // bisection search
                        vs = (lower+upper)/2;
                        if (ra1 > lam*accum_kInow[vs]) lower = vs+1;
                        else upper=vs;
                    }

                    // vs is the tentative starting infectious node
                    kInow = (vs>=1)? accum_kInow[vs]-accum_kInow[vs-1] : accum_kInow[0];
                    ve = E[accum_k[Inow[vs]] - genrand_int32()%kInow - 1];
                    // nonoptimal. In fact, it is better not to generate a random number here.

                    if (st[ve]==0) { // S -> I allowed, infection: vs -> ve
                        Inow[nI]=ve;
                        accum_kInow[nI] = (ve==0)? accum_kInow[nI-1]+accum_k[0] : accum_kInow[nI-1]+accum_k[ve]-accum_k[ve-1];
                        nI++;
                        st[ve]=1;
                    }
                } // a single transition done
            } // one trial completed
        } // all trials completed

        // cauculate and output summary statistics
        p_I[0] = 1.0; // fraction of infectious nodes. Initially, all nodes are assumed to be infectious in each run.
        for (i=1 ; i<(int)(t_max/skip) ; i++) // i indexes the time
            p_I[i] /= nV * trials;

        i = 0;
        while (p_I[i] > 1e-4 && i < (int)(t_max/skip)) // we exclude large time points at which p_I[i] is not more than 10^{-4}
            i++;
        nX = i; // # samples used for linear regression to calculate the decay rate
        if (nX < (int)(t_max/skip))
            cerr << ", new nX = " << nX; // number of samples fed to linear regression

        // calculate the decay rate as a Pearson correlation coefficient (= linear regression)
        for (i=0 ; i<5 ; i++) tmp[i] = 0.0; // working variables
        for (i=0 ; i<nX ; i++) {
            tmp[0] += i; // X
        	tmp[1] += log(p_I[i]); // Y
        	tmp[2] += i * i;
        	tmp[3] += log(p_I[i]) * log(p_I[i]);
        	tmp[4] += i * log(p_I[i]);
        }
        corr= (tmp[4]-tmp[0]*tmp[1]/nX)/
        sqrt((tmp[2]-tmp[0]*tmp[0]/nX) * (tmp[3]-tmp[1]*tmp[1]/nX));
        slope=(tmp[4]-tmp[0]*tmp[1]/nX)/
        (tmp[2]-tmp[0]*tmp[0]/nX);
        intercept = tmp[1]/nX - slope * tmp[0]/nX;

        cerr << ", decay rate = " << - slope * skip << ", corr = " << corr << endl;
        cout << lam << " " << - slope * skip << " " << corr << " " << nX << endl;
    } // one value of lambda completed

    free(accum_k); // free memory
    free(E);
    free(p_I);
    return 0;
}