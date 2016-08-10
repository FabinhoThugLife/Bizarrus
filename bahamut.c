// gcc tese.c -o tese -std=c99 -L/usr/local/lib  -lgsl -lgslcblas  -lm -lmatheval

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <matheval.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


#define N_ELEMENTS 8 //the number of elements in each node of the PPT
#define N_POPULATION 2000 //the number of functions in each population
#define N_GENERATIONS 300 // number of generations
#define BUFFER_SIZE 1024 // size of the buffer for the functions
#define N_DATA 64 // number of observations
#define MAXHEIGHT 4// the max height of the tree representing a function
#define MAXNODES  17 // max number of nodes used
#define LEARNING_RATE 0.1 //the learning rate
#define EPS_PIPE 0.1 // epsilon for the learning factor
#define MUTATION_PARAMETER (1e-6) // mutation parameter
#define MUTATION_RATE 0.1 //mutation rate

// codes for the elements in each ppt node	
#define C_TIMES 0
#define C_DIVIDED 1
#define C_PLUS 2
#define C_MINUS 3
#define C_EXP 4
#define C_LOG 5
#define C_X 6 
#define C_NUMBER 7


struct List{ //this stores the elements of F
	char element[20]; //each string stores a function "writen" as in a text (e.g. exp, log, +, etc)
	struct List *next;
}*elements;
typedef struct List list;

struct PPT_{ //this stores the probabilities in the PPT in a binary tree
	float x[N_ELEMENTS]; // the probabilities of choosing a certain element from the list
	//int key; // this identifies the node
	struct PPT_ *left;
	struct PPT_ *right;
}*PPT;
typedef struct PPT_ ppt;

struct Candidate{ // each of the candidate  functions
  char function[BUFFER_SIZE]; //the function written as a string
  int nodes[MAXHEIGHT]; // this keeps track of which nodes of the PPT were used in which order
};
typedef struct Candidate candidate;

// list functions
void printList(){
	if (elements==NULL){
		printf("nothing!\n");
	} else {
		list *l = elements;
		while(l!=NULL){
			printf("%s\n",l->element);
			l = l->next;
		}
	}
}

void add(char *s){ // adds a string to the list of elements of F
	if (elements==NULL){
		elements = (list *) malloc(sizeof(list));
		strcpy(elements->element,s);
		elements->next = NULL;
	} else {
		list *l = elements;
		list *lastLink = (list*)malloc(sizeof(list));
		strcpy(lastLink->element,s);
		while(l->next!=NULL){
			l = l->next;
		}
		l->next = lastLink;
	}	
}

void get(int i, char* c){ // gets the ith element of the elements of F
  int actual = 0;
  list *l = elements;
  if (elements==NULL){
    printf("Nothing to get from an empty list!\n");
  } else {
    while (actual<i){
      if(l->next!=NULL){
        l = l->next;
        actual++;
      } else {
        printf("No element in the position %d.\n",i);
        exit(-1);
        break;
      }
    }
    strcpy(c,l->element);
  }
}

// chooses an element from a list based on a vector of probabilities
// returns the position of the element (from 0)
int sampleElement(float *probs,char *s,gsl_rng * rng){
  float r = (float) gsl_rng_uniform (rng);
  float acc = 0.0; // accumulator
  int i = -1;
  while(acc<=r && acc<0.99){
    i++;
    acc = acc + probs[i];
  }
  get(i,s);
  return i;
}

// end of list functions

// PPT functions

// recursivelly populates the PPT
// must allocate memory for PPT before using this
void createPPT(ppt* temp,float* x,int height){ //creates the PPT with probabilities given by x
  for(int i = 0;i<N_ELEMENTS;i++){
    (temp->x)[i] = x[i];
  }
  
  // if the tree is not tall enough, go on and add the children
  if(height>0){
  // allocate memory for the children
    temp->left = (ppt*) malloc(sizeof(ppt));
    temp->right = (ppt*) malloc(sizeof(ppt));
    createPPT(temp->left,x,height-1);
    createPPT(temp->right,x,height-1);
  }
  
  // when the height reaches zero you're left with NULL pointers in the children
  // then recursion ends
}

// helps to print a vector of floats
void myprint (float *x){
  for (int i=0;i<N_ELEMENTS;i++){
    printf("%.4f ",x[i]);
  }
  printf(";\n");
}

// prints the probabilities in the PPT using inorder
void printPPT(ppt *p){
  if (p != NULL){
    printf("Probs: ");
    myprint(p->x);
    if(p->right!=NULL)
      printPPT(p->right);
    if(p->left!=NULL)
      printPPT(p->left);
  }
}

// end of PPT functions

// generates a candidate from the PPT and stores on c
// nodesVisited keeps track of how many nodes have been visited and is passed to recursive calls
// each node is registered in nodesUsed as soon as it is selected.
void generateFunction(ppt* p, char c[BUFFER_SIZE], int nodesUsed[MAXNODES], int *nodesVisited, int height, gsl_rng * rng){
	 // height indicates the maximum possible height of the tree function
  	 
	 float r; // random number
	 char aux[5000], auxl[5000], auxr[5000], auxn[100]; // auxiliar string
	 
    // if height is one, then generate a literal or a number
    if(height==0){
    	float p_number, p_literal; // probabilities of choosing a literal or a number
    	float sum = ((p->x)[N_ELEMENTS-1]+(p->x)[N_ELEMENTS-2]);
    	p_number = (p->x)[N_ELEMENTS-1]/sum;
    	p_literal = (p->x)[N_ELEMENTS-2]/sum;
    	
    	r = (float) gsl_rng_uniform(rng);
    	
    	if(r<p_number){ // choose a number
  			r = (float) gsl_rng_uniform(rng); // generates a number
	  		snprintf(auxn,10,"%f",r); // changes the float to string
    		strcpy(c,auxn); 
    		nodesUsed[*nodesVisited] = C_NUMBER;
    		(*nodesVisited)++;
    	} else { // choose an x
    		strcpy(c,"x");
    		nodesUsed[*nodesVisited] = C_X;
    		(*nodesVisited)++;
    	}

    } else { // the height is not one, we will use recursion
    	switch(sampleElement(p->x,aux,rng)){
    		case 0: //*
//    			printf("0\n");
    			nodesUsed[*nodesVisited] = C_TIMES;
    			(*nodesVisited)++;
    			generateFunction(p->left,auxl,nodesUsed,nodesVisited,height-1,rng);
    			strcpy(c,"(");
    			strcat(c,auxl);
    			strcat(c,")*(");
    			generateFunction(p->right,auxr,nodesUsed,nodesVisited,height-1,rng);
    			strcat(c,auxr);
    			strcat(c,")");
    			break;
    		case 1: // /
//    			printf("1\n");
          nodesUsed[*nodesVisited] = C_DIVIDED;
    			(*nodesVisited)++;
    			generateFunction(p->left,auxl,nodesUsed,nodesVisited,height-1,rng);
    			strcpy(c,"(");
    			strcat(c,auxl);
    			strcat(c,")/(");
    			generateFunction(p->right,auxr,nodesUsed,nodesVisited,height-1,rng);
    			strcat(c,auxr);
    			strcat(c,")");
    			break;
    		case 2: //+
//    			printf("2\n");
          nodesUsed[*nodesVisited] = C_PLUS;
    			(*nodesVisited)++;
    			generateFunction(p->left,auxl,nodesUsed,nodesVisited,height-1,rng);
    			strcpy(c,"(");
    			strcat(c,auxl);
    			strcat(c,")+(");
    			generateFunction(p->right,auxr,nodesUsed,nodesVisited,height-1,rng);
    			strcat(c,auxr);
    			strcat(c,")");
    			break;
    		case 3: // -
//    			printf("3\n");
          nodesUsed[*nodesVisited] = C_MINUS;
    			(*nodesVisited)++;
    			generateFunction(p->left,auxl,nodesUsed,nodesVisited,height-1,rng);
    			strcpy(c,"(");
    			strcat(c,auxl);
    			strcat(c,")-(");
     			generateFunction(p->right,auxr,nodesUsed,nodesVisited,height-1,rng);
    			strcat(c,auxr);
    			strcat(c,")");
    			break;
    		case 4: //exp
//    			printf("4\n");
          nodesUsed[*nodesVisited] = C_EXP;
    			(*nodesVisited)++;
    			generateFunction(p->left,auxl,nodesUsed,nodesVisited,height-1,rng);
    			strcpy(c,"exp(");
    			strcat(c,auxl);
    			strcat(c,")");
    			break;
    		case 5: //log
//    			printf("5\n");
          nodesUsed[*nodesVisited] = C_LOG;
    			(*nodesVisited)++;
    			generateFunction(p->left,auxl,nodesUsed,nodesVisited,height-1,rng);
    			strcpy(c,"log(");
    			strcat(c,auxl);
    			strcat(c,")");
    			break;
    		case 6: // x
//    			printf("6\n");
          nodesUsed[*nodesVisited] = C_X;
    			(*nodesVisited)++;
    			strcpy(c,"x");
    			break;
    		case 7: // numero
 //   			printf("0\n");
           nodesUsed[*nodesVisited] = C_NUMBER;
    			(*nodesVisited)++;
    			r = (float) gsl_rng_uniform(rng); // generates a number
				snprintf(auxn,10,"%f",r); // changes the float to string
	    		strcpy(c,auxn); 
    			break;
    	} // end of switch-case 
    }//end of if-else
}

// generates n functions with a given height
void generateFunctions(ppt *p, int n, char functions[N_POPULATION][BUFFER_SIZE], int nodesUsed[N_POPULATION][MAXNODES], int height, gsl_rng *rng){
	for(int i = 0;i<n;i++){
		int nodesVisited = 0;
		generateFunction(p,functions[i],nodesUsed[i],&nodesVisited,height,rng);
	}
}
	
// converts	several matheval functions to gsl functions
void matheval2gsl(gsl_function F[N_POPULATION], char functions[N_POPULATION][BUFFER_SIZE]){
	void *f; //for the matheval to create a function
	for(int i =0;i<N_POPULATION;i++){
		f = evaluator_create (functions[i]);
	   assert (f); // functions exists now
		
		// create a function
		double function (double x, void *f){
			return(evaluator_evaluate_x(f,x));
		}
		// and place t in a gsl_function structure
		F[i].function = function;
		F[i].params = f;
	}
}

// calculates the fit of a candidate function 
float getFit(gsl_function F, float *data, float *dataF){
  float candidateF[N_DATA];
  float fit =0;
  for(int i=0;i<N_DATA;i++){
    candidateF[i]=GSL_FN_EVAL(&F,data[i]);
    fit = fit + (candidateF[i]-dataF[i])*(candidateF[i]-dataF[i]);
  }
  if(!gsl_isnan(fit)){
    return fit;
  }else{
    return GSL_POSINF;
  }
}

// binds the functions as strings, the functions as gsl_functions and the list of nodes used 
// as an array of cadidate structs
void bindFunctionAndNodes(candidate candidates[N_POPULATION], char c[N_POPULATION][BUFFER_SIZE], 
                          gsl_function F[N_POPULATION], int nodesUsed[N_POPULATION][MAXNODES])
{
  for (int i = 0; i<N_POPULATION; i++)  {
    //TODO
  }
}

int myCompare(const void *a, const void *b){
  return (int)(*(float*)a - *(float*)b);
}

void evaluateFit(float fit[N_POPULATION], gsl_function F[N_POPULATION], float data[N_DATA], float dataF[N_DATA]){
	for(int i =0;i<N_POPULATION;i++)
	  fit[i] = getFit(F[i],data,dataF);
 }
 
// finds the index of the best fit
int findBestFit (float fit[N_POPULATION]){
  int bIndex = 0;
  for(int i=1; i<N_POPULATION; i++){
    if(fit[i]<fit[bIndex])
      bIndex=i;
  }
  return bIndex;
}


// calculates the probability of a program
float getProbProgram(ppt* p, int nodes[MAXNODES],int* currentIndex){
	float prob; // the probability
	prob = (p->x)[nodes[*currentIndex]];
	// if i'm beyond the last node, stop.
	if(*currentIndex==MAXNODES-1)
		return(prob);
	
	// if the current function from this node was  *,/,+ or - then get the prob of the right and left functions
	if(nodes[*currentIndex]<4){
		*currentIndex = *currentIndex + 1;
		prob = prob*getProbProgram(p->left,nodes,currentIndex);
		*currentIndex = *currentIndex + 1;
		prob = prob*getProbProgram(p->right,nodes,currentIndex);
	} else if(nodes[*currentIndex]==4||nodes[*currentIndex]==5){//if i got an exp or a log
		// then get the probability of the only child
		*currentIndex = *currentIndex + 1;
		prob = prob*getProbProgram(p->left,nodes,currentIndex);
	}
	return prob;	
}

// adapts the ppt to raise the probability of getting the best program
void adapt_PPT_towards2(ppt* p, int prog_b_nodes[MAXNODES], float p_prog_b, float p_target, int* currentIndex){
  // raise the probability for this node's selected element
	(p->x)[prog_b_nodes[*currentIndex]] += LEARNING_RATE*(1-(p->x)[prog_b_nodes[*currentIndex]]);
	// if i'm beyond the last node, stop.
	if(*currentIndex==MAXNODES)
		return;
	
	// if the current function from this node was  *,/,+ or - then get the prob of the right and left functions
	if(prog_b_nodes[*currentIndex]<4){
		*currentIndex = *currentIndex + 1;
		adapt_PPT_towards2(p->left,prog_b_nodes,p_prog_b,p_target,currentIndex);
		*currentIndex = *currentIndex + 1;
		adapt_PPT_towards2(p->right,prog_b_nodes,p_prog_b,p_target,currentIndex);
	} else if(prog_b_nodes[*currentIndex]==4||prog_b_nodes[*currentIndex]==5){//if i got an exp or a log
		// then get the probability of the only child
		*currentIndex = *currentIndex + 1;
		adapt_PPT_towards2(p->left,prog_b_nodes,p_prog_b,p_target,currentIndex);
	}
  return;
}

// adapts the ppt to raise the probability of getting the best program
void adapt_PPT_towards(int prog_b_nodes[MAXNODES], float p_prog_b, float p_target){
  int currentIndex;
  while(p_prog_b<p_target){
    currentIndex=0;
    adapt_PPT_towards2(PPT,prog_b_nodes,p_prog_b,p_target,&currentIndex);
    currentIndex=0;
    p_prog_b = getProbProgram(PPT, prog_b_nodes,&currentIndex);
  }
}

// mutates the PPT
void mutate_PPT(ppt* p, float p_prog_b, gsl_rng *rng){
  float mutation_prob = MUTATION_PARAMETER/sqrt(p_prog_b);
  float r;
  for (int i=0; i<N_ELEMENTS;i++){
    r = gsl_rng_uniform(rng);
    if(r<mutation_prob){
      (p->x)[i] += MUTATION_RATE*(1-(p->x)[i]);
    }
  }
  if(p->left!=NULL) 
    mutate_PPT(p->left,p_prog_b,rng);
  
  if(p->right!=NULL)
    mutate_PPT(p->right,p_prog_b,rng);
  
  return;
}

// keeps the sum of the probabilities equal to 1
void normalize_PPT(ppt* p){
  float sum=0.0;
  for (int i=0; i<N_ELEMENTS;i++){
    sum+=(p->x)[i];
  }
  for (int i=0; i<N_ELEMENTS;i++){
    (p->x)[i] /= sum;
  }
  if(p->left!=NULL) 
    normalize_PPT(p->left);
  if(p->right!=NULL)
    normalize_PPT(p->right);
  return;
}

int main (void){
  // setting up the random number generator
	const gsl_rng_type * T;
  gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
  r = gsl_rng_alloc (T);
	gsl_rng_set (r, 24895);
	
	// the probabilities for selecting each node
  float x[8]={0.1,0.2,0.1,0.2,0.1,0.1,0.1,0.1};

	// the elements of the list of node labels
	add("*");
	add("/");
	add("+");
	add("-");
	add("exp");
	add("log");
	add("x");
	add("number");

	// It is necessary to allocate space for the PPT before creating it
	PPT = (ppt*)malloc(sizeof(ppt));
	
	// creating the PPT
	createPPT(PPT,x,MAXHEIGHT);

  // the data, sorted
	float data[N_DATA] = { 7, 8,   8,   8,   8 ,  8,   9 ,  9,  10,  10,  10,  10,  11,  11,  12,  14,  15,  16,  17,  17,  17,  18,  18,  18,  18, 21,  25,  25,  26,  27,  27,  28,  28,  29,  29,  34,  35,  36,  36,  37,  40,  42 , 47,  51,  54,  55,  56,  60 , 60,  61, 61,  68,  69,  73,  77,  82,  83,  83,  87,  89,  91,  95, 146, 169};
	// the empirical distribution of the data
	float dataF[N_DATA] = {0.015625,0.093750,0.093750,0.093750,0.093750,0.093750,0.125000,0.125000,0.187500,0.187500,0.187500,0.187500,0.218750,0.218750,0.234375,0.250000,0.265625,0.281250,0.328125,0.328125,0.328125,0.390625,0.390625,0.390625,0.390625,0.406250,0.437500,0.437500,0.453125,0.484375,0.484375,0.515625,0.515625,
0.546875,0.546875,0.562500,0.578125,0.609375,0.609375,0.625000,0.640625,0.656250,0.671875,0.687500,0.703125,0.718750,0.734375,0.765625,0.765625,0.796875,0.796875,0.812500,0.828125,0.843750,0.859375,0.875000,0.906250,0.906250,0.921875,0.937500,0.953125,0.968750,0.984375,1.000000};
	
	// these variables represent the best program across generations (the elite)
	char functionEL[BUFFER_SIZE]; // its expression as a string
	int nodesEL[MAXNODES]; // which elements from the nodes were used to create it
	gsl_function FEL; // its expression as a gsl function
	float fitEL = GSL_POSINF; // its fit value
	float p_prog_el =0.0;

	// repeat for every generation
	for(int generation=0; generation<N_GENERATIONS; generation++){
	  // setup the describers of the programs in this generation

  	// functions generate by the PPT for this generation
  	char functions[N_POPULATION][BUFFER_SIZE]; 
  	
  	// the nodes used in each function of this generation
  	int nodesUsed[N_POPULATION][MAXNODES];
    // -1 in a node slot means it was not used.
  	for(int i = 0; i< N_POPULATION;i++){
  	  for (int j = 0; j<MAXNODES; j++)
  	    nodesUsed[i][j] = -1;
  	}
  	
    // functions created before but now as gsl_function structures
  	gsl_function F[N_POPULATION];
    
    // the fit values for this generation
    float fit[N_POPULATION];

  	// the index of the best program at this generation
  	int bIndex;
	
  	// end of setting up the describers for the programs 	  

  	// generate the mathematical expressions for the functions as strings
  	generateFunctions(PPT,N_POPULATION,functions,nodesUsed,MAXHEIGHT,r);
  	
  	// generate the gsl_functions from the strings
	  matheval2gsl(F,functions);
	
	  // evaluate the fitness of each individual in this generation
	  evaluateFit(fit, F, data, dataF);
	
	  // chose the most fit in this generation
	  bIndex = findBestFit(fit);
	  
	  // this is used to keep track of which node we are visiting during the adaptation phase
	  int currentIndex =0;
    
    // get the probability of getting the best program from the PPT
    float p_prog_b = getProbProgram(PPT, nodesUsed[bIndex], &currentIndex);
    
    // if a better than the elite program is found in this generation, store it
    if(fit[bIndex]<fitEL){
      p_prog_el = p_prog_b;
      strcpy(functionEL,functions[bIndex]);
      for(int i=0; i<MAXNODES;i++)
        nodesEL[i] = nodesUsed[bIndex][i];
      FEL = F[bIndex];
      fitEL = fit[bIndex];
    }
    float p_target = p_prog_b + (1-p_prog_b)*LEARNING_RATE*((EPS_PIPE+p_prog_el)/(EPS_PIPE+p_prog_b));
	  
  	// adapts the ppt to raise the probability of getting the best program
    adapt_PPT_towards(nodesUsed[bIndex], p_prog_b, p_target);

    // mutate the PPT
    mutate_PPT(PPT,p_prog_b,r);

    // normalize the PPT
    normalize_PPT(PPT);	
    
    printf(".");
	} // end of generation loop
		
	void* f = evaluator_create(functionEL);
  void* f_prim = evaluator_derivative_x (f);
	printf("\nThe best fit was for the program ");
	for(int i = 0; i< MAXNODES;i++)
	  printf("%d",nodesEL[i]);
	printf(" and it was %.4f.\n",fitEL);
	printf(" %s with probability %.12f.\n",functionEL,p_prog_el);
	printf(" Its density is %s.\n", evaluator_get_string (f_prim));
	
		
	gsl_rng_free(r);
}
