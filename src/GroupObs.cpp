/*includes*/
#include <mex.h>
#include <vector>
#include <iostream>
using namespace std;
#include <cmath>
#include <cassert>

/*MACROS {{{*/
/* 
 * 
 *    J    j
 *    ^    ^
 *    |    | +--------+--------+
 *    |    | |        |        |
 * 1X |    | |   2    |   3    |
 *    |    | |        |        |
 *    |    | +--------+--------+
 *    |    | |        |        |
 * 0X |    | |   0    |   1    |
 *    |    | |        |        |
 *    |    | +--------+--------+
 *    |    +-----------------------> i
 *    |         
 *    |----------------------------> I
 *              X0        X1  
 *
 * box 0 -> I=0 J=0 IJ=00  = 0
 * box 1 -> I=1 J=0 IJ=01  = 1
 * box 2 -> I=0 J=1 IJ=10  = 2
 * box 3 -> I=1 J=1 IJ=11  = 3
 */
//IJ(i,j,l) returns the box number of i and j with respect to l
//if !j&l and !i&l -> 0 (box zero: lower left )
//if !j&l and  i&l -> 1 (box one:  lower right)
//if  j&l and !i&l -> 2 (box two:  upper left )
//if  j&l and  i&l -> 3 (box three:upper right)
#define IJ(i,j,l)  ((j&l) ? ((i&l) ? 3:2 ) :((i&l) ? 1:0 ))
/*}}}*/
/*Inputs{{{*/
#define DATAX   (mxArray*)prhs[0]
#define DATAY   (mxArray*)prhs[1]
#define DATA    (mxArray*)prhs[2]
#define MINSPACING (mxArray*)prhs[3]
/*}}}*/
/*Outputs{{{*/
#define INDICES (mxArray**)&plhs[0]
/*}}}*/
/*Prototypes{{{*/
void  FetchData(double* pdouble,const mxArray* dataref);
void  FetchVectorPointer(double** pvector,int *pN,const mxArray* dataref);
void  WriteVector(mxArray** pdataref,double* vector,int N);
class Object{/*{{{*/
	public: 
		virtual       ~Object() {};
};/*}}}*/
class DataSet{/*{{{*/
	private:
		std::vector<Object*> objects;
	public:
		DataSet();
		DataSet(double* observations_list,double* x,double* y,int n);
		~DataSet();
		int AddObject(Object* obj);
		int   Size();
		Object* GetObjectByOffset(int offset);
};
/*}}}*/
class Observation: public Object{/*{{{*/
	public:
		double x,y;
		int    xi,yi;
		int    index;
		int    sindex;
		double value;
		Observation(double x_in,double y_in,int xi_in,int yi_in,int index_in,int sindex_in,double value_in);
		~Observation();
};
class Observations;
/*}}}*/
class Quadtree{/*{{{*/
	private:
		class QuadtreeBox: public Object{ 
			public:
				int    nbitems; // number of current vertices in the box
				double xcenter; // x position of the center (double)
				double ycenter; // x position of the center (double)
				double length;  // width of the box
				union{
					QuadtreeBox *box[4];
					Observation *obs[4];
				};
				int  IsWithinRange(double  x,double y,double range);
				void RangeSearch(int *indices,int *pnobs,double x,double y,double range);
				void WriteObservations(int *indices,int *pnobs);
		};
		DataSet* boxcontainer;
	public:
		int          MaxDepth;          // maximum number of subdivision
		QuadtreeBox *root;              // main box
		int          NbQuadtreeBox;     // total number of boxes
		int          NbObs;             // number of points
		Quadtree();
		Quadtree(double xmin,double xmax,double ymin,double ymax,int maxdepth_in);
		~Quadtree();
		void         Add(Observation *observation);
		void         IntergerCoordinates(int *xi,int *yi,double x,double y);
		QuadtreeBox *NewQuadtreeBox(QuadtreeBox *master,int index);
		QuadtreeBox *NewQuadtreeBox(double xcenter,double ycenter,double length);
		void         QuadtreeDepth2(int *A,int xi,int yi);
		void         RangeSearch(int **pindices,int *pnobs,double x,double y,double range);
};
/*}}}*/
class Observations: public DataSet{/*{{{*/
	private:
		Quadtree* quadtree;
	public:
		Observations();
		Observations(double* observations_list,double* x,double* y,int n,double minspacing);
		~Observations();
		void GetIndices(double** pindices,int *pnum);
};
/*}}}*/
/*}}}*/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){

	double *datax   = NULL;
	double *datay   = NULL;
	double *data    = NULL;
	double  minspacing;
	int     dataN,outN;
	double *indices = NULL;
	int     test,num;

	/*Check arguments to avoid crash*/
	if(nlhs>1 || nrhs!=4) mexErrMsgTxt("Wrong usage (indices=GroupObs(x,y,data,spacing);)");

	/*Get variables from matlab to C*/
	FetchVectorPointer(&datax,&dataN,DATAX);
	FetchVectorPointer(&datay,&test, DATAY); if(test!=dataN)   mexErrMsgTxt("x and y do not have the same size");
	FetchVectorPointer(&data, &test, DATA);  if(test!=dataN)   mexErrMsgTxt("x and data do not have the same size");
	FetchData(&minspacing,MINSPACING);       if(minspacing<=0) mexErrMsgTxt("minspacing must be positive");

	/*Process observation dataset*/
	Observations *observations=new Observations(data,datax,datay,dataN,minspacing);
	observations->GetIndices(&indices,&num);

	/*Write output vector*/
	WriteVector(INDICES,indices,num);

	/*Clean-up and return*/
	delete observations;
	/*Do not erase pointers!*/
	return;
}

/*FUNCTION DataSet::DataSet(){{{*/
DataSet::DataSet(){
	return;
}
/*}}}*/
/*FUNCTION DataSet::~DataSet{{{*/
DataSet::~DataSet(){

	/*  use reverse_iterator for efficiency in matlab memory manager
		 (keeping old code in case it needs to revert back)  */

	//	vector<Object*>::iterator object;
	vector<Object*>::reverse_iterator object;

	//	for ( object=objects.begin() ; object < objects.end(); object++ ){
	//		delete (*object);
	//	}
	for ( object=objects.rbegin() ; object < objects.rend(); object++ ){
		delete (*object);
	}
}
/*}}}*/
/*FUNCTION DataSet::AddObject{{{*/
int  DataSet::AddObject(Object* obj){

	objects.push_back(obj);

	return 1;
}
/*}}}*/
/*FUNCTION DataSet::GetObjectByOffset{{{*/
Object* DataSet::GetObjectByOffset(int offset){

	/*Check index in debugging mode*/
	assert(this!=NULL);
	assert(offset<this->Size());

	return objects[offset];

}
/*}}}*/
/*FUNCTION DataSet::Size{{{*/
int  DataSet::Size(void){
	assert(this!=NULL);
	return objects.size();
}
/*}}}*/
/*FUNCTION Observations::Observations{{{*/
Observations::Observations(double* observations_list,double* x,double* y,int n,double minspacing){

	/*Intermediaries*/
	int          i,j,maxdepth,level,counter;
	int          xi,yi;
	int          num;
	double       xmin,xmax,ymin,ymax;
	double       offset,minlength;
	int         *indices     = NULL;
	Observation *observation = NULL;

	/*Get extrema*/
	if(!x || !y) mexErrMsgTxt("x and y are empty");
	xmin=x[0]; ymin=y[0];
	xmax=x[0]; ymax=y[0];
	for(i=1;i<n;i++){
		if(x[i]<xmin) xmin=x[i];
		if(y[i]<ymin) ymin=y[i];
		if(x[i]>xmax) xmax=x[i];
		if(y[i]>ymax) ymax=y[i];
	}
	offset=0.05*(xmax-xmin); xmin-=offset; xmax+=offset;
	offset=0.05*(ymax-ymin); ymin-=offset; ymax+=offset;

	/*Get Minimum box size*/
	maxdepth = 20;
	minlength=(xmax-xmin)/double((1L<<maxdepth)-1);

	/*Initialize Quadtree*/
	cout << "Generating quadtree with a maximum box size " << minlength << " (depth=" << maxdepth << ")... " << endl;
	this->quadtree = new Quadtree(xmin,xmax,ymin,ymax,maxdepth);

	/*Add observations one by one*/
	counter = 0;
	for(i=0;i<n;i++){

		/*check that this observation is not too close from another one*/
		this->quadtree->RangeSearch(&indices,&num,x[i],y[i],minspacing);
		if(!num){
			this->quadtree->IntergerCoordinates(&xi,&yi,x[i],y[i]);
			this->quadtree->QuadtreeDepth2(&level,xi,yi);
			if((int)level <= maxdepth){
				observation = new Observation(x[i],y[i],xi,yi,counter++,i+1,observations_list[i]);
				this->quadtree->Add(observation);
				this->AddObject(observation);
			}
			else{
				printf("obs number %i : %i > %i (minspacing = %g)\n",i+1,level,maxdepth,minspacing);
				mexErrMsgTxt("stop (inconsistent depths)");
			}
		}
		delete [] indices;
	}
	cout << "done" << endl;
	cout << "Initial number of observations: " << n << endl;
	cout << "  Final number of observations: " << this->quadtree->NbObs << endl;
}
/*}}}*/
/*FUNCTION Observations::~Observations(){{{*/
Observations::~Observations(){
	delete quadtree;
	return;
}
/*}}}*/
/*FUNCTION Observations::GetIndices{{{*/
void Observations::GetIndices(double** pindices,int *pnum){

	/*Allocate output*/
	double *indices=(double*)mxMalloc(this->Size()*sizeof(double));

	for(int i=0;i<this->Size();i++){
		Observation *observation=(Observation*)this->GetObjectByOffset(i);
		indices[i]=observation->sindex;
	}

	*pindices=indices;
	*pnum = this->Size();
}
/*}}}*/
/*FUNCTION Observation::Observation{{{*/
Observation::Observation(double x_in,double y_in,int xi_in,int yi_in,int index_in,int sindex_in,double value_in){

	this->x      = x_in;
	this->y      = y_in;
	this->xi     = xi_in;
	this->yi     = yi_in;
	this->index  = index_in;
	this->sindex  = sindex_in;
	this->value  = value_in;

}
/*}}}*/
/*FUNCTION Observation::~Observation(){{{*/
Observation::~Observation(){
	return;
}
/*}}}*/
/*FUNCTION Quadtree::Quadtree{{{*/
Quadtree::Quadtree(double xmin,double xmax,double ymin,double ymax,int maxdepth){

	/*Intermediaries*/
	double length;

	/*Initialize fields*/
	this->MaxDepth=maxdepth;
	this->NbQuadtreeBox=0;
	this->NbObs=0;

	/*Create container*/
	this->boxcontainer=new DataSet();

	/*Create Root, pointer toward the main box*/
	length=max(xmax-xmin,ymax-ymin);
	this->root=NewQuadtreeBox(xmin+length/2,ymin+length/2,length);
}
/*}}}*/
/*FUNCTION Quadtree::~Quadtree(){{{*/
Quadtree::~Quadtree(){

	delete boxcontainer;
	root=NULL;

}
/*}}}*/
/*FUNCTION Quadtree::Add{{{*/
void  Quadtree::Add(Observation* observation){

	/*Intermediaries*/
	int          xi,yi,ij,level,levelbin;
	QuadtreeBox **pbox    = NULL; // pointer toward current box b
	QuadtreeBox **pmaster = NULL; // pointer toward master of b
	QuadtreeBox  *box     = NULL; // current box b
	QuadtreeBox  *slave   = NULL; // suslaveox of b (if necessary)
	Observation  *obs[4];

	/*Get integer coodinates*/
	xi = observation->xi;
	yi = observation->yi;

	/*Initialize levels*/
	level    = 0;
	levelbin = (1L<<this->MaxDepth);// = 2^30

	/*Get inital box (the largest)*/
	pmaster = &root;
	pbox    = &root;

	/*Find the smallest box where the observation is located*/
	while((box=*pbox) && (box->nbitems<0)){ 

		/*Go down one level (levelbin = 00100 -> 00010)*/
		levelbin>>=1; level+=1; assert(level<this->MaxDepth);

		/*Get next box according to the bit value (levelbin)*/
		pmaster = pbox;
		pbox    = &box->box[IJ(xi,yi,levelbin)];
	}
	assert(levelbin>0);

	/*Now, try to add the vertex, if the box is full (nbitems=4), we have to divide it in 4 new boxes*/
	while((box=*pbox) && (box->nbitems==4)){

		/*Copy the 4 observation in the current Quadtreebox*/
		obs[0] = box->obs[0];
		obs[1] = box->obs[1];
		obs[2] = box->obs[2];
		obs[3] = box->obs[3];

		/*set nbitems as -1 (now holding boxes instead of observations)*/
		box->nbitems = -1;
		box->box[0]  = NULL;
		box->box[1]  = NULL;
		box->box[2]  = NULL;
		box->box[3]  = NULL;

		/*Go down one level (levelbin = 00010 -> 00001)*/
		levelbin>>=1; level+=1; assert(level<this->MaxDepth);

		/*Put the four observations in the new boxes*/
		for (int k=0;k<4;k++){

			/*Get box for observation number k*/
			ij    = IJ(obs[k]->xi,obs[k]->yi,levelbin);
			slave = box->box[ij];
			if(!slave){
				box->box[ij] = NewQuadtreeBox(box,ij);
				slave        = box->box[ij];
			}
			slave->obs[slave->nbitems++] = obs[k];
		}

		/*Get the suslaveox where the current observation is located*/
		ij      = IJ(xi,yi,levelbin);
		pmaster = pbox;
		pbox    = &box->box[ij];
	}

	/*alloc the QuadtreeBox if necessary and add current observation*/
	box = *pbox;
	if(!box){
		ij  = IJ(xi,yi,levelbin);
		box = *pbox = NewQuadtreeBox(*pmaster,ij);
	}
	box->obs[box->nbitems++]=observation;
	NbObs++;

}/*}}}*/
/*FUNCTION Quadtree::IntergerCoordinates{{{*/
void  Quadtree::IntergerCoordinates(int *xi,int *yi,double x,double y){

	/*Intermediaries*/
	double coefficient;
	double xmin,ymin;

	/*Checks in debugging mode*/
	assert(xi && yi);
	assert(this->root);

	/*coeffIcoor is the coefficient used for integer coordinates:
	 *                (x-xmin)
	 * xi = (2^30 -1) --------- 
	 *                 length
	 * coefficient = (2^30 -1)/length
	 */
	coefficient = double((1L<<this->MaxDepth)-1)/(this->root->length);
	xmin        = this->root->xcenter - this->root->length/2;
	ymin        = this->root->ycenter - this->root->length/2;

	*xi=int(coefficient*(x - xmin));
	*yi=int(coefficient*(y - ymin));
}/*}}}*/
/*FUNCTION Quadtree::NewQuadtreeBox(QuadtreeBox* master,int index) {{{*/
Quadtree::QuadtreeBox* Quadtree::NewQuadtreeBox(QuadtreeBox* master,int index){

	/*Output*/
	QuadtreeBox* newbox=NULL;

	/*Checks in debugging mode*/
	assert(master);

	/*Create and initialize a new box*/
	newbox=new QuadtreeBox();
	newbox->nbitems=0;
	newbox->box[0]=NULL;
	newbox->box[1]=NULL;
	newbox->box[2]=NULL;
	newbox->box[3]=NULL;
	switch(index){
		case 0:
			newbox->xcenter=master->xcenter - master->length/4;
			newbox->ycenter=master->ycenter - master->length/4;
			break;
		case 1:
			newbox->xcenter=master->xcenter + master->length/4;
			newbox->ycenter=master->ycenter - master->length/4;
			break;
		case 2:
			newbox->xcenter=master->xcenter - master->length/4;
			newbox->ycenter=master->ycenter + master->length/4;
			break;
		case 3:
			newbox->xcenter=master->xcenter + master->length/4;
			newbox->ycenter=master->ycenter + master->length/4;
			break;
		default:
			mexErrMsgTxt("Case not supported");
	}
	newbox->length=master->length/2;

	/*Add to container*/
	this->boxcontainer->AddObject(newbox);
	NbQuadtreeBox++;

	/*currentbox now points toward next quadtree box*/
	return newbox;
}/*}}}*/
/*FUNCTION Quadtree::NewQuadtreeBox(double xcenter,double ycenter,double length){{{*/
Quadtree::QuadtreeBox* Quadtree::NewQuadtreeBox(double xcenter,double ycenter,double length){

	/*Output*/
	QuadtreeBox* newbox=NULL;

	/*Create and initialize a new box*/
	newbox=new QuadtreeBox();
	newbox->nbitems=0;
	newbox->xcenter=xcenter;
	newbox->ycenter=ycenter;
	newbox->length=length;
	newbox->box[0]=NULL;
	newbox->box[1]=NULL;
	newbox->box[2]=NULL;
	newbox->box[3]=NULL;

	/*Add to container*/
	this->boxcontainer->AddObject(newbox);
	NbQuadtreeBox++;

	/*currentbox now points toward next quadtree box*/
	return newbox;
}/*}}}*/
/*FUNCTION Quadtree::QuadtreeDepth2{{{*/
void Quadtree::QuadtreeDepth2(int* A,int xi,int yi){

	QuadtreeBox **pbox = NULL;
	QuadtreeBox  *box  = NULL;
	int           level,levelbin;

	/*Initialize levels*/
	level    = 0;
	levelbin = (1L<<this->MaxDepth);// = 2^30

	/*Get inital box (the largest)*/
	pbox=&root;

	/*Find the smallest box where this point is located*/
	while((box=*pbox) && (box->nbitems<0)){ 

		levelbin>>=1; level+=1; 

		pbox = &box->box[IJ(xi,yi,levelbin)];
	}
	if(box && box->nbitems>0){
		/*This box is not empty, add one level*/
		level+=1;
	}

	/*If we were to add the vertex, get level*/
	if(box && box->nbitems==4){
		int ij;
		bool flag=true;
		while(flag){

			levelbin>>=1; level+=1;
			if(level>this->MaxDepth){
				level+=1;
				break;
			}

			/*loop over the four observations*/
			ij=IJ(box->obs[0]->xi,box->obs[0]->yi,levelbin);
			for (int k=1;k<4;k++){
				if(IJ(box->obs[k]->xi,box->obs[k]->yi,levelbin) != ij){
					flag = false;
				}
			}
			if(IJ(xi,yi,levelbin)!=ij){
				flag = false;
			}
		}
	}

	*A=level;
}/*}}}*/
/*FUNCTION Quadtree::RangeSearch{{{*/
void Quadtree::RangeSearch(int **pindices,int *pnobs,double x,double y,double range){

	/*Intermediaries*/
	int  nobs;
	int *indices = NULL;

	/*Allocate indices (maximum by default*/
	if(this->NbObs) indices = new int[this->NbObs];
	nobs = 0;

	if(this->root) this->root->RangeSearch(indices,&nobs,x,y,range);

	/*Clean-up and return*/
	*pnobs=nobs;
	*pindices=indices;

}/*}}}*/
/*FUNCTION QuadtreeBox::IsWithinRange{{{*/
int Quadtree::QuadtreeBox::IsWithinRange(double x,double y,double range){

	/*Return 0 if the 2 boxes do not overlap*/
	if(this->xcenter+this->length/2 < x-range) return 0;
	if(this->xcenter-this->length/2 > x+range) return 0;
	if(this->ycenter+this->length/2 < y-range) return 0;
	if(this->ycenter-this->length/2 > y+range) return 0;

	/*Return 2 if the this box is included in the range*/
	if(this->xcenter+this->length/2 <= x+range &&
				this->ycenter+this->length/2 <= y+range &&
				this->xcenter-this->length/2 >= x-range &&
				this->ycenter-this->length/2 >= y-range) return 2;

	/*This is a simple overlap*/
	return 1;

}/*}}}*/
/*FUNCTION QuadtreeBox::RangeSearch{{{*/
void Quadtree::QuadtreeBox::RangeSearch(int* indices,int *pnobs,double x,double y,double range){

	/*Intermediaries*/
	int i,nobs;

	/*Recover current number of observations*/
	nobs = *pnobs;

	switch(this->IsWithinRange(x,y,range)){
		case 0:
			/*If this box is not within range, return*/
			break;
		case 2:
			/*This box is included in range*/
			this->WriteObservations(indices,&nobs);
			break;
		case 1:
			/*This box is partly included*/
			if(this->nbitems>0){
				/*If this box has only observations, add indices that are within range*/
				for(i=0;i<this->nbitems;i++){
					if(fabs(this->obs[i]->x-x) <= range && fabs(this->obs[i]->y-y) <= range){
						indices[nobs++]=this->obs[i]->index;
					}
				}
			}
			else{
				/*This box points toward boxes*/
				if(this->box[0]) this->box[0]->RangeSearch(indices,&nobs,x,y,range);
				if(this->box[1]) this->box[1]->RangeSearch(indices,&nobs,x,y,range);
				if(this->box[2]) this->box[2]->RangeSearch(indices,&nobs,x,y,range);
				if(this->box[3]) this->box[3]->RangeSearch(indices,&nobs,x,y,range);
			}
			break;
		default:
			cout << "Case " << this->IsWithinRange(x,y,range) << " not supported";
			mexErrMsgTxt("");
	}

	/*Assign output pointers: */
	*pnobs=nobs;
}/*}}}*/
/*FUNCTION QuadtreeBox::WriteObservations{{{*/
void Quadtree::QuadtreeBox::WriteObservations(int* indices,int *pnobs){

	/*Intermediaries*/
	int i,nobs;

	/*Recover current number of observations*/
	nobs = *pnobs;

	if(this->nbitems>0){
		/*If this box has only observations, add all indices*/
		for(i=0;i<this->nbitems;i++){
			indices[nobs++]=this->obs[i]->index;
		}
	}
	else{
		/*This box points toward boxes, */
		if(this->box[0]) this->box[0]->WriteObservations(indices,&nobs);
		if(this->box[1]) this->box[1]->WriteObservations(indices,&nobs);
		if(this->box[2]) this->box[2]->WriteObservations(indices,&nobs);
		if(this->box[3]) this->box[3]->WriteObservations(indices,&nobs);
	}

	/*Assign output pointers: */
	*pnobs=nobs;
}/*}}}*/
/*FUNCTION FetchData(double* pscalar,const mxArray* dataref){{{*/
void FetchData(double* pscalar,const mxArray* dataref){

	double scalar;

	if (!mxIsClass(dataref,"double")){
		mexErrMsgTxt("input data_type is not a double!");
	}
	else{
		/*Recover the double: */
		scalar=mxGetScalar(dataref);
	}

	/*Assign output pointers:*/
	*pscalar=scalar;
}
/*}}}*/
/*FetchVectorPointer {{{*/
void FetchVectorPointer(double** pvector,int *pN,const mxArray* dataref){

	double *vector=NULL;
	double *values=NULL;
	int     N;

	if(mxIsEmpty(dataref) ){
		N=0;
		vector=NULL;
	}
	else if (mxIsDouble(dataref) ){
		if(mxGetM(dataref)!=1 && mxGetN(dataref)!=1){
			mexErrMsgTxt("input is a matrix and not a vector");
		}
		N=mxGetN(dataref)*mxGetM(dataref);
		vector=(double*)mxGetPr(dataref);
	}
	else{
		mexErrMsgTxt("vector type not supported");
	}

	*pvector=vector;
	if (pN)*pN=N;
}/*}}}*/
/*WriteVector {{{*/
void WriteVector(mxArray** pdataref,double* vector,int N){

	mxArray* dataref=NULL;

	if(vector){
		/*data is a double* pointer. Copy into a vector: */
		dataref = mxCreateDoubleMatrix(0,0,mxREAL);
		mxSetM(dataref,(mwSize)N);
		mxSetN(dataref,(mwSize)1);
		mxSetPr(dataref,(double*)vector);
	}
	else{
		dataref = mxCreateDoubleScalar(0.0);
	}
	*pdataref=dataref;
}
/*}}}*/
