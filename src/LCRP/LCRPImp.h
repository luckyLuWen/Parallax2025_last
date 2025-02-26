

#pragma once
#include "LCRP.h"
#include <map>
#include <vector>
#include <set>
using namespace std;
using namespace Eigen;

//use sba_crsm structure from SBA (http://www.ics.forth.gr/~lourakis/sba/) to store sparse matrix
struct sba_crsm	
{
    int nr, nc;   //稀疏矩阵的行列
    int nnz;      //非零元素个数
    int *val;     //存储非零元素
    int *colidx;  //非零元素的列号
    int *rowptr;  //指向每行第一个非零元素的索引数组 (size: nr+1)
};

class CLCRP : public ILCRP
{
public:
	CLCRP(void);
	~CLCRP(void);
	virtual bool CLCRP::lcrp_run( bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL, 
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6 );

	virtual	bool lcrp_run( int argc, char** argv );

	virtual bool lcrp_initialize( char* szCamera, char* szFeature, char* szCalib = NULL, char* szXYZ = NULL );

	virtual bool lcrp_motstr_levmar( );

	virtual bool lcrp_motstr_gn(  FixType ft = LCRP_FixDefault  );

private:
	//initialize feature. You can provide xyz or system also provide them by itself;
	bool    lcrp_initializeMainArchor( double* imgpts, double* camera,double* K,double* feature, int nP, int FID, double* KR );
	bool    lcrp_initializeAssoArchor( double* imgpts, int* photo, double* camera,double* K,double* feature,int nMI, int nAI, int FID, bool bLast );
	bool	lcrp_initializeOtheArchors( double* imgpts, int* photo, double* camera,double* K,double* feature,int* archorSort,int nfeacout, int nOI, int FID );
	


	//update KR KdA KdB KdG
	void	lcrp_updateKR( double *R, double *KR, double *KdA, double *KdB, double *KdG, double *K, double *p );


	//compute reprojection error
	void	lcrp_cost(double *p, double *hx, int* archor );

	//compute Jacobian, not save Jp, Jc etc
	void	lcrp_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );
	void	lcrp_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );

	//construct S matrix, S = U - W*V^-1*W^T
	void	lcrp_constructSLM( double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij, double mu );
	void	lcrp_constructSGN( double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij );

	//prepare for ordering of CHOLMOD
	void	lcrp_constructAuxCSSLM( int *Ap, int *Aii );
	void	lcrp_constructAuxCSSGN( int *Ap, int *Aii );

	//set CSS format
	void	lcrp_constructCSSLM( int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init);
	void	lcrp_constructCSSGN( int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init, int nft);

	//solve sparse system using CHOLMOD
	bool	lcrp_solveCholmodLM( int* Ap, int* Aii, bool init, bool ordering);
	bool	lcrp_solveCholmodGN( int* Ap, int* Aii, bool init, bool ordering);

	//compute feature parameters
	void	lcrp_solveFeatures( double *W, double *IV, double *ea, double *eb, double *dpa, double *dpb);

private:
	//compute smask of S matrix
	int		lcrp_ConstructSmask( sba_crsm& Sidxij, sba_crsm& Uidxij );

	//compute inverse V matrix
	void	lcrp_inverseVLM( double* V, double* IV, sba_crsm& Uidxij, double mu );
	void	lcrp_inverseVGN( double* U, double* V, sba_crsm& Uidxij );

	//compute initial mu of LM
	double	lcrp_computeInitialmu( double* U, double* V, sba_crsm& Uidxij, double tau, int nvars );

	//transform angle into XYZ
	int		lcrp_angle2xytGN( double *p );
	int		lcrp_angle2xytLM( double *p );
	void	lcrp_saveXYZ( const char* camera, const char* sz3Dpt, double *p, bool gn = true );


	bool	lcrp_parseArgs( int argc, char* argv[] );
	void	lcrp_printHelp();

private:
	void	lcrp_readAndInitialize( char *camsfname, char *ptsfname, int *ncams, int *n3Dpts, int *n2Dprojs,double **motstruct, 
		double** imgpts, double** lc, int** archor, char** vmask, char** umask, int** nphoto, int** nfeature, int** archorSort);

	void	lcrp_readCameraPoseration( char *fname, double* ical );

	void	lcrp_readProjectionAndInitilizeFeature(	FILE *fp, double *params, double *projs,double *lc, char *vmask, int ncams, 
		int *archor,char* umask,int* nphoto, int* nfeature, int* archorSort );

	//Add By Zuo
	bool lcrp_initializeOtheArchors_Mindw(double* imgpts, int* photo, double* camera, double* K, double* feature, int* archorSort, int nfeacout, int nOI, int FID);
	void lcrp_readProjectionAndTriangulateFeature(FILE* fp, double* projs, int ncams);
	void lcrp_saveTriangulatedxyz(char* sz3Dpt, double* p);
	void lcrp_saveInitialXYZ(char* sz3Dpt, double* p);
	void lcrp_saveInitialParallax(char* sz3Dpt, double* p);
	void lcrp_constructP(double* P, double* K, double* p);
	double* lcrp_angle2xyz(double* p);
	int zu = 0;

	//Add bu Lu
	int savePara = 0;

	void	lcrp_readCameraPose(FILE *fp, double *params);

	void	lcrp_readCablibration(FILE* fp, double *K );

private:
	int		findNcameras(FILE *fp);
	int		countNDoubles(FILE *fp);
	int		readNInts(FILE *fp, int *vals, int nvals);
	int		readNDoubles(FILE *fp, double *vals, int nvals);
	int		skipNDoubles(FILE *fp, int nvals);	
	void	readNpointsAndNprojections(FILE *fp, int *n3Dpts, int pnp, int *nprojs, int mnp);
	double	nrmL2xmy(double *const e, const double *const x, const double *const y, const int n);
	void	readNpointsAndNprojectionsFromProj(FILE *fp, int &n3Dpts, int &nprojs);
	void	readPointProjections(FILE *fp,double *imgpts, int *photo,int* imgptsSum, int n3Dpts, int n2Dprojs );
	void    readImagePts( const char* szProj, double **imgpts, int **photo,int** imgptsSum, int &n3Dpts, int &n2Dprojs );

private:
	int		m_ncams, m_n3Dpts, m_n2Dprojs, m_nS;  //number of camera, 3D points, 2D projection points, non-zero element of S matrix
	int		*m_archor;							  
	int		*m_photo, *m_feature;	  			  
	double* m_motstruct, * m_imgpts, * m_lc;			  //6 camera pose and 3 feature parameters/PBA parameter,
	double	*m_XYZ;								  //initial XYZ provided 	
	double  *m_K;								  //calibration parameters

	char	*m_vmask, *m_umask, *m_smask;		  
	int		*m_imgptsSum, *m_struct, *m_pnt2main, *m_archorSort;
	double  * m_KR, * m_KdA, * m_KdB, * m_KdG, * m_P, * m_R; 

	bool    m_bProvideXYZ, m_bFocal;	
	char*	m_szCameraInit;
	char*	m_szFeatures;
	char*   m_szCalibration;
	char*	m_szXYZ;
	char*   m_szCamePose;
	char*   m_sz3Dpts;
	char*   m_szReport;
	int		m_nMaxIter;
	double  m_Tau, m_e1, m_e2, m_e3, m_e4;
	bool	m_bRobustKernel;
	bool	m_bsolverLM;
	bool	m_bsolverGN;
	int		m_nRobustType;
	double  m_delt;

private:
	//Solve Sparse Matrix using CHOLMOD (http://www.cise.ufl.edu/research/sparse/SuiteSparse/) 
	cholmod_sparse *m_cholSparseS;				
	cholmod_factor *m_cholFactorS; 
	cholmod_common m_cS; 
	cholmod_dense  *m_cholSparseR, *m_cholSparseE;
};

