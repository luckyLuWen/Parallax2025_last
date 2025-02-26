

#include "stdafx.h"
#include "LCRPImp.h"

#define MAXSTRLEN  2048 /* 2K */
#define SKIP_LINE(f){                                                       \
	char buf[MAXSTRLEN];                                                        \
	while(!feof(f))                                                           \
	if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}

//cite from SBA code in order to search sub matrix of S according to (i,j)
static void sba_crsm_alloc(struct sba_crsm *sm, int nr, int nc, int nnz)
{
	int msz;
	sm->nr=nr;
	sm->nc=nc;
	sm->nnz=nnz;
	msz=2*nnz+nr+1;
	sm->val=(int *)malloc(msz*sizeof(int));  /* required memory is allocated in a single step */
	if(!sm->val){
		fprintf(stderr, "memory allocation request failed in sba_crsm_alloc() [nr=%d, nc=%d, nnz=%d]\n", nr, nc, nnz);
		exit(1);
	}
	sm->colidx=sm->val+nnz;
	sm->rowptr=sm->colidx+nnz;
}

static void sba_crsm_free(struct sba_crsm *sm)
{
	 sm->nr=sm->nc=sm->nnz=-1;
	free(sm->val);
	sm->val=sm->colidx=sm->rowptr=NULL;
}


/* returns the index of the (i, j) element. No bounds checking! */
static int sba_crsm_elmidx(struct sba_crsm *sm, int i, int j)//����i��j��Ԫ�ص�������
{
	register int low, high, mid, diff;

	low=sm->rowptr[i];
	high=sm->rowptr[i+1]-1;

	/* binary search for finding the element at column j */
	while(low<=high)
	{
		mid=(low+high)>>1; //(low+high)/2;
		diff=j-sm->colidx[mid];
		if(diff<0)//��j�Ҳ�
			 high=mid-1;
		else if(diff>0)//��j���
			low=mid+1;
		else
		return mid;
	}

	return -1; /* not found */
}


//compute reprojection image coordinates of each image points
static void lcrp_reprojectEachPts( double *R, double *KR, double* pa, double* pb, int nM, int nN, int nP, double lc[1] )
{
	//double *posei, *posek;
	//double pti2k[3];
	//double ptXUnit[3];
	//double dDot, dDisi2k, dW2;
	//double ptXk[3];
	//double *posel;
	//double pti2l[3];
	double ptXj[3], n[3];
	double* pKR, * pR;
	
	//pKR = KR + nP * 9;//nP�����ͼ��KR����
	pR = R + nP * 9;//nP�����ͼ��R����

	ptXj[0] = pb[0] - pa[0];
	ptXj[1] = pb[1] - pa[1];
	ptXj[2] = pb[2] - pa[2];

	n[0] = pR[0] * ptXj[0] + pR[1] * ptXj[1] + pR[2] * ptXj[2];
	n[1] = pR[3] * ptXj[0] + pR[4] * ptXj[1] + pR[5] * ptXj[2];
	n[2] = pR[6] * ptXj[0] + pR[7] * ptXj[1] + pR[8] * ptXj[2];

	lc[0] = atan(sqrt(n[0] * n[0] + n[1] * n[1]) / n[2]);

	//n[0] = (pKR[0] * ptXj[0] + pKR[1] * ptXj[1] + pKR[2] * ptXj[2]) /
	//	(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

	//n[1] = (pKR[3] * ptXj[0] + pKR[4] * ptXj[1] + pKR[5] * ptXj[2]) /
	//	(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

	//if ( nP == nM )	//Main archor
	//{
	//	ptXj[0] = sin( pb[0] ) * cos( pb[1] );
	//	ptXj[1] = sin( pb[1] );
	//	ptXj[2] = cos( pb[0] ) * cos( pb[1] );

	//	pKR = KR + nP*9;
	//	n[0] = (pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
	//		(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

	//	n[1] = (pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
	//		(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);
	//}else 
	//	if ( nP == nN )	//Associate archor
	//	{
	//		posei = pa+6*nM+3;
	//		posek = pa+6*nN+3;

	//		//��ê�㵽��ê���ƽ������
	//		pti2k[0] = posek[0]-posei[0];	pti2k[1] = posek[1]-posei[1];	pti2k[2] = posek[2]-posei[2];	

	//		//��ê�㵽������ĵ�λ����
	//		ptXUnit[0] = sin( pb[0] ) * cos( pb[1] );
	//		ptXUnit[1] = sin( pb[1] );
	//		ptXUnit[2] = cos( pb[0] ) * cos( pb[1] );

	//		//compute angle w2
	//		dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1] + ptXUnit[2]*pti2k[2];
	//		dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
	//		if (dDot/dDisi2k > 1)
	//			dW2 = 0;
	//		if ( dDot/dDisi2k < -1)
	//			dW2 = PI;
	//		else
	//			dW2  = acos( dDot/dDisi2k );

	//		//compute Xk vector according sin theory
	//		ptXk[0] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[0] - sin( pb[2] ) * pti2k[0];
	//		ptXk[1] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[1] - sin( pb[2] ) * pti2k[1];
	//		ptXk[2] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[2] - sin( pb[2] ) * pti2k[2];

	//		pKR = KR + nN*9;		
	//		n[0] = ( pKR[0]*ptXk[0]	+ pKR[1]*ptXk[1] + pKR[2]*ptXk[2])/
	//			( pKR[6]*ptXk[0] + pKR[7]*ptXk[1] + pKR[8]*ptXk[2]);

	//		n[1] = ( pKR[3]*ptXk[0]	+ pKR[4]*ptXk[1] + pKR[5]*ptXk[2])/
	//			(pKR[6]*ptXk[0]	+ pKR[7]*ptXk[1] + pKR[8]*ptXk[2]);
	//	}
	//	else
	//	{
	//		posel = pa + nP*6 + 3;
	//		posek = pa + nN *6 + 3;
	//		posei = pa + nM *6 + 3;

	//		//��ê�㵽��ê���ƽ������
	//		pti2k[0] = posek[0] - posei[0];		pti2k[1] = posek[1] - posei[1];		pti2k[2] = posek[2] - posei[2];
	//		//��ê�㵽����ê���ƽ������
	//		pti2l[0] = posel[0] - posei[0];		pti2l[1] = posel[1] - posei[1];		pti2l[2] = posel[2] - posei[2];
	//		
	//		//XUnit ��ê�㵽������ĵ�λ����
	//		ptXUnit[0] = sin( pb[0] ) * cos( pb[1] );
	//		ptXUnit[1] = sin( pb[1] );
	//		ptXUnit[2] = cos( pb[0] ) * cos( pb[1] );

	//		//compute angle w2
	//		dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1]+ ptXUnit[2]*pti2k[2];
	//		dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
	//		//dW2  = acos( dDot/dDisi2k );
	//		if (dDot/dDisi2k > 1)
	//			dW2 = 0;
	//		if ( dDot/dDisi2k < -1)
	//			dW2 = PI;
	//		else
	//			dW2  = acos( dDot/dDisi2k );

	//		//compute Xl vector according sin theory
	//		ptXk[0] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[0] - sin( pb[2] ) * pti2l[0];
	//		ptXk[1] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[1] - sin( pb[2] ) * pti2l[1];
	//		ptXk[2] = dDisi2k * sin( dW2+pb[2] ) * ptXUnit[2] - sin( pb[2] ) * pti2l[2];

	//		pKR = KR + nP*9;			
	//		n[0] = (pKR[0]*ptXk[0] + pKR[1]*ptXk[1] + pKR[2]*ptXk[2])/
	//			( pKR[6]*ptXk[0] + pKR[7]*ptXk[1] + pKR[8]*ptXk[2]);

	//		n[1] = (pKR[3]*ptXk[0] + pKR[4]*ptXk[1] + pKR[5]*ptXk[2])/
	//			( pKR[6]*ptXk[0] + pKR[7]*ptXk[1] + pKR[8]*ptXk[2]);
	//	}
}
//#pragma optimize("", off) // �ر��Ż�
//compute Jacobian for each image point
//pPA pPB,uv�������������ά���һ�׵�
//pAM uv����ê�����λ�Ʋ�����һ�׵�
//pAA uv�Ը�ê�����λ�Ʋ�����һ�׵�
static void lcrp_jacobianEachPts(double *R, double* KR, double* KdA, double *KdB, double* KdG, double* pa, double* ppt, int nM, int nN, int nP,double* pAM, double* pAA, double* pPA, double* pPB)
{
	double matXj[3];
	double matxyt[3];
	double matDlcDxyt[3];
	double matDxytDRA[3], matDxytDRB[3], matDxytDRG[3], matDxytDXc[3], matDxytDYc[3], matDxytDZc[3];
	double matDlcDRA[1], matDlcDRB[1], matDlcDRG[1], matDlcDXc[1], matDlcDYc[1], matDlcDZc[1], matDlcDX[1], matDlcDY[1], matDlcDZ[1];
	//double matDXjDDH[6];
	//double matDxytDDH[6];
	//double matDuvDDH[4];
	//double *ptAngle = NULL;

	//double *posei, *posek, *posel;
	//double matPosei2k[3], matPosei2l[3];
	//double	dDisTi2k, dDot, dArcCosIn, dW;

	//double matXk[3],matXl[3];
	//double matDDotDDH[2];
	//double matDArcCosInDDH[2];
	//double dDWDArcCosIn;
	//double matDsinwWDDH[2];
	//double matDXkDpa[3];
	//double matDuvDpa[2];
	//double tmp1[6];
	//double tmp2[6];
	//
	//double matDdisDTi[3];
	//double matDDotDTi[3];
	//double matDArcCosInDTi[3];
	//double matDsinwWDTi[3];
	//double matDTi2kDTi[9], matDTi2kDTk[9];
	//double matDXkDTi[9], matDXkDTk[9], matDXlDTl[9];
	//double matDuvDTi[6], matDuvDTk[6], matDuvDTl[6];
	//double matDdisDTk[3];
	//double matDDotDTk[3], matDArcCosInDTk[3], matDsinwWDTk[3];
	double *pR, *pKR, *pKdA, *pKdB, *pKdG;

	double* ppa = pa + nP * 5 + 2;//��ǰ��ͼλ��ָ��
	//printf("%f %f %f %f %f %f\n", ppt[0], ppt[1], ppt[2], ppa[0], ppa[1], ppa[2]);
	matXj[0] = ppt[0] - ppa[0];
	matXj[1] = ppt[1] - ppa[1];
	matXj[2] = ppt[2] - ppa[2];
	pR = R + nP * 9;
	matxyt[0] = pR[0] * matXj[0] + pR[1] * matXj[1] + pR[2] * matXj[2];
	matxyt[1] = pR[3] * matXj[0] + pR[4] * matXj[1] + pR[5] * matXj[2];
	matxyt[2] = pR[6] * matXj[0] + pR[7] * matXj[1] + pR[8] * matXj[2];

	double r = sqrt(matxyt[0] * matxyt[0] + matxyt[1] * matxyt[1]);
	double T = r / matxyt[2];
	double T1 = 1 / (1 + T * T);

	matDlcDxyt[0] = T1 * matxyt[0] / (matxyt[2] * r);
	matDlcDxyt[1] = T1 * matxyt[1] / (matxyt[2] * r);
	matDlcDxyt[2] = -T1 * r / (matxyt[2] * matxyt[2]);
	//matDuvDxyt[3] = 0;
	//matDuvDxyt[4] = 1 / matxyt[2];
	//matDuvDxyt[5] = -matxyt[1] / (matxyt[2] * matxyt[2]);

	matDxytDXc[0] = -pR[0];
	matDxytDXc[1] = -pR[3];
	matDxytDXc[2] = -pR[6];

	matDxytDYc[0] = -pR[1];
	matDxytDYc[1] = -pR[4];
	matDxytDYc[2] = -pR[7];

	matDxytDZc[0] = -pR[2];
	matDxytDZc[1] = -pR[5];
	matDxytDZc[2] = -pR[8];

	matDlcDXc[0] = matDlcDxyt[0] * matDxytDXc[0] + matDlcDxyt[1] * matDxytDXc[1] + matDlcDxyt[2] * matDxytDXc[2];//lc����Xc��һ�׵�
	matDlcDYc[0] = matDlcDxyt[0] * matDxytDYc[0] + matDlcDxyt[1] * matDxytDYc[1] + matDlcDxyt[2] * matDxytDYc[2];//lc����Yc��һ�׵�
	matDlcDZc[0] = matDlcDxyt[0] * matDxytDZc[0] + matDlcDxyt[1] * matDxytDZc[1] + matDlcDxyt[2] * matDxytDZc[2];//lc����Zc��һ�׵�

	matDlcDX[0] = -matDlcDXc[0];
	matDlcDY[0] = -matDlcDYc[0];
	matDlcDZ[0] = -matDlcDZc[0];

	//pKdG:pKR�������omega��һ�׵�
	//pKdB:pKR�������phi��һ�׵�
	//camera angles
	pKdG = KdG + nP * 9;
	matDxytDRG[0] = pKdG[0] * matXj[0] + pKdG[1] * matXj[1] + pKdG[2] * matXj[2];//xyt����omega��һ�׵�
	matDxytDRG[1] = pKdG[3] * matXj[0] + pKdG[4] * matXj[1] + pKdG[5] * matXj[2];
	matDxytDRG[2] = pKdG[6] * matXj[0] + pKdG[7] * matXj[1] + pKdG[8] * matXj[2];

	pKdB = KdB + nP * 9;
	matDxytDRB[0] = pKdB[0] * matXj[0] + pKdB[1] * matXj[1] + pKdB[2] * matXj[2];//xyt����phi��һ�׵�
	matDxytDRB[1] = pKdB[3] * matXj[0] + pKdB[4] * matXj[1] + pKdB[5] * matXj[2];
	matDxytDRB[2] = pKdB[6] * matXj[0] + pKdB[7] * matXj[1] + pKdB[8] * matXj[2];

	matDlcDRB[0] = matDlcDxyt[0] * matDxytDRB[0] + matDlcDxyt[1] * matDxytDRB[1] + matDlcDxyt[2] * matDxytDRB[2];//lc����phi��һ�׵�
	matDlcDRG[0] = matDlcDxyt[0] * matDxytDRG[0] + matDlcDxyt[1] * matDxytDRG[1] + matDlcDxyt[2] * matDxytDRG[2];//lc����omega��һ�׵�

	pPA[0] = matDlcDRB[0];			pPA[1] = matDlcDRG[0];//lc����ε�һ�׵�
	pPA[2] = matDlcDXc[0];			pPA[3] = matDlcDYc[0];			pPA[4] = matDlcDZc[0];

	pPB[0] = matDlcDX[0];		pPB[1] = matDlcDY[0];		pPB[2] = matDlcDZ[0];//lc����ά���һ�׵�
	//pPB[3] = matDuvDX[1];		pPB[4] = matDuvDY[1];		pPB[5] = matDuvDZ[1];

//	if ( nP == nM )	//main archor
//	{
//		matXj[0] = sin( ppt[0] ) * cos( ppt[1] );
//		matXj[1] = sin( ppt[1] );
//		matXj[2] = cos( ppt[0] ) * cos( ppt[1] );
//
//		pKR = KR + nP*9;
//		matxyt[0] = pKR[0]*matXj[0] +pKR[1]*matXj[1] +pKR[2]*matXj[2];
//		matxyt[1] = pKR[3]*matXj[0] +pKR[4]*matXj[1] +pKR[5]*matXj[2];
//		matxyt[2] = pKR[6]*matXj[0] +pKR[7]*matXj[1] +pKR[8]*matXj[2];
//
//		matDuvDxyt[0] = 1/matxyt[2];	
//		matDuvDxyt[1] = 0;
//		matDuvDxyt[2] = -matxyt[0]/(matxyt[2]*matxyt[2]);
//		matDuvDxyt[3] = 0;	
//		matDuvDxyt[4] = 1/matxyt[2];		
//		matDuvDxyt[5] = -matxyt[1]/(matxyt[2]*matxyt[2]);
//		
//		//pKdG:pKR�������omega��һ�׵�
//		//pKdB:pKR�������phi��һ�׵�
//		//pKdA:pKR�������kappa��һ�׵�
//		//camera angles
//		pKdG = KdG + nP*9;
//		matDxytDRG[0] =  pKdG[0]*matXj[0] + pKdG[1]*matXj[1] + pKdG[2]*matXj[2];//xyt����omega��һ�׵�
//		matDxytDRG[1] =  pKdG[3]*matXj[0] + pKdG[4]*matXj[1] + pKdG[5]*matXj[2];
//		matDxytDRG[2] =  pKdG[6]*matXj[0] + pKdG[7]*matXj[1] + pKdG[8]*matXj[2];
//
//		pKdB = KdB + nP*9;
//		matDxytDRB[0] =  pKdB[0]*matXj[0] + pKdB[1]*matXj[1] + pKdB[2]*matXj[2];//xyt����phi��һ�׵�
//		matDxytDRB[1] =  pKdB[3]*matXj[0] + pKdB[4]*matXj[1] + pKdB[5]*matXj[2];
//		matDxytDRB[2] =  pKdB[6]*matXj[0] + pKdB[7]*matXj[1] + pKdB[8]*matXj[2];
//
//		pKdA = KdA + nP*9;
//		matDxytDRA[0] =  pKdA[0]*matXj[0] + pKdA[1]*matXj[1] + pKdA[2]*matXj[2];//xyt����kappa��һ�׵�
//		matDxytDRA[1] =  pKdA[3]*matXj[0] + pKdA[4]*matXj[1] + pKdA[5]*matXj[2];
//		matDxytDRA[2] =  pKdA[6]*matXj[0] + pKdA[7]*matXj[1] + pKdA[8]*matXj[2];
//
//		matDuvDRA[0] = matDuvDxyt[0]*matDxytDRA[0] + matDuvDxyt[1]*matDxytDRA[1] + matDuvDxyt[2]*matDxytDRA[2];//uv����kappa��һ�׵�
//		matDuvDRA[1] = matDuvDxyt[3]*matDxytDRA[0] + matDuvDxyt[4]*matDxytDRA[1] + matDuvDxyt[5]*matDxytDRA[2];
//
//		matDuvDRB[0] = matDuvDxyt[0]*matDxytDRB[0] + matDuvDxyt[1]*matDxytDRB[1] + matDuvDxyt[2]*matDxytDRB[2];//uv����phi��һ�׵�
//		matDuvDRB[1] = matDuvDxyt[3]*matDxytDRB[0] + matDuvDxyt[4]*matDxytDRB[1] + matDuvDxyt[5]*matDxytDRB[2];
//
//		matDuvDRG[0] = matDuvDxyt[0]*matDxytDRG[0] + matDuvDxyt[1]*matDxytDRG[1] + matDuvDxyt[2]*matDxytDRG[2];//uv����omega��һ�׵�
//		matDuvDRG[1] = matDuvDxyt[3]*matDxytDRG[0] + matDuvDxyt[4]*matDxytDRG[1] + matDuvDxyt[5]*matDxytDRG[2];
//
//		pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv����ê����ת������һ�׵�
//		pPA[3] = 0;						pPA[4] = 0;						pPA[5] = 0;
//		pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
//		pPA[9] = 0;						pPA[10] = 0;					pPA[11] = 0;
//		
//		//matXj[0] = sin(ppt[0]) * cos(ppt[1]);
//		//matXj[1] = sin(ppt[1]);
//		//matXj[2] = cos(ppt[0]) * cos(ppt[1]);
//		//matxyt[0] = pKR[0] * matXj[0] + pKR[1] * matXj[1] + pKR[2] * matXj[2];
//		//matxyt[1] = pKR[3] * matXj[0] + pKR[4] * matXj[1] + pKR[5] * matXj[2];
//		//matxyt[2] = pKR[6] * matXj[0] + pKR[7] * matXj[1] + pKR[8] * matXj[2];
//		//azimuth and elevation angle
//		matDXjDDH[0] = cos(ppt[0])*cos(ppt[1]);			matDXjDDH[2] = 0;
//		matDXjDDH[4] = -sin(ppt[0])*cos(ppt[1]);		matDXjDDH[1] = -sin(ppt[0])*sin(ppt[1]);
//		matDXjDDH[3] = cos(ppt[1]);						matDXjDDH[5] = -cos(ppt[0])*sin(ppt[1]);
//
//		matDxytDDH[0] = pKR[0]*matDXjDDH[0]	+ pKR[1]*matDXjDDH[2] + pKR[2]*matDXjDDH[4];//x����azimuth�ĵ���
//		matDxytDDH[1] = pKR[0]*matDXjDDH[1] + pKR[1]*matDXjDDH[3] + pKR[2]*matDXjDDH[5];//x����elevation�ĵ���
//
//		matDxytDDH[2] = pKR[3]*matDXjDDH[0]	+ pKR[4]*matDXjDDH[2] + pKR[5]*matDXjDDH[4];//y����azimuth�ĵ���
//		matDxytDDH[3] = pKR[3]*matDXjDDH[1] + pKR[4]*matDXjDDH[3] + pKR[5]*matDXjDDH[5];//y����elevation�ĵ���
//
//		matDxytDDH[4] = pKR[6]*matDXjDDH[0]	+ pKR[7]*matDXjDDH[2] + pKR[8]*matDXjDDH[4];//t����azimuth�ĵ���
//		matDxytDDH[5] = pKR[6]*matDXjDDH[1] + pKR[7]*matDXjDDH[3] + pKR[8]*matDXjDDH[5];//t����elevation�ĵ���
//
//		matDuvDDH[0] = matDuvDxyt[0]*matDxytDDH[0] + matDuvDxyt[1]*matDxytDDH[2] + matDuvDxyt[2]*matDxytDDH[4];//uv����azimuth�ĵ���
//		matDuvDDH[1] = matDuvDxyt[0]*matDxytDDH[1] + matDuvDxyt[1]*matDxytDDH[3] + matDuvDxyt[2]*matDxytDDH[5];
//		matDuvDDH[2] = matDuvDxyt[3]*matDxytDDH[0] + matDuvDxyt[4]*matDxytDDH[2] + matDuvDxyt[5]*matDxytDDH[4];//uv����elevation�ĵ���
//		matDuvDDH[3] = matDuvDxyt[3]*matDxytDDH[1] + matDuvDxyt[4]*matDxytDDH[3] + matDuvDxyt[5]*matDxytDDH[5];
//
//		pPB[0] = matDuvDDH[0];		pPB[1] = matDuvDDH[1];		pPB[2] = 0;//uv��azimuth��elevation��һ�׵�
//		pPB[3] = matDuvDDH[2];		pPB[4] = matDuvDDH[3];		pPB[5] = 0;
//		
//	}else 
//		if ( nP == nN )	//associate archor
//		{			
//			ptAngle = pa + nN*6;
//			matXj[0] = sin( ppt[0] ) * cos( ppt[1] );
//			matXj[1] = sin( ppt[1] );
//			matXj[2] = cos( ppt[0] ) * cos( ppt[1] );
//			
//			pKR = KR + nN*9;
//			pKdA= KdA+ nN*9;
//			pKdB= KdB+ nN*9;
//			pKdG= KdG+ nN*9;
//			
//			posei = pa + nM*6 + 3;//��ê��Xc
//			posek = pa + nN*6 + 3;//��ê��Xc
//			matPosei2k[0] = posek[0]-posei[0];		matPosei2k[1] = posek[1]-posei[1];		matPosei2k[2] = posek[2]-posei[2];//��ê�㵽��ê�������
//
//			dDisTi2k = sqrt( matPosei2k[0]*matPosei2k[0] + matPosei2k[1]*matPosei2k[1] + matPosei2k[2]*matPosei2k[2] );//��ê�㵽��ê���������ģ��
//			dDot     = matXj[0]*matPosei2k[0] + matXj[1]*matPosei2k[1] + matXj[2]*matPosei2k[2];
//			dArcCosIn= dDot / dDisTi2k;
//			dW       = acos( dArcCosIn );
//
//			matXk[0] = dDisTi2k * sin( dW + ppt[2] ) * matXj[0] - sin(ppt[2])*matPosei2k[0];
//			matXk[1] = dDisTi2k * sin( dW + ppt[2] ) * matXj[1] - sin(ppt[2])*matPosei2k[1];
//			matXk[2] = dDisTi2k * sin( dW + ppt[2] ) * matXj[2] - sin(ppt[2])*matPosei2k[2];
//
//			matxyt[0] = pKR[0]*matXk[0] + pKR[1]*matXk[1] + pKR[2]*matXk[2];
//			matxyt[1] = pKR[3]*matXk[0] + pKR[4]*matXk[1] + pKR[5]*matXk[2];
//			matxyt[2] = pKR[6]*matXk[0] + pKR[7]*matXk[1] + pKR[8]*matXk[2];
//
//			matDuvDxyt[0] = 1/matxyt[2];	
//			matDuvDxyt[1] = 0;
//			matDuvDxyt[2] = -matxyt[0]/(matxyt[2]*matxyt[2]);
//			matDuvDxyt[3] = 0;	
//			matDuvDxyt[4] = 1/matxyt[2];		
//			matDuvDxyt[5] = -matxyt[1]/(matxyt[2]*matxyt[2]);
//			
////pKdG:pKR�������omega��һ�׵�
////pKdB:pKR�������phi��һ�׵�
////pKdA:pKR�������kappa��һ�׵�
//			//camera angles
//			matDxytDRG[0] =  pKdG[0]*matXk[0] + pKdG[1]*matXk[1] + pKdG[2]*matXk[2];
//			matDxytDRG[1] =  pKdG[3]*matXk[0] + pKdG[4]*matXk[1] + pKdG[5]*matXk[2];
//			matDxytDRG[2] =  pKdG[6]*matXk[0] + pKdG[7]*matXk[1] + pKdG[8]*matXk[2];
//
//			
//			matDxytDRB[0] =  pKdB[0]*matXk[0] + pKdB[1]*matXk[1] + pKdB[2]*matXk[2];
//			matDxytDRB[1] =  pKdB[3]*matXk[0] + pKdB[4]*matXk[1] + pKdB[5]*matXk[2];
//			matDxytDRB[2] =  pKdB[6]*matXk[0] + pKdB[7]*matXk[1] + pKdB[8]*matXk[2];
//
//			matDxytDRA[0] =  pKdA[0]*matXk[0] + pKdA[1]*matXk[1] + pKdA[2]*matXk[2];
//			matDxytDRA[1] =  pKdA[3]*matXk[0] + pKdA[4]*matXk[1] + pKdA[5]*matXk[2];
//			matDxytDRA[2] =  pKdA[6]*matXk[0] + pKdA[7]*matXk[1] + pKdA[8]*matXk[2];
//
//			matDuvDRA[0] = matDuvDxyt[0]*matDxytDRA[0] + matDuvDxyt[1]*matDxytDRA[1] + matDuvDxyt[2]*matDxytDRA[2];//uv����kappa��һ�׵�
//			matDuvDRA[1] = matDuvDxyt[3]*matDxytDRA[0] + matDuvDxyt[4]*matDxytDRA[1] + matDuvDxyt[5]*matDxytDRA[2];
//
//			matDuvDRB[0] = matDuvDxyt[0]*matDxytDRB[0] + matDuvDxyt[1]*matDxytDRB[1] + matDuvDxyt[2]*matDxytDRB[2];//uv����phi��һ�׵�
//			matDuvDRB[1] = matDuvDxyt[3]*matDxytDRB[0] + matDuvDxyt[4]*matDxytDRB[1] + matDuvDxyt[5]*matDxytDRB[2];
//
//			matDuvDRG[0] = matDuvDxyt[0]*matDxytDRG[0] + matDuvDxyt[1]*matDxytDRG[1] + matDuvDxyt[2]*matDxytDRG[2];//uv����omega��һ�׵�
//			matDuvDRG[1] = matDuvDxyt[3]*matDxytDRG[0] + matDuvDxyt[4]*matDxytDRG[1] + matDuvDxyt[5]*matDxytDRG[2];
//			
//			//Xj����azimuth��elevation��һ�׵�
//			//azimuth and elevation angle
//			matDXjDDH[0] = cos(ppt[0])*cos(ppt[1]);			matDXjDDH[2] = 0;
//			matDXjDDH[4] = -sin(ppt[0])*cos(ppt[1]);		matDXjDDH[1] = -sin(ppt[0])*sin(ppt[1]);
//			matDXjDDH[3] = cos(ppt[1]);						matDXjDDH[5] = -cos(ppt[0])*sin(ppt[1]);
//
////matPosei2k[0] = posek[0] - posei[0];		matPosei2k[1] = posek[1] - posei[1];		matPosei2k[2] = posek[2] - posei[2];//��ê�㵽��ê�������
//
////dDisTi2k = sqrt(matPosei2k[0] * matPosei2k[0] + matPosei2k[1] * matPosei2k[1] + matPosei2k[2] * matPosei2k[2]);//��ê�㵽��ê���������ģ��
////dDot = matXj[0] * matPosei2k[0] + matXj[1] * matPosei2k[1] + matXj[2] * matPosei2k[2];
////dArcCosIn = dDot / dDisTi2k;
////dW = acos(dArcCosIn);
//
//			matDDotDDH[0] = cos(ppt[0])*cos(ppt[1])*matPosei2k[0] - sin(ppt[0])*cos(ppt[1])*matPosei2k[2];//dDot����azimuth��һ�׵�
//			matDDotDDH[1] = -sin(ppt[0])*sin(ppt[1])*matPosei2k[0] + cos(ppt[1])*matPosei2k[1]-cos(ppt[0])*sin(ppt[1])*matPosei2k[2];//dDot����elevation��һ�׵�
//
//			matDArcCosInDDH[0] = matDDotDDH[0]/dDisTi2k;//ArcCosIn����azimuth��һ�׵�
//			matDArcCosInDDH[1] = matDDotDDH[1]/dDisTi2k;//ArcCosIn����elevation��һ�׵�
//
//			dDWDArcCosIn = -1/sqrt(1-dArcCosIn*dArcCosIn);//dW����ArccosIn��һ�׵�
//
////matXk[0] = dDisTi2k * sin(dW + ppt[2]) * matXj[0] - sin(ppt[2]) * matPosei2k[0];
////matXk[1] = dDisTi2k * sin(dW + ppt[2]) * matXj[1] - sin(ppt[2]) * matPosei2k[1];
////matXk[2] = dDisTi2k * sin(dW + ppt[2]) * matXj[2] - sin(ppt[2]) * matPosei2k[2];
//
//			matDsinwWDDH[0] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[0];//sin(dW + ppt[2])����azimuth��һ�׵�
//			matDsinwWDDH[1] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[1];//sin(dW + ppt[2])����elevation��һ�׵�
//
//			tmp2[0] =  dDisTi2k*matXj[0]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[0];//matXk[0]����azimuth��һ�׵�
//			tmp2[1] =  dDisTi2k*matXj[0]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[1];//matXk[0]����elevation��һ�׵�
//			tmp2[2] =  dDisTi2k*matXj[1]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[2];//matXk[1]����azimuth��һ�׵�
//			tmp2[3] =  dDisTi2k*matXj[1]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[3];//matXk[1]����elevation��һ�׵�
//			tmp2[4] =  dDisTi2k*matXj[2]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[4];//matXk[2]����azimuth��һ�׵�
//			tmp2[5] =  dDisTi2k*matXj[2]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[5];//matXk[2]����elevation��һ�׵�
//
//			matDxytDDH[0] = pKR[0]*tmp2[0] + pKR[1]*tmp2[2] + pKR[2]*tmp2[4];//x����azimuth��һ�׵�
//			matDxytDDH[1] = pKR[0]*tmp2[1] + pKR[1]*tmp2[3] + pKR[2]*tmp2[5];//x����elevation��һ�׵�
//			matDxytDDH[2] = pKR[3]*tmp2[0] + pKR[4]*tmp2[2] + pKR[5]*tmp2[4];//y����azimuth��һ�׵�
//			matDxytDDH[3] = pKR[3]*tmp2[1] + pKR[4]*tmp2[3] + pKR[5]*tmp2[5];//y����elevation��һ�׵�
//			matDxytDDH[4] = pKR[6]*tmp2[0] + pKR[7]*tmp2[2] + pKR[8]*tmp2[4];//t����azimuth��һ�׵�
//			matDxytDDH[5] = pKR[6]*tmp2[1] + pKR[7]*tmp2[3] + pKR[8]*tmp2[5];//t����elevation��һ�׵�
//
//			matDuvDDH[0] = matDuvDxyt[0]*matDxytDDH[0] + matDuvDxyt[1]*matDxytDDH[2] + matDuvDxyt[2]*matDxytDDH[4];//u����azimuth��һ�׵�
//			matDuvDDH[1] = matDuvDxyt[0]*matDxytDDH[1] + matDuvDxyt[1]*matDxytDDH[3] + matDuvDxyt[2]*matDxytDDH[5];//u����elevation��һ�׵�
//			matDuvDDH[2] = matDuvDxyt[3]*matDxytDDH[0] + matDuvDxyt[4]*matDxytDDH[2] + matDuvDxyt[5]*matDxytDDH[4];//v����azimuth��һ�׵�
//			matDuvDDH[3] = matDuvDxyt[3]*matDxytDDH[1] + matDuvDxyt[4]*matDxytDDH[3] + matDuvDxyt[5]*matDxytDDH[5];//v����elevation��һ�׵�
//			
//			//parallax angle
//			matDXkDpa[0] = dDisTi2k*cos(dW+ppt[2])*matXj[0] - cos(ppt[2])*matPosei2k[0];//matXk[0]����parallax��һ�׵�
//			matDXkDpa[1] = dDisTi2k*cos(dW+ppt[2])*matXj[1] - cos(ppt[2])*matPosei2k[1];//matXk[1]����parallax��һ�׵�
//			matDXkDpa[2] = dDisTi2k*cos(dW+ppt[2])*matXj[2] - cos(ppt[2])*matPosei2k[2];//matXk[2]����parallax��һ�׵�
//
//			tmp1[0] = pKR[0]*matDuvDxyt[0] + pKR[3]*matDuvDxyt[1] + pKR[6]*matDuvDxyt[2];
//			tmp1[1] = pKR[1]*matDuvDxyt[0] + pKR[4]*matDuvDxyt[1] + pKR[7]*matDuvDxyt[2];
//			tmp1[2] = pKR[2]*matDuvDxyt[0] + pKR[5]*matDuvDxyt[1] + pKR[8]*matDuvDxyt[2];
//			tmp1[3] = pKR[0]*matDuvDxyt[3] + pKR[3]*matDuvDxyt[4] + pKR[6]*matDuvDxyt[5];
//			tmp1[4] = pKR[1]*matDuvDxyt[3] + pKR[4]*matDuvDxyt[4] + pKR[7]*matDuvDxyt[5];
//			tmp1[5] = pKR[2]*matDuvDxyt[3] + pKR[5]*matDuvDxyt[4] + pKR[8]*matDuvDxyt[5];
//
//			matDuvDpa[0] = tmp1[0]*matDXkDpa[0] + tmp1[1]*matDXkDpa[1] + tmp1[2]*matDXkDpa[2];//u����parallax��һ�׵�
//			matDuvDpa[1] = tmp1[3]*matDXkDpa[0] + tmp1[4]*matDXkDpa[1] + tmp1[5]*matDXkDpa[2];//v����parallax��һ�׵�
//
//			//uv����azimuth��elevation��parallax��һ�׵�
//			pPB[0] = matDuvDDH[0];		pPB[1] = matDuvDDH[1];		pPB[2] = matDuvDpa[0];
//			pPB[3] = matDuvDDH[2];		pPB[4] = matDuvDDH[3];		pPB[5] = matDuvDpa[1];
//
////posei = pa + nM * 6 + 3;//��ê��Xc
////posek = pa + nN * 6 + 3;//��ê��Xc
////matPosei2k[0] = posek[0] - posei[0];		matPosei2k[1] = posek[1] - posei[1];		matPosei2k[2] = posek[2] - posei[2];//��ê�㵽��ê�������
////dDisTi2k = sqrt(matPosei2k[0] * matPosei2k[0] + matPosei2k[1] * matPosei2k[1] + matPosei2k[2] * matPosei2k[2]);//��ê�㵽��ê���������ģ��
////dDot = matXj[0] * matPosei2k[0] + matXj[1] * matPosei2k[1] + matXj[2] * matPosei2k[2];
////dArcCosIn = dDot / dDisTi2k;
////dW = acos(dArcCosIn);
////matXk[0] = dDisTi2k * sin(dW + ppt[2]) * matXj[0] - sin(ppt[2]) * matPosei2k[0];
////matXk[1] = dDisTi2k * sin(dW + ppt[2]) * matXj[1] - sin(ppt[2]) * matPosei2k[1];
////matXk[2] = dDisTi2k * sin(dW + ppt[2]) * matXj[2] - sin(ppt[2]) * matPosei2k[2];
//// 
//			//Ti
//			matDdisDTi[0] = -matPosei2k[0]/dDisTi2k;//dDisTi2k����ê������Ti��һ�׵�
//			matDdisDTi[1] = -matPosei2k[1]/dDisTi2k;
//			matDdisDTi[2] = -matPosei2k[2]/dDisTi2k;
//
//
//			matDDotDTi[0] = -sin(ppt[0])*cos(ppt[1]);//dDot����ê������Ti��һ�׵�
//			matDDotDTi[1] = -sin(ppt[1]);
//			matDDotDTi[2] = -cos(ppt[0])*cos(ppt[1]);
//
//			
//			matDArcCosInDTi[0] = (dDisTi2k*matDDotDTi[0]-dDot*matDdisDTi[0])/(dDisTi2k*dDisTi2k);//dArcCosIn����ê������Ti��һ�׵�
//			matDArcCosInDTi[1] = (dDisTi2k*matDDotDTi[1]-dDot*matDdisDTi[1])/(dDisTi2k*dDisTi2k);
//			matDArcCosInDTi[2] = (dDisTi2k*matDDotDTi[2]-dDot*matDdisDTi[2])/(dDisTi2k*dDisTi2k);
//
//			matDsinwWDTi[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[0];//sin(dW + ppt[2])����ê������Ti��һ�׵�
//			matDsinwWDTi[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[1];
//			matDsinwWDTi[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[2];
//
//			matDTi2kDTi[0] = matDTi2kDTi[4] = matDTi2kDTi[8] = -1;//matPosei2k����ê������Ti��һ�׵�
//			matDTi2kDTi[1] = matDTi2kDTi[2] = matDTi2kDTi[3] = 0;
//			matDTi2kDTi[5] = matDTi2kDTi[6] = matDTi2kDTi[7] = 0;
//
//			matDXkDTi[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[0] + dDisTi2k*matXj[0]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[0];//matXk����ê������Ti��һ�׵�
//			matDXkDTi[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[1] + dDisTi2k*matXj[0]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[1];
//			matDXkDTi[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[2] + dDisTi2k*matXj[0]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[2];
//			matDXkDTi[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[0] + dDisTi2k*matXj[1]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[3];
//			matDXkDTi[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[1] + dDisTi2k*matXj[1]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[4];
//			matDXkDTi[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[2] + dDisTi2k*matXj[1]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[5];
//			matDXkDTi[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[0] + dDisTi2k*matXj[2]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[6];
//			matDXkDTi[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[1] + dDisTi2k*matXj[2]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[7];
//			matDXkDTi[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[2] + dDisTi2k*matXj[2]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[8];
//
//			matDuvDTi[0] = tmp1[0]*matDXkDTi[0] + tmp1[1]*matDXkDTi[3] + tmp1[2]*matDXkDTi[6];//uv����ê������Ti��һ�׵�
//			matDuvDTi[1] = tmp1[0]*matDXkDTi[1] + tmp1[1]*matDXkDTi[4] + tmp1[2]*matDXkDTi[7];
//			matDuvDTi[2] = tmp1[0]*matDXkDTi[2] + tmp1[1]*matDXkDTi[5] + tmp1[2]*matDXkDTi[8];
//			matDuvDTi[3] = tmp1[3]*matDXkDTi[0] + tmp1[4]*matDXkDTi[3] + tmp1[5]*matDXkDTi[6];
//			matDuvDTi[4] = tmp1[3]*matDXkDTi[1] + tmp1[4]*matDXkDTi[4] + tmp1[5]*matDXkDTi[7];
//			matDuvDTi[5] = tmp1[3]*matDXkDTi[2] + tmp1[4]*matDXkDTi[5] + tmp1[5]*matDXkDTi[8];
//
//			pAM[0] = matDuvDTi[0];		pAM[1] = matDuvDTi[1];		pAM[2] = matDuvDTi[2];
//			pAM[3] = matDuvDTi[3];		pAM[4] = matDuvDTi[4];		pAM[5] = matDuvDTi[5];				 
//			
//			//Tk
//			matDdisDTk[0] = matPosei2k[0]/dDisTi2k;
//			matDdisDTk[1] = matPosei2k[1]/dDisTi2k;
//			matDdisDTk[2] = matPosei2k[2]/dDisTi2k;
//
//			matDDotDTk[0] = sin(ppt[0])*cos(ppt[1]);
//			matDDotDTk[1] = sin(ppt[1]);
//			matDDotDTk[2] = cos(ppt[0])*cos(ppt[1]);
//
//			matDArcCosInDTk[0] = (dDisTi2k*matDDotDTk[0] - dDot*matDdisDTk[0])/(dDisTi2k*dDisTi2k);
//			matDArcCosInDTk[1] = (dDisTi2k*matDDotDTk[1] - dDot*matDdisDTk[1])/(dDisTi2k*dDisTi2k);
//			matDArcCosInDTk[2] = (dDisTi2k*matDDotDTk[2] - dDot*matDdisDTk[2])/(dDisTi2k*dDisTi2k);
//
//			matDsinwWDTk[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[0];
//			matDsinwWDTk[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[1];
//			matDsinwWDTk[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[2];
//
//			matDTi2kDTk[0] = matDTi2kDTk[4] = matDTi2kDTk[8] = 1;
//			matDTi2kDTk[1] = matDTi2kDTk[2] = matDTi2kDTk[3] = 0;
//			matDTi2kDTk[5] = matDTi2kDTk[6] = matDTi2kDTk[7] = 0;
//
//
//			matDXkDTk[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[0] + dDisTi2k*matXj[0]*matDsinwWDTk[0] - sin(ppt[2])*matDTi2kDTk[0];
//			matDXkDTk[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[1] + dDisTi2k*matXj[0]*matDsinwWDTk[1] - sin(ppt[2])*matDTi2kDTk[1];
//			matDXkDTk[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[2] + dDisTi2k*matXj[0]*matDsinwWDTk[2] - sin(ppt[2])*matDTi2kDTk[2];
//			matDXkDTk[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[0] + dDisTi2k*matXj[1]*matDsinwWDTk[0] - sin(ppt[2])*matDTi2kDTk[3];
//			matDXkDTk[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[1] + dDisTi2k*matXj[1]*matDsinwWDTk[1] - sin(ppt[2])*matDTi2kDTk[4];
//			matDXkDTk[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[2] + dDisTi2k*matXj[1]*matDsinwWDTk[2] - sin(ppt[2])*matDTi2kDTk[5];
//			matDXkDTk[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[0] + dDisTi2k*matXj[2]*matDsinwWDTk[0] - sin(ppt[2])*matDTi2kDTk[6];
//			matDXkDTk[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[1] + dDisTi2k*matXj[2]*matDsinwWDTk[1] - sin(ppt[2])*matDTi2kDTk[7];
//			matDXkDTk[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[2] + dDisTi2k*matXj[2]*matDsinwWDTk[2] - sin(ppt[2])*matDTi2kDTk[8];
//
//			matDuvDTk[0] = tmp1[0]*matDXkDTk[0] + tmp1[1]*matDXkDTk[3] + tmp1[2]*matDXkDTk[6];//uv�Ը�ê������Tk��һ�׵�
//			matDuvDTk[1] = tmp1[0]*matDXkDTk[1] + tmp1[1]*matDXkDTk[4] + tmp1[2]*matDXkDTk[7];
//			matDuvDTk[2] = tmp1[0]*matDXkDTk[2] + tmp1[1]*matDXkDTk[5] + tmp1[2]*matDXkDTk[8];
//			matDuvDTk[3] = tmp1[3]*matDXkDTk[0] + tmp1[4]*matDXkDTk[3] + tmp1[5]*matDXkDTk[6];
//			matDuvDTk[4] = tmp1[3]*matDXkDTk[1] + tmp1[4]*matDXkDTk[4] + tmp1[5]*matDXkDTk[7];
//			matDuvDTk[5] = tmp1[3]*matDXkDTk[2] + tmp1[4]*matDXkDTk[5] + tmp1[5]*matDXkDTk[8];
//			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//			pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv�Ը�ê�������һ�׵�
//			pPA[3] = matDuvDTk[0];			pPA[4] = matDuvDTk[1];			pPA[5] = matDuvDTk[2];
//			pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
//			pPA[9] = matDuvDTk[3];			pPA[10] = matDuvDTk[4];			pPA[11] = matDuvDTk[5];
//		}
//		else
//		{
//			ptAngle = pa + nP*6;
//			matXj[0] = sin( ppt[0] ) * cos( ppt[1] );
//			matXj[1] = sin( ppt[1] );
//			matXj[2] = cos( ppt[0] ) * cos( ppt[1] );
//
//			pKR = KR + nP*9;
//			pKdA= KdA+ nP*9;
//			pKdB= KdB+ nP*9;
//			pKdG= KdG+ nP*9;
//			
//			posei = pa + nM*6 + 3;//��ê��
//			posek = pa + nN*6 + 3;//��ê��
//			posel = pa + nP*6  + 3;//����ê��
//
//			matPosei2k[0] = posek[0]-posei[0];		matPosei2k[1] = posek[1]-posei[1];		matPosei2k[2] = posek[2]-posei[2];//��ê�㵽��ê��
//			matPosei2l[0] = posel[0]-posei[0];		matPosei2l[1] = posel[1]-posei[1];		matPosei2l[2] = posel[2]-posei[2];//��ê�㵽����ê��
//
//			dDisTi2k = sqrt( matPosei2k[0]*matPosei2k[0] + matPosei2k[1]*matPosei2k[1] + matPosei2k[2]*matPosei2k[2] );
//			dDot     = matXj[0]*matPosei2k[0] + matXj[1]*matPosei2k[1] + matXj[2]*matPosei2k[2];
//			dArcCosIn= dDot / dDisTi2k;
//			dW       = acos( dArcCosIn );
//
//			matXl[0] = dDisTi2k * sin( dW + ppt[2] ) * matXj[0] - sin(ppt[2])*matPosei2l[0];
//			matXl[1] = dDisTi2k * sin( dW + ppt[2] ) * matXj[1] - sin(ppt[2])*matPosei2l[1];
//			matXl[2] = dDisTi2k * sin( dW + ppt[2] ) * matXj[2] - sin(ppt[2])*matPosei2l[2];
//
//			matxyt[0] = pKR[0]*matXl[0] + pKR[1]*matXl[1] + pKR[2]*matXl[2];
//			matxyt[1] = pKR[3]*matXl[0] + pKR[4]*matXl[1] + pKR[5]*matXl[2];
//			matxyt[2] = pKR[6]*matXl[0] + pKR[7]*matXl[1] + pKR[8]*matXl[2];
//
//			matDuvDxyt[0] = 1/matxyt[2];	
//			matDuvDxyt[1] = 0;
//			matDuvDxyt[2] = -matxyt[0]/(matxyt[2]*matxyt[2]);
//			matDuvDxyt[3] = 0;	
//			matDuvDxyt[4] = 1/matxyt[2];		
//			matDuvDxyt[5] = -matxyt[1]/(matxyt[2]*matxyt[2]);	
//		
//			//camera angle
//			matDxytDRG[0] =  pKdG[0]*matXl[0] + pKdG[1]*matXl[1] + pKdG[2]*matXl[2];
//			matDxytDRG[1] =  pKdG[3]*matXl[0] + pKdG[4]*matXl[1] + pKdG[5]*matXl[2];
//			matDxytDRG[2] =  pKdG[6]*matXl[0] + pKdG[7]*matXl[1] + pKdG[8]*matXl[2];
//
//			matDxytDRB[0] =  pKdB[0]*matXl[0] + pKdB[1]*matXl[1] + pKdB[2]*matXl[2];
//			matDxytDRB[1] =  pKdB[3]*matXl[0] + pKdB[4]*matXl[1] + pKdB[5]*matXl[2];
//			matDxytDRB[2] =  pKdB[6]*matXl[0] + pKdB[7]*matXl[1] + pKdB[8]*matXl[2];
//
//			matDxytDRA[0] =  pKdA[0]*matXl[0] + pKdA[1]*matXl[1] + pKdA[2]*matXl[2];
//			matDxytDRA[1] =  pKdA[3]*matXl[0] + pKdA[4]*matXl[1] + pKdA[5]*matXl[2];
//			matDxytDRA[2] =  pKdA[6]*matXl[0] + pKdA[7]*matXl[1] + pKdA[8]*matXl[2];			
//
//			matDuvDRA[0] = matDuvDxyt[0]*matDxytDRA[0] + matDuvDxyt[1]*matDxytDRA[1] + matDuvDxyt[2]*matDxytDRA[2];//uv������ê��kappa������һ�׵�
//			matDuvDRA[1] = matDuvDxyt[3]*matDxytDRA[0] + matDuvDxyt[4]*matDxytDRA[1] + matDuvDxyt[5]*matDxytDRA[2];
//
//			matDuvDRB[0] = matDuvDxyt[0]*matDxytDRB[0] + matDuvDxyt[1]*matDxytDRB[1] + matDuvDxyt[2]*matDxytDRB[2];//uv������ê��phi������һ�׵�
//			matDuvDRB[1] = matDuvDxyt[3]*matDxytDRB[0] + matDuvDxyt[4]*matDxytDRB[1] + matDuvDxyt[5]*matDxytDRB[2];
//
//			matDuvDRG[0] = matDuvDxyt[0]*matDxytDRG[0] + matDuvDxyt[1]*matDxytDRG[1] + matDuvDxyt[2]*matDxytDRG[2];//uv������ê��omega������һ�׵�
//			matDuvDRG[1] = matDuvDxyt[3]*matDxytDRG[0] + matDuvDxyt[4]*matDxytDRG[1] + matDuvDxyt[5]*matDxytDRG[2];
//		
//			//azimuth and elevation angle
//			matDXjDDH[0] = cos(ppt[0])*cos(ppt[1]);			matDXjDDH[2] = 0;
//			matDXjDDH[4] = -sin(ppt[0])*cos(ppt[1]);		matDXjDDH[1] = -sin(ppt[0])*sin(ppt[1]);
//			matDXjDDH[3] = cos(ppt[1]);						matDXjDDH[5] = -cos(ppt[0])*sin(ppt[1]);
//
//			matDDotDDH[0] = cos(ppt[0])*cos(ppt[1])*matPosei2k[0] - sin(ppt[0])*cos(ppt[1])*matPosei2k[2];
//			matDDotDDH[1] = -sin(ppt[0])*sin(ppt[1])*matPosei2k[0] + cos(ppt[1])*matPosei2k[1]-cos(ppt[0])*sin(ppt[1])*matPosei2k[2];
//
//			matDArcCosInDDH[0] = matDDotDDH[0]/dDisTi2k;
//			matDArcCosInDDH[1] = matDDotDDH[1]/dDisTi2k;
//
//			dDWDArcCosIn = -1/sqrt(1-dArcCosIn*dArcCosIn);
//
//			matDsinwWDDH[0] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[0];
//			matDsinwWDDH[1] = cos(ppt[2]+dW)*dDWDArcCosIn*matDArcCosInDDH[1];
//
//			tmp2[0] =  dDisTi2k*matXj[0]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[0];
//			tmp2[1] =  dDisTi2k*matXj[0]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[1];
//			tmp2[2] =  dDisTi2k*matXj[1]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[2];
//			tmp2[3] =  dDisTi2k*matXj[1]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[3];
//			tmp2[4] =  dDisTi2k*matXj[2]*matDsinwWDDH[0] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[4];
//			tmp2[5] =  dDisTi2k*matXj[2]*matDsinwWDDH[1] + dDisTi2k*sin(dW+ppt[2])*matDXjDDH[5];
//
//			matDxytDDH[0] = pKR[0]*tmp2[0] + pKR[1]*tmp2[2] + pKR[2]*tmp2[4];
//			matDxytDDH[1] = pKR[0]*tmp2[1] + pKR[1]*tmp2[3] + pKR[2]*tmp2[5];
//			matDxytDDH[2] = pKR[3]*tmp2[0] + pKR[4]*tmp2[2] + pKR[5]*tmp2[4];
//			matDxytDDH[3] = pKR[3]*tmp2[1] + pKR[4]*tmp2[3] + pKR[5]*tmp2[5];
//			matDxytDDH[4] = pKR[6]*tmp2[0] + pKR[7]*tmp2[2] + pKR[8]*tmp2[4];
//			matDxytDDH[5] = pKR[6]*tmp2[1] + pKR[7]*tmp2[3] + pKR[8]*tmp2[5];
//
//			matDuvDDH[0] = matDuvDxyt[0]*matDxytDDH[0] + matDuvDxyt[1]*matDxytDDH[2] + matDuvDxyt[2]*matDxytDDH[4];//u��azimuth��һ�׵�
//			matDuvDDH[1] = matDuvDxyt[0]*matDxytDDH[1] + matDuvDxyt[1]*matDxytDDH[3] + matDuvDxyt[2]*matDxytDDH[5];//u��elevation��һ�׵�
//			matDuvDDH[2] = matDuvDxyt[3]*matDxytDDH[0] + matDuvDxyt[4]*matDxytDDH[2] + matDuvDxyt[5]*matDxytDDH[4];//v��azimuth��һ�׵�
//			matDuvDDH[3] = matDuvDxyt[3]*matDxytDDH[1] + matDuvDxyt[4]*matDxytDDH[3] + matDuvDxyt[5]*matDxytDDH[5];//v��elevation��һ�׵�
//		
//			//parallax angle
//			matDXkDpa[0] = dDisTi2k*cos(dW+ppt[2])*matXj[0] - cos(ppt[2])*matPosei2l[0];
//			matDXkDpa[1] = dDisTi2k*cos(dW+ppt[2])*matXj[1] - cos(ppt[2])*matPosei2l[1];
//			matDXkDpa[2] = dDisTi2k*cos(dW+ppt[2])*matXj[2] - cos(ppt[2])*matPosei2l[2];	
//
//			tmp1[0] = pKR[0]*matDuvDxyt[0] + pKR[3]*matDuvDxyt[1] + pKR[6]*matDuvDxyt[2];
//			tmp1[1] = pKR[1]*matDuvDxyt[0] + pKR[4]*matDuvDxyt[1] + pKR[7]*matDuvDxyt[2];
//			tmp1[2] = pKR[2]*matDuvDxyt[0] + pKR[5]*matDuvDxyt[1] + pKR[8]*matDuvDxyt[2];
//			tmp1[3] = pKR[0]*matDuvDxyt[3] + pKR[3]*matDuvDxyt[4] + pKR[6]*matDuvDxyt[5];
//			tmp1[4] = pKR[1]*matDuvDxyt[3] + pKR[4]*matDuvDxyt[4] + pKR[7]*matDuvDxyt[5];
//			tmp1[5] = pKR[2]*matDuvDxyt[3] + pKR[5]*matDuvDxyt[4] + pKR[8]*matDuvDxyt[5];
//
//			matDuvDpa[0] = tmp1[0]*matDXkDpa[0] + tmp1[1]*matDXkDpa[1] + tmp1[2]*matDXkDpa[2];//u��parallax��һ�׵�
//			matDuvDpa[1] = tmp1[3]*matDXkDpa[0] + tmp1[4]*matDXkDpa[1] + tmp1[5]*matDXkDpa[2];//v��parallax��һ�׵�
//
//			pPB[0] = matDuvDDH[0];		pPB[1] = matDuvDDH[1];		pPB[2] = matDuvDpa[0];
//			pPB[3] = matDuvDDH[2];		pPB[4] = matDuvDDH[3];		pPB[5] = matDuvDpa[1];
//		
//			//Ti
//			matDdisDTi[0] = -matPosei2k[0]/dDisTi2k;
//			matDdisDTi[1] = -matPosei2k[1]/dDisTi2k;
//			matDdisDTi[2] = -matPosei2k[2]/dDisTi2k;
//
//			matDDotDTi[0] = -sin(ppt[0])*cos(ppt[1]);
//			matDDotDTi[1] = -sin(ppt[1]);
//			matDDotDTi[2] = -cos(ppt[0])*cos(ppt[1]);
//
//			matDArcCosInDTi[0] = (dDisTi2k*matDDotDTi[0]-dDot*matDdisDTi[0])/(dDisTi2k*dDisTi2k);
//			matDArcCosInDTi[1] = (dDisTi2k*matDDotDTi[1]-dDot*matDdisDTi[1])/(dDisTi2k*dDisTi2k);
//			matDArcCosInDTi[2] = (dDisTi2k*matDDotDTi[2]-dDot*matDdisDTi[2])/(dDisTi2k*dDisTi2k);
//
//			matDsinwWDTi[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[0];
//			matDsinwWDTi[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[1];
//			matDsinwWDTi[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTi[2];
//
//			matDTi2kDTi[0] = matDTi2kDTi[4] = matDTi2kDTi[8] = -1;
//			matDTi2kDTi[1] = matDTi2kDTi[2] = matDTi2kDTi[3] = 0;
//			matDTi2kDTi[5] = matDTi2kDTi[6] = matDTi2kDTi[7] = 0;
//
//			matDXkDTi[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[0] + dDisTi2k*matXj[0]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[0];
//			matDXkDTi[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[1] + dDisTi2k*matXj[0]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[1];
//			matDXkDTi[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTi[2] + dDisTi2k*matXj[0]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[2];
//			matDXkDTi[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[0] + dDisTi2k*matXj[1]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[3];
//			matDXkDTi[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[1] + dDisTi2k*matXj[1]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[4];
//			matDXkDTi[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTi[2] + dDisTi2k*matXj[1]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[5];
//			matDXkDTi[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[0] + dDisTi2k*matXj[2]*matDsinwWDTi[0] - sin(ppt[2])*matDTi2kDTi[6];
//			matDXkDTi[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[1] + dDisTi2k*matXj[2]*matDsinwWDTi[1] - sin(ppt[2])*matDTi2kDTi[7];
//			matDXkDTi[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTi[2] + dDisTi2k*matXj[2]*matDsinwWDTi[2] - sin(ppt[2])*matDTi2kDTi[8];
//
//			matDuvDTi[0] = tmp1[0]*matDXkDTi[0] + tmp1[1]*matDXkDTi[3] + tmp1[2]*matDXkDTi[6];
//			matDuvDTi[1] = tmp1[0]*matDXkDTi[1] + tmp1[1]*matDXkDTi[4] + tmp1[2]*matDXkDTi[7];
//			matDuvDTi[2] = tmp1[0]*matDXkDTi[2] + tmp1[1]*matDXkDTi[5] + tmp1[2]*matDXkDTi[8];
//			matDuvDTi[3] = tmp1[3]*matDXkDTi[0] + tmp1[4]*matDXkDTi[3] + tmp1[5]*matDXkDTi[6];
//			matDuvDTi[4] = tmp1[3]*matDXkDTi[1] + tmp1[4]*matDXkDTi[4] + tmp1[5]*matDXkDTi[7];
//			matDuvDTi[5] = tmp1[3]*matDXkDTi[2] + tmp1[4]*matDXkDTi[5] + tmp1[5]*matDXkDTi[8];
//
//			pAM[0] = matDuvDTi[0];		pAM[1] = matDuvDTi[1];		pAM[2] = matDuvDTi[2];//uv����ê��ƽ��������һ�׵�
//			pAM[3] = matDuvDTi[3];		pAM[4] = matDuvDTi[4];		pAM[5] = matDuvDTi[5];	
//		
//			//Tk
//			matDdisDTk[0] = matPosei2k[0]/dDisTi2k;
//			matDdisDTk[1] = matPosei2k[1]/dDisTi2k;
//			matDdisDTk[2] = matPosei2k[2]/dDisTi2k;
//
//			matDDotDTk[0] = sin(ppt[0])*cos(ppt[1]);
//			matDDotDTk[1] = sin(ppt[1]);
//			matDDotDTk[2] = cos(ppt[0])*cos(ppt[1]);
//
//			matDArcCosInDTk[0] = (dDisTi2k*matDDotDTk[0] - dDot*matDdisDTk[0])/(dDisTi2k*dDisTi2k);
//			matDArcCosInDTk[1] = (dDisTi2k*matDDotDTk[1] - dDot*matDdisDTk[1])/(dDisTi2k*dDisTi2k);
//			matDArcCosInDTk[2] = (dDisTi2k*matDDotDTk[2] - dDot*matDdisDTk[2])/(dDisTi2k*dDisTi2k);
//
//			matDsinwWDTk[0] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[0];
//			matDsinwWDTk[1] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[1];
//			matDsinwWDTk[2] = cos(dW+ppt[2])*dDWDArcCosIn*matDArcCosInDTk[2];
//
//			matDTi2kDTk[0] = matDTi2kDTk[4] = matDTi2kDTk[8] = 1;
//			matDTi2kDTk[1] = matDTi2kDTk[2] = matDTi2kDTk[3] = 0;
//			matDTi2kDTk[5] = matDTi2kDTk[6] = matDTi2kDTk[7] = 0;
//
//			matDXkDTk[0] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[0] + dDisTi2k*matXj[0]*matDsinwWDTk[0];
//			matDXkDTk[1] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[1] + dDisTi2k*matXj[0]*matDsinwWDTk[1];
//			matDXkDTk[2] = sin(dW+ppt[2])*matXj[0]*matDdisDTk[2] + dDisTi2k*matXj[0]*matDsinwWDTk[2];
//			matDXkDTk[3] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[0] + dDisTi2k*matXj[1]*matDsinwWDTk[0];
//			matDXkDTk[4] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[1] + dDisTi2k*matXj[1]*matDsinwWDTk[1];
//			matDXkDTk[5] = sin(dW+ppt[2])*matXj[1]*matDdisDTk[2] + dDisTi2k*matXj[1]*matDsinwWDTk[2];
//			matDXkDTk[6] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[0] + dDisTi2k*matXj[2]*matDsinwWDTk[0];
//			matDXkDTk[7] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[1] + dDisTi2k*matXj[2]*matDsinwWDTk[1];
//			matDXkDTk[8] = sin(dW+ppt[2])*matXj[2]*matDdisDTk[2] + dDisTi2k*matXj[2]*matDsinwWDTk[2];
//
//			matDuvDTk[0] = tmp1[0]*matDXkDTk[0] + tmp1[1]*matDXkDTk[3] + tmp1[2]*matDXkDTk[6];
//			matDuvDTk[1] = tmp1[0]*matDXkDTk[1] + tmp1[1]*matDXkDTk[4] + tmp1[2]*matDXkDTk[7];
//			matDuvDTk[2] = tmp1[0]*matDXkDTk[2] + tmp1[1]*matDXkDTk[5] + tmp1[2]*matDXkDTk[8];
//			matDuvDTk[3] = tmp1[3]*matDXkDTk[0] + tmp1[4]*matDXkDTk[3] + tmp1[5]*matDXkDTk[6];
//			matDuvDTk[4] = tmp1[3]*matDXkDTk[1] + tmp1[4]*matDXkDTk[4] + tmp1[5]*matDXkDTk[7];
//			matDuvDTk[5] = tmp1[3]*matDXkDTk[2] + tmp1[4]*matDXkDTk[5] + tmp1[5]*matDXkDTk[8];
//			
//			pAA[0] = matDuvDTk[0];		pAA[1] = matDuvDTk[1];		pAA[2] = matDuvDTk[2];//uv�Ը�ê��ƽ��������һ�׵�
//			pAA[3] = matDuvDTk[3];		pAA[4] = matDuvDTk[4];		pAA[5] = matDuvDTk[5];
//		
//			//Tl
//			matDXlDTl[0] = matDXlDTl[4] = matDXlDTl[8] = -sin(ppt[2]);
//			matDXlDTl[1] = matDXlDTl[2] = matDXlDTl[3] = 0;
//			matDXlDTl[5] = matDXlDTl[6] = matDXlDTl[7] = 0;	
//
//			matDuvDTl[0] = tmp1[0]*matDXlDTl[0] + tmp1[1]*matDXlDTl[3] + tmp1[2]*matDXlDTl[6];
//			matDuvDTl[1] = tmp1[0]*matDXlDTl[1] + tmp1[1]*matDXlDTl[4] + tmp1[2]*matDXlDTl[7];
//			matDuvDTl[2] = tmp1[0]*matDXlDTl[2] + tmp1[1]*matDXlDTl[5] + tmp1[2]*matDXlDTl[8];
//			matDuvDTl[3] = tmp1[3]*matDXlDTl[0] + tmp1[4]*matDXlDTl[3] + tmp1[5]*matDXlDTl[6];
//			matDuvDTl[4] = tmp1[3]*matDXlDTl[1] + tmp1[4]*matDXlDTl[4] + tmp1[5]*matDXlDTl[7];
//			matDuvDTl[5] = tmp1[3]*matDXlDTl[2] + tmp1[4]*matDXlDTl[5] + tmp1[5]*matDXlDTl[8];
//		
//			pPA[0] = matDuvDRA[0];			pPA[1] = matDuvDRB[0];			pPA[2] = matDuvDRG[0];//uv������ê�������һ�׵�
//			pPA[3] = matDuvDTl[0];			pPA[4] = matDuvDTl[1];			pPA[5] = matDuvDTl[2];
//			pPA[6] = matDuvDRA[1];			pPA[7] = matDuvDRB[1];			pPA[8] = matDuvDRG[1];
//			pPA[9] = matDuvDTl[3];			pPA[10] = matDuvDTl[4];			pPA[11] = matDuvDTl[5];
//		}
}

CLCRP::CLCRP(void)
{
	m_e1 = m_e2 = m_e3 = 1E-10;
	m_e4 = 0;	
	m_bProvideXYZ = false;
	m_bFocal = false;
	m_szCameraInit = m_szFeatures = m_szCalibration = m_szXYZ = m_sz3Dpts = m_szCamePose = m_szReport = NULL;

	m_nMaxIter = 100;
	m_Tau      = 1E-6;	
	m_bRobustKernel = false;
	m_nRobustType = 1;
	m_bsolverLM = true;
	m_bsolverGN = false;
	m_delt      = 1;
}


CLCRP::~CLCRP(void)
{
	if ( m_szCameraInit!= NULL)	free(m_szCameraInit);
	if ( m_szFeatures!= NULL)	free(m_szFeatures);
	if ( m_szCalibration!= NULL)	free(m_szCalibration);
	if ( m_szXYZ!= NULL)		free(m_szXYZ);
	if ( m_szCamePose!= NULL)	free(m_szCamePose);
	if ( m_sz3Dpts!= NULL)		free(m_sz3Dpts);
	if ( m_szReport!= NULL)		free(m_szReport);
}

bool CLCRP::lcrp_run( bool bRobust,
	bool bLM,
	int nMaxIter,
	char* szCam,
	char* szFea,
	char* szXYZ,
	char* szCalib,
	char* szReport,
	char* szPose,
	char* sz3D,
	double Tau )
{
	m_szCameraInit = szCam;                                                       //����ⷽλԪ��
	m_szFeatures   = szFea;                                                       //ƥ���
	m_szCalibration= szCalib;                                                     //����ڷ�λԪ��
	m_szXYZ        = szXYZ;                                                       //
	m_bRobustKernel= bRobust;                                                     //
	m_bsolverGN    = !bLM;                                                        //
	m_nMaxIter     = nMaxIter;                                                    //
	m_szCamePose   = szPose;                                                      //
	m_sz3Dpts      = sz3D;                                                        //
	m_szReport     = szReport;                                                    //
	m_Tau          = Tau;                                                         //

	lcrp_initialize( m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ );     //������ֵ

	if (m_bsolverGN)
		lcrp_motstr_gn( );
	else
	{
		//double t1 = clock();
		lcrp_motstr_levmar();
		//double t2 = clock();
		//double t_diff = (t2 - t1) / CLOCKS_PER_SEC;
		//printf("%s %f\n", "��������ʼ���ڴ��ƽ��ʱ�䣺", t_diff);
	}
	return true;
}

bool CLCRP::lcrp_run( int argc, char** argv )
{
	bool bTrue = lcrp_parseArgs( argc, argv );
	if (!bTrue)
	{
		fprintf( stderr, "LCRP: Input wrong commands, please check them!\n");
		return false;
	}
	
	lcrp_initialize( m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ );

	if (m_bsolverGN)		
		lcrp_motstr_gn( );
	else					
	{
		lcrp_motstr_levmar();
	}
	return true;
}


bool CLCRP::lcrp_initialize( char* szCamera, char* szFeature,  char* szCalib, char* szXYZ )
{
	printf("LCRP: Light Cone Reprojection Principle Version 1.0\n");
	FILE* fp;

	//must input initial initial camera pose file and projection image points file
	fp = fopen( szCamera, "r" );
	if ( fp == NULL )
	{
		fprintf( stderr, "LCRP: Missing initial camera poses file! \n");
		exit(1);
	}
	else
		fclose(fp);

	fp = fopen( szFeature, "r" );
	if ( fp == NULL )
	{	
		fprintf( stderr, "LCRP: Missing feature projection points file! \n");
		exit(1); 
	}
	else
		fclose(fp);

	if ( szCalib != NULL )
	{
		m_bFocal = false;																											
		m_K = (double*)malloc(9*sizeof(double));                           //�������ڷ�λԪ��
		lcrp_readCameraPoseration(szCalib, m_K);                            //
		///*-------------------------------------------------������---------------------------------------------------*/
		//for (int i = 0;i < 9;i++)
		//{
		//	printf("%f ", m_K[i]);
		//	if ((i + 1) % 3 == 0)
		//		printf("\n");
		//}
		///*-------------------------------------------------������---------------------------------------------------*/
	}	

	if ( szXYZ != NULL )
		m_bProvideXYZ = true;	
	zu = 0;
	savePara = 0;//Add by Lu
	//read camera pose & features images projs, and initialize features points( three kinds of angle )
	lcrp_readAndInitialize( szCamera, szFeature, &m_ncams, &m_n3Dpts, &m_n2Dprojs,&m_motstruct,//number of camera, 3D points, 2D projection points,6 camera pose and 3 feature parameters
		&m_imgpts, &m_lc, &m_archor, &m_vmask, &m_umask, &m_photo, &m_feature, &m_archorSort);
	
	//for (int i = 0; i < 10; i++)
	//{
	//	printf("%f %f %f\n", m_imgpts[2 * i], m_imgpts[2 * i + 1], m_lc[i]);
	//}
	if(savePara == 1)//Add by Lu
		lcrp_saveInitialParallax("C:/txtdata/Village/InitialParallax.txt", m_motstruct);//Add by Lu
	//Modify the directory to export point file
	if (zu == 1)
		lcrp_saveTriangulatedxyz("C:/txtdata/Village/Triangulatedxyz.txt", m_motstruct);
	else if (zu == 2)//Add by Lu
	{
		lcrp_saveInitialXYZ("C:/txtdata/Village/xyz_Maxdw.txt", lcrp_angle2xyz(m_motstruct));//ת��Ϊxyz
	}
	else//Add by Lu
		;//Other codes

	printf( "Number of cameras: %d\n", m_ncams );
	printf( "Number of points: %d\n", m_n3Dpts );
	printf( "Number of projections: %d\n", m_n2Dprojs );
	return true;
}
double *CLCRP::lcrp_angle2xyz(double* p)
{
	double *n3DptsT = (double*)malloc((m_n3Dpts * 3) * sizeof(double));
	static int i, j;
	double* pAngle;
	double xj[3], xk[3];
	double Tik[3];
	int nM, nA;
	double Dik;
	double w, w2;
	int cnp = 6, pnp = 3;
	for (i = 0; i < m_n3Dpts; i++)
	{
		pAngle = p + m_ncams * 6 + i * 3;//azimuth angle
		w = *(p + m_ncams * 6 + i * 3 + 2);//parallax angle

		xj[0] = sin(*(pAngle)) * cos(*(pAngle + 1));//x
		xj[1] = sin(*(pAngle + 1));//y
		xj[2] = cos(*(pAngle)) * cos(*(pAngle + 1));//z

		nM = *(m_archor + i * 3 + 1);
		nA = *(m_archor + i * 3 + 2);

		Tik[0] = -*(p + nM * 6 + 3) + *(p + nA * 6 + 3);
		Tik[1] = -*(p + nM * 6 + 4) + *(p + nA * 6 + 4);
		Tik[2] = -*(p + nM * 6 + 5) + *(p + nA * 6 + 5);

		Dik = sqrt(Tik[0] * Tik[0] + Tik[1] * Tik[1] + Tik[2] * Tik[2]);

		w2 = acos((xj[0] * Tik[0] + xj[1] * Tik[1] + xj[2] * Tik[2]) / Dik);

		//�������������ê���xyz����
		xk[0] = (Dik * sin(w2 + w) * xj[0]) / sin(w);
		xk[1] = (Dik * sin(w2 + w) * xj[1]) / sin(w);
		xk[2] = (Dik * sin(w2 + w) * xj[2]) / sin(w);

		//�������ȫ��xyz����
		*(n3DptsT + i * 3) = *(p + nM * 6 + 3) + xk[0];
		*(n3DptsT + i * 3 + 1) = *(p + nM * 6 + 4) + xk[1];
		*(n3DptsT + i * 3 + 2) = *(p + nM * 6 + 5) + xk[2];
		if (i ==2)
		{
			printf("%f %f %f\n", *(n3DptsT + i * 3), *(n3DptsT + i * 3 + 1), *(n3DptsT + i * 3 + 2));
		}
	}

	return n3DptsT;
}
void CLCRP::lcrp_saveInitialXYZ(char* sz3Dpt, double* p)
{
	static int i = 0;
	double dx, dy, dz;
	FILE* fp = NULL, * fpc = NULL;
	int nM, nA;
	//save features xyz
	if (sz3Dpt != NULL)
	{
		fp = fopen(sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			dx = *(p + i * 3);
			dy = *(p + i * 3 + 1);
			dz = *(p + i * 3 + 2);
			nM = *(m_archor + i * 3 + 1);
			nA = *(m_archor + i * 3 + 2);
			//fprintf(fp, "%0.5lf %0.5lf %0.5lf %d %d\n", dx, dy, dz,nM, nA);
			//fprintf(fp, "%d %0.5lf %0.5lf %0.5lf\n", i+1,dx, dy, dz);
			fprintf(fp, "%0.5lf %0.5lf %0.5lf\n", dx, dy, dz);
			if (i == 2)
			{
				printf("%0.5lf %0.5lf %0.5lf\n", dx, dy, dz);
				printf("%d\n", m_n3Dpts);
			}
				
		}
		printf("%d\n", i);
		fclose(fp);
	}
	free(p);
}
void CLCRP::lcrp_saveInitialParallax(char* sz3Dpt, double* p)
{
	//double* n3DptsT = (double*)malloc((m_n3Dpts * 3) * sizeof(double));

	static int i = 0;
	double azimuth, elevation, parallax;
	FILE* fp = NULL, * fpc = NULL;
	double* pAngle;
	//save features xyz
	if (sz3Dpt != NULL)
	{
		fp = fopen(sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			pAngle = p + m_ncams * 6 + i * 3;

			azimuth = *pAngle;
			elevation = *(pAngle + 1);
			parallax = *(pAngle + 2);
			fprintf(fp, "%0.5lf     %0.5lf     %0.5lf\n", azimuth, elevation, parallax);
			//fprintf(fp, "%0.5lf\n", parallax*180/PI);
		}
		fclose(fp);
	}
	free(p);
}
void CLCRP::lcrp_saveTriangulatedxyz(char* sz3Dpt, double* p)
{
	static int i = 0;
	double x, y, z;
	FILE* fp = NULL;
	//save features xyz
	if (sz3Dpt != NULL)
	{
		fp = fopen(sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			x = *(p + m_ncams * 6 + i * 3);
			y = *(p + m_ncams * 6 + i * 3 + 1);
			z = *(p + m_ncams * 6 + i * 3 + 2);
			//fprintf(fp, "%d %0.5lf %0.5lf %0.5lf\n", i+1, x, y, z);
			fprintf(fp, "%0.5lf %0.5lf %0.5lf\n", x, y, z);
			//fprintf(fp, "%0.5lf\n", parallax * 180 / PI);
		}
		fclose(fp);
	}
	free(p);
}

bool CLCRP::lcrp_motstr_levmar( )
{
	if ( m_motstruct == NULL || m_imgpts == NULL || m_K == NULL )
	{	
		fprintf( stderr, "LCRP: Missing related source file, please input them first\n");	
		return false;
	}
	
	FILE *fpRe = NULL;
	if (m_szReport!= NULL)
	{
		fpRe = fopen(m_szReport, "w ");
		fprintf(fpRe, "%d  poses, %d 3D features, %d projection\n", m_ncams, m_n3Dpts, m_n2Dprojs );
		fprintf(fpRe, "Levenberg-Marquardt is used\n" );
	}

	int i, ii;
	int n = m_n3Dpts, m = m_ncams, cnp = 6, pnp = 3, mnp = 2;
	int nvis, nuis, itno, issolved, nobs, nvars, nMaxS, nu=2, stop=0;

	double *p = m_motstruct, *x = m_imgpts;	//p pointer refers to unknown parameters, x pointer refers to image coordinate	
	double *U, *V, *W, *e, *eab, *E, *S, *dp, *IV; // pointers into U V W E S IV
	double *pa = NULL, *pb = NULL, *ea, *eb, *dpa, *dpb, *hx, *Ex, *rx, *pdp; // pointers into p, jac, eab and dp respectively 	double *Ex, *rx; 
	double initialerror = 0, error = 0;
	int    nIter = 0, nLinear = 0;
	int Usz, Vsz, Wsz, easz, esz, ebsz, Sblsz, Sdim; 
	static double mu, tmp, p_eL2, eab_inf, pdp_eL2, x2, delt2, init_p_eL2 = 0; 
	double p_L2, dp_L2=DBL_MAX, dF, dL;
	double tau=fabs(m_Tau), eps1=fabs(m_e1),eps3_sq=m_e3*m_e3;
	bool init = false, ordering = true;
	double tStart, tEnd, tTimeUse, t0, t1, t11, t2, t3, t4, t5;
	struct sba_crsm Sidxij, Uidxij;		// S mask and U mask
	double tCons = 0, tSolve = 0, tCost = 0, tTimeIter = 0, tTimeCal = 0;
	int nLmIterations = 0;

	Usz=cnp*cnp;												//6*6
	Vsz=pnp * pnp;												//3*3
	Wsz=cnp * pnp;												//6*3
	esz=mnp; easz=cnp; ebsz=pnp;								//2 �� 6 �� 3
	Sblsz=cnp*cnp;												//6*6
	Sdim=m * cnp;												//m_ncams�������Ŀ��*6
	nvis = m_n2Dprojs;											//��ά���Ӧ��2DͶӰ������
	mu=eab_inf=0.0;												//0.0

	nuis = lcrp_ConstructSmask( Sidxij, Uidxij );
	nobs=nvis*mnp;
	nvars=m*cnp + n*pnp;
			
	S	=	(double *)malloc(m_nS*Sblsz*sizeof(double));
	W	=	(double *)malloc(nvis*Wsz*sizeof(double));
	U	=	(double *)malloc(nuis*Usz*sizeof(double));
	V	=	(double *)malloc(n*Vsz*sizeof(double));
	IV	=	(double *)malloc(n*Vsz*sizeof(double));
	e	=	(double *)malloc(nobs*sizeof(double));
	eab	=	(double *)malloc(nvars*sizeof(double));
	E	=	(double *)malloc(m*cnp*sizeof(double));		
	dp	=	(double *)malloc(nvars*sizeof(double));
	hx	=	(double *)malloc(nobs*sizeof(double));	
	pdp	=	(double *)malloc(nvars*sizeof(double));	
	pa=p; pb=p+m*cnp; ea=eab; eb=eab+m*cnp;	dpa=dp; dpb=dp+m*cnp;	

	cholmod_start (&m_cS) ;  	
	//m_cS.print_function = NULL;	
	int *Ap  = (int*)malloc((m+1)*sizeof(int));
	int * Aii = (int*)malloc(m_nS*sizeof(int));
	lcrp_constructAuxCSSLM( Ap, Aii );

	m_cholSparseE = cholmod_zeros( cnp*m, 1, CHOLMOD_REAL, &m_cS);
	Ex = (double*)m_cholSparseE->x;
	nMaxS = (m_nS-m_ncams)*36+m_ncams*21;	//maximum non-zero element in S matrix 

	m_cholSparseS = cholmod_allocate_sparse(m_ncams*6,m_ncams*6,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
	int *Sp, *Si;
	double* Sx = NULL;
	Sp = (int*)m_cholSparseS->p;		//column pointer
	Si = (int*)m_cholSparseS->i;		//row pointer
	
	//Compute initial error and initial reprojection error 
	tStart = clock();
	lcrp_cost(p, hx, m_archor ); 

	//for (int i = 0; i < nobs; i++)
	//{
	//	if(x[i]== 2027.86||x[i] == 456.817)
	//		printf("%f  %f\n", x[i], hx[i]);
	//}
	

	p_eL2 = nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */

	initialerror = p_eL2/nvis;
	printf("Initial Error %0.1lf [%0.8lf]\n", p_eL2, initialerror);

	if( m_szReport != NULL )
		fprintf( fpRe, "Initial Error  %lf\n", initialerror );

	init_p_eL2 = p_eL2;

	//Iteration 
	for(itno=0; itno<m_nMaxIter && !stop; ++itno)
	{
		//Setup S matrix, include two step
		memset( U, 0, nuis*Usz*sizeof(double) );
		memset( ea, 0,m*easz*sizeof(double) );
		memset( V, 0, n*Vsz*sizeof(double));
		memset( IV, 0, n*Vsz*sizeof(double));
		memset( eb, 0, n*ebsz*sizeof(double));
		memset( W, 0, nvis*Wsz*sizeof(double));	
						
		//Step one; compute W V U directly, don't save each projection image Jacobian 
		t0 = clock();
		//��һ��.1���Ǽ����ſ˱Ⱦ���
		//��������Ÿ������� U��V��W ������������������������� 3D ��Ĳ��ֵ�����
		if ( m_bRobustKernel)
			lcrp_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		else
			lcrp_jacobian(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		t1 = clock();
		//��һ��.2����ȷ�������������������������Ʋ�����
		if( itno == 0)//iter No����һ�ξ���0�Ϳ϶���ִ��
			mu = lcrp_computeInitialmu( U, V, Uidxij, tau, nvars );

		nLmIterations = 0;//������0��ʼ
		while(1) //determine increment using adaptive damping ������������pba_computeInitialmu()����
		{//�Ÿ��������ʾ��Ͷ�������ڲ����ĵ�����
			nLmIterations++;
			t11 = clock();
			//��V�������V^-1
			lcrp_inverseVLM( V, IV, Uidxij, mu ); //compute inverse matrix of V
			
			//Step two: construct S matrix using U V W, S = U - W*V^-1*W^T
			memset( E, 0, m*easz*sizeof(double));
			memset( S, 0, m_nS*Sblsz*sizeof(double) );
			lcrp_constructSLM( S, E, U, IV, W, ea, eb, Sidxij, mu );
			t2 = clock();
			
			//Solve linear equation
			//set CSS format using S matrix
			lcrp_constructCSSLM( Si, Sp, Sx, S, m_cholSparseS, Sidxij, init ); 
			for ( ii = 0; ii < cnp*m; ii++  )
				Ex[ii] = E[ii];	
			
			lcrp_solveCholmodLM( Ap, Aii, init, ordering);
			nLinear++;
			
			init = true;
			rx = (double*)m_cholSparseR->x;

			if (m_cS.status != CHOLMOD_NOT_POSDEF )
			{
				for ( ii = 0; ii < cnp*m; ii++ )
					dpa[ii] = rx[ii];
				issolved = 1;
			}
			else
				issolved = 1;
			
			t3 = clock();
						
			if(issolved)
			{
				//Solve features
				lcrp_solveFeatures( W, IV, ea, eb, dpa, dpb );
			
				// Compute ||J^T e||_inf and ||p||^2 
				for(i=0, p_L2=eab_inf=0.0; i<nvars; ++i)
				{
					if(eab_inf < (tmp=fabs(eab[i])))
						eab_inf=tmp;
					
					p_L2+=p[i]*p[i];
				}								
	
				//update
				for(i=0, dp_L2=0.0; i<nvars; ++i)// compute p's new estimate and ||dp||^2 
				{
					pdp[i]=p[i] + (tmp=dp[i]);
					dp_L2+=tmp*tmp;
				}

				//m_e1 = 0;
				if (sqrt(dp_L2)<=m_e1*(sqrt(p_L2)+m_e1))
				{	stop = 1;	break;	}

				lcrp_updateKR(m_R, m_KR, m_KdA, m_KdB, m_KdG, m_K, pdp );
				t4 = clock();
				
				lcrp_cost(pdp, hx, m_archor );
				pdp_eL2=nrmL2xmy(hx, x, hx, nobs); 	
				error = pdp_eL2/nvis;
				t5 = clock();
								
				if ( m_bRobustKernel )
				{
					pdp_eL2 = 0;
					delt2 = m_delt*m_delt;
					if ( m_nRobustType==1)						//Cauchy Kernel Function
					{
						for ( i = 0; i < m_n2Dprojs; i++ )
						{					
							x2 = hx[i*2]*hx[i*2]+hx[i*2+1]*hx[i*2+1];
							x2 = delt2*log( x2/delt2 + 1 );
							pdp_eL2 += x2;
						}
					}
					else										//Huber Kernel Function
					{
						for ( i = 0; i < m_n2Dprojs; i++ )
						{					
							x2 = hx[i*2]*hx[i*2]+hx[i*2+1]*hx[i*2+1];
							
							if (x2 <= delt2)  // inlier
								x2 = x2;
							else  // outliers
								x2 = 2*sqrt(x2)*m_delt - delt2;
								
							pdp_eL2 += x2;
						}
					}
					error = pdp_eL2/nvis;
				}				

				for(i=0, dL=0.0; i<nvars; ++i)
					dL+=dp[i]*(mu*dp[i]+eab[i]);  //low  
				dF=p_eL2-pdp_eL2;	

				if((dF/dL)>0.0)
				{ 
					if((sqrt(p_eL2)-sqrt(pdp_eL2))<m_e3*sqrt(p_eL2)) 
					{	stop=2;		break;	}

					for(i=0; i<nvars; ++i) 
						p[i]=pdp[i];
					for(i=0; i<nobs; ++i) 
						e[i]=hx[i];		
			

					p_eL2=pdp_eL2;
					if((eab_inf <= eps1))
					{	dp_L2=0.0; 		stop=4;		break;	}

					tmp=(2.0*dF/dL-1.0);
					tmp=1.0-tmp*tmp*tmp;
					mu=mu*( (tmp>=1.0/3.0)? tmp : 1.0/3.0 );
					nu=2;

					tTimeIter = (t5 - t0)*0.001;
					tTimeCal += tTimeIter;

					printf( "Iteration=%d  MSE=%0.8lf   LmIters=%d  Pertime=%0.2lf TotalTime=%0.2lf\n", itno, pdp_eL2/nvis, nLmIterations, tTimeIter, tTimeCal ); 					
					if( m_szReport!= NULL )
						fprintf( fpRe, "Iteration %d  Error  %0.8lf\n", itno, pdp_eL2/nvis );
					nIter++;

					break;
				}
				else
				{
					mu*=nu;
					nu *= 2;
				}
			} 			
		}
		
		if(p_eL2<=eps3_sq) stop=5; 
	}
	if(itno>=m_nMaxIter)
		stop=3;	
	
	//clear memory and print
iterstop:
	sba_crsm_free(&Uidxij);
	sba_crsm_free(&Sidxij);

	cholmod_free_factor(&m_cholFactorS, &m_cS) ;              
	cholmod_l_free_dense(&m_cholSparseE, &m_cS);
	cholmod_l_free_dense(&m_cholSparseR, &m_cS);
	cholmod_finish (&m_cS) ;  
	free(Ap);
	free(Aii);
	cholmod_free_sparse(&m_cholSparseS, &m_cS) ;
	
	tEnd  = clock();
	tTimeUse = tEnd - tStart;

	//save optimal camera pose and feature
	//lcrp_saveXYZ( m_szCamePose, m_sz3Dpts, p, false );

	free(S);	free(W);	free(U);	free(V);	free(IV);
	free(e);	free(eab);	free(E);   	free(dp);	free(hx);	free(pdp);	
	free(m_KR); free(m_KdA);free(m_KdB);free(m_KdG);	

	printf( "%d parameters, %d observations, Levenberg-Marquardt, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
		nvars, m_n2Dprojs*2, stop, p_eL2/nvis, initialerror, nIter+1, nLinear, tTimeUse / CLOCKS_PER_SEC );
	printf( "LCRP reasons is listed as following:\n" );
	printf( "reason 1: relative change of state vector is small\n" );
	printf( "reason 2: relative change of projection error is small\n" );
	printf( "reason 3: maximum iteration\n");
	printf( "reason 4: maximum value of b is small\n");
	printf( "reason 5: total reprojection error is small\n");

	if( m_szReport!=NULL )
	{
		fprintf( fpRe, "%d parameters, %d observations, Levenberg-Marquardt, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
			nvars, m_n2Dprojs*2, stop, error, initialerror, nIter+1, nLinear, tTimeUse*0.001 );

		fprintf( fpRe, "LCRP reasons is listed as following:\n" );
		fprintf( fpRe, "reason 1: relative change of state vector is small\n" );
		fprintf( fpRe, "reason 2: relative change of projection error is small\n" );
		fprintf( fpRe, "reason 3: maximum iteration\n");
		fprintf( fpRe, "reason 4: maximum value of b is small\n");
		fprintf( fpRe, "reason 5: total reprojection error is small\n");

		fclose(fpRe);
	}	
	
	return true;
}

//#pragma optimize("", off) // �ر��Ż�
bool CLCRP::lcrp_motstr_gn( FixType ft )
{	
	if ( m_motstruct == NULL || m_imgpts == NULL || m_K == NULL )
	{	
		printf( "LCRP: Missing related file, please input them first\n");	
		return false;
	}

	FILE *fpRe = NULL;
	if( m_szReport != NULL )
	{
		fpRe = fopen(m_szReport, "w ");
		fprintf(fpRe, "%d  poses, %d 3D features, %d projection\n", m_ncams, m_n3Dpts, m_n2Dprojs );
		fprintf(fpRe, "Gauss-Newton is used\n" );
	}

	int i, ii, nft;
	int n = m_n3Dpts, m = m_ncams, cnp = 5, pnp = 3, mnp = 1;
	int nvis, nuis, itno, nobs, nvars, nMaxS;

	double *p = m_motstruct, *x = m_lc;	//p pointer refers to unknown parameters, x pointer refers to image coordinate	
	double *U, *V, *W, *e, *eab, *E, *S, *dp; // pointers into U V W E S IV
	double *pa, *pb, *ea, *eb, *dpa, *dpb, *hx, *Ex, *rx; // pointers into p, jac, eab and dp respectively 	double *Ex, *rx; 
	double initialerror = 0, error = 0, lasterror = 0;
	int Usz, Vsz, Wsz, esz, easz, ebsz, Sblsz, Sdim; 

	double tmp, tmp1, init_p_eL2, p_eL2, eab_inf, pdp_eL2, x2, delt2; 
	double dp_L2=DBL_MAX;
	int nu=2, stop=0;
	eab_inf=0.0;		

	bool init = false;
	bool ordering = true;
	double tStart, tEnd, tTimeUse, tPerTime, t0, t1, t2, t3, t4, t5, tTotalTime = 0;
	struct sba_crsm Sidxij, Uidxij;		// S mask and U mask          �洢ϡ�����
	nuis = lcrp_ConstructSmask( Sidxij, Uidxij );   //����U�������Ԫ�ظ���

	Usz=cnp*cnp; Vsz=pnp * pnp; Wsz=cnp * pnp; 	
	esz=mnp; easz=cnp; ebsz=pnp; Sblsz=cnp*cnp;	Sdim=m * cnp;	
	nvis = m_n2Dprojs;                                            //������
	nobs=nvis*mnp;                                                //�������������۲ⷽ�̸���
	nvars=m*cnp + n*pnp;                                          //δ֪��������

	S	=	(double *)malloc(m_nS*Sblsz*sizeof(double));          //S���󣺷�0Ԫ�ظ���*5*5
	W	=	(double *)malloc(nvis*Wsz*sizeof(double));            //W���������*5*3
	U	=	(double *)malloc(nuis*Usz*sizeof(double));            //U���󣺷�0Ԫ�ظ�����5*5
	V	=	(double *)malloc(n*Vsz*sizeof(double));               //V������ά������3*3
	e	=	(double *)malloc(nobs*sizeof(double));                //����� 
	eab	=	(double *)malloc(nvars*sizeof(double));               //δ֪����
	E	=	(double *)malloc(m*cnp*sizeof(double));	              //�ⷽλԪ��	
	dp	=	(double *)malloc(nvars*sizeof(double));               //δ֪����
	hx	=	(double *)malloc(nobs*sizeof(double));	              //�۲ⷽ��
	pa=p; //paָ��δ֪����
	pb=p+m*cnp; //pbָ����ά�����

	ea=eab;//eaָ��
	eb=eab+m*cnp;

	dpa=dp; 
	dpb=dp+m*cnp;
	
	//Select fix axis for Gauss-Newton
	//���� BA ֻ��������͵� ���λ�ã�������̶�ĳЩ�������Ż������п��ܻᷢ�� ����Ư�� (Drift)����������������������ƽ�ơ���ת�����ţ�������ı�ͶӰ�����½ⲻΨһ�����̶���һ֡�����һ֡���᣿
	if ( ft == LCRP_FixDefault )
	{
		//double test1 = *(m_motstruct + 7);
		//double test2 = *(m_motstruct + 2);
		//double test3 = *(m_motstruct + 8);
		//double test4 = *(m_motstruct + 3);
		//double test5 = *(m_motstruct + 9);
		//double test6 = *(m_motstruct + 4);
		tmp = abs(*(m_motstruct+7)-*(m_motstruct+2));//��2�������Xc-��1�������Xc
		ft = LCRP_FixX;									
		tmp1 = abs(*(m_motstruct+8)-*(m_motstruct+3));//��2�������Yc-��1�������Yc
		if ( tmp1 > tmp )
		{ft= LCRP_FixY; tmp = tmp1;	}
		tmp1 = abs(*(m_motstruct+9)-*(m_motstruct+4));//��2�������Zc-��1�������Zc
		if ( tmp1 > tmp )
		{ft= LCRP_FixZ; }
	}
	nft = ft;
	
	// ʹ�� CHOLMOD ����ϡ��������ӷֽ⣨Cholesky �ֽ⣩��׼������
	cholmod_start (&m_cS) ;  	//��ʼ�� CHOLMOD �ڲ������ݽṹ����Ҫ�Ƿ��乤���ռ䣬������Ĭ�ϲ���
	//m_cS.print_function = NULL;	
	//�ڽ��㷽�̵�ʱ���õ�
	int *Ap  = (int*)malloc((m)*sizeof(int));
	int * Aii = (int*)malloc(m_nS*sizeof(int));//m_nS��S�����з���Ԫ������
	lcrp_constructAuxCSSGN( Ap, Aii );//����ϡ�����������ṹ����1��1�п�ʼ���洢ÿ�е���ʼ����,�洢����Ԫ�ص�������

	//�� CHOLMOD �m_cholSparseE ��һ�� cholmod_dense �ṹ��
	//typedef struct cholmod_dense
	//{
	//	size_t nrow;  // ����
	//	size_t ncol;  // ����
	//	size_t nzmax; // ������Ԫ������ͨ���� nrow * ncol��
	//	void* x;      // ָ�����ݵ�ָ��
	//} cholmod_dense;

	m_cholSparseE = cholmod_zeros( cnp*m-6, 1, CHOLMOD_REAL, &m_cS);//�̶���1֡����κ͵�2֡��Xc�ᣬʵ��δ֪����������5*m-6
	Ex = (double*)m_cholSparseE->x;
	
	//(m_nS-m_ncams)��ʾS�����зǶԽǿ飨���֮���Լ���������������е�ÿ���ǶԽǿ�ռ��5*5������Ԫ�أ����Գƣ�
	//m_ncams��ʾS�����жԽǿ飨������������Լ����������������ÿ���Խǿ�ռ��15������Ԫ�أ��Գƣ�ֻ�洢�����ǲ��֣�
	nMaxS = (m_nS-m_ncams)*25+m_ncams*15;	//maximum non-zero element in S matrix 

	//����ϡ�����S����4������Ϊtrue��ʾ�����������򣬵�5������Ϊtrue������ܴ洢����6������Ϊ1�������Ĵ洢����Ϊʵ����
	//���ܴ洢��Packed Storage��ָ����ֻ�洢�����з���Ԫ�أ�����������������������ڴ��У���Ԥ������ռ䡣�������Լ����ڴ�ռ�ã���߼���Ч�ʡ�
	m_cholSparseS = cholmod_allocate_sparse(m_ncams*5-6,m_ncams*5-6,nMaxS,true,true,1,CHOLMOD_REAL,&m_cS);
	int *Sp, *Si;
	double* Sx = NULL;
	Sp = (int*)m_cholSparseS->p;		//column pointer
	Si = (int*)m_cholSparseS->i;		//row pointer
	
	//Compute initial error and initial reprojection error 
	tStart = clock();
	//dpa[0] = dpa[1] = dpa[2] = dpa[3] = dpa[4] = dpa[5] = 0;//δ�õ�
	//dpa[9+nft] = 0;
	lcrp_cost(p, hx, m_archor );
	//double sum = 0;
	//for (int i = 0; i < nvis; i ++)
	//{
	//	printf("%f %f\n", hx[i], x[i]);
	//	sum += pow(hx[i] - x[i], 2);
	//}
	p_eL2=nrmL2xmy(e, x, hx, nobs); /* e=x-hx, p_eL2=||e|| */
	//printf("%f %f\n", sum, p_eL2);
	initialerror = p_eL2/nvis;
	lasterror = initialerror;//Add by zuo
	printf("Initial Error %0.1lf [%0.8lf]\n", p_eL2, initialerror);

	if( m_szReport!= NULL )
		fprintf( fpRe, "Initial Error  %0.8lf\n", initialerror );
	init_p_eL2=p_eL2;	
	
	//Iteration 
	for(itno=0; itno<m_nMaxIter && !stop; ++itno)
	{
		//Setup S matrix, include two step
		memset( U, 0, nuis*Usz*sizeof(double) );//nuis��U����ķ���Ԫ�ظ�����Usz=5*5���洢U����ķ���Ԫ��
		memset( ea, 0,m*easz*sizeof(double) );//m�����������easz=5��ea�洢�������
		memset( V, 0, n*Vsz*sizeof(double));//n����ά��������Vsz=3*3��
		memset( eb, 0, n*ebsz*sizeof(double));//n����ά��������ebsz=3��eb�洢��ά�����
		memset( W, 0, nvis*Wsz*sizeof(double));	//nvis���������Wsz=5*3

		//Step one; compute W V U directly, don't save each projection image Jacobian 
		t0 = clock();
		if ( m_bRobustKernel)
			lcrp_jacobian_RobustKernel(p, m_archor, &Uidxij,e,U,ea,V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature );
		else
			lcrp_jacobian(p, m_archor, &Uidxij, e, U, ea, V, eb, W, n, m, m_n2Dprojs, m_photo, m_feature);


		t1 = clock();

		lcrp_inverseVGN( U, V, Uidxij ); //compute inverse matrix of V

		//Step two: construct S matrix using U V W, S = U - W*V^-1*W^T
		memset( E, 0, m*easz*sizeof(double));
		memset( S, 0, m_nS*Sblsz*sizeof(double) );
		lcrp_constructSGN( S, E, U, V, W, ea, eb, Sidxij );//SΪ������࣬EΪ�����Ҳ�
		t2 = clock();
	
		//Solve equation
		lcrp_constructCSSGN( Si, Sp, Sx, S, m_cholSparseS, Sidxij, init, nft ); //set CSS format using S matrix//break
		for ( ii = 0; ii < 2+nft; ii++ )
			Ex[ii] = E[5+ii];//��2֡����̬��

		for ( ii = 0; ii < 5*m-(8+nft); ii++  )
			Ex[2+nft+ii] = E[8+nft+ii];//������2֡��Xc��ֻ�洢Yc��Zc

		//SxΪ��࣬ExΪ�Ҳ�

		lcrp_solveCholmodGN( Ap, Aii, init, ordering);//break
		init = true;
		rx = (double*)m_cholSparseR->x;//������������

		if (m_cS.status == CHOLMOD_NOT_POSDEF)
		{
			printf ( "Cholesky failure, writing debug.txt (Hessian loadable by Octave)" );
			goto iterstop;
		}
		else
		{
			for ( ii = 0; ii < 2+nft; ii++ )
				dpa[5+ii] = rx[ii];//�洢��2֡��2����̬
			for ( ii = 0; ii < 5*m-(8+nft); ii++ )
				dpa[ii+8+nft] = rx[2+nft+ii];//������2֡��Xc���洢������������
		}			
		//for (int i = 0; i < m; i++)
		//	printf("%f %f %f %f %f\n", dpa[5 * i], dpa[5 * i + 1], dpa[5 * i + 2], dpa[5 * i + 3], dpa[5 * i + 4]);
		t3 = clock();
		
		//Solve features
		lcrp_solveFeatures( W, V, ea, eb, dpa, dpb );//break
		//for (int i = 0; i < n; i++)
		//	printf("%f %f %f\n", dpb[i * 3], dpb[i * 3 + 1], dpb[i * 3 + 2]);
		
		//update
		for(i=0, dp_L2=0.0; i<nvars; ++i)//nvars��δ֪��������
		{
			//printf("%f ", dp[i]);
			p[i]=p[i] + (tmp=dp[i]);//dpaָ��dp�е���θ�������dpbָ��dp�е���ά���������pָ��δ֪����
			dp_L2+=tmp*tmp;//�洢��������ƽ����
		}

		//dp_L2=sqrt(dp_L2);
		if (dp_L2<1E-8*1E-8)//break
		{	stop = 1;	goto iterstop;	}	//�����ı仯�Ѿ���С						

		lcrp_updateKR( m_R,m_KR, m_KdA, m_KdB, m_KdG, m_K, p );
		t4 = clock();
		
		lcrp_cost(p, hx, m_archor ); 
		pdp_eL2=nrmL2xmy(e, x, hx, nobs); 				
		error = pdp_eL2/nvis;//break ����mse
		t5 = clock();
				
		if ( m_bRobustKernel )
		{
			pdp_eL2 = 0;
			delt2 = m_delt*m_delt;
			if ( m_nRobustType==1)						//Cauchy Kernel Function
			{
				for ( i = 0; i < m_n2Dprojs; i++ )
				{					
					x2 = e[i*2]*e[i*2]+e[i*2+1]*e[i*2+1];
					x2 = delt2*log( x2/delt2 + 1 );
					pdp_eL2 += x2;
				}
			}
			else										//Huber Kernel Function
			{
				for ( i = 0; i < m_n2Dprojs; i++ )
				{					
					x2 = e[i*2]*e[i*2]+e[i*2+1]*e[i*2+1];

					if (x2 <= delt2)  // inlier
						x2 = x2;
					else  // outlier
						x2 = 2*sqrt(x2)*m_delt - delt2;

					pdp_eL2 += x2;
				}
			}
			error = pdp_eL2/nvis;
		}
		
		if (abs(error-lasterror) < 1E-9)//change by zuo
		{	stop = 2;	goto iterstop;		}//���ϴ�mse��ȣ��Ƿ���������
		
		tPerTime = (t5-t0)*0.001;//0.001�����صļ�ʱֵת������
		tTotalTime += tPerTime;//���ε�������

		printf( "Iteration=%d  MSE=%0.8lf Pertime=%0.2lf TotalTime=%0.2lf\n", itno, pdp_eL2/nvis, tPerTime, tTotalTime );
		//printf("setup: %0.2lf  solve: %0.2lf  update: %0.2lf cost: %0.2lf totoal: %0.2lf\n",
		//(t2-t0)/1000.0,(t3-t2)/1000.0, (t4-t3)/1000.0, (t5-t4)/1000.0,(t5-t0)/1000.0);
		//printf("\n");

		if( m_szReport != NULL )
			fprintf( fpRe, "Iteration %d  Error  %0.8lf\n", itno, pdp_eL2/nvis );

		lasterror = error;		//change by zuo
	}
	if(itno>=m_nMaxIter) stop=3;
	
	//clear memory and print
iterstop:
	cholmod_finish (&m_cS) ;  
	sba_crsm_free(&Uidxij);
	sba_crsm_free(&Sidxij);

	cholmod_free_factor(&m_cholFactorS, &m_cS) ;              
	cholmod_l_free_dense(&m_cholSparseE, &m_cS);
	cholmod_l_free_dense(&m_cholSparseR, &m_cS);
	free(Ap);
	free(Aii);
	cholmod_free_sparse(&m_cholSparseS, &m_cS) ;

	tEnd  = clock();
	tTimeUse = tEnd - tStart;
	
	//save optimal camera pose and feature
	//sba_saveXYZ( m_szCamePose, m_sz3Dpts, p );

	free(S);	free(W);	free(U);	free(V);	
	free(e);	free(eab);	free(E);   	free(dp);	free(hx);	
	free(m_KR); free(m_KdA);free(m_KdB);free(m_KdG);

	printf( "%d parameters, %d observations, Gauss-Newton, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
		nvars, m_n2Dprojs*2, stop, error, initialerror, itno, itno, tTimeUse*0.001 );
	printf( "LCRP reasons is listed as following:\n " );
	printf( "reason 1: relative change of state vector is small\n" );
	printf( "reason 2: relative change of projection error is small\n " );
	printf( "reason 3: maximum iteration\n");

	if( m_szReport != NULL )
	{
		fprintf( fpRe, "%d parameters, %d observations, Gauss-Newton, reason %d, error %lf [initial %g], %d iterations [%d linear equations], time %lf sec.\n",
			nvars, m_n2Dprojs*2, stop, error, initialerror, itno, itno, tTimeUse*0.001 );
		fprintf( fpRe, "LCRP reasons is listed as following:\n" );
		fprintf( fpRe, "reason 1: relative change of state vector is small\n" );
		fprintf( fpRe, "reason 2: relative change of projection error is small\n" );
		fprintf( fpRe, "reason 3: maximum iteration\n");
		fclose(fpRe);
	}
	
	return true;
}

double CLCRP::lcrp_computeInitialmu(double* U, double* V, sba_crsm& Uidxij, double tau, int nvars )
{
	int i, j;
	int pos, m = m_ncams,n = m_n3Dpts, cnp = 6, pnp = 3, Usz = 36, Vsz = 9;
	double tmp = 0;
	double *ptr1, *ptr2;
	double mu;

	double *diagUV = (double *)malloc(nvars*sizeof(double));

	double *diagU = diagUV; 
	double *diagV = diagUV + m*cnp;

	for(j=0; j<m; ++j)
	{
		pos  = sba_crsm_elmidx(&Uidxij, j, j);
		ptr1 = U + pos*Usz;
		ptr2 = diagU + j*cnp; 
		for(i=0; i<cnp; ++i)
			ptr2[i]=ptr1[i*cnp+i];
	}
	for(i=0; i<n; ++i)
	{
		ptr1=V + i*Vsz; // set ptr1 to point to V_i
		ptr2=diagV + i*pnp; // set ptr2 to point to diagV_i
		for(j=0; j<pnp; ++j)
			ptr2[j]=ptr1[j*pnp+j];
	}	

	/* find max diagonal element */
	for(i=0, tmp=DBL_MIN; i<m*cnp; ++i)
		if(diagUV[i]>tmp) 
			tmp=diagUV[i];
	for(i=m*cnp; i<nvars; ++i) /* tmp is not re-initialized! */
		if(diagUV[i]>tmp) 
			tmp=diagUV[i];

	mu=m_Tau*tmp;
	
	free(diagUV);
	diagU = diagV = NULL;

	return mu;
}

void CLCRP::lcrp_inverseVLM( double* V, double* IV, sba_crsm& Uidxij, double mu )
{
	int i, j;
	int m = m_ncams, n = m_n3Dpts;
	int Usz = 36, Vsz = 9, pnp = 3, cnp = 6;
	double *ptr1, *ptr2;
	Matrix3d MatInv;
	
	//IV save inverse V matrix, V must unchange for the next step
	memcpy( IV, V, n*Vsz*sizeof(double) );	
	for(i=0; i<n; ++i)
	{
		ptr1=V + i*Vsz; 
		ptr2=IV+ i*Vsz;

		for(j=0; j<pnp; ++j)
			ptr2[j*pnp+j] += mu;

		Eigen::Matrix3d matV(ptr2);
		MatInv = matV.inverse();
		ptr2[0] = MatInv(0,0);
		ptr2[4] = MatInv(1,1);
		ptr2[8] = MatInv(2,2);
		ptr2[1] = ptr2[3] = MatInv(0,1);
		ptr2[2] = ptr2[6] = MatInv(0,2);
		ptr2[5] = ptr2[7] = MatInv(1,2);
	} 
}

void CLCRP::lcrp_inverseVGN( double* U, double* V, sba_crsm& Uidxij )
{
	int i;
	int m = m_ncams, n = m_n3Dpts;
	int Usz = 25, Vsz = 9, pnp = 3, cnp = 5;
	double *ptr1;
	Matrix3d matV, MatInv;
	
	//compute V inverse matrix using Eigen that has better performance than Lapack
	for(i=0; i<n; ++i)
	{
		ptr1=V + i*Vsz; // set ptr1 to point to V_i
		matV << ptr1[0],ptr1[1],ptr1[2],ptr1[3],ptr1[4],ptr1[5],ptr1[6],ptr1[7],ptr1[8];		
		MatInv = matV.inverse();
		ptr1[0] = MatInv(0,0);
		ptr1[4] = MatInv(1,1);
		ptr1[8] = MatInv(2,2);
		ptr1[1] = ptr1[3] = MatInv(0,1);
		ptr1[2] = ptr1[6] = MatInv(0,2);
		ptr1[5] = ptr1[7] = MatInv(1,2);
	}	
}

int CLCRP::lcrp_ConstructSmask( sba_crsm& Sidxij, sba_crsm& Uidxij )//�ֱ���S�����U�����ϡ������
{
	int i, j, k, ii, jj;
	int nuis, m = m_ncams;//nuis��U����ķ���Ԫ�ظ���
	//compute total smask
	for ( i = 0; i < m; i++ ) for ( j = 0; j < m; j++ )
	{
		//��δ���δִ��
		if ( m_umask[i*m+j] == 1 && m_smask[i*m+j] == 0 )//m_umask[i*m+j] == 1 ��ʾ U �����ڸ�λ���Ƿ���ģ���� S ��û��ǣ��Ͱ�������ϡ�
		{
			m_smask[i*m+j] = 1;//S �����ϡ������ (Sidxij �Ľṹ),m_smask[i*m+j] ��һ�� m��m ά�Ķ�ά���飨�� 1D �洢������ʾ S ������ (i,j) λ���Ƿ��з���Ԫ�ء�
			m_nS += 1;
		}
	}	
	//���� S ����Ĵ洢�ռ�,Sidxij��S ����Ľṹ��,m, mΪ�����С
	sba_crsm_alloc(&Sidxij, m, m, m_nS);//m_nS��S����ķ���Ԫ�ظ���
	for(i=k=0; i<m; ++i)
	{
		Sidxij.rowptr[i]=k;// ��¼ S ����� i �е���ʼ����
		ii=i*m;
		for(j=0; j<m; ++j)
			if(m_smask[ii+j])// �����λ���Ƿ���Ԫ��
			{
				Sidxij.val[k]=k;// ������ֵ������ k
				Sidxij.colidx[k++]=j; // ��¼�÷���Ԫ�ص�������
			}
	}
	Sidxij.rowptr[m]=m_nS;// ���һ��ָ��ĩβ

	for(i=nuis=0, jj=m*m; i<jj; ++i)
		nuis+=(m_umask[i]!=0);//U �����ϡ������ (Uidxij �Ľṹ)

	sba_crsm_alloc(&Uidxij, m, m, nuis);
	for(i=k=0; i<m; ++i)
	{
		Uidxij.rowptr[i]=k;
		ii=i*m;
		for(j=0; j<m; ++j)
			if(m_umask[ii+j])
			{
				Uidxij.val[k]=k;
				Uidxij.colidx[k++]=j;
			}
	}
	Uidxij.rowptr[m]=nuis;

	//����
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			printf("%d ", m_umask[i * m + j]);
		}
		printf("\n");
	}
	printf("\n\n\n");
	for (i = 0; i <= m; ++i)
		printf("%d ", Uidxij.rowptr[i]);
	printf("\n\n\n");
	for (i = 0; i < nuis; ++i)
		printf("%d ", Uidxij.val[i]);
	printf("\n\n\n");
	for (i = 0; i < nuis; ++i)
		printf("%d ", Uidxij.colidx[i]);
	printf("\n\n\n");

	//����
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			printf("%d ", m_smask[i * m + j]);
		}
		printf("\n");
	}
	printf("\n\n\n");
	for (i = 0; i <= m; ++i)
		printf("%d ", Sidxij.rowptr[i]);
	printf("\n\n\n");
	for (i = 0; i < m_nS; ++i)
		printf("%d ", Sidxij.val[i]);
	printf("\n\n\n");
	for (i = 0; i < m_nS; ++i)
		printf("%d ", Sidxij.colidx[i]);
	printf("\n\n\n");


	return nuis;

}

void CLCRP::lcrp_constructAuxCSSLM( int *Ap, int *Aii )
{
	int* Cp = Ap;
	int* Ci = Aii;
	int ii, jj;
	int m = m_ncams, nZ = 0;
	for ( ii = 0; ii < m; ii++ ) 
	{
		*Cp = nZ;
		for( jj=0; jj<=ii; jj++ )
		{
			if (m_smask[jj*m+ii]==1)
			{
				*Ci++ = jj;
				nZ++;
			}
		}
		Cp++;
	}
	*Cp=nZ;
}
//����ϡ�����������ṹ����1��1�п�ʼ���洢ÿ�е���ʼ����,�洢����Ԫ�ص�������
void CLCRP::lcrp_constructAuxCSSGN( int *Ap, int *Aii )
{
	//��0��0�п�ʼ����Ϊ�̶���һ֡
	int* Cp = Ap;//�洢ÿ�е���ʼ����
	int* Ci = Aii;//�洢����Ԫ�ص�������
	int ii, jj;
	int m = m_ncams, nZ = 0;
	//printf("\n\n\n");
	for ( ii = 1; ii < m; ii++ ) //��
	{
		*Cp = nZ;
		for( jj=1; jj<=ii; jj++ )//��
		{
			if (m_smask[jj*m+ii]==1)//���б���S���󣨶Գƾ��󣩵�������
			{
				*Ci++ = jj-1;//���ܴ�0�п�ʼ
				//printf("%d ", jj-1);
				nZ++;
			}
		}
		Cp++;
	}
	*Cp=nZ;
	
	//����
	Cp = Ap;//��ָ���Ƶ��ļ�ͷ
	Ci = Aii;
	printf("\n\n\n");
	for (int i = 0; i < m; ++i)
		printf("%d ", *Cp++);
	printf("\n\n\n");
	for (int i = 0; i < m_nS; ++i)
		printf("%d ", *Ci++);
	printf("\n\n\n");
}

void CLCRP::lcrp_constructSLM( double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij, double mu )
{
	int i, j, ii, jj, k, l;
	int m = m_ncams, cnp = 6, pnp = 3, Usz = 36, ebsz=3;
	int pos, pos1, numfea;
	int nF1, nP1, nP2;
	double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5, *ptrS, *ptrE;
	double WV[6*3], sum;

	//Copy U matrix to S matrix 
	pos = 0;
	for ( i = 0; i < m; i++ ) for( j = 0; j < m; j++ )
	{
		if ( m_umask[i*m+j] == 1)// save upper triangle for diagonal element S
		{
			pos1 = sba_crsm_elmidx( &Sidxij, i, j);
			ptr2 = S + pos1*36;
			if ( i == j )
			{					
				ptr1 = U + pos * Usz;
				for(ii=0; ii<cnp; ++ii, ptr2+=6)
				{
					ptr2[ii]= ptr1[ii*cnp+ii] + mu;
					for(jj=ii+1; jj<cnp; ++jj)
						ptr2[jj]= ptr1[ii*cnp+jj];
				}
				pos++;
			}
			else
			{					
				ptr1 = U + pos * Usz;
				for(ii=0; ii<cnp; ++ii, ptr2+=6)
					for(jj=0; jj<cnp; ++jj)
						ptr2[jj]= ptr1[ii*cnp+jj];
				pos++;
			}
		}
	}

	for ( i = 0; i < m*cnp; i++ )
		E[i] = ea[i];

	//Create integrated S matrix, S = U - W(V^-1)W^T
	pos = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];
		for ( j = 0; j < numfea; j++ )
		{
			nF1 = m_feature[pos];
			nP1 = m_photo[pos];
			memset( WV, 0, sizeof(double)*cnp*3 );

			ptr1 = W + pos*cnp*3;	
			ptr2 = V + nF1*3*3;	
			ptrE = E + nP1*cnp;


			//WV
			for(ii=0; ii<cnp; ++ii)
			{
				ptr3=ptr1+ii*pnp;
				for(jj=0; jj<pnp; ++jj)
				{
					for(k=0, sum=0.0; k<=jj; ++k)
						sum+=ptr3[k]*ptr2[jj*pnp+k]; 
					for( ; k<pnp; ++k)
						sum+=ptr3[k]*ptr2[k*pnp+jj]; 
					for(k=0, sum= 0.0; k<pnp; k++ )
						sum+=ptr3[k]*ptr2[jj*pnp+k];
					WV[ii*pnp+jj]=sum;
				}
			}

			for ( k = j; k < numfea; k++ )
			{
				nP2 = m_photo[pos+(k-j)];

				//W(V^-1)W^T
				ptr3 = W + (pos+(k-j))*cnp*3;
				//ptrS = S + (nP1*m*36) + nP2*cnp;

				if ( nP1 == nP2 )
				{
					pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
					ptrS = S + pos1*36;
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;									
						for(jj=ii; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; 
						}
					}
				}else 
				{
					pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
					ptrS = S + pos1*36;
					for(ii=0; ii<cnp; ++ii,ptrS+=6)
					{
						ptr5=WV+ii*pnp;									
						for(jj=0; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; 
						}
					}
				}
			}
			//-W^tb
			ptr5 = eb + nF1*ebsz;
			for(ii=0; ii<cnp; ++ii)
			{
				ptr4=WV+ii*pnp;
				for(jj=0, sum=0.0; jj<pnp; ++jj)
					sum+=ptr4[jj]*ptr5[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];
				ptrE[ii]-= sum;
			}
			pos++;					
		}
	}	
}

void CLCRP::lcrp_constructSGN( double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij )
{
	int i, j, ii, jj, k, l;
	int m = m_ncams, cnp = 5, pnp = 3, Usz = 25, ebsz=3;
	int pos, pos1, numfea;
	int nF1, nP1, nP2;
	double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5, *ptrS, *ptrE;
	double WV[5*3], sum;

	//Copy U matrix to S matrix 
	pos = 0;
	for ( i = 0; i < m; i++ ) for( j = 0; j < m; j++ )
	{
		if ( m_umask[i*m+j] == 1)// save upper triangle for diagonal element S
		{
			pos1 = sba_crsm_elmidx( &Sidxij, i, j);
			//printf("%d ", pos1);
			ptr2 = S + pos1*25;
			if ( i == j )
			{					
				ptr1 = U + pos * Usz;
				for(ii=0; ii<cnp; ++ii, ptr2+=5)
				{
					ptr2[ii]= ptr1[ii*cnp+ii];
					for(jj=ii+1; jj<cnp; ++jj)
						ptr2[jj]= ptr1[ii*cnp+jj];
				}
				pos++;
			}
			else
			{					
				ptr1 = U + pos * Usz;
				for(ii=0; ii<cnp; ++ii, ptr2+=5)
					for(jj=0; jj<cnp; ++jj)
						ptr2[jj]= ptr1[ii*cnp+jj];
				pos++;
			}
		}
	}

	for ( i = 0; i < m*cnp; i++ )
		E[i] = ea[i];//�ⷽλԪ��

	//Create integrated S matrix, S = U - W(V^-1)W^T
	pos = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];
		for ( j = 0; j < numfea; j++ )
		{
			nF1 = m_feature[pos];
			nP1 = m_photo[pos];
			memset( WV, 0, sizeof(double)*cnp*3 );

			ptr1 = W + pos*cnp*3;	
			ptr2 = V + nF1*3*3;	
			ptrE = E + nP1*cnp;

			//W(V^-1)��V�Ѿ������V^-1
			for(ii=0; ii<cnp; ++ii)
			{
				ptr3=ptr1+ii*pnp;//W�����ii��
				for(jj=0; jj<pnp; ++jj)
				{
					for(k=0, sum=0.0; k<=jj; ++k)
						sum+=ptr3[k]*ptr2[jj*pnp+k]; 
					for( ; k<pnp; ++k)
						sum+=ptr3[k]*ptr2[k*pnp+jj]; 
					for(k=0, sum= 0.0; k<pnp; k++ )//(V^-1)�ǶԳƾ���
						sum+=ptr3[k]*ptr2[jj*pnp+k];//�˴���sum�������sum���
					WV[ii*pnp+jj]=sum;
				}
			}

			for ( k = j; k < numfea; k++ )
			{
				nP2 = m_photo[pos+(k-j)];

				//W(V^-1)W^T
				ptr3 = W + (pos+(k-j))*cnp*3;
				pos1 = sba_crsm_elmidx( &Sidxij, nP1, nP2 );
				ptrS = S + pos1*25;
				
				if ( nP1 == nP2 )
				{
					for(ii=0; ii<cnp; ++ii,ptrS+=5)
					{
						ptr5=WV+ii*pnp;	//WV�����ii��								
						for(jj=ii; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;//W�����jj��

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; //-W(V^-1)W^T

							ptrS[jj]-=sum; //S-W(V^-1)W^T
						}
					}
				}else 
				{
					for(ii=0; ii<cnp; ++ii,ptrS+=5)
					{
						ptr5=WV+ii*pnp;									
						for(jj=0; jj<cnp; ++jj)
						{
							ptr4=ptr3+jj*pnp;

							for(l=0, sum=0.0; l<pnp; ++l)
								sum+=ptr5[l]*ptr4[l]; 

							ptrS[jj]-=sum; //�������
						}
					}
				}
			}
			//-W^tb �����Ҳ�
			ptr5 = eb + nF1*ebsz;//��ά��
			for(ii=0; ii<cnp; ++ii)
			{
				ptr4=WV+ii*pnp;//W(V^-1)��ii��
				for(jj=0, sum=0.0; jj<pnp; ++jj)
					sum+=ptr4[jj]*ptr5[jj]; //ptr2[ii*pnp+jj]*ptr3[jj];   W(V^-1)eb
				ptrE[ii]-= sum; //ea-W(V^-1)eb
			}
			pos++;					
		}
	}	
}

void CLCRP::lcrp_solveFeatures( double *W, double *IV, double *ea, double *eb, double *dpa, double *dpb)
{
	int i, j, ii, jj, pos, numfea;
	int nP1, cnp = 5, pnp = 3;
	double *ptr1, *ptr2, *ptr3, *ptr4, *ptr5;
	double sum, eb2[5];
	pos = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		ptr1 = eb  + i*3;//ָ����ά���������
		ptr2 = IV   + i*3*3;//ָ��V����
		ptr5 = dpb + i*3;//ָ����ά��δ֪����
		memset( eb2, 0, sizeof(double)*cnp );
		numfea = m_archor[i*3];

		for ( j = 0; j < numfea; j++ )
		{
			nP1  = m_photo[pos];//��ͼ���
			ptr3 = W + pos*cnp*3;//ָ��W����
			ptr4 = dpa + nP1*cnp;//ָ�����δ֪����
			//Wta
			for(ii=0; ii<pnp; ++ii)
			{
				for(jj=0, sum=0; jj < cnp; ++jj )
					sum += ptr3[jj*3+ii]*ptr4[jj];
				eb2[ii] += sum;
			}
			pos++;
		}

		//V*(eb-Wta)
		for(ii=0; ii<pnp; ++ii)
		{
			for(jj=0, sum = 0; jj < pnp; jj++ )
				sum += ptr2[ii*3+jj]*(ptr1[jj]-eb2[jj]);
			ptr5[ii] = sum;
		}
	}
}

void CLCRP::lcrp_constructCSSLM( int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init)
{
	int ii, jj, jjj, k;
	int pos1, m = m_ncams;
	//Copy S matrix and E matrix to specific format structure for Cholmod 
	double *ptr5;
	int nZ = 0;

	Sx = (double*)m_cholSparseS->x;
	if ( !init)
	{
		for ( ii = 0; ii < m; ii++ )  //colum
		{
			for ( k = 0; k < 6; k++ )
			{
				*Sp = nZ;
				for ( jj = 0; jj <= ii; jj++ )	//row
				{
					if (m_smask[jj*m+ii]==1)
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*36;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
							{
								*Si++ = jj*6 + jjj;
								*Sx++ = ptr5[jjj*6+k];
								nZ++;
							}
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
							{
								*Si++ = jj*6 + jjj;
								*Sx++ = ptr5[jjj*6+k];
								nZ++;
							}
						}
					}
				}
				Sp++;
			}
		}
		*Sp=nZ;
	}
	else
	{
		for ( ii = 0; ii < m; ii++ )  //colum
		{
			for ( k = 0; k < 6; k++ )
			{
				for ( jj = 0; jj <= ii; jj++ )	//row
				{
					if (m_smask[jj*m+ii]==1)
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*36;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
								*Sx++ = ptr5[jjj*6+k];
						}
						else
						{
							for ( jjj = 0; jjj < 6; jjj++)
								*Sx++ = ptr5[jjj*6+k];
						}
					}
				}
			}
		}
	}
}

void CLCRP::lcrp_constructCSSGN( int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init, int nft)
{
	int ii, jj, jjj, k;
	int pos1, m = m_ncams;
	//Copy S matrix and E matrix to specific format structure for Cholmod 
	double *ptr5;
	int nZ = 0;
	Sx = (double*)m_cholSparseS->x;//�洢ϡ����������з���Ԫ�ص�ʵ��ֵ
	//printf("\n\n\n\n");
	if ( !init)
	{
		for ( ii = 1; ii < m; ii++ )  //column����0����ͼ�ǲο�֡��ͨ�����̶�������Ϊ�Ż�����
		{
			for ( k = 0; k < 5; k++ )
			{
				*Sp = nZ;//��ʾϡ�������ÿһ�е���ʼλ�ã�����Ԫ�ص�������,0,1,3,6,
				//printf("%d ", *Sp);
				if ((ii*5+k)==(7+nft))//k=2����Xc��3����Yc��4����Zc��nft=0����̶�Xc��1����̶�Yc��2����̶�Zc
					continue;//�̶���2֡��Xc��Ϊ�߶�Լ��

				for ( jj = 1; jj <= ii; jj++ )	//row
				{
					if ((m_smask[jj*m+ii]==1))
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						//printf("%d ", pos1);
						ptr5 = S + pos1*25;//ָ��S�����Ӧ��jj,ii�����Ӿ�������
						
						if( ii == jj )//�Խ��߿�
						{
							for ( jjj = 0; jjj <= k; jjj++)
							{
								if ( (jj*5+jjj) != (7+nft))//����jj=1��jjj=3�������2�������Xc�̶����⣬
								{
									if (jj * 5 + jjj < 7 + nft)//��2���������̬��
										*Si++ = jj * 5 + jjj - 5;
									else//��2�������Yc��Zc
										*Si++ = jj * 5 + jjj - 6;
									
									*Sx++ = ptr5[jjj*5+k];
									//printf("%d ", ptr5[jjj * 6 + k]);
									nZ++;
								}	
							}
						}
						else
						{
							for ( jjj = 0; jjj < 5; jjj++)
							{
								if ((jj*5+jjj) != (7+nft) )
								{
									if (jj * 5 + jjj < 7 + nft)
										*Si++ = jj * 5 + jjj - 5;
									else
										*Si++ = jj * 5 + jjj - 6;
																			
									*Sx++ = ptr5[jjj*5+k];
									nZ++;
								}
							}
						}
					}
				}
				Sp++;
			}
		}
		*Sp=nZ;
	}
	else
	{
		for ( ii = 1; ii < m; ii++ )  //column
		{
			for ( k = 0; k < 5; k++ )
			{
				if ((ii*5+k)==(7+nft))
					continue;

				for ( jj = 1; jj <= ii; jj++ )	//row
				{
					if ((m_smask[jj*m+ii]==1))
					{
						pos1 = sba_crsm_elmidx( &Sidxij, jj, ii );
						ptr5 = S + pos1*25;

						if( ii == jj )
						{
							for ( jjj = 0; jjj <= k; jjj++)
							{
								if ( (jj*5+jjj) != (7+nft))
									*Sx++ = ptr5[jjj*5+k];
							}
						}
						else
						{
							for ( jjj = 0; jjj < 5; jjj++)
							{
								if ((jj*5+jjj) != (7+nft) )
									*Sx++ = ptr5[jjj*5+k];
							}
						}
					}
				}
			}
		}
	}
}

//learning this skill from G2O
bool CLCRP::lcrp_solveCholmodLM( int* Ap, int* Aii, bool init, bool ordering)
{
	int i, j;
	VectorXi scalarPermutation, blockPermutation;

	ordering = true;
	if (!init)
	{
		if (!ordering)
			m_cholFactorS = cholmod_analyze(m_cholSparseS, &m_cS); // symbolic factorization
		else
		{
			// get the ordering for the block matrix
			if (blockPermutation.size() == 0)
				blockPermutation.resize(m_ncams);
			if (blockPermutation.size() < m_ncams) // double space if resizing
				blockPermutation.resize(2*m_ncams);

			// prepare AMD call via CHOLMOD
			cholmod_sparse auxCholmodSparse;
			auxCholmodSparse.nzmax = m_nS;
			auxCholmodSparse.nrow = auxCholmodSparse.ncol = m_ncams;
			auxCholmodSparse.p = Ap;
			auxCholmodSparse.i = Aii;
			auxCholmodSparse.nz = 0;
			auxCholmodSparse.x = 0;
			auxCholmodSparse.z = 0;
			auxCholmodSparse.stype = 1;
			auxCholmodSparse.xtype = CHOLMOD_PATTERN;
			auxCholmodSparse.itype = CHOLMOD_INT;
			auxCholmodSparse.dtype = CHOLMOD_DOUBLE;
			auxCholmodSparse.sorted = 1;
			auxCholmodSparse.packed = 1;
			int amdStatus = cholmod_amd(&auxCholmodSparse, NULL, 0, blockPermutation.data(), &m_cS);
			if (! amdStatus) 
				return false;
			

			// blow up the permutation to the scalar matrix
			if (scalarPermutation.size() == 0)
				scalarPermutation.resize(m_cholSparseS->ncol);
			if (scalarPermutation.size() < (int)m_cholSparseS->ncol)
				scalarPermutation.resize(2*m_cholSparseS->ncol);
			size_t scalarIdx = 0;

			for ( i = 0; i < m_ncams; ++i)
			{
				const int &pp = blockPermutation(i);
				int base = (pp==0) ? 0 : pp*6;
				int nCols= 6;

				for ( j = 0; j < nCols; ++j)
					scalarPermutation(scalarIdx++) = base++;

			}
			assert(scalarIdx == m_cholSparseS->ncol);

			// apply the ordering
			m_cS.nmethods = 1 ;
			m_cS.method[0].ordering = CHOLMOD_GIVEN;
			m_cholFactorS = cholmod_analyze_p(m_cholSparseS, scalarPermutation.data(), NULL, 0, &m_cS);
		}
	}

	cholmod_factorize(m_cholSparseS, m_cholFactorS, &m_cS); 
	m_cholSparseR = cholmod_solve (CHOLMOD_A, m_cholFactorS, m_cholSparseE, &m_cS) ;
			
	return true;
}
//��Ч���GN�г��ֵ�ϡ������ϵͳ
//ʹ��Cholmod����о���Cholesky���ֽ�����
//Ap����ָ�룬Aii������������CSC��ʽ
bool CLCRP::lcrp_solveCholmodGN( int* Ap, int* Aii, bool init, bool ordering)
{
	int i, j;
	int m = m_ncams;
	VectorXi scalarPermutation, blockPermutation;

	ordering = true;
	if (!init)//�״γ�ʼ��
	{
		if (!ordering)//����Ҫ�ֶ�����
		{
			m_cS.nmethods = 1;
			m_cS.method[0].ordering = CHOLMOD_AMD; //���򷽷�����ΪApproximately Minimum Degree��������С������
			m_cholFactorS = cholmod_analyze(m_cholSparseS, &m_cS); // symbolic factorization�����ɷֽ�������ڲ����ݽṹ
		}
		else
		{
			// get the ordering for the block matrix
			if (blockPermutation.size() == 0)
				blockPermutation.resize(m_ncams-1);//���ڴ洢ÿ���������˳��

			// prepare AMD call via CHOLMOD
			cholmod_sparse auxCholmodSparse;//����cholmod_sparse���͵ĸ�������auxcholmodSparse����ϡ�����Ԫ���ݣ���Ap��Aii��ӳ�䵽Cholmod��ʽ
			auxCholmodSparse.nzmax = m_nS;
			auxCholmodSparse.nrow = auxCholmodSparse.ncol = m-1;//������1֡
			auxCholmodSparse.p = Ap;
			auxCholmodSparse.i = Aii;
			auxCholmodSparse.nz = 0;
			auxCholmodSparse.x = 0;
			auxCholmodSparse.z = 0;
			auxCholmodSparse.stype = 1;
			auxCholmodSparse.xtype = CHOLMOD_PATTERN;
			auxCholmodSparse.itype = CHOLMOD_INT;
			auxCholmodSparse.dtype = CHOLMOD_DOUBLE;
			auxCholmodSparse.sorted = 1;
			auxCholmodSparse.packed = 1;
			//AMD�������ڼ���ϡ�����ֽ��е����
			int amdStatus = cholmod_amd(&auxCholmodSparse, NULL, 0, blockPermutation.data(), &m_cS);//ִ��AMD���򣬲�������洢��blockPermutation
			if (! amdStatus) {
				return false;
			}

			// blow up the permutation to the scalar matrix
			//����������չΪ��������
			if (scalarPermutation.size() == 0)
				scalarPermutation.resize(m_cholSparseS->ncol);
			size_t scalarIdx = 0;

			int a = 0;
			for ( i = 0; i < m_ncams-1; ++i)
			{
				const int &pp = blockPermutation(i);
				int base = (pp==0) ? 0 : pp*5-1;
				int nCols= (pp==0) ? 4 : 5;
				for (j = 0; j < nCols; ++j)
					scalarPermutation(scalarIdx++) = base++;
			}
			assert(scalarIdx == m_cholSparseS->ncol);//m_cholSparseS->ncol=9*5-6=39

			// apply the ordering
			m_cS.nmethods = 1 ;
			m_cS.method[0].ordering = CHOLMOD_GIVEN;//����CHOLMODʹ�ø���������˳��
			m_cholFactorS = cholmod_analyze_p(m_cholSparseS, scalarPermutation.data(), NULL, 0, &m_cS);//����cholmod_analyze_p���з��ŷ���
		}
		init = true;//��ɳ�ʼ��
	}
		
	//Cholmod package for solving sparse linear equation              
	cholmod_factorize(m_cholSparseS, m_cholFactorS, &m_cS); //��ϡ����������ֵ�ֽ�
	m_cholSparseR = cholmod_solve (CHOLMOD_A, m_cholFactorS, m_cholSparseE, &m_cS) ;//������Է�����

	return true;
}
//--------------------------------------------------------

//�����ļ��з�ע���е�����
int CLCRP::findNcameras(FILE *fp)
{
	int lineno, ncams, ch;

	lineno=ncams=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#'){ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if(feof(fp)) break;

		ungetc(ch, fp);

		SKIP_LINE(fp);
		++lineno;
		if(ferror(fp))
		{
			fprintf(stderr, "findNcameras(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		++ncams;
	}
	return ncams;
}

int CLCRP::countNDoubles(FILE *fp)
{
	int lineno, ch, np, i;
	char buf[MAXSTRLEN], *s;
	double dummy;

	lineno=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if(feof(fp)) return 0;

		ungetc(ch, fp);
		++lineno;
		if(!fgets(buf, MAXSTRLEN-1, fp)){ /* read the line found... */
			fprintf(stderr, "countNDoubles(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		/* ...and count the number of doubles it has */
		for(np=i=0, s=buf; 1 ; ++np, s+=i){
			ch=sscanf(s, "%lf%n", &dummy, &i);
			if(ch==0 || ch==EOF) break;
		}

		rewind(fp);
		return np;
	}
	return 0; // should not reach this point
}

int CLCRP::skipNDoubles(FILE *fp, int nvals)
{
	register int i;
	int j;

	for(i=0; i<nvals; ++i)
	{
		j=fscanf(fp, "%*f");
		if(j==EOF) return EOF;

		if(ferror(fp)) return EOF-1;
	}

	return nvals;
}


void CLCRP::readNpointsAndNprojections(FILE *fp, int *n3Dpts, int pnp, int *nprojs, int mnp)
{
	int nfirst, lineno, npts, nframes, ch, n;

	/* #parameters for the first line */
	nfirst=countNDoubles(fp);

	*n3Dpts=*nprojs=lineno=npts=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if(feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n=readNInts(fp, &nframes, 1);
		if(n!=1)
			exit(1);		

		SKIP_LINE(fp);
		*nprojs+=nframes;
		++npts;
	}

	*n3Dpts=npts;
}


void CLCRP::lcrp_readCablibration(FILE* fp, double *K )
{
	int n = fscanf( fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &K[0], &K[3], &K[6], &K[1], &K[4], &K[7], &K[2], &K[5], &K[8] );
	
	if ( n!= 9 )
	{
		fprintf( stderr, "SparseBA error: Format of Calibaration is wrong\n" );
		exit(1);
	}
}


void CLCRP::lcrp_readCameraPose( FILE *fp, double *params )
{
	int n, num, lineno = 0 ;
	double *tofilter;
	double * pPrams = params;

	//the number of element per line is 8, it represents that focal length vary, or it is constant
	num = countNDoubles(fp);
	if ( num==8 )
	{	
		m_bFocal = true;
		m_K = (double*)malloc(m_ncams*2*sizeof(double));
		tofilter=(double *)malloc(8*sizeof(double));
	}
	else
		tofilter=(double *)malloc(6*sizeof(double));

	while(!feof(fp))
	{
		if ( num == 6 )
		{
			n	=	readNDoubles(fp, tofilter, 6);
			pPrams[0] = tofilter[1];	pPrams[1] = tofilter[2];	pPrams[2] = tofilter[3]; 
			pPrams[3] = tofilter[4];	pPrams[4] = tofilter[5];	//pPrams[5] = tofilter[5]; 
			/*------------------------------------������-------------------------------------------*/
			//for (int i = 0;i < 5;i++)
			//	printf("%f ", pPrams[i]);
			//printf("\n");
			/*-------------------------------------������------------------------------------------*/
		}
		if ( num == 8 )
		{
			n	=	readNDoubles(fp, tofilter, 8);
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2]; 
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5]; 

			m_K[lineno*2] = tofilter[6];
			m_K[lineno*2+1] = tofilter[7];
		}		

		if( n==-1 ) 
			break;
		pPrams += 5;
		++lineno;
	}
}


int CLCRP::readNInts(FILE *fp, int *vals, int nvals)
{
	register int i;
	int n, j;

	for(i=n=0; i<nvals; ++i){//��һ�Σ�nvals����1
		j=fscanf(fp, "%d", vals+i);//��һ�Σ���ȡ��һ����Ҳ����1����ά�������2DͶӰ����
		if(j==EOF) return EOF;//��һ�Σ�һ�㲻���ǣ�����j=1����ʾ�ɹ�

		if(j!=1 || ferror(fp)) return EOF-1;//��һ�Σ�һ�㶼�ǣ�j=1����ʾ�ɹ�

		n+=j;//��һ����Ϊn=0�����Ծ���j��Ҳ���ǳɹ���־j=1
	}//֮��Ļ���fscanfÿ��һ����ô�ͻ��ƶ�ָ�룬Ȼ��%d��������ո񡢻��з��������հ��ַ���ֱ��������%d����������
	//����������EOF����j!=0��
	return n;
}


int CLCRP::readNDoubles(FILE* fp, double* vals, int nvals)
{
	register int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i)//nvals��2����ȡxy����
	{
		j = fscanf(fp, "%lf", vals + i);//���￪ʼ������float��Ҳ��������ֵ
		if (j == EOF) return EOF;

		if (j != 1 || ferror(fp)) return EOF - 1;

		n += j;//����֮���2��Ϊ����2��
	}

	return n;//����2
}

//#pragma optimize("", off) // �ر��Ż�
void CLCRP::lcrp_readProjectionAndInitilizeFeature(FILE *fp,
	double *params, double *projs,double *lcpt, char *vmask, int ncams, 
	int *archor,char* umask,int* nphoto, int* nfeature, int* archorSort )
{
	int n;
	int nframes;
	int ptno = 0;// , cur;

	int nproj2D = 0;
	
	//int count = 0;

	int frameno;
	int feastart = 0;

	int nP, nP2;
	
	double* ptr1 = projs;//��¼ÿ��ͶӰ���xy��������ŵģ�2��2��������¼�ʹ洢��
	double* ptr_lc = lcpt;

	int i, j; 
	int  sum, cnp = 5;

	//int nFlag;
	
	int *ptr2;
	//bool bAdjust;	

	double f = m_K[0];
	double cx = m_K[6];
	double cy = m_K[7];
	//printf("%f %f %f\n", f, cx, cy);

	m_smask = (char*)malloc(m_ncams*m_ncams*sizeof(char));
	memset( m_smask, 0, m_ncams*m_ncams*sizeof(char) );

	//int* archorEx = new int[m_n3Dpts*2];

	//bool bM, bN;
	//int max_nframes = -1;
	//read all projection point, initialize three feature angle at the same time
	//��ȡ����ͶӰ�㣬ͬʱ��ʼ������������***�ص�***
	while(!feof(fp))//һ��һ�ж�Match��
	{
		//����
		//if (m_bProvideXYZ)//����������3D������֪
		//{
		//	printf("%d  %f  %f  %f\n", ptno+1, m_XYZ[ptno * 3], m_XYZ[ptno * 3 + 1], m_XYZ[ptno * 3 + 2]);
		//}

		//nFlag = 0;
		n = readNInts( fp, &nframes, 1 );  //��ȡMatch-FeaturePoint�ļ�ÿһ�еĵ�һ�У���ͬ����ĸ���Ϊnframes/3D���2DͶӰ����
		if( n!=1 )
			break;//һ�㶼��1����ʾ�ɹ�
		
		//������Ҫ����ê��
		
		archor[ptno*3] = nframes;//��archor[0]��ʼ����¼3D���2DͶӰ����nframes
		//cur = 0;

		//
		//bM = bN = false;//�������ж���ģ�Ϊʵ�����Ķ�������ص㣺
		for( i=0, sum = 0; i<nframes; ++i )//һ��ͶӰ��һ��ͶӰ��Ķ�ȡ��������nframes
		{
			n = readNInts( fp, &frameno, 1 ); //�ڶ��Σ���ȡ��һ��ͶӰ������-image index

			nphoto[nproj2D] = frameno;//����һ��ͶӰ���id������nphoto��nproj2D��0��ʼ�ģ�
			nfeature[nproj2D] = ptno;//ptnoҲ�Ǵ�0��ʼ�ģ�Ҳ����nfeature[0]=0?
			nproj2D++;

			if(frameno>=ncams)//һ�㲻�ᷢ��
			{
				fprintf(stderr, "LCRP: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles( fp, ptr1, 2 ); //�ڶ��Σ���ȡ��һ��ͶӰ���������� ptr1[0] ptr1[1]
			//��uv����任Ϊ��׶
			double r = sqrt(pow(ptr1[0] - cx, 2) + pow(ptr1[1] - cy, 2));
			*ptr_lc = atan(r / f);
			//printf("%f %f %f\n", ptr1[0], ptr1[1], lcpt[0]);
			ptr1+=2;//һ�μ�¼��������+2
			ptr_lc+=1;
		//	//Ҳ����1��id+2������
			if(n!=3)//һ�㲻�ᷢ��
			{
				fprintf(stderr, "LCRP:reading image projections wrong!\n");
				return;
			}


		//	//�����Σ�ִ����if ( bM && !bN )��ִ������
		//	if ( bM && bN )//����������
		//	{
		//		ptr2 = archorSort+ptno*2;
		//		//bAdjust = pba_initializeOtheArchors_Mindw( //_Mindw
		//		bAdjust = sba_initializeOtheArchors( //Maxdw
		//			projs+feastart*2,
		//			nphoto+feastart,
		//			m_motstruct,
		//			m_K,
		//			m_motstruct + m_ncams*cnp + ptno * 3, 
		//			ptr2,
		//			sum,
		//			i,
		//			ptno );
		//		if ( bAdjust )
		//		{
		//			archor[ptno*3+1] = *(nphoto+feastart+ptr2[0]);
		//			archor[ptno*3+2] = *(nphoto+feastart+ptr2[1]);

		//			archorEx[ptno*2] = ptr2[0];
		//			archorEx[ptno*2+1] = ptr2[1];
		//		}
		//		sum++;
		//	}

		//	//�ڶ��Σ�ִ����if ( !bM )��ִ������
		//	if ( bM && !bN )
		//	{	//bLast����ʱ����
		//		bool bLast = (i == nframes-1);//��һ�α���3��i=0������bLastҲ����false���������һ��
		//		bool bT = sba_initializeAssoArchor( 
		//			projs+feastart*2,	//��¼ÿ��ͶӰ���xy��2��������¼&�洢, feastartһ��ʼ��0�����Ϊ2
		//			nphoto+feastart,	//��¼ÿ�������ţ�feastartһ��ʼ��0�����Ϊ1
		//			m_motstruct,
		//			m_K,
		//			m_motstruct+m_ncams*cnp+ptno*3,
		//			0,
		//			1,
		//			ptno,
		//			bLast );
		//		///*������*/
		//		//if (bT == true)
		//		//	printf("%s\n", "true");
		//		//else
		//		//	printf("%s\n", "false");
		//		//printf("%f %f %f\n", (m_motstruct + m_ncams*cnp + ptno * 3)[0], (m_motstruct + m_ncams*cnp + ptno * 3)[1], (m_motstruct + m_ncams*cnp + ptno * 3)[2]);
		//		///*������*/
		//		if (bT)
		//		{
		//			archorSort[ptno*2+1] = i;//Ϊê�����򣬱���archorSort[0]=2,���3��ͶӰ�����ڵ���ƬΪ��ê�㣬archorSort[1]=0��ʾ��1��ͶӰ��Ϊ��ê��
		//			//archorSort[2]=1��ʾ��2��ͶӰ��Ϊ����ê��
		//			archor[ptno*3+2] = nphoto[count];
		//			sum++;

		//			archorEx[ptno*2+1] = i;

		//			bN = true;
		//		}
		//	}


		//	//��һ�Σ�ǰ���bM = bN = false;�Ļ�������ִ�������
		//	if ( !bM )
		//	{//bLast����ʱ����,���nframes��2�Ļ���ô���������ê����
		//		bool bLast = (i == nframes-2);
		//		bool bT = sba_initializeMainArchor( 
		//			projs+feastart*2,	//��¼ÿ��ͶӰ���xy��2��������¼&�洢, feastartһ��ʼ��0�����Ϊ2
		//			m_motstruct,		//��������������� + ��ά����������
		//			m_K,				//�ڲΣ���֪����
		//			m_motstruct+m_ncams*cnp+ptno*3,//m_ncams�Ƕ��ٸ����/��Ƭ,cnp��6��ʼҲ������εļ����ptno��0��ʼ��*3ΪPBA���������ļ��
		//			nphoto[count],		//nphoto����ȫ��2DͶӰ������
		//			ptno,				//ptno����3D��ĵ�ǰ���
		//			m_KR );				//KR���߷��̱任����
		//		/*�ص�ע�⣡��
		//		* ����*1Ϊ��Ŵ洢
		//		* ����*2Ϊ2D��洢
		//		* ����*3Ϊ3D��洢
		//		* ����������Ϊ��Σ��ڲεȵȴ洢
		//		* ������Ϊ��int/float����doubleΪ��ǰ���
		//		*/

		//		///*������*/
		//		//if (bT == true)
		//		//	printf("%s\n", "true");
		//		//else
		//		//	printf("%s\n", "false");
		//		//printf("%f %f %f\n", (m_motstruct + m_ncams*cnp + ptno * 3)[0], (m_motstruct + m_ncams*cnp + ptno * 3)[1], (m_motstruct + m_ncams*cnp + ptno * 3)[2]);
		//		///*������*/
		//		archorSort[ptno*2] = i;//��i��ͶӰ�㼴��i��anchor
		//		archor[ptno*3+1] = nphoto[count];
		//		sum++;

		//		archorEx[ptno*2] = i;
		//		bM = true;//�ڶ��Σ�ִ������
		//	}	
		//	count++;//�����iһ���ѣ�ÿѭ��һ�ξ͡�+1������ֻ�����һ��������		
		}
		
		//set masks for U and S matrix              umask�洢U����m_smask�洢S����
		for( i = 0; i < nframes; i++ )
		{
			nP = nphoto[feastart+i];                         //��i���۲����ͼ��
			//int nM_ = archor[ptno * 3 + 1];
			//int nA_ = archor[ptno * 3 + 2];
			//int tmp3 = archor[ptno * 3 + 3];
			////umask ������� U �����ϡ��ṹ��ֻ�� ��ê�㡢��ê�㡢��ǰ�۲���� ֮���¼ 1����Ϊ��Щ�������ֱ��Ӱ��õ�Ĺ۲ⷽ�̡�
			//if (nM_<nP)                         
			//	umask[nM_*(ncams)+nP] = 1;//(��ê�㣬��ǰê��)
			//else                                             
			//	umask[nP*(ncams)+nM_] = 1;//(��ǰê�㣬��ê��)

			//
			//if (nA_<nP)                         
			//	umask[nA_*(ncams)+nP] = 1;//����ê�㣬��ǰê�㣩
			//else                                             
			//	umask[nP*(ncams)+nA_] = 1;//����ǰê�㣬��ê�㣩

			umask[nP*ncams+nP] = 1;//����ǰê�㣬��ǰê�㣩

			for ( j = i; j < nframes; j++  )
			{
				nP2 = nphoto[feastart+j];//��j���۲����ͼ��                     

				if ( nP == nP2 )                              //
					m_smask[nP*m_ncams+nP2] = 1;//����i���۲����ͼ�ţ���j���۲����ͼ�ţ�
				else if ( nP < nP2 )
					m_smask[nP*m_ncams+nP2] = 1;//����i���۲����ͼ�ţ���j���۲����ͼ�ţ�
				else
				{
					m_smask[nP2*m_ncams + nP] = 1;//����j���۲����ͼ�ţ���i���۲����ͼ�ţ�δִ�й�
				}
					
			}
		}					
		/*����umask*/
		//for (int k1 = 0; k1 < m_ncams; k1++)
		//{
		//	for (int k2 = 0; k2 < m_ncams; k2++)
		//	{
		//		printf("%d", m_smask[k1*m_ncams + k2]);
		//	}
		//	printf("\n");
		//}
		//max_nframes = nframes > max_nframes ? nframes : max_nframes;
		/*����umask*/
		feastart += nframes;
		ptno++;//3D���������ţ�point NO. OR photo NO.
	}
	/*����umask*/
	//for (int k1 = 0; k1 < m_ncams; k1++)
	//{
	//	for (int k2 = 0; k2 < m_ncams; k2++)
	//	{
	//		printf("%d", umask[k1 * m_ncams + k2]);// -umask[k1 * m_ncams + k2]);
	//	}
	//	printf("\n");
	//}
	//printf("\n\n\n\n\n");
	/*����umask*/
	//count number of non-zero element in S matrix
	m_nS = 0;
	for ( i = 0; i < m_ncams; i++ ) 
	{
		for (j = 0; j < m_ncams; j++)
		{
			//printf("%d", m_smask[i*m_ncams + j]);
			if (m_smask[i*m_ncams + j] == 1)
			{
				m_nS++;
			}
		}
		//printf("\n");
	}
	//printf("\n\n");
	//for (i = 0; i < m_ncams; i++)
	//{
	//	for (j = 0; j < m_ncams; j++)
	//	{
	//		printf("%d", m_smask[i*m_ncams + j]);
	//	}
	//	printf("\n");
	//}
}
void CLCRP::lcrp_readProjectionAndTriangulateFeature(FILE* fp, double* projs, int ncams)
{
	int n;
	int nframes;
	int ptno = 0;
	int frameno;
	double* ptr1 = projs;//��¼ÿ��ͶӰ���xy��������ŵģ�2��2��������¼�ʹ洢��

	int i, j;
	//read all projection point, initialize three feature angle at the same time
	//��ȡ����ͶӰ�㣬ͬʱ��ʼ������������***�ص�***
	while (!feof(fp))//һ��һ�ж�Match��
	{
		n = readNInts(fp, &nframes, 1);  //��ȡMatch-FeaturePoint�ļ�ÿһ�еĵ�һ�У���ͬ����ĸ���Ϊnframes/3D���2DͶӰ����
		if (n != 1)
			break;//һ�㶼��1����ʾ�ɹ�
		
		Eigen::MatrixXd A(2 * nframes, 4);
		//if (nframes > 3)
		//	nframes = 3;


		for (i = 0; i < nframes; ++i)//һ��ͶӰ��һ��ͶӰ��Ķ�ȡ��������nframes
		{
			n = readNInts(fp, &frameno, 1); //�ڶ��Σ���ȡ��һ��ͶӰ������-image index

			if (frameno >= ncams)//һ�㲻�ᷢ��
			{
				fprintf(stderr, "SparseBA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles(fp, ptr1, 2); //�ڶ��Σ���ȡ��һ��ͶӰ���������� ptr1[0] ptr1[1]
			
			//Ҳ����1��id+2������
			if (n != 3)//һ�㲻�ᷢ��
			{
				fprintf(stderr, "SparseBA:reading image projections wrong!\n");
				return;
			}

			const Eigen::Vector2d pt(ptr1[0], ptr1[1]);//��������
			double* ptr = m_P + frameno * 12;
			const Eigen::Matrix<double, 3, 4> Pmat = (Eigen::Matrix<double, 3, 4>() <<
				ptr[0], ptr[1], ptr[2], ptr[3],
				ptr[4], ptr[5], ptr[6], ptr[7],
				ptr[8], ptr[9], ptr[10], ptr[11]).finished();
			A.row(2 * i) = pt(0) * Pmat.row(2) - Pmat.row(0);
			A.row(2 * i + 1) = pt(1) * Pmat.row(2) - Pmat.row(1);
		}

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
		Eigen::Vector4d X = svd.matrixV().col(3);
		X /= X(3);

		(m_motstruct + m_ncams * 6 + ptno * 3)[0] = X[0];
		(m_motstruct + m_ncams * 6 + ptno * 3)[1] = X[1];
		(m_motstruct + m_ncams * 6 + ptno * 3)[2] = X[2];
		ptno++;//3D���������ţ�point NO. OR photo NO.
	}
}
//#pragma optimize("", off) // �ر��Ż�
void CLCRP::lcrp_readAndInitialize(char *camsfname, char *ptsfname, int *ncams,
	int *n3Dpts, int *n2Dprojs,
	double **motstruct, double **imgpts,double **lcpt,
	int **archor, char **vmask,
	char **umask, int **nphoto,
	int** nfeature, int** archorSort)
{
	FILE *fpc, *fpp, *fpXYZ;
	int i, tmp1, tmp2;
	double ptMain[3], ptA[3];
	double dW1, dW2;	

	//calculate number of cameras, 3D points and projection points
	fpc		=	fopen( camsfname, "r" );//���ⷽλԪ���ļ�
	*ncams	=	findNcameras( fpc );//��ȡ��Ƭ������
	m_ncams =	*ncams;

	fpp		=	fopen( ptsfname, "r" );//��ƥ���ͬ�����ļ�
	readNpointsAndNprojections( fpp, n3Dpts, 3, n2Dprojs, 2 );//��ȡ3D������ĸ�����ͶӰ�����

	*motstruct	=	(double *)malloc( (*ncams*5 + *n3Dpts*3)*sizeof(double) );//���ڴ洢���е��ⷽλԪ�� �� �����������3D����
	if(	*motstruct==NULL )
	{
		fprintf(stderr, "LCRP error: Memory allocation for 'motstruct' failed \n");
		exit(1);
	}

	*imgpts	=	(double *)malloc(*n2Dprojs*2*sizeof(double));//���ڴ洢����ͶӰ���2D����
	*lcpt = (double*)malloc(*n2Dprojs * 1 * sizeof(double));
	if (*imgpts == NULL || *lcpt == NULL)
	{
		fprintf(stderr, "LCRP error: Memory allocation for 'imgpts' or 'lcpt' failed\n");
		exit(1);
	}
	//���Ҫ���ļ���ͷ���������ָ���ƶ����ļ�ͷ�����øú���
	rewind(fpc);//ÿ��ȡһ���ַ����ļ��ڲ�λ��ָ������ƶ�һ���ֽڣ���ȡ��ϣ���ָ����ָ���ļ�ĩβ��
	rewind(fpp);

	//allocate indicator of U
	*umask = (char*)malloc(*ncams * *ncams );
	memset(*umask, 0, *ncams * *ncams * sizeof(char));//��1������һ��Ҫ��һ����֪�ģ��Ѿ��������ڴ�ĵ�ַ����3������һ��Ҫʹ��sizeof��������
	//��һ���Ѿ������ַ���ڴ���г�ʼ��������ͨ����ʼ��0�����ַ�'\0'

	//allocate main and associate anchors
	*archor = (int*)malloc(*n3Dpts*3*sizeof(int));//
	memset( *archor, -1, *n3Dpts*3*sizeof(int) ); //��������ê���λ�ó�ʼ��Ϊ-1

	*nphoto		= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*nfeature	= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*archorSort = (int*)malloc(*n3Dpts*3*sizeof(int));

	lcrp_readCameraPose(fpc, *motstruct);
	///*-------------------------------������-----------------------------------------*/
	//for (int i = 0;i < *ncams * 6;i++)
	//{
	//	printf("%f ", motstruct[0][i]);
	//	if ((i + 1) % 6 == 0)
	//		printf("\n");
	//}
	///*-------------------------------������-----------------------------------------*/
	fclose(fpc);//�ѹر�
	
	//Update KR
	m_R = (double*)malloc(m_ncams * 9 * sizeof(double));
	m_KR  = (double*)malloc(m_ncams*9*sizeof(double));//����ڲξ��� �� ��ת����ĳ˻�
	m_KdA = (double*)malloc(m_ncams*9*sizeof(double));//��ת����M_x (��) M_y (��) M_z (��)��һ�׵�
	m_KdB = (double*)malloc(m_ncams*9*sizeof(double));//��ת����M_x (��) M_y (��) M_z (��)��һ�׵�
	m_KdG = (double*)malloc(m_ncams*9*sizeof(double));//��ת����M_x (��) M_y (��) M_z (��)��һ�׵�
	
	if (zu == 1)
	{
		//����P����
		m_P = (double*)malloc(m_ncams * 12 * sizeof(double));
		lcrp_constructP(m_P, m_K, *motstruct);
	}
	else
		lcrp_updateKR(m_R, m_KR, m_KdA, m_KdB, m_KdG, m_K, *motstruct);//��5������Ϊ����ڲΣ���6������Ϊ������


	//test
	//for (int i = 0; i < m_ncams; i++)
	//{
	//	printf("%f %f %f  %f %f %f\n", m_KR[i * 9], m_KR[i * 9 + 1], m_KR[i * 9 + 2], m_P[i * 12], m_P[i * 12 + 1], m_P[i * 12 + 2]);
	//	printf("%f %f %f  %f %f %f\n", m_KR[i * 9 + 3], m_KR[i * 9 + 4], m_KR[i * 9 + 5], m_P[i * 12 + 4], m_P[i * 12 + 5], m_P[i * 12 + 6]);
	//	printf("%f %f %f  %f %f %f\n", m_KR[i * 9 + 6], m_KR[i * 9 + 7], m_KR[i * 9 + 8], m_P[i * 12 + 8], m_P[i * 12 + 9], m_P[i * 12 + 10]);
	//}


	//if XYZ are provided, we can use them as feature initialization.
	if (m_bProvideXYZ)
	{
		fpXYZ = fopen( m_szXYZ, "r");
		m_XYZ = (double*)malloc(m_n3Dpts*3*sizeof(double));
		for (i = 0; i < m_n3Dpts; i++)
		{
			fscanf(fpXYZ, "%lf  %lf  %lf", m_XYZ + i * 3, m_XYZ + i * 3 + 1, m_XYZ + i * 3 + 2);
			(*motstruct)[m_ncams * 5 + i * 3 + 0] = *(m_XYZ + i * 3 + 0);
			(*motstruct)[m_ncams * 5 + i * 3 + 1] = *(m_XYZ + i * 3 + 1);
			(*motstruct)[m_ncams * 5 + i * 3 + 2] = *(m_XYZ + i * 3 + 2);
		}
		fclose(fpXYZ);
	}
	//����
	//for (int i = 0; i < m_ncams; i++)
	//{
	//	printf("%f %f %f %f %f\n", (*motstruct)[i * 5 + 0], (*motstruct)[i * 5 + 1], (*motstruct)[i * 5 + 2], (*motstruct)[i * 5 + 3], (*motstruct)[i * 5 + 4]);
	//}
	//for (int i = 0; i < m_n3Dpts; i++)
	//{
	//	printf("%f %f %f\n", (*motstruct)[m_ncams * 5 + i * 3 + 0], (*motstruct)[m_ncams * 5 + i * 3 + 1], (*motstruct)[m_ncams * 5 + i * 3 + 2]);
	//}
	
	
	if(zu==1)
		lcrp_readProjectionAndTriangulateFeature(fpp, *imgpts, *ncams);
	else
	{
		//ֻ��������umask��smask
		lcrp_readProjectionAndInitilizeFeature(fpp,
			*motstruct + *ncams * 6,
			*imgpts,
			*lcpt,
			*vmask,
			*ncams,
			*archor,//ɾ
			*umask,
			*nphoto,
			*nfeature,//ɾ
			*archorSort);//ɾ
	}

	fclose(fpp);//���������ر�

	
	//int nCount = 0;
	//double pti2k[3];
	//int cur = 0;
	//if ( m_bProvideXYZ )
	//{
	//	for ( i = 0; i < m_n3Dpts; i++ )
	//	{
			//int nM = m_archor[i*3+1];
			//int nN = m_archor[i*3+2];

			////��ê��Position
			//ptMain[0] = *(*motstruct + nM*6 + 3);
			//ptMain[1] = *(*motstruct + nM*6 + 4);
			//ptMain[2] = *(*motstruct + nM*6 + 5);

			////��ê��Position
			//ptA[0] = *(*motstruct + nN*6 + 3);
			//ptA[1] = *(*motstruct + nN*6 + 4);
			//ptA[2] = *(*motstruct + nN*6 + 5);

			////��ê�㵽��ê�������
			//pti2k[0] = ptA[0] - ptMain[0];
			//pti2k[1] = ptA[1] - ptMain[1];
			//pti2k[2] = ptA[2] - ptMain[2];

			//double dispti2k;
			//dispti2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1]  + pti2k[2]*pti2k[2] );//ģ��

			////��ê�㵽��ά�������
			//ptMain[0] = m_XYZ[i*3] - ptMain[0];
			//ptMain[1] = m_XYZ[i*3+1] - ptMain[1];
			//ptMain[2] = m_XYZ[i*3+2] - ptMain[2];	
			////��ê�㵽��ά�������
			//ptA[0] = m_XYZ[i*3] - ptA[0];
			//ptA[1] = m_XYZ[i*3+1] - ptA[1];
			//ptA[2] = m_XYZ[i*3+2] - ptA[2];

			//dW1 = ptMain[0]*ptMain[0] + ptMain[1]*ptMain[1] + ptMain[2]*ptMain[2];
			//dW2 = ptA[0]*ptA[0] + ptA[1]*ptA[1] + ptA[2]*ptA[2];

			//double disDot2;
			//disDot2 = ptMain[0]*pti2k[0] + ptMain[1]*pti2k[1] + ptMain[2]*pti2k[2]; //���
			//double dww = disDot2/(dispti2k * sqrt(dW1));//�������н�

			//double* pKR = m_KR + nM*9;//��ê���KR����
			//double n[2], n2[2], ptXj[3];

			////��ά������ê���uv����
			//ptXj[0] = ptMain[0];	ptXj[1] = ptMain[1];	ptXj[2] = ptMain[2];
			//n[0] = (pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
			//	(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

			//n[1] = (pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
			//	(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

			//pKR = m_KR + nN*9;//��ê���KR����
			//ptXj[0] = ptA[0];	ptXj[1] = ptA[1];	ptXj[2] = ptA[2];
			////��ά���ڸ�ê���uv����
			//n2[0] = (pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
			//	(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

			//n2[1] = (pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
			//	(pKR[6]*ptXj[0] + pKR[7]*ptXj[1] + pKR[8]*ptXj[2]);

			//
			//int id1 = cur + m_archorSort[i*2];//��ê����
			//int id2 = cur + m_archorSort[i*2+1];//��ê����
			//double err1 = (m_imgpts[id1*2]-n[0])*(m_imgpts[id1*2]-n[0])+(m_imgpts[id1*2+1]-n[1])*(m_imgpts[id1*2+1]-n[1]);//��ά������ê�����ͶӰ�в�ƽ��
			//double err2 = (m_imgpts[id2*2]-n2[0])*(m_imgpts[id2*2]-n2[0])+(m_imgpts[id2*2+1]-n2[1])*(m_imgpts[id2*2+1]-n2[1]);//��ά���ڸ�ê�����ͶӰ�в�ƽ��

			//cur += m_archor[i*3];//��i����ά��Ĺ��Ӵ���
			////������ģ���븱����ģ���ı�ֵ>30
			//if ( (sqrt(dW1)/sqrt(dW2)>30) || (err1>err2)&&(m_archor[3*i]==2) )//��ê��в�>��ê��в� �� ���Ӵ���Ϊ2
			//{
			//	nCount++;
			//	//����ê�����
			//	m_archor[i*3+1] = nN;
			//	m_archor[i*3+2] = nM;

			//	tmp1 = m_archorSort[i*2] ;
			//	tmp2 = m_archorSort[i*2+1] ;

			//	m_archorSort[i*2] = tmp2;
			//	m_archorSort[i*2+1] = tmp1;

			//	//��������ڸ�ê��ķ�λ�ǡ��߳̽�
			//	double dDAngle = atan2( ptA[0], ptA[2] );
			//	double dHAngle = atan2( ptA[1], sqrt(ptA[0]*ptA[0]+ ptA[2]*ptA[2]) );

			//	(*motstruct)[m_ncams*6+i*3] = dDAngle;
			//	(*motstruct)[m_ncams*6+i*3+1] = dHAngle;

			//	double dwwDot = ptMain[0]*ptA[0] + ptMain[1]*ptA[1] + ptMain[2]*ptA[2];				
			//}	 
	//	}
	//}

	if(m_bProvideXYZ)
		free(m_XYZ);

	fpXYZ = NULL;
}


void CLCRP::lcrp_readCameraPoseration(char *fname, double* ical )
{
	FILE *fp;
	int  ch=EOF;

	if((fp=fopen(fname, "r"))==NULL)
	{
		fprintf(stderr, "LCRP: Cannot open calbration file %s, exiting\n", fname);
		return;
	}
	//���ж�ȡ1-2-3 || 4-5-6 || 7-8-9
	int num = fscanf( fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &ical[0], &ical[3], &ical[6], &ical[1], &ical[4], &ical[7], &ical[2], &ical[5], &ical[8] );	
	if ( num != 9 )
	{
		fprintf(stderr, "LCRP error: Format of Calibration file is wrong");
		return;
	}

	fclose(fp);
}

void CLCRP::lcrp_updateKR( double *R, double *KR, double *KdA, double *KdB, double *KdG, double *K, double *p )
{
	if ( !m_bFocal )
	{
		int i = 0;
		double *ptAngle;
		double *pR, *pKR, *pKdA, *pKdB, *pKdG;
		double matR[9];//�����ת����
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];

		for ( i = 0; i < m_ncams; i++ )
		{
			ptAngle = p + i*5;//ָ�������ξ���
			/*kappa phi omegaϵͳ*/
			//matR=matRG*matRB*matRA
			//ptAngle=[kappa,phi,omega]
			matR[0] = cos(ptAngle[0]);
			matR[1] = 0;
			matR[2] = -sin(ptAngle[0]);
			matR[3] = sin(ptAngle[1]) * sin(ptAngle[0]);
			matR[4] = cos(ptAngle[1]);
			matR[5] = sin(ptAngle[1])*cos(ptAngle[0]);
			matR[6] = cos(ptAngle[1])*sin(ptAngle[0]);
			matR[7] = -sin(ptAngle[1]);
			matR[8] = cos(ptAngle[1])*cos(ptAngle[0]);
			
			//matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			//matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			//matR[2] = -sin(ptAngle[1]);
			//matR[3] = sin(ptAngle[2])*sin(ptAngle[1])*cos(ptAngle[0])-cos(ptAngle[2])*sin(ptAngle[0]);
			//matR[4] = sin(ptAngle[2])*sin(ptAngle[1])*sin(ptAngle[0])+cos(ptAngle[2])*cos(ptAngle[0]);
			//matR[5] = sin(ptAngle[2])*cos(ptAngle[1]);
			//matR[6] = cos(ptAngle[2])*sin(ptAngle[1])*cos(ptAngle[0]) + sin(ptAngle[2])*sin(ptAngle[0]);
			//matR[7] = cos(ptAngle[2])*sin(ptAngle[1])*sin(ptAngle[0]) - sin(ptAngle[2])*cos(ptAngle[0]);
			//matR[8] = cos(ptAngle[2])*cos(ptAngle[1]);

			//omega��ת����
			matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
			matRG[3] = 0;		matRG[4] = cos(ptAngle[1]);	matRG[5] = sin(ptAngle[1]);
			matRG[6] = 0;		matRG[7] = -sin(ptAngle[1]);	matRG[8] = cos(ptAngle[1]);

			//phi��ת����
			matRB[0] = cos(ptAngle[0]);		matRB[1] = 0;		matRB[2] = -sin(ptAngle[0]);
			matRB[3] = 0;					matRB[4] = 1;		matRB[5] = 0;
			matRB[6] = sin(ptAngle[0]);		matRB[7] = 0;		matRB[8] = cos(ptAngle[0]);

			//kappa��ת����
			//matRA[0] = cos(ptAngle[0]);		matRA[1] = sin(ptAngle[0]);			matRA[2] = 0;
			//matRA[3] = -sin(ptAngle[0]);	matRA[4] = cos(ptAngle[0]);			matRA[5] = 0;
			//matRA[6] = 0;					matRA[7] = 0;						matRA[8] = 1;

			//matRG�������omega��һ�׵�
			matDRG[0] = 0;		matDRG[1] = 0;			matDRG[2] = 0;
			matDRG[3] = 0;		matDRG[4] = -sin(ptAngle[1]);	matDRG[5] = cos(ptAngle[1]);
			matDRG[6] = 0;		matDRG[7] = -cos(ptAngle[1]);	matDRG[8] = -sin(ptAngle[1]);

			//matRB�������phi��һ�׵�
			matDRB[0] = -sin(ptAngle[0]);		matDRB[1] = 0;		matDRB[2] = -cos(ptAngle[0]);
			matDRB[3] = 0;						matDRB[4] = 0;		matDRB[5] = 0;
			matDRB[6] = cos(ptAngle[0]);		matDRB[7] = 0;		matDRB[8] = -sin(ptAngle[0]);

			//matRA�������kappa��һ�׵�
			//matDRA[0] = -sin(ptAngle[0]);		matDRA[1] = cos(ptAngle[0]);		matDRA[2] = 0;
			//matDRA[3] = -cos(ptAngle[0]);		matDRA[4] = -sin(ptAngle[0]);		matDRA[5] = 0;
			//matDRA[6] = 0;						matDRA[7] = 0;						matDRA[8] = 0;

			//pKR=KR*matR
			pKR = KR + i*9;
			pR = R + i * 9;
			pR[0] = matR[0];
			pR[1] = matR[1];
			pR[2] = matR[2];
			pR[3] = matR[3];
			pR[4] = matR[4];
			pR[5] = matR[5];
			pR[6] = matR[6];
			pR[7] = matR[7];
			pR[8] = matR[8];
			//for (int k1 = 0; k1 < 3; k1++)
			//{
			//	for (int k2 = 0; k2 < 3; k2++)
			//	{
			//		printf("%f ", K[k1 * 3 + k2]);
			//	}
			//	printf("\n");
			//}
			pKR[0] = K[0]*matR[0]+K[3]*matR[3]+K[6]*matR[6];
			pKR[1] = K[0]*matR[1]+K[3]*matR[4]+K[6]*matR[7];
			pKR[2] = K[0]*matR[2]+K[3]*matR[5]+K[6]*matR[8];
			pKR[3] = K[1]*matR[0]+K[4]*matR[3]+K[7]*matR[6];
			pKR[4] = K[1]*matR[1]+K[4]*matR[4]+K[7]*matR[7];
			pKR[5] = K[1]*matR[2]+K[4]*matR[5]+K[7]*matR[8];
			pKR[6] = K[2]*matR[0]+K[5]*matR[3]+K[8]*matR[6];
			pKR[7] = K[2]*matR[1]+K[5]*matR[4]+K[8]*matR[7];
			pKR[8] = K[2]*matR[2]+K[5]*matR[5]+K[8]*matR[8];		

			//pR�������omega��һ�׵�
			pKdG = KdG + i*9;
			//tmp1[0] = matDRG[0] + matDRG[3] + matDRG[6];
			//tmp1[1] = matDRG[0] + matDRG[3] + matDRG[6];
			//tmp1[2] = matDRG[0] + matDRG[3] + matDRG[6];
			//tmp1[3] = matDRG[1] + matDRG[4] + matDRG[7];
			//tmp1[4] = matDRG[1] + matDRG[4] + matDRG[7];
			//tmp1[5] = matDRG[1] + matDRG[4] + matDRG[7];
			//tmp1[6] = matDRG[2] + matDRG[5] + matDRG[8];
			//tmp1[7] = matDRG[2] + matDRG[5] + matDRG[8];
			//tmp1[8] = matDRG[2] + matDRG[5] + matDRG[8];

			//tmp2[0] = tmp1[0]*matRB[0]+tmp1[3]*matRB[3]+tmp1[6]*matRB[6];
			//tmp2[1] = tmp1[1]*matRB[0]+tmp1[4]*matRB[3]+tmp1[7]*matRB[6];
			//tmp2[2] = tmp1[2]*matRB[0]+tmp1[5]*matRB[3]+tmp1[8]*matRB[6];
			//tmp2[3] = tmp1[0]*matRB[1]+tmp1[3]*matRB[4]+tmp1[6]*matRB[7];
			//tmp2[4] = tmp1[1]*matRB[1]+tmp1[4]*matRB[4]+tmp1[7]*matRB[7];
			//tmp2[5] = tmp1[2]*matRB[1]+tmp1[5]*matRB[4]+tmp1[8]*matRB[7];
			//tmp2[6] = tmp1[0]*matRB[2]+tmp1[3]*matRB[5]+tmp1[6]*matRB[8];
			//tmp2[7] = tmp1[1]*matRB[2]+tmp1[4]*matRB[5]+tmp1[7]*matRB[8];
			//tmp2[8] = tmp1[2]*matRB[2]+tmp1[5]*matRB[5]+tmp1[8]*matRB[8];

			pKdG[0] = matDRG[0] * matRB[0] + matDRG[1] * matRB[3] + matDRG[2] * matRB[6];
			pKdG[1] = matDRG[0] * matRB[1] + matDRG[1] * matRB[4] + matDRG[2] * matRB[7];
			pKdG[2] = matDRG[0] * matRB[2] + matDRG[1] * matRB[5] + matDRG[2] * matRB[8];
			pKdG[3] = matDRG[3] * matRB[0] + matDRG[4] * matRB[3] + matDRG[5] * matRB[6];
			pKdG[4] = matDRG[3] * matRB[1] + matDRG[4] * matRB[4] + matDRG[5] * matRB[7];
			pKdG[5] = matDRG[3] * matRB[2] + matDRG[4] * matRB[5] + matDRG[5] * matRB[8];
			pKdG[6] = matDRG[6] * matRB[0] + matDRG[7] * matRB[3] + matDRG[8] * matRB[6];
			pKdG[7] = matDRG[6] * matRB[1] + matDRG[7] * matRB[4] + matDRG[8] * matRB[7];
			pKdG[8] = matDRG[6] * matRB[2] + matDRG[7] * matRB[5] + matDRG[8] * matRB[8];

			//pKR�������phi��һ�׵�
			pKdB = KdB + i*9;
			//tmp1[0] = K[0]*matRG[0] + K[3]*matRG[3] + K[6]*matRG[6];
			//tmp1[1] = K[1]*matRG[0] + K[4]*matRG[3] + K[7]*matRG[6];
			//tmp1[2] = K[2]*matRG[0] + K[5]*matRG[3] + K[8]*matRG[6];
			//tmp1[3] = K[0]*matRG[1] + K[3]*matRG[4] + K[6]*matRG[7];
			//tmp1[4] = K[1]*matRG[1] + K[4]*matRG[4] + K[7]*matRG[7];
			//tmp1[5] = K[2]*matRG[1] + K[5]*matRG[4] + K[8]*matRG[7];
			//tmp1[6] = K[0]*matRG[2] + K[3]*matRG[5] + K[6]*matRG[8];
			//tmp1[7] = K[1]*matRG[2] + K[4]*matRG[5] + K[7]*matRG[8];
			//tmp1[8] = K[2]*matRG[2] + K[5]*matRG[5] + K[8]*matRG[8];

			//tmp2[0] = tmp1[0]*matDRB[0]+tmp1[3]*matDRB[3]+tmp1[6]*matDRB[6];
			//tmp2[1] = tmp1[1]*matDRB[0]+tmp1[4]*matDRB[3]+tmp1[7]*matDRB[6];
			//tmp2[2] = tmp1[2]*matDRB[0]+tmp1[5]*matDRB[3]+tmp1[8]*matDRB[6];
			//tmp2[3] = tmp1[0]*matDRB[1]+tmp1[3]*matDRB[4]+tmp1[6]*matDRB[7];
			//tmp2[4] = tmp1[1]*matDRB[1]+tmp1[4]*matDRB[4]+tmp1[7]*matDRB[7];
			//tmp2[5] = tmp1[2]*matDRB[1]+tmp1[5]*matDRB[4]+tmp1[8]*matDRB[7];
			//tmp2[6] = tmp1[0]*matDRB[2]+tmp1[3]*matDRB[5]+tmp1[6]*matDRB[8];
			//tmp2[7] = tmp1[1]*matDRB[2]+tmp1[4]*matDRB[5]+tmp1[7]*matDRB[8];
			//tmp2[8] = tmp1[2]*matDRB[2]+tmp1[5]*matDRB[5]+tmp1[8]*matDRB[8];

			pKdB[0] = matRG[0] * matDRB[0] + matRG[1] * matDRB[3] + matRG[2] * matDRB[6];
			pKdB[1] = matRG[0] * matDRB[1] + matRG[1] * matDRB[4] + matRG[2] * matDRB[7];
			pKdB[2] = matRG[0] * matDRB[2] + matRG[1] * matDRB[5] + matRG[2] * matDRB[8];
			pKdB[3] = matRG[3] * matDRB[0] + matRG[4] * matDRB[3] + matRG[5] * matDRB[6];
			pKdB[4] = matRG[3] * matDRB[1] + matRG[4] * matDRB[4] + matRG[5] * matDRB[7];
			pKdB[5] = matRG[3] * matDRB[2] + matRG[4] * matDRB[5] + matRG[5] * matDRB[8];
			pKdB[6] = matRG[6] * matDRB[0] + matRG[7] * matDRB[3] + matRG[8] * matDRB[6];
			pKdB[7] = matRG[6] * matDRB[1] + matRG[7] * matDRB[4] + matRG[8] * matDRB[7];
			pKdB[8] = matRG[6] * matDRB[2] + matRG[7] * matDRB[5] + matRG[8] * matDRB[8];


			//pKR�������kappa��һ�׵�
			//pKdA = KdA + i*9;
			//tmp2[0] = tmp1[0]*matRB[0]+tmp1[3]*matRB[3]+tmp1[6]*matRB[6];
			//tmp2[1] = tmp1[1]*matRB[0]+tmp1[4]*matRB[3]+tmp1[7]*matRB[6];
			//tmp2[2] = tmp1[2]*matRB[0]+tmp1[5]*matRB[3]+tmp1[8]*matRB[6];
			//tmp2[3] = tmp1[0]*matRB[1]+tmp1[3]*matRB[4]+tmp1[6]*matRB[7];
			//tmp2[4] = tmp1[1]*matRB[1]+tmp1[4]*matRB[4]+tmp1[7]*matRB[7];
			//tmp2[5] = tmp1[2]*matRB[1]+tmp1[5]*matRB[4]+tmp1[8]*matRB[7];
			//tmp2[6] = tmp1[0]*matRB[2]+tmp1[3]*matRB[5]+tmp1[6]*matRB[8];
			//tmp2[7] = tmp1[1]*matRB[2]+tmp1[4]*matRB[5]+tmp1[7]*matRB[8];
			//tmp2[8] = tmp1[2]*matRB[2]+tmp1[5]*matRB[5]+tmp1[8]*matRB[8];

			//pKdA[0] = tmp2[0]*matDRA[0]+tmp2[3]*matDRA[3]+tmp2[6]*matDRA[6];
			//pKdA[3] = tmp2[1]*matDRA[0]+tmp2[4]*matDRA[3]+tmp2[7]*matDRA[6];
			//pKdA[6] = tmp2[2]*matDRA[0]+tmp2[5]*matDRA[3]+tmp2[8]*matDRA[6];
			//pKdA[1] = tmp2[0]*matDRA[1]+tmp2[3]*matDRA[4]+tmp2[6]*matDRA[7];
			//pKdA[4] = tmp2[1]*matDRA[1]+tmp2[4]*matDRA[4]+tmp2[7]*matDRA[7];
			//pKdA[7] = tmp2[2]*matDRA[1]+tmp2[5]*matDRA[4]+tmp2[8]*matDRA[7];
			//pKdA[2] = tmp2[0]*matDRA[2]+tmp2[3]*matDRA[5]+tmp2[6]*matDRA[8];
			//pKdA[5] = tmp2[1]*matDRA[2]+tmp2[4]*matDRA[5]+tmp2[7]*matDRA[8];
			//pKdA[8] = tmp2[2]*matDRA[2]+tmp2[5]*matDRA[5]+tmp2[8]*matDRA[8];

		}
	}
	else
	{
		int i = 0;
		double *ptAngle;
		double *pKR, *pKdA, *pKdB, *pKdG;
		double matR[9];
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];
		double K[9];
		memset( K, 0, 9*sizeof(double));
		K[8] = 1;
		for ( i = 0; i < m_ncams; i++ )
		{
			ptAngle = p + i*6;

			matR[0] = cos(ptAngle[1]);
			matR[1] = 0;
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]);
			matR[4] = cos(ptAngle[2]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]);
			matR[7] = -sin(ptAngle[2]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);
			//matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			//matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			//matR[2] = -sin(ptAngle[1]);
			//matR[3] = sin(ptAngle[2])*sin(ptAngle[1])*cos(ptAngle[0])-cos(ptAngle[2])*sin(ptAngle[0]);
			//matR[4] = sin(ptAngle[2])*sin(ptAngle[1])*sin(ptAngle[0])+cos(ptAngle[2])*cos(ptAngle[0]);
			//matR[5] = sin(ptAngle[2])*cos(ptAngle[1]);
			//matR[6] = cos(ptAngle[2])*sin(ptAngle[1])*cos(ptAngle[0]) + sin(ptAngle[2])*sin(ptAngle[0]);
			//matR[7] = cos(ptAngle[2])*sin(ptAngle[1])*sin(ptAngle[0]) - sin(ptAngle[2])*cos(ptAngle[0]);
			//matR[8] = cos(ptAngle[2])*cos(ptAngle[1]);

			matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
			matRG[3] = 0;		matRG[4] = cos(ptAngle[2]);	matRG[5] = sin(ptAngle[2]);
			matRG[6] = 0;		matRG[7] = -sin(ptAngle[2]);	matRG[8] = cos(ptAngle[2]);

			matRB[0] = cos(ptAngle[1]);		matRB[1] = 0;		matRB[2] = -sin(ptAngle[1]);
			matRB[3] = 0;					matRB[4] = 1;		matRB[5] = 0;
			matRB[6] = sin(ptAngle[1]);		matRB[7] = 0;		matRB[8] = cos(ptAngle[1]);

			matRA[0] = cos(ptAngle[0]);		matRA[1] = sin(ptAngle[0]);			matRA[2] = 0;
			matRA[3] = -sin(ptAngle[0]);	matRA[4] = cos(ptAngle[0]);			matRA[5] = 0;
			matRA[6] = 0;					matRA[7] = 0;						matRA[8] = 1;

			matDRG[0] = 0;		matDRG[1] = 0;			matDRG[2] = 0;
			matDRG[3] = 0;		matDRG[4] = -sin(ptAngle[2]);	matDRG[5] = cos(ptAngle[2]);
			matDRG[6] = 0;		matDRG[7] = -cos(ptAngle[2]);	matDRG[8] = -sin(ptAngle[2]);

			matDRB[0] = -sin(ptAngle[1]);		matDRB[1] = 0;		matDRB[2] = -cos(ptAngle[1]);
			matDRB[3] = 0;						matDRB[4] = 0;		matDRB[5] = 0;
			matDRB[6] = cos(ptAngle[1]);		matDRB[7] = 0;		matDRB[8] = -sin(ptAngle[1]);

			//matDRA[0] = -sin(ptAngle[0]);		matDRA[1] = cos(ptAngle[0]);		matDRA[2] = 0;
			//matDRA[3] = -cos(ptAngle[0]);		matDRA[4] = -sin(ptAngle[0]);		matDRA[5] = 0;
			//matDRA[6] = 0;						matDRA[7] = 0;						matDRA[8] = 0;

			//KR

			K[0] = m_K[i*2];
			K[4] = m_K[i*2+1];

			pKR = KR + i*9;
			pKR[0] = K[0]*matR[0]+K[3]*matR[3]+K[6]*matR[6];
			pKR[1] = K[0]*matR[1]+K[3]*matR[4]+K[6]*matR[7];
			pKR[2] = K[0]*matR[2]+K[3]*matR[5]+K[6]*matR[8];
			pKR[3] = K[1]*matR[0]+K[4]*matR[3]+K[7]*matR[6];
			pKR[4] = K[1]*matR[1]+K[4]*matR[4]+K[7]*matR[7];
			pKR[5] = K[1]*matR[2]+K[4]*matR[5]+K[7]*matR[8];
			pKR[6] = K[2]*matR[0]+K[5]*matR[3]+K[8]*matR[6];
			pKR[7] = K[2]*matR[1]+K[5]*matR[4]+K[8]*matR[7];
			pKR[8] = K[2]*matR[2]+K[5]*matR[5]+K[8]*matR[8];		

			//KdG
			pKdG = KdG + i*9;
			tmp1[0] = K[0]*matDRG[0] + K[3]*matDRG[3] + K[6]*matDRG[6];
			tmp1[1] = K[1]*matDRG[0] + K[4]*matDRG[3] + K[7]*matDRG[6];
			tmp1[2] = K[2]*matDRG[0] + K[5]*matDRG[3] + K[8]*matDRG[6];
			tmp1[3] = K[0]*matDRG[1] + K[3]*matDRG[4] + K[6]*matDRG[7];
			tmp1[4] = K[1]*matDRG[1] + K[4]*matDRG[4] + K[7]*matDRG[7];
			tmp1[5] = K[2]*matDRG[1] + K[5]*matDRG[4] + K[8]*matDRG[7];
			tmp1[6] = K[0]*matDRG[2] + K[3]*matDRG[5] + K[6]*matDRG[8];
			tmp1[7] = K[1]*matDRG[2] + K[4]*matDRG[5] + K[7]*matDRG[8];
			tmp1[8] = K[2]*matDRG[2] + K[5]*matDRG[5] + K[8]*matDRG[8];

			tmp2[0] = tmp1[0]*matRB[0]+tmp1[3]*matRB[3]+tmp1[6]*matRB[6];
			tmp2[1] = tmp1[1]*matRB[0]+tmp1[4]*matRB[3]+tmp1[7]*matRB[6];
			tmp2[2] = tmp1[2]*matRB[0]+tmp1[5]*matRB[3]+tmp1[8]*matRB[6];
			tmp2[3] = tmp1[0]*matRB[1]+tmp1[3]*matRB[4]+tmp1[6]*matRB[7];
			tmp2[4] = tmp1[1]*matRB[1]+tmp1[4]*matRB[4]+tmp1[7]*matRB[7];
			tmp2[5] = tmp1[2]*matRB[1]+tmp1[5]*matRB[4]+tmp1[8]*matRB[7];
			tmp2[6] = tmp1[0]*matRB[2]+tmp1[3]*matRB[5]+tmp1[6]*matRB[8];
			tmp2[7] = tmp1[1]*matRB[2]+tmp1[4]*matRB[5]+tmp1[7]*matRB[8];
			tmp2[8] = tmp1[2]*matRB[2]+tmp1[5]*matRB[5]+tmp1[8]*matRB[8];

			pKdG[0] = tmp2[0]*matRA[0]+tmp2[3]*matRA[3]+tmp2[6]*matRA[6];
			pKdG[3] = tmp2[1]*matRA[0]+tmp2[4]*matRA[3]+tmp2[7]*matRA[6];
			pKdG[6] = tmp2[2]*matRA[0]+tmp2[5]*matRA[3]+tmp2[8]*matRA[6];
			pKdG[1] = tmp2[0]*matRA[1]+tmp2[3]*matRA[4]+tmp2[6]*matRA[7];
			pKdG[4] = tmp2[1]*matRA[1]+tmp2[4]*matRA[4]+tmp2[7]*matRA[7];
			pKdG[7] = tmp2[2]*matRA[1]+tmp2[5]*matRA[4]+tmp2[8]*matRA[7];
			pKdG[2] = tmp2[0]*matRA[2]+tmp2[3]*matRA[5]+tmp2[6]*matRA[8];
			pKdG[5] = tmp2[1]*matRA[2]+tmp2[4]*matRA[5]+tmp2[7]*matRA[8];
			pKdG[8] = tmp2[2]*matRA[2]+tmp2[5]*matRA[5]+tmp2[8]*matRA[8];

			//KdB
			pKdB = KdB + i*9;
			tmp1[0] = K[0]*matRG[0] + K[3]*matRG[3] + K[6]*matRG[6];
			tmp1[1] = K[1]*matRG[0] + K[4]*matRG[3] + K[7]*matRG[6];
			tmp1[2] = K[2]*matRG[0] + K[5]*matRG[3] + K[8]*matRG[6];
			tmp1[3] = K[0]*matRG[1] + K[3]*matRG[4] + K[6]*matRG[7];
			tmp1[4] = K[1]*matRG[1] + K[4]*matRG[4] + K[7]*matRG[7];
			tmp1[5] = K[2]*matRG[1] + K[5]*matRG[4] + K[8]*matRG[7];
			tmp1[6] = K[0]*matRG[2] + K[3]*matRG[5] + K[6]*matRG[8];
			tmp1[7] = K[1]*matRG[2] + K[4]*matRG[5] + K[7]*matRG[8];
			tmp1[8] = K[2]*matRG[2] + K[5]*matRG[5] + K[8]*matRG[8];

			tmp2[0] = tmp1[0]*matDRB[0]+tmp1[3]*matDRB[3]+tmp1[6]*matDRB[6];
			tmp2[1] = tmp1[1]*matDRB[0]+tmp1[4]*matDRB[3]+tmp1[7]*matDRB[6];
			tmp2[2] = tmp1[2]*matDRB[0]+tmp1[5]*matDRB[3]+tmp1[8]*matDRB[6];
			tmp2[3] = tmp1[0]*matDRB[1]+tmp1[3]*matDRB[4]+tmp1[6]*matDRB[7];
			tmp2[4] = tmp1[1]*matDRB[1]+tmp1[4]*matDRB[4]+tmp1[7]*matDRB[7];
			tmp2[5] = tmp1[2]*matDRB[1]+tmp1[5]*matDRB[4]+tmp1[8]*matDRB[7];
			tmp2[6] = tmp1[0]*matDRB[2]+tmp1[3]*matDRB[5]+tmp1[6]*matDRB[8];
			tmp2[7] = tmp1[1]*matDRB[2]+tmp1[4]*matDRB[5]+tmp1[7]*matDRB[8];
			tmp2[8] = tmp1[2]*matDRB[2]+tmp1[5]*matDRB[5]+tmp1[8]*matDRB[8];

			pKdB[0] = tmp2[0]*matRA[0]+tmp2[3]*matRA[3]+tmp2[6]*matRA[6];
			pKdB[3] = tmp2[1]*matRA[0]+tmp2[4]*matRA[3]+tmp2[7]*matRA[6];
			pKdB[6] = tmp2[2]*matRA[0]+tmp2[5]*matRA[3]+tmp2[8]*matRA[6];
			pKdB[1] = tmp2[0]*matRA[1]+tmp2[3]*matRA[4]+tmp2[6]*matRA[7];
			pKdB[4] = tmp2[1]*matRA[1]+tmp2[4]*matRA[4]+tmp2[7]*matRA[7];
			pKdB[7] = tmp2[2]*matRA[1]+tmp2[5]*matRA[4]+tmp2[8]*matRA[7];
			pKdB[2] = tmp2[0]*matRA[2]+tmp2[3]*matRA[5]+tmp2[6]*matRA[8];
			pKdB[5] = tmp2[1]*matRA[2]+tmp2[4]*matRA[5]+tmp2[7]*matRA[8];
			pKdB[8] = tmp2[2]*matRA[2]+tmp2[5]*matRA[5]+tmp2[8]*matRA[8];

			//KdA
			//pKdA = KdA + i*9;
			//tmp2[0] = tmp1[0]*matRB[0]+tmp1[3]*matRB[3]+tmp1[6]*matRB[6];
			//tmp2[1] = tmp1[1]*matRB[0]+tmp1[4]*matRB[3]+tmp1[7]*matRB[6];
			//tmp2[2] = tmp1[2]*matRB[0]+tmp1[5]*matRB[3]+tmp1[8]*matRB[6];
			//tmp2[3] = tmp1[0]*matRB[1]+tmp1[3]*matRB[4]+tmp1[6]*matRB[7];
			//tmp2[4] = tmp1[1]*matRB[1]+tmp1[4]*matRB[4]+tmp1[7]*matRB[7];
			//tmp2[5] = tmp1[2]*matRB[1]+tmp1[5]*matRB[4]+tmp1[8]*matRB[7];
			//tmp2[6] = tmp1[0]*matRB[2]+tmp1[3]*matRB[5]+tmp1[6]*matRB[8];
			//tmp2[7] = tmp1[1]*matRB[2]+tmp1[4]*matRB[5]+tmp1[7]*matRB[8];
			//tmp2[8] = tmp1[2]*matRB[2]+tmp1[5]*matRB[5]+tmp1[8]*matRB[8];

			//pKdA[0] = tmp2[0]*matDRA[0]+tmp2[3]*matDRA[3]+tmp2[6]*matDRA[6];
			//pKdA[3] = tmp2[1]*matDRA[0]+tmp2[4]*matDRA[3]+tmp2[7]*matDRA[6];
			//pKdA[6] = tmp2[2]*matDRA[0]+tmp2[5]*matDRA[3]+tmp2[8]*matDRA[6];
			//pKdA[1] = tmp2[0]*matDRA[1]+tmp2[3]*matDRA[4]+tmp2[6]*matDRA[7];
			//pKdA[4] = tmp2[1]*matDRA[1]+tmp2[4]*matDRA[4]+tmp2[7]*matDRA[7];
			//pKdA[7] = tmp2[2]*matDRA[1]+tmp2[5]*matDRA[4]+tmp2[8]*matDRA[7];
			//pKdA[2] = tmp2[0]*matDRA[2]+tmp2[3]*matDRA[5]+tmp2[6]*matDRA[8];
			//pKdA[5] = tmp2[1]*matDRA[2]+tmp2[4]*matDRA[5]+tmp2[7]*matDRA[8];
			//pKdA[8] = tmp2[2]*matDRA[2]+tmp2[5]*matDRA[5]+tmp2[8]*matDRA[8];

		}
	}	
}

//By Zuo
void CLCRP::lcrp_constructP(double* P, double* K, double* p)
{
	if (!m_bFocal)
	{
		int i = 0;
		double* ptAngle;
		double* pP;
		double matR[9];
		double matT[3];

		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;//ָ�������ξ���
			/*kappa phi omegaϵͳ*/
			//matR=matRG*matRB*matRA
			//ptAngle=[kappa,phi,omega]
			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			matT[0] = -matR[0] * ptAngle[3] - matR[1] * ptAngle[4] - matR[2] * ptAngle[5];
			matT[1] = -matR[3] * ptAngle[3] - matR[4] * ptAngle[4] - matR[5] * ptAngle[5];
			matT[2] = -matR[6] * ptAngle[3] - matR[7] * ptAngle[4] - matR[8] * ptAngle[5];

			//pP=K*[matR matT]
			pP = P + i * 12;
			pP[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pP[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pP[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pP[3] = K[0] * matT[0] + K[3] * matT[1] + K[6] * matT[2];
			pP[4] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pP[5] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pP[6] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pP[7] = K[1] * matT[0] + K[4] * matT[1] + K[7] * matT[2];
			pP[8] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pP[9] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pP[10] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];
			pP[11] = K[2] * matT[0] + K[5] * matT[1] + K[8] * matT[2];

			//test
			//double t1 = pP[0] * (2 - ptAngle[3]) + pP[1] * (2 - ptAngle[4]) + pP[2] * (2 - ptAngle[5]);
			//double t2 = pP[4] * (2 - ptAngle[3]) + pP[5] * (2 - ptAngle[4]) + pP[6] * (2 - ptAngle[5]);
			//double t3 = pP[8] * (2 - ptAngle[3]) + pP[9] * (2 - ptAngle[4]) + pP[10] * (2 - ptAngle[5]);

			//double u1 = t1 / t3;
			//double v1 = t2 / t3;

			//t1 = pP[0] * 2 + pP[1] * 2 + pP[2] * 2 + pP[3];
			//t2 = pP[4] * 2 + pP[5] * 2 + pP[6] * 2 + pP[7];
			//t3 = pP[8] * 2 + pP[9] * 2 + pP[10] * 2 + pP[11];
			//double u2 = t1 / t3;
			//double v2 = t2 / t3;
		}
	}
	else
	{
		int i = 0;
		double* ptAngle;
		double* pP;
		double matR[9],matT[3];
		double K[9];
		memset(K, 0, 9 * sizeof(double));
		K[8] = 1;
		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;

			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			matT[0] = -matR[0] * ptAngle[3] - matR[1] * ptAngle[4] - matR[2] * ptAngle[5];
			matT[1] = -matR[3] * ptAngle[3] - matR[4] * ptAngle[4] - matR[5] * ptAngle[5];
			matT[2] = -matR[6] * ptAngle[3] - matR[7] * ptAngle[4] - matR[8] * ptAngle[5];

			//KR

			K[0] = m_K[i * 2];
			K[4] = m_K[i * 2 + 1];

			//pP=K*[matR matT]
			pP = P + i * 12;
			pP[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pP[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pP[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pP[3] = K[0] * matT[0] + K[3] * matT[1] + K[6] * matT[2];
			pP[4] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pP[5] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pP[6] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pP[7] = K[1] * matT[0] + K[4] * matT[1] + K[7] * matT[2];
			pP[8] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pP[9] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pP[10] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];
			pP[11] = K[2] * matT[0] + K[5] * matT[1] + K[8] * matT[2];
		}
	}
}

#pragma optimize("", off) // �ر��Ż�
void CLCRP::lcrp_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
{
	int nM, nN;	  
	register int i, j, ii, jj, k;
	int cnp, pnp, mnp;
	double *pa, *pb, *ppt;	
	double *ppUpa, *pea, *peb, *pe, *pV, *pW;
	int pos, pos2, nP, nF, numfea, cur;
	double sum;
	register double tmp = 0;

	double pAM[3], pAA[3], pPA[5], pPB[3];

	cnp = 5; pnp = 3; mnp = 1;
	pa=p; pb=p+m*cnp;

	pos2 = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];//��ά��Ĺ���ͼ����
		cur = 0;
		nF = i;
		//printf("%d %d\n", m_archorSort[i * 2], m_archorSort[i * 2 + 1]);
		
		for ( j = 0; j < numfea; j++ )
		{
			
			nP = nphoto[pos2];//����ͼ���
			ppt=pb + nF*pnp;//��ά������  
			pe= e + pos2*1;		//error

			//nM = archor[nF*3+1];//��ê��
			//nN = archor[nF*3+2];//��ê��	

			//paָ���������
			//pptָ����ά�����
			//nP��ǰ��ͼ���
			//nM��ê����
			//nA��ê����
			//pPAָ��uv�Ե�ǰ��ͼ������һ�׵�
			//pPBָ��uv����ά���һ�׵�
			//pAMָ��uv������ê����ͼ��uv������ê����ͼλ�Ʋ�����һ�׵�
			//pAAָ��uv��������ê����ͼ��uv���Ը�ê����ͼλ�Ʋ�����һ�׵�

			lcrp_jacobianEachPts(m_R, m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nM, nN, nP, pAM, pAA, pPA, pPB );
		

			//����
			//double matXj[3], matxyt[3], n1[1], n2[1];
			//double *ppa = pa + nP * 5 + 2;//��ǰ��ͼλ��ָ��
			//matXj[0] = ppt[0] - ppa[0];
			//matXj[1] = ppt[1] - ppa[1];
			//matXj[2] = ppt[2] - ppa[2];
			//double* pR = m_R + nP * 9;
			//matxyt[0] = pR[0] * matXj[0] + pR[1] * matXj[1] + pR[2] * matXj[2];
			//matxyt[1] = pR[3] * matXj[0] + pR[4] * matXj[1] + pR[5] * matXj[2];
			//matxyt[2] = pR[6] * matXj[0] + pR[7] * matXj[1] + pR[8] * matXj[2];
			//n1[0] = atan(sqrt(matxyt[0] * matxyt[0] + matxyt[1] * matxyt[1]) / matxyt[2]);

			//double sigma = 1e-6;
			////double* ppa_new = pa + nP * 5;
			//double* ppa_new = pa + nP * 5;
			//double t[5];
			//for (int i = 0; i < 5; i++)
			//{
			//	ppa_new[i] += sigma;
			//	//for (int hh = 0; hh < 90; hh++)
			//	//{
			//	//	printf("%f %f %f %f %f %f\n", pa[hh * 6], pa[hh * 6 + 1], pa[hh * 6 + 2], pa[hh * 6 + 3], pa[hh * 6 + 4], pa[hh * 6 + 5]);
			//	//}
			//	lcrp_updateKR(m_R, m_KR, m_KdA, m_KdB, m_KdG, m_K, pa);
			//	//printf("%f %f %f\n", ppt[0], ppt[1], ppt[2]);
			//	matXj[0] = ppt[0] - ppa[0];
			//	matXj[1] = ppt[1] - ppa[1];
			//	matXj[2] = ppt[2] - ppa[2];
			//	matxyt[0] = pR[0] * matXj[0] + pR[1] * matXj[1] + pR[2] * matXj[2];
			//	matxyt[1] = pR[3] * matXj[0] + pR[4] * matXj[1] + pR[5] * matXj[2];
			//	matxyt[2] = pR[6] * matXj[0] + pR[7] * matXj[1] + pR[8] * matXj[2];
			//	n2[0] = atan(sqrt(matxyt[0] * matxyt[0] + matxyt[1] * matxyt[1]) / matxyt[2]);
			//	t[i] = (n2[0] - n1[0]) / sigma;
			//	ppa_new[i] -= sigma;
			//}

			//printf("%f %f %f %f %f    %f %f %f %f %f\n", pPA[0], pPA[1], pPA[2], pPA[3], pPA[4], t[0], t[1], t[2], t[3], t[4]);

			//U
			pos = sba_crsm_elmidx(Uidxij, nP, nP);//����U������nP��nP�е�Ԫ����ϡ��ṹ�е�������
			//printf("%d  ", pos);
			ppUpa = U + pos*(5*5);//������䣬�᲻�����ͬһ�鱻�ظ��������⣿ͬһ���ǲ�ͬ�۲�ƫ������ӹ��̣����ᱻ�ظ����

			//printf("\n\n\n\n");
			for ( ii = 0; ii < 5; ii++ ) for( jj = ii; jj < 5; jj++ )	//�����ǶԳƾ���ֻ��������
			{
				for ( k = 0, sum = 0; k < 1; k++ )
					sum += pPA[k*5+ii]*pPA[k*5+jj];//U^T*U
				ppUpa[ii*5+jj] += sum;
				//printf("%d ", ppUpa[ii * 6 + jj]);
			}

			//ea
			pea = ea + nP*5;
			for ( ii = 0; ii < 5; ii++ )
			{
				for( jj = 0, sum = 0; jj < 1; jj++ )
					sum += pPA[jj*5+ii] * pe[jj];
				pea[ii] += sum;
			}
			//V
			pV = V + nF*3*3;
			for ( ii = 0; ii < 3; ii++ ) for ( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 1; k++ )
					sum += pPB[k*3+ii]*pPB[k*3+jj];
				pV[ii*3+jj] += sum;

			}
			//eb
			peb = eb + nF*3;
			for ( ii = 0; ii < 3; ii++ )
			{
				for( k = 0, sum = 0; k < 1; k++ )
					sum += pPB[k*3+ii] * pe[k];
				peb[ii] += sum;
			}
			//W
			pW = W + pos2 * 3*5;
			for ( ii = 0; ii < 5; ii++ )	for( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 1; k++ )
					sum += pPA[k*5+ii] * pPB[k*3+jj];
				pW[ii*3+jj] += sum;
			}
			//if ( nP == nM )
			//{
			//	cur++;
			//	pos2++;
			//	continue;
			//}
			//else
			//if( nP == nN )
			//{
			//	//U		  
			//	pos = sba_crsm_elmidx(Uidxij, nM, nM );		//main archor * main archor
			//	ppUpa = U + pos*(6*6);
			//	for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
			//	{
			//		for ( k = 0, sum = 0; k < 2; k++ )
			//			sum += pAM[k*3+ii] * pAM[k*3+jj];
			//		ppUpa[(ii+3)*6+3+jj] += sum;
			//	}
			//	if( nM < nP )
			//	{
			//		pos = sba_crsm_elmidx(Uidxij, nM, nP );		//main archor * associate archor
			//		ppUpa = U + pos*(6*6);
			//		for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
			//		{
			//			for ( k = 0, sum = 0; k < 2; k++ )
			//				sum += pAM[k*3+ii] * pPA[k*6+jj];
			//			ppUpa[(ii+3)*6+jj] += sum;
			//		}
			//	}
			//	else
			//	{
			//		pos = sba_crsm_elmidx(Uidxij, nP, nM );		
			//		ppUpa = U + pos*(6*6);
			//		for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
			//		{
			//			for ( k = 0, sum = 0; k < 2; k++ )
			//				sum +=  pPA[k*6+ii]*pAM[k*3+jj];
			//			ppUpa[ii*6+jj+3] += sum;
			//		}
			//	}

			//	//ea
			//	pea = ea + nM*6;
			//	for ( ii = 0; ii < 3; ii++ )
			//	{
			//		for( k = 0, sum = 0; k < 2; k++ )
			//			sum += pAM[k*3+ii] * pe[k];//������ͼ������ê�㣩�������ê��λ�Ʋ�����ƫ�� * ������ͼ�����Ԥ��в�
			//		pea[ii+3] += sum;
			//	}
			//	//W
			//	pos = pos2 + (m_archorSort[i*2]-j);//ѭ��β��pos2++����������ط�-j
			//	pW = W + pos*6*3;
			//	for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
			//	{
			//		for ( k = 0, sum = 0; k < 2; k++ )
			//			sum += pAM[k*3+ii] * pPB[k*3+jj];
			//		pW[(ii+3)*3+jj] += sum;
			//	}
			//}
			//else
			//{
			//	//U
			//	pos = sba_crsm_elmidx(Uidxij, nM, nM );		
			//	ppUpa = U + pos*(6*6);
			//	for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
			//	{
			//		for ( k = 0, sum = 0; k < 2; k++ )
			//		{
			//			sum += pAM[k*3+ii] * pAM[k*3+jj];
			//		}
			//		
			//		ppUpa[(ii+3)*6+3+jj] += sum;
			//	}

			//	pos = sba_crsm_elmidx(Uidxij, nN, nN );		
			//	ppUpa = U + pos*(6*6);
			//	for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
			//	{
			//		for ( k = 0, sum = 0; k < 2; k++ )
			//		{
			//			sum += pAA[k*3+ii] * pAA[k*3+jj];
			//		}
			//		ppUpa[(ii+3)*6+3+jj] += sum;
			//	}
			//	if ( nM < nN )
			//	{
			//		pos = sba_crsm_elmidx(Uidxij, nM, nN );		
			//		ppUpa = U + pos*(6*6);
			//		for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
			//		{
			//			for ( k = 0, sum = 0; k < 2; k++ )
			//			{
			//				sum += pAM[k*3+ii] * pAA[k*3+jj];
			//			}
			//			ppUpa[(ii+3)*6+3+jj] += sum;
			//		}
			//	}	
			//	else
			//	{
			//		pos = sba_crsm_elmidx(Uidxij, nN, nM );		
			//		ppUpa = U + pos*(6*6);
			//		for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
			//		{
			//			for ( k = 0, sum = 0; k < 2; k++ )
			//			{
			//				sum += pAA[k*3+ii] * pAM[k*3+jj];
			//			}
			//			ppUpa[(ii+3)*6+3+jj] += sum;
			//		}
			//	}
			//	if ( nM < nP )
			//	{
			//		pos = sba_crsm_elmidx(Uidxij, nM, nP );		
			//		ppUpa = U + pos*(6*6);
			//		for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
			//		{
			//			for ( k = 0, sum = 0; k < 2; k++ )
			//				sum += pAM[k*3+ii] * pPA[k*6+jj];
			//			ppUpa[(ii+3)*6+jj] += sum;
			//		}
			//	}
			//	else
			//	{
			//		pos = sba_crsm_elmidx(Uidxij, nP, nM );		
			//		ppUpa = U + pos*(6*6);
			//		for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
			//		{
			//			for ( k = 0, sum = 0; k < 2; k++ )
			//				sum +=  pPA[k*6+ii]*pAM[k*3+jj];
			//			ppUpa[ii*6+jj+3] += sum;
			//		}
			//	}
			//	
			//	if ( nN < nP )
			//	{
			//		pos = sba_crsm_elmidx(Uidxij, nN, nP );		
			//		ppUpa = U + pos*(6*6);
			//		for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
			//		{
			//			for ( k = 0, sum = 0; k < 2; k++ )
			//				sum += pAA[k*3+ii] * pPA[k*6+jj];
			//			ppUpa[(ii+3)*6+jj] += sum;
			//		}
			//	}	
			//	else
			//	{
			//		pos = sba_crsm_elmidx(Uidxij, nP, nN );		
			//		ppUpa = U + pos*(6*6);
			//		for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
			//		{
			//			for ( k = 0, sum = 0; k < 2; k++ )
			//				sum +=  pPA[k*6+ii]*pAA[k*3+jj];
			//			ppUpa[ii*6+jj+3] += sum;
			//		}
			//	}
			//	//ea
			//	pea = ea + nM*6;
			//	for ( ii = 0; ii < 3; ii++ )
			//	{
			//		for( k = 0, sum = 0; k < 2; k++ )
			//			sum += pAM[k*3+ii] * pe[k];
			//		pea[ii+3] += sum;
			//	}
			//	pea = ea + nN*6;
			//	for ( ii = 0; ii < 3; ii++ )
			//	{
			//		for( k = 0, sum = 0; k < 2; k++ )
			//			sum += pAA[k*3+ii] * pe[k];
			//		pea[ii+3] += sum;
			//	}  
			//	//W
			//	pos = pos2 + (m_archorSort[i*2]-j);
			//	pW = W + pos*6*3;
			//	for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
			//	{
			//		for ( k = 0, sum = 0; k < 2; k++ )
			//			sum += pAM[k*3+ii] * pPB[k*3+jj];
			//		pW[(ii+3)*3+jj] += sum;
			//	}
			//	pos = pos2 + (m_archorSort[i*2+1]-j);
			//	pW = W + pos*6*3;
			//	for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
			//	{
			//		for ( k = 0, sum = 0; k < 2; k++ )
			//			sum += pAA[k*3+ii] * pPB[k*3+jj];
			//		pW[(ii+3)*3+jj] += sum;
			//	}
			//}
			pos2++;
			cur++;
		}
	}
}

void CLCRP::lcrp_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
{
	int nM, nN;	  
	register int i, j, ii, jj, k;
	int cnp, pnp, mnp;
	double *pa, *pb, *ppt;	
	double *ppUpa, *pea, *peb, *pe, *pV, *pW;
	int pos, pos2, nP, nF, numfea, cur;
	double sum, x2;
	register double tmp = 0;

	double pAM[6], pAA[6], pPA[12], pPB[6];

	cnp = 6; pnp = 3; mnp = 2;
	pa=p; pb=p+m*cnp;

	pos2 = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];
		cur = 0;
		nF = i;
		for ( j = 0; j < numfea; j++ )
		{
			nP = nphoto[pos2];
			ppt=pb + nF*pnp;      
			pe= e + pos2*2;		//error

			nM = archor[nF*3+1];
			nN = archor[nF*3+2];

			double delt2 = m_delt*m_delt;
			if( m_nRobustType == 1)
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];
				x2 = 1.0/( x2/delt2+1 );
			}
			else
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];

				if (sqrt(x2) < m_delt) //inlier
					x2 = 1;
				else // outliers 
					x2 = m_delt/sqrt(x2);
			}			

			lcrp_jacobianEachPts(m_R, m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nM, nN, nP, pAM, pAA, pPA, pPB );

			//U
			pos = sba_crsm_elmidx(Uidxij, nP, nP);		
			ppUpa = U + pos*(6*6);

			for ( ii = 0; ii < 6; ii++ ) for( jj = ii; jj < 6; jj++ )	//diag
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii]*pPA[k*6+jj];
				//ppUpa[ii*6+jj] += sum;
				ppUpa[ii*6+jj] += sum*x2;
			}

			//ea
			pea = ea + nP*6;
			for ( ii = 0; ii < 6; ii++ )
			{
				for( jj = 0, sum = 0; jj < 2; jj++ )
					sum += pPA[jj*6+ii] * pe[jj];
				//pea[ii] += sum;
				pea[ii] += sum*x2;
			}
			//V
			pV = V + nF*3*3;
			for ( ii = 0; ii < 3; ii++ ) for ( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii]*pPB[k*3+jj];
				//pV[ii*3+jj] += sum;
				pV[ii*3+jj] += sum*x2;

			}
			//eb
			peb = eb + nF*3;
			for ( ii = 0; ii < 3; ii++ )
			{
				for( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii] * pe[k];
				//peb[ii] += sum;
				peb[ii] += sum*x2;
			}
			//W
			pW = W + pos2 * 3*6;
			for ( ii = 0; ii < 6; ii++ )	for( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii] * pPB[k*3+jj];
				//pW[ii*3+jj] += sum;
				pW[ii*3+jj] += sum*x2;
			}

			if ( nP == nM )
			{
				cur++;
				pos2++;
				continue;
			}
			else
				if( nP == nN )
				{
					//U		  
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		//main archor * main archor
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						//ppUpa[(ii+3)*6+3+jj] += sum;
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		//main archor * associate archor
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							//ppUpa[(ii+3)*6+jj] += sum;
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							//ppUpa[ii*6+jj+3] += sum;
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						//pea[ii+3] += sum;
						pea[ii+3] += sum*x2;
					}
					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						//pW[(ii+3)*3+jj] += sum;
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}
				else
				{
					//U
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						}
						//ppUpa[(ii+3)*6+3+jj] += sum;
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					pos = sba_crsm_elmidx(Uidxij, nN, nN );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAA[k*3+ii] * pAA[k*3+jj];
						}
						//ppUpa[(ii+3)*6+3+jj] += sum;
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if ( nM < nN )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAM[k*3+ii] * pAA[k*3+jj];
							}
							//ppUpa[(ii+3)*6+3+jj] += sum;
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAA[k*3+ii] * pAM[k*3+jj];
							}
							//ppUpa[(ii+3)*6+3+jj] += sum;
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}

					if ( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							//ppUpa[(ii+3)*6+jj] += sum;
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							//ppUpa[ii*6+jj+3] += sum;
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					if ( nN < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAA[k*3+ii] * pPA[k*6+jj];
							//ppUpa[(ii+3)*6+jj] += sum;
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAA[k*3+jj];
							//ppUpa[ii*6+jj+3] += sum;
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}


					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						//pea[ii+3] += sum;
						pea[ii+3] += sum*x2;
					}

					pea = ea + nN*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pe[k];
						//pea[ii+3] += sum;
						pea[ii+3] += sum*x2;
					}  

					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						//pW[(ii+3)*3+jj] += sum;
						pW[(ii+3)*3+jj] += sum*x2;
					}

					pos = pos2 + (m_archorSort[i*2+1]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pPB[k*3+jj];
						//pW[(ii+3)*3+jj] += sum;
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}

				pos2++;
				cur++;

		}
	}
}
/*void CParallaxBA::pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature )
{
	int nM, nN;	  
	register int i, j, ii, jj, k;
	int cnp, pnp, mnp;
	double *pa, *pb, *ppt;	
	double *ppUpa, *pea, *peb, *pe, *pV, *pW;
	int pos, pos2, nP, nF, numfea, cur;
	static double sum, x2, delt2;
	register double tmp = 0;

	double pAM[6], pAA[6], pPA[12], pPB[6];

	cnp = 6; pnp = 3; mnp = 2;
	pa=p; pb=p+m*cnp;
	delt2 = m_delt*m_delt;

	pos2 = 0;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		numfea = m_archor[i*3];
		cur = 0;
		nF = i;
		for ( j = 0; j < numfea; j++ )
		{
			nP = nphoto[pos2];
			ppt=pb + nF*pnp;      
			pe= e + pos2*2;		//error

			nM = archor[nF*3+1];
			nN = archor[nF*3+2];	

			if( m_nRobustType == 1)
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];
				x2 = 1.0/( x2/delt2+1 );
			}
			else
			{
				x2 = pe[0]*pe[0]+pe[1]*pe[1];

				if (sqrt(x2) < m_delt) //inlier
					x2 = 1;
				else // outliers 
					x2 = m_delt/sqrt(x2);
			}			

			pba_jacobianEachPts( m_KR, m_KdA, m_KdB, m_KdG, pa, ppt, nM, nN, nP, pAM, pAA, pPA, pPB );

			//U
			pos = sba_crsm_elmidx(Uidxij, nP, nP);		
			ppUpa = U + pos*(6*6);

			for ( ii = 0; ii < 6; ii++ ) for( jj = ii; jj < 6; jj++ )	//diag
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii]*pPA[k*6+jj];
				ppUpa[ii*6+jj] += sum*x2;
			}

			//ea
			pea = ea + nP*6;
			for ( ii = 0; ii < 6; ii++ )
			{
				for( jj = 0, sum = 0; jj < 2; jj++ )
					sum += pPA[jj*6+ii] * pe[jj];
				pea[ii] += sum*x2;
			}
			//V
			pV = V + nF*3*3;
			for ( ii = 0; ii < 3; ii++ ) for ( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii]*pPB[k*3+jj];
				pV[ii*3+jj] += sum*x2;

			}
			//eb
			peb = eb + nF*3;
			for ( ii = 0; ii < 3; ii++ )
			{
				for( k = 0, sum = 0; k < 2; k++ )
					sum += pPB[k*3+ii] * pe[k];
				peb[ii] += sum*x2;
			}
			//W
			pW = W + pos2 * 3*6;
			for ( ii = 0; ii < 6; ii++ )	for( jj = 0; jj < 3; jj++ )
			{
				for ( k = 0, sum = 0; k < 2; k++ )
					sum += pPA[k*6+ii] * pPB[k*3+jj];
				pW[ii*3+jj] += sum*x2;
			}

			if ( nP == nM )
			{
				cur++;
				pos2++;
				continue;
			}
			else
				if( nP == nN )
				{
					//U		  
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		//main archor * main archor
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		//main archor * associate archor
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						pea[ii+3] += sum*x2;
					}
					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}
				else
				{
					//U
					pos = sba_crsm_elmidx(Uidxij, nM, nM );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAM[k*3+ii] * pAM[k*3+jj];
						}
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					pos = sba_crsm_elmidx(Uidxij, nN, nN );		
					ppUpa = U + pos*(6*6);
					for ( ii = 0; ii < 3; ii++ ) for( jj = ii; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
						{
							sum += pAA[k*3+ii] * pAA[k*3+jj];
						}
						ppUpa[(ii+3)*6+3+jj] += sum*x2;
					}

					if ( nM < nN )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAM[k*3+ii] * pAA[k*3+jj];
							}
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
							{
								sum += pAA[k*3+ii] * pAM[k*3+jj];
							}
							ppUpa[(ii+3)*6+3+jj] += sum*x2;
						}
					}

					if ( nM < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nM, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAM[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nM );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAM[k*3+jj];
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}

					if ( nN < nP )
					{
						pos = sba_crsm_elmidx(Uidxij, nN, nP );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 3; ii++ ) for( jj = 0; jj < 6; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum += pAA[k*3+ii] * pPA[k*6+jj];
							ppUpa[(ii+3)*6+jj] += sum*x2;
						}
					}	
					else
					{
						pos = sba_crsm_elmidx(Uidxij, nP, nN );		
						ppUpa = U + pos*(6*6);
						for ( ii = 0; ii < 6; ii++ ) for( jj = 0; jj < 3; jj++ )
						{
							for ( k = 0, sum = 0; k < 2; k++ )
								sum +=  pPA[k*6+ii]*pAA[k*3+jj];
							ppUpa[ii*6+jj+3] += sum*x2;
						}
					}


					//ea
					pea = ea + nM*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pe[k];
						pea[ii+3] += sum*x2;
					}

					pea = ea + nN*6;
					for ( ii = 0; ii < 3; ii++ )
					{
						for( k = 0, sum = 0; k < 2; k++ )
						pea[ii+3] += sum*x2;
					}  

					//W
					pos = pos2 + (m_archorSort[i*2]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAM[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum*x2;
					}

					pos = pos2 + (m_archorSort[i*2+1]-j);
					pW = W + pos*6*3;
					for ( ii = 0; ii < 3; ii++ )	for( jj = 0; jj < 3; jj++ )
					{
						for ( k = 0, sum = 0; k < 2; k++ )
							sum += pAA[k*3+ii] * pPB[k*3+jj];
						pW[(ii+3)*3+jj] += sum*x2;
					}
				}

				pos2++;
				cur++;
		}
	}
}
*/

double CLCRP::nrmL2xmy(double *const e, const double *const x, const double *const y, const int n)
{
	const int blocksize=8, bpwr=3; /* 8=2^3 */
	register int i;
	int j1, j2, j3, j4, j5, j6, j7;
	int blockn;
	register double sum0=0.0, sum1=0.0, sum2=0.0, sum3=0.0;

	/* n may not be divisible by blocksize, 
	* go as near as we can first, then tidy up.
	*/
	//printf("nobs%d\n", n);//Add by Lu
	blockn = (n>>bpwr)<<bpwr; /* (n / blocksize) * blocksize; */
	//printf("cycle_num%d, blocksize%d", blockn - 1, blocksize);//Add by Lu

	/* unroll the loop in blocks of `blocksize'; looping downwards gains some more speed */
	for(i=blockn-1; i>0; i-=blocksize){
		e[i ]=x[i ]-y[i ]; sum0+=e[i ]*e[i ];
		j1=i-1; e[j1]=x[j1]-y[j1]; sum1+=e[j1]*e[j1];
		j2=i-2; e[j2]=x[j2]-y[j2]; sum2+=e[j2]*e[j2];
		j3=i-3; e[j3]=x[j3]-y[j3]; sum3+=e[j3]*e[j3];
		j4=i-4; e[j4]=x[j4]-y[j4]; sum0+=e[j4]*e[j4];
		j5=i-5; e[j5]=x[j5]-y[j5]; sum1+=e[j5]*e[j5];
		j6=i-6; e[j6]=x[j6]-y[j6]; sum2+=e[j6]*e[j6];
		j7=i-7; e[j7]=x[j7]-y[j7]; sum3+=e[j7]*e[j7];
	}
	
	i=blockn;
	if(i<n){ 
		switch(n - i){ 
		case 7 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 6 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 5 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 4 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 3 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 2 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		case 1 : e[i]=x[i]-y[i]; sum0+=e[i]*e[i]; ++i;
		}
	}

	return sum0+sum1+sum2+sum3;
}

void CLCRP::lcrp_cost(double* p, double* hx, int* archor)
{
	register int i;
	int cnp, pnp, mnp;
	double* pa, * pb, * ppt, * pmeas, * ppa;
	int nF, nP, nM, nN;
	int m;

	cnp = 5, pnp = 3, mnp = 1;
	m = m_ncams;
	pa = p; pb = p + m * cnp;

	for (i = 0; i < m_n2Dprojs; i++)
	{
		nF = m_feature[i];
		nP = m_photo[i];

		ppt = pb + nF * pnp;//��ά������
		ppa = pa + nP * cnp + 2;//ê��λ��ָ��
		pmeas = hx + i * mnp; // set pmeas to point to hx_ij �洢��ͶӰ��׶

		//printf("%d %d\n", nF,nP);
		//nM = archor[nF * 3 + 1];
		//nN = archor[nF * 3 + 2];
		
		//KR�������ָ�룬��ά��ָ�룬��ê�㣬��ê�㣬��ǰê�㣬��ͶӰuv
		lcrp_reprojectEachPts(m_R,m_KR, ppa, ppt, nM, nN, nP, pmeas);
		//printf("%f %f\n", pmeas[0], pmeas[1]);
	}
}

void CLCRP::readNpointsAndNprojectionsFromProj(FILE *fp, int &n3Dpts, int &nprojs)
{
	int nfirst, lineno, npts, nframes, ch, n;
	nprojs = 0;
	n3Dpts = 0;
	npts = 0;

	/* #parameters for the first line */
	nfirst=countNDoubles(fp);

	//*n3Dpts=*nprojs=lineno=npts=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if(feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n=readNInts(fp, &nframes, 1);
		if(n!=1)
		{
			fprintf(stderr, "readNpointsAndNprojections(): error reading input file, line %d: "
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		SKIP_LINE(fp);
		nprojs+=nframes;
		++npts;
	}

	n3Dpts=npts;
}

void CLCRP::readPointProjections(FILE *fp,double *imgpts, int *photo,int* imgptsSum, int n3Dpts, int n2Dprojs )
{
	int nframes, ch, lineno, ptno, frameno, n;
	int i;
	int nproj2D = 0;

	lineno=ptno=0;
	while(!feof(fp))
	{
		if((ch=fgetc(fp))=='#')
		{ /* skip comments */
			SKIP_LINE(fp);
			lineno++;

			continue;
		}

		if(feof(fp)) break;

		ungetc(ch, fp);

		n=readNInts(fp, &nframes, 1);  /* read in number of image projections */
		if(n!=1)
		{
			fprintf(stderr, "lcrp_readProjectionAndInitilizeFeature(): error reading input file, line %d:\n"
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		imgptsSum[ptno] = nframes;

		for(i=0; i<nframes; ++i)
		{
			n=readNInts(fp, &frameno, 1); /* read in frame number... */

			photo[nproj2D] = frameno;

			n+=readNDoubles(fp, imgpts+nproj2D*2, 2); /* ...and image projection */

			nproj2D++;
		}
		fscanf(fp, "\n"); // consume trailing newline

		lineno++;
		ptno++;
	}
}

void CLCRP::readImagePts( const char* szProj, double **imgpts, int **photo,int** imgptsSum, int &n3Dpts, int &n2Dprojs )
{
	FILE  *fpp;
	if((fpp=fopen(szProj, "r"))==NULL){
		fprintf(stderr, "cannot open file %s, exiting\n", szProj);
		exit(1);
	}
	readNpointsAndNprojectionsFromProj(fpp, n3Dpts, n2Dprojs);

	*imgpts = (double*)malloc(n2Dprojs*2*sizeof(double));
	if(*imgpts==NULL){
		fprintf(stderr, "memory allocation for 'imgpts' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*photo = (int*)malloc(n2Dprojs*sizeof(int)); 
	if (*photo==NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*imgptsSum = (int*)malloc(n3Dpts*sizeof(int));
	if (*imgptsSum==NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	rewind(fpp);
	readPointProjections(fpp, *imgpts,*photo,*imgptsSum, n3Dpts, n2Dprojs );

	fclose(fpp);
}

/*��ʼ����ê��*/
bool CLCRP::lcrp_initializeMainArchor( 
	double* imgpts,	// ͶӰ��� 2D ͼ�����꣨x, y����
	double* camera,	// �����Σ�����λ�ú���ת����
	double* K,		// ����ڲξ��󣨽��ࡢ����λ�õȣ���
	double* feature,// �������������Ϣ����λ��Azimuth��������Elevation��3D ��ʼ������
	int nP,			// ��ǰê���Ӧ��ͼ���ţ����ID����
	int FID,		// ������ı�š�
	double* KR )	// �ڲ�����ת����ĳ˻���K * R����
{
	/*
	bool bT = pba_initializeMainArchor( 
					projs+feastart*2,	//��¼ÿ��ͶӰ���xy��2��������¼&�洢, feastartһ��ʼ��0�����Ϊ2
					m_motstruct,		//��������������� + ��ά����������
					m_K,				//�ڲΣ���֪����
					m_motstruct+m_ncams*cnp+ptno*3,//m_ncams�Ƕ��ٸ����/��Ƭ,cnp��6��ʼҲ������εļ����ptno��0��ʼ��*3Ϊ��ά��ļ��
					nphoto[count],		//nphoto����ȫ��2DͶӰ������
					ptno,				//ptno����3D��ĵ�ǰ���
					m_KR );				//KR���߷��̱任����
	*/
	//solve  KRX = x
	Vector3d x;//�����������/ͶӰ����ê���XYZ����
	//������ά��϶��ǲ�֪���ġ�
	if (m_bProvideXYZ)//����������3D������֪
	{
		x(0) = m_XYZ[FID*3] - *(camera + nP*6+3);
		x(1) = m_XYZ[FID*3+1] - *(camera + nP*6+4);
		x(2) = m_XYZ[FID*3+2] - *(camera + nP*6+5);
	}
	
	else//����������3D����δ֪��������������3D���꣬ע���������Ƕȱ�ʾ�������3D����
	{
		double *ptr = m_KR + nP*9;	//m_KR���ڲξ���x��ת����np�ǵ�ǰ3D�����ţ�*9�Ǽ��
		Matrix3d  A;
		A << ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7], ptr[8];//�����ʼ����KR����
		//�� KR ����ȡ 3x3 ����A��
		double matx[3];				//��ε��������
		matx[0] = imgpts[0];
		matx[1] = imgpts[1];
		matx[2] = 1;
		//Matrix A��KR���󣬶������������飩��2DͶӰ����
		Vector3d  b(matx);
		//Matrix A��KR���󣬶�Matrix b���������
		//���¾�������Ͳ�˵�ˡ�����2DͶӰ��ͨ��KR���󣬿������3D�����꣬��ʱ��û�м������أ�����XYZ����
		x = A.colPivHouseholderQr().solve(b);//���Ax=b  ������������������ê���XYZ����
		///*------------------------------------������*/
		//printf("%f %f %f\n", x(0), x(1), x(2));
		//printf("%f %f %f\n", ptr[6], ptr[7], ptr[8]);
		///*------------------------------------������*/
		//x�Ѿ���IJRR�Ĺ�ʽ��15��
	}

	double* pKR = KR + nP*9;//KR��ʾ�ڲξ�������ת����ĳ˻�����m_KR��һ�������KR�Ǵ������ģ���m_KR�Ǽ�¼��
	double t = pKR[6]*x(0) + pKR[7]*x(1) + pKR[8]*x(2);	//������
	///*------------------------------------������----------------------------------*/
	//printf("%f %f %f\n", pKR[6], pKR[7], pKR[8]);
	///*------------------------------------������----------------------------------*/
	//compute azimuth and elevation angle
	double dDAngle = atan2( x(0), x(2) );				//\Psi��Ҳ����Azimuth angle-��ê��
	double dHAngle = atan2( x(1), sqrt(x(0)*x(0)+ x(2)*x(2)) );//\theta��Ҳ����Elevation angle
	//feature�����Ӳ�ǵ�����
	feature[0] = dDAngle;
	feature[1] = dHAngle;
	feature[2] = 0;

	if ( t < 0 )//��
		return true;
	else
		return false;
}

/*��ʼ����ê��*/
bool CLCRP::lcrp_initializeAssoArchor( 
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int nMI,
	int nAI,
	int FID,
	bool bLast )
{
	/*
	bool bT = pba_initializeAssoArchor( 
					projs+feastart*2,	//��¼ÿ��ͶӰ���xy��2��������¼&�洢, feastartһ��ʼ��0�����Ϊ2
					nphoto+feastart,	//��¼ÿ�������ţ�feastartһ��ʼ��0�����Ϊ1
					m_motstruct,
					m_K,
					m_motstruct+m_ncams*cnp+ptno*3,
					0,
					1,
					ptno,
					bLast );
	*/
	int nM = photo[nMI];                              //��ê���
	int nA = photo[nAI];                              //��ê���

	Vector3d  xM, xA;

	if (m_bProvideXYZ)
	{
		xM[0] = m_XYZ[FID*3]   - *(camera + nM*6+3);
		xM[1] = m_XYZ[FID*3+1] - *(camera + nM*6+4);
		xM[2] = m_XYZ[FID*3+2] - *(camera + nM*6+5);

		xA[0] = m_XYZ[FID*3]   - *(camera + nA*6+3);
		xA[1] = m_XYZ[FID*3+1] - *(camera + nA*6+4);
		xA[2] = m_XYZ[FID*3+2] - *(camera + nA*6+5);
	}
	else
	{
		//Main anchor ray
		double *ptr1 = m_KR + nM*9;
		Matrix3d  AM;	
		AM << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

		double matxM[3];
		matxM[0] = *(imgpts+2*nMI);
		matxM[1] = *(imgpts+2*nMI+1);
		matxM[2] = 1;

		Vector3d  bM(matxM);
		xM = AM.colPivHouseholderQr().solve(bM);			//�������������ê���XYZ����

		//Associate archor ray
		double *ptr2 = m_KR + nA*9;
		Matrix3d  AA;	
		AA << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

		double matxA[3];
		matxA[0] = *(imgpts+2*nAI);
		matxA[1] = *(imgpts+2*nAI+1);
		matxA[2] = 1;

		Vector3d  bA(matxA);
		xA = AA.colPivHouseholderQr().solve(bA);			//����������ڸ�ê���XYZ����
	}
		
	//Parallax Angle
	double dDot = xM(0)*xA(0) + xM(1)*xA(1) + xM(2)*xA(2);			//IJRR�Ĺ�ʽ��3.13���ķ���

	double dDisM = sqrt( xM(0)*xM(0)+xM(1)*xM(1)+xM(2)*xM(2) );		//IJRR�Ĺ�ʽ��3.13���ķ�ĸ
	double dDisA = sqrt( xA(0)*xA(0)+xA(1)*xA(1)+xA(2)*xA(2) );		//IJRR�Ĺ�ʽ��3.13���ķ�ĸ

	if (dDot/(dDisM*dDisA)>1)			//�������
		feature[2] = 0;
	else if (dDot/(dDisM*dDisA)<-1)		//�������
		feature[2] = PI;
	else
	{
		double dw = acos( dDot/(dDisM*dDisA) );
		feature[2] = dw;
	}
	//�������/omega��

	
	double pti2k[3];
	//+3��X��Ҳ���Ǹ�ê��-��ê�㣬Ҳ����t^m_a��Ҳ���Ǵ���ê��m����ê���������X����XA-XM
	pti2k[0] = *(camera + nA*6+3) - *(camera + nM*6+3);    
	//YA-YM��Y����
	pti2k[1] = *(camera + nA*6+4) - *(camera + nM*6+4);    
	//ZA-ZM��Z����
	pti2k[2] = *(camera + nA*6+5) - *(camera + nM*6+5);    
	//��������pti2k��������t_m->a

	//dDot1������ê�㵽������F_j��x����xM��t_m->a�ĵ��
	double dDot1 = xM[0]*pti2k[0] + xM[1]*pti2k[1] + xM[2]*pti2k[2];
	double dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
	double tmp = dDot1/(dDisM*dDisi2k);//��ʵ����IJRRͼ2��/varphi��cos
	double dW2;
	if (tmp > 1)						//�������
		dW2 = 0;
	if ( tmp < -1)						//�������
		dW2 = PI;
	else
		dW2  = acos( tmp );//��ʵ����IJRRͼ2��/varphi

	return true;//dW2��/varphi_j;dw=feature[2]��/omega_j;�����ƺ�����ҪdW2������˵��tmp�м�ֵ��IJRR��BA��ʽ�У�
}

/*��ʼ������ê��*/
bool CLCRP::lcrp_initializeOtheArchors( 
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int* archorSort,
	int nfeacout,
	int nOI,
	int FID )
{
	/*
	* bAdjust = pba_initializeOtheArchors( 
					projs+feastart*2,
					nphoto+feastart,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams*cnp + ptno * 3,
					ptr2,
					sum,
					i,
					ptno );
	*/
	static int i = 0;
	double dw = feature[2];                   //�Ӳ��
	double dwNew;
	double dmaxw = dw;
	int   nNewI = 0;
	bool bAdjust = false;
	double dDot,dDisM,dDisA;

	if ( dw < MAXARCHOR   )
	{
		//current archor vector 
		int nO = photo[nOI];

		Vector3d  xO;
		if ( m_bProvideXYZ)
		{
			xO(0) = m_XYZ[FID*3] - *(camera + nO*6 + 3);
			xO(1) = m_XYZ[FID*3+1]-*(camera + nO*6 + 4);
			xO(2) = m_XYZ[FID*3+2]-*(camera + nO*6 + 5);
		}
		else
		{
			double *ptr1 = m_KR + nO*9;
			Matrix3d  AO;	
			AO << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

			double matxO[3];
			matxO[0] = *(imgpts+nOI*2);
			matxO[1] = *(imgpts+nOI*2+1);
			matxO[2] = 1;
			///*-------------------------------������--------------------------------*/
			//printf("%f %f %f\n", matxO[0], matxO[1], matxO[2]);
			///*-------------------------------������--------------------------------*/
			Vector3d  bO(matxO);
			xO = AO.colPivHouseholderQr().solve(bO);
		}

		double dDAngle = atan2( xO(0), xO(2) );
		double dHAngle = atan2( xO(1), sqrt(xO(0)*xO(0)+ xO(2)*xO(2)) );

		for ( i = 0; i < nfeacout; i++ )
		{
			//Main Archor Vector
			int nM = photo[i];                             //��Ƭ��
			Vector3d  xM;

			if (m_bProvideXYZ)
			{
				xM(0) = m_XYZ[FID*3]  -*(camera + nM*6+3);
				xM(1) = m_XYZ[FID*3+1]-*(camera + nM*6+4);
				xM(2) = m_XYZ[FID*3+2]-*(camera + nM*6+5);
			}
			else
			{
				double *ptr2 = m_KR + nM*9;
				Matrix3d  AM;	
				AM << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

				double matxM[3];
				matxM[0] = *(imgpts+i*2);
				matxM[1] = *(imgpts+i*2+1);
				matxM[2] = 1;

				Vector3d  bM(matxM);
				xM = AM.colPivHouseholderQr().solve(bM);
			}

			//Parallax angle between current archor and main archor
			dDot = xM(0)*xO(0) + xM(1)*xO(1) + xM(2)*xO(2);
			dDisM = sqrt( xM(0)*xM(0)+xM(1)*xM(1)+xM(2)*xM(2) );
			dDisA = sqrt( xO(0)*xO(0)+xO(1)*xO(1)+xO(2)*xO(2) );

			if( dDot/(dDisM*dDisA) > 1 )
				dwNew = 0;
			else if(dDot/(dDisM*dDisA)<-1)
				dwNew = PI;
			else
				dwNew = acos( dDot/(dDisM*dDisA) );

			if ( dwNew > dmaxw )
			{
				dmaxw = dwNew;
				archorSort[0] = nOI;   //�����ê���
				archorSort[1] = i;     //��Ÿ�ê���
				feature[0] = dDAngle;
				feature[1] = dHAngle;
				feature[2] = dmaxw;
				bAdjust = true;
			}
		}
	}
	return bAdjust;
}


bool CLCRP::lcrp_initializeOtheArchors_Mindw(
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int* archorSort,
	int nfeacout,
	int nOI,
	int FID)
{
	/*
	* bAdjust = pba_initializeOtheArchors(
					projs+feastart*2,
					nphoto+feastart,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams*cnp + ptno * 3,
					ptr2,
					sum,
					i,
					ptno );
	*/
	static int i = 0;
	double dw = feature[2];                   //�Ӳ��
	double dwNew;
	double dminw = dw;
	int   nNewI = 0;
	bool bAdjust = false;
	double dDot, dDisM, dDisA;

	if (dw < MAXARCHOR)
	{
		//current archor vector 
		int nO = photo[nOI];

		Vector3d  xO;
		if (m_bProvideXYZ)
		{
			xO(0) = m_XYZ[FID * 3] - *(camera + nO * 6 + 3);
			xO(1) = m_XYZ[FID * 3 + 1] - *(camera + nO * 6 + 4);
			xO(2) = m_XYZ[FID * 3 + 2] - *(camera + nO * 6 + 5);
		}
		else
		{
			double* ptr1 = m_KR + nO * 9;
			Matrix3d  AO;
			AO << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

			double matxO[3];
			matxO[0] = *(imgpts + nOI * 2);
			matxO[1] = *(imgpts + nOI * 2 + 1);
			matxO[2] = 1;
			///*-------------------------------������--------------------------------*/
			//printf("%f %f %f\n", matxO[0], matxO[1], matxO[2]);
			///*-------------------------------������--------------------------------*/
			Vector3d  bO(matxO);
			xO = AO.colPivHouseholderQr().solve(bO);
		}

		double dDAngle = atan2(xO(0), xO(2));
		double dHAngle = atan2(xO(1), sqrt(xO(0) * xO(0) + xO(2) * xO(2)));

		for (i = 0; i < nfeacout; i++)
		{
			//Main Archor Vector
			int nM = photo[i];                             //��Ƭ��
			Vector3d  xM;

			if (m_bProvideXYZ)
			{
				xM(0) = m_XYZ[FID * 3] - *(camera + nM * 6 + 3);
				xM(1) = m_XYZ[FID * 3 + 1] - *(camera + nM * 6 + 4);
				xM(2) = m_XYZ[FID * 3 + 2] - *(camera + nM * 6 + 5);
			}
			else
			{
				double* ptr2 = m_KR + nM * 9;
				Matrix3d  AM;
				AM << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

				double matxM[3];
				matxM[0] = *(imgpts + i * 2);
				matxM[1] = *(imgpts + i * 2 + 1);
				matxM[2] = 1;

				Vector3d  bM(matxM);
				xM = AM.colPivHouseholderQr().solve(bM);
			}

			//Parallax angle between current archor and main archor
			dDot = xM(0) * xO(0) + xM(1) * xO(1) + xM(2) * xO(2);
			dDisM = sqrt(xM(0) * xM(0) + xM(1) * xM(1) + xM(2) * xM(2));
			dDisA = sqrt(xO(0) * xO(0) + xO(1) * xO(1) + xO(2) * xO(2));

			if (dDot / (dDisM * dDisA) > 1)
				dwNew = 0;
			else if (dDot / (dDisM * dDisA) < -1)
				dwNew = PI;
			else
				dwNew = acos(dDot / (dDisM * dDisA));

			if (dwNew < dminw)
			{
				dminw = dwNew;
				archorSort[0] = nOI;   //�����ê���
				archorSort[1] = i;     //��Ÿ�ê���
				feature[0] = dDAngle;
				feature[1] = dHAngle;
				feature[2] = dminw;
				bAdjust = true;
			}
		}
	}
	return bAdjust;
}

//--------------------------------------------------------
int CLCRP::lcrp_angle2xytGN( double *p )
{
	static int i, j;
	double* pAngle;
	double xj[3], xk[3];
	double Tik[3];
	int nM, nA;
	double Dik;
	double w, w2;
	static double aa;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		pAngle = p + m_ncams*6 + i*3;
		w = *(p + m_ncams*6 + i*3 + 2);

		if ( i == 408741 )
		{
			aa = w;
		}

		xj[0] = sin(*(pAngle)) * cos(*(pAngle+1));
		xj[1] = sin(*(pAngle+1));
		xj[2] = cos(*(pAngle)) * cos(*(pAngle+1));

		nM = *(m_archor+i*3+1);
		nA = *(m_archor+i*3+2);

		Tik[0] = -*(p+nM*6+3)+*(p+nA*6+3);
		Tik[1] = -*(p+nM*6+4)+*(p+nA*6+4);
		Tik[2] = -*(p+nM*6+5)+*(p+nA*6+5);
		
		Dik = sqrt( Tik[0]*Tik[0] + Tik[1]*Tik[1] + Tik[2]*Tik[2] );

		w2 = acos( (xj[0]*Tik[0]+xj[1]*Tik[1]+xj[2]*Tik[2])/Dik );

		xk[0] = ( Dik * sin(w2+w) * xj[0] )/sin(w);
		xk[1] = ( Dik * sin(w2+w) * xj[1] )/sin(w);
		xk[2] = ( Dik * sin(w2+w) * xj[2] )/sin(w);

		*(p+m_ncams*6+i*3)   = *(p+nM*6+3) + xk[0];
		*(p+m_ncams*6+i*3+1) = *(p+nM*6+4) + xk[1];
		*(p+m_ncams*6+i*3+2) = *(p+nM*6+5) + xk[2];
	}

	return 1;
}

void CLCRP::lcrp_saveXYZ( const char* camera, const char* sz3Dpt, double *p,bool gn  )
{
	static int i = 0;
	double dx, dy, dz;
	FILE *fp = NULL, *fpc = NULL;
	
	//transform from angle to xyz
	if( gn )
		lcrp_angle2xytGN(p);
	else
		lcrp_angle2xytLM(p);
	
	//save camera poss
	if( camera != NULL )
	{
		fpc = fopen( camera, "w" );

		for( i = 0; i < m_ncams; i++ )
		{
			fprintf( fpc, "%0.5lf     %0.5lf      %0.5lf     %0.5lf       %0.5lf     %0.5lf\n",
				*(p+i*6), *(p+i*6+1), *(p+i*6+2), *(p+i*6+3), *(p+i*6+4), *(p+i*6+5) );
		}
		fclose(fpc);
	}	

	//save features xyz
	if ( sz3Dpt!= NULL )
	{
		fp	= fopen( sz3Dpt, "w" );
		for ( i = 0; i < m_n3Dpts; i++ )
		{
			dx = *(p+m_ncams*6+i*3);
			dy = *(p+m_ncams*6+i*3+1);
			dz = *(p+m_ncams*6+i*3+2);
			fprintf( fp, "%0.5lf     %0.5lf     %0.5lf\n", dx, dy, dz );
		}
		fclose(fp);
	}

}

int	CLCRP::lcrp_angle2xytLM( double *p )
{
	static int i, j;
	double* pAngle;
	double xj[3], xk[3];
	double Tik[3];
	int nM, nA;
	double Dik;
	double w, w2;
	int cnp = 6, pnp = 3;
	for ( i = 0; i < m_n3Dpts; i++ )
	{
		pAngle = p + m_ncams*6 + i*3;
		w = *(p + m_ncams*6 + i*3 + 2);

		xj[0] = sin(*(pAngle)) * cos(*(pAngle+1));//x
		xj[1] = sin(*(pAngle+1));//y
		xj[2] = cos(*(pAngle)) * cos(*(pAngle+1));//z

		nM = *(m_archor+i*3+1);
		nA = *(m_archor+i*3+2);

		Tik[0] = -*(p+nM*6+3)+*(p+nA*6+3);
		Tik[1] = -*(p+nM*6+4)+*(p+nA*6+4);
		Tik[2] = -*(p+nM*6+5)+*(p+nA*6+5);

		Dik = sqrt( Tik[0]*Tik[0] + Tik[1]*Tik[1] + Tik[2]*Tik[2] );

		w2 = acos( (xj[0]*Tik[0]+xj[1]*Tik[1]+xj[2]*Tik[2])/Dik );

		//�������������ê���xyz����
		xk[0] = ( Dik * sin(w2+w) * xj[0] )/sin(w);
		xk[1] = ( Dik * sin(w2+w) * xj[1] )/sin(w);
		xk[2] = ( Dik * sin(w2+w) * xj[2] )/sin(w);

		//�������ȫ��xyz����
		*(p+m_ncams*6+i*3)   = *(p+nM*6+3) + xk[0];
		*(p+m_ncams*6+i*3+1) = *(p+nM*6+4) + xk[1];
		*(p+m_ncams*6+i*3+2) = *(p+nM*6+5) + xk[2];
	}

	return 1;
}

bool CLCRP::lcrp_parseArgs( int argc, char* argv[] )
{
	int i;
	string param;
	bool bSuccess, bRKF;

	for (i = 1; i < argc; i++) 
	{
		bSuccess = false;
		string name = argv[i];

		if (name[0] != '-') { // each param has to start with at least one dash
			return false;
		}
		
		string::size_type dashPos = name.find_first_not_of('-');
		if (dashPos != string::npos)
			name = name.substr(dashPos);

		if ( strcmp(name.c_str(), "help") == 0 )
		{
			lcrp_printHelp();
			return false;
		}
		
		if ( strcmp( name.c_str(), "cam") == 0 )
		{
			i++;
			param = argv[i];
			m_szCameraInit = (char*)malloc(param.length());
			strcpy( m_szCameraInit, param.c_str() );
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "fea") == 0 )
		{
			i++;
			param = argv[i];
			m_szFeatures = (char*)malloc(param.length());
			strcpy( m_szFeatures, param.c_str() );
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "calib") == 0 )
		{
			i++;
			param = argv[i];
			m_szCalibration = (char*)malloc(param.length());
			strcpy( m_szCalibration, param.c_str() );			
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "pose") == 0 )
		{
			i++;
			param = argv[i];
			m_szCamePose = (char*)malloc(param.length());
			strcpy( m_szCamePose, param.c_str() );			
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "3D") == 0 )
		{
			i++;
			param = argv[i];
			m_sz3Dpts = (char*)malloc(param.length());
			strcpy( m_sz3Dpts, param.c_str() );			
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "report") == 0 )
		{
			i++;
			param = argv[i];
			m_szReport = (char*)malloc(param.length());
			strcpy( m_szReport, param.c_str() );
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "xyz") == 0 )
		{
			i++;
			param = argv[i];
			m_szXYZ = (char*)malloc(param.length());
			strcpy( m_szXYZ, param.c_str() );
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "i") == 0 )
		{
			i++;
			param = argv[i];
			m_nMaxIter = atoi(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "robustKernel") == 0 )
		{
			m_bRobustKernel = true;
			bRKF = false;
			i++;
			param = argv[i];
			if ( strcmp(param.c_str(), "Huber") == 0 )
			{
				m_nRobustType = 2;
				bRKF = true;
			}

			if ( strcmp(param.c_str(), "Cauchy") == 0 )
			{
				m_nRobustType = 1;
				bRKF = true;
			}
			
			bSuccess = true;

			if ( !bRKF )
			{
				printf( "LCRP: Must input right robust kernel function!\n" );
				return false;
			}
		}

		if ( strcmp( name.c_str(), "solve") == 0 )
		{
			i++;
			param = argv[i];
			if( strcmp(param.c_str(), "LM") == 0)
				m_bsolverLM = true;

			if( strcmp(param.c_str(), "GN") == 0)
				m_bsolverGN = true;

			bSuccess = true;
		}		

		if ( strcmp( name.c_str(), "t") == 0 )
		{
			i++;
			param = argv[i];
			
			m_Tau = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "e1") == 0 )
		{
			i++;
			param = argv[i];

			m_e1 = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "e2") == 0 )
		{
			i++;
			param = argv[i];

			m_e2 = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "e3") == 0 )
		{
			i++;
			param = argv[i];

			m_e3 = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "e4") == 0 )
		{
			i++;
			param = argv[i];

			m_e4 = atof(param.c_str());
			bSuccess = true;
		}

		if ( strcmp( name.c_str(), "robustKernelWidth") == 0 )
		{
			i++;
			param = argv[i];

			m_delt = atof(param.c_str());
			bSuccess = true;
		}

		if (!bSuccess)
		{
			printf( "LCRP error: %s command is wrong!\n", name.c_str() );
			return false;
		}
		
	}

	return true;
			
}

void CLCRP::lcrp_printHelp()
{
	printf( "Light Cone Reprojection Principle General Options\n" );
	printf( "\n" );

	printf( "-cam			Provide initial camera pose.\n" );
	printf( "-fea			Provide features.\n" );
	printf( "-calib			Provide calibration.\n" );
	printf( "-xyz			Provide initial XYZ.\n" );
	printf( "-pose			Output optimal camera pose.\n" );
	printf( "-3D			Output optimal 3D point cloud.\n" );
	printf( "-report			Output report.\n" );
	printf( "-solve			Solve method including LevenbergMarquart(LM) and Gauss-Newton(GN).\n" );
	printf( "-i			Max Iteration.\n" );
	printf( "-t			LevenbergMarquart parameters.\n" );
	printf( "-robustKernel		use Cauchy Robust Kernel Function.\n" );
	printf( "-robustKernelWidth		width for the robust Kernel.\n" );	
}
