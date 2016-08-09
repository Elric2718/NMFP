#include <R.h>
#include <Rinternals.h>
#include <math.h>



SEXP Cpenalty(SEXP RW, SEXP RH, SEXP RV,SEXP RCset, SEXP Ralpha){
  const int alpha= asReal(Ralpha);  	//character	-->char
	
	//matrix	--> array[][]
	SEXP dimnamesW = getAttrib(RW, R_DimNamesSymbol);
	SEXP rownamesW = VECTOR_ELT(dimnamesW, 0);
	SEXP colnamesW = VECTOR_ELT(dimnamesW, 1);
	SEXP dimnamesH = getAttrib(RH, R_DimNamesSymbol);
    //SEXP rownamesH = VECTOR_ELT(dimnamesH, 0);
	SEXP colnamesH = VECTOR_ELT(dimnamesH, 1);
	const int nCol = length(colnamesH);
	const int nRow = length(rownamesW);
	const int rank	= length(colnamesW);
	
	//transform RW into W
	double W [nRow][rank];
	for(int m=0; m<nRow; m++){
		for(int n=0; n<rank; n++){
            W[m][n]=REAL(RW)[m+n*nRow];
		}
	}

    //transform RH into H
    double H [rank][nCol];
        for(int m=0; m<rank; m++){
        for(int n=0; n<nCol; n++){
            H[m][n]=REAL(RH)[m+n*rank];
        }
    }

    //transform RV into V
    double V [nRow][nCol];
    for(int m=0; m<nRow; m++){
        for(int n=0; n<nCol; n++){
            V[m][n]=REAL(RV)[m+n*nRow];
        }
    }

    //transform RCset into Cset
    int Cset [nRow][nRow];
    for(int m=0; m<nRow; m++){
        for(int n=0; n<nRow; n++){
            Cset[m][n]=REAL(RCset)[m+n*nRow];
        }
    }

    //U=W%*%t(W)
    double U [nRow][nRow];
    for(int i=0;i<nRow;i++){
        for(int j=0;j<nRow;j++){
            U[i][j]=0;
            for(int k=0;k<rank;k++){
                U[i][j]=U[i][j]+W[i][k]*W[j][k];
            }
        }
    }

    //U=W%*%t(W)
    double WH [nRow][nCol];
    for(int i=0;i<nRow;i++){
        for(int j=0;j<nCol;j++){
            WH[i][j]=0;
            for(int k=0;k<rank;k++){
                WH[i][j]=WH[i][j]+W[i][k]*H[k][j];
            }
        }
    }


    double Cpenalty=0;
    double Psum=0;
    for(int i=0;i<nRow;i++){
        Psum=0;
        for(int j=0;j<nCol;j++){
            if(V[i][j]>0&&WH[i][j]>0){
                Psum=Psum+V[i][j]*log(V[i][j]/WH[i][j])-V[i][j]+WH[i][j];
            }else{
                Psum=Psum+abs(WH[i][j]-V[i][j]);
            }
            if(Cset[i][j]==1){
                Psum=Psum+alpha*U[i][j];
            }
        }
        Cpenalty=Cpenalty+Psum;
    }

    SEXP ans;
    ans=PROTECT(allocVector(REALSXP,1));
    REAL(ans)[0]=Cpenalty;
    UNPROTECT(1);
    return(ans);


			

}
