#include <R.h>
#include <Rinternals.h>
#include <math.h>



SEXP LogProb(SEXP RG, SEXP RGvisual, SEXP RNormCoef,SEXP REindex, SEXP RBLen, SEXP RIsoLen, SEXP RLenRead){
    const int LenRead= asInteger(RLenRead);    //character	-->char
	
	//matrix	--> array[][]
    SEXP dimnamesG = getAttrib(RG, R_DimNamesSymbol);
    SEXP rownamesG = VECTOR_ELT(dimnamesG, 0);
    SEXP colnamesG = VECTOR_ELT(dimnamesG, 1);
    SEXP dimnamesNc = getAttrib(RNormCoef, R_DimNamesSymbol);
    //SEXP rownamesNc = VECTOR_ELT(dimnamesNc, 0);
    SEXP colnamesNc = VECTOR_ELT(dimnamesNc, 1);
    SEXP dimnamesGv = getAttrib(RGvisual, R_DimNamesSymbol);
    SEXP rownamesGv = VECTOR_ELT(dimnamesGv, 0);
    //SEXP colnamesGvisual = VECTOR_ELT(dimnamesGvisual, 1);
    const int nCol = length(colnamesNc);
    const int nRow = length(rownamesG);
    const int rank	= length(colnamesG);
    const int nExon= length(rownamesGv);
    
	
    //transform RG into G
    double G [nRow][rank];
	for(int m=0; m<nRow; m++){
		for(int n=0; n<rank; n++){
            G[m][n]=REAL(RG)[m+n*nRow];
		}
	}

    //transform RGvisual into Gvisual
    double Gvisual [nExon][rank];
    for(int m=0; m<nExon; m++){
        for(int n=0; n<rank; n++){
            Gvisual[m][n]=REAL(RGvisual)[m+n*nExon];
        }
    }

    //transform RNormCoef into NormCoef
    double NormCoef [rank][nCol];
    for(int m=0; m<rank; m++){
        for(int n=0; n<nCol; n++){
            NormCoef[m][n]=REAL(RNormCoef)[m+n*rank];
        }
    }

    //transform REindex into Eindex
    int Eindex [nRow][2];
 
    for(int m=0; m<nRow; m++){
        for(int n=0; n<2; n++){
        	Eindex[m][n]=INTEGER(REindex)[m+n*nRow];       	
        }
    }
    

    //transform RBLen into BLen
    double BLen [nRow];
    for(int m=0; m<nRow; m++){
        BLen[m]=REAL(RBLen)[m];
    }

    //transform RIsoLen into IsoLen
    double IsoLen [rank];
    for(int m=0; m<rank; m++){
        IsoLen[m]=REAL(RIsoLen)[m];
    }


    double logProb=0;
    double SumBin=0;
    int eindex=1;



    for(int j=0;j<nCol;j++){
        for(int m=0;m<nRow;m++){
            SumBin=0;
            for(int n=0;n<rank;n++){
                if(G[m][n]==1){
                    eindex=1;
                    for(int k=0;k<2;k++){
                        if(Eindex[m][k]>0){
                        eindex=eindex*Gvisual[Eindex[m][k]-1][n];
                        }else{break;}                       
                    }           
                 if(eindex==1){
                        SumBin=SumBin+NormCoef[n][j]*BLen[m]/(IsoLen[n]-LenRead+1);
                	 }
                }
            }
            if(SumBin>0){
                logProb=logProb+log(SumBin);
            }          	
        }        
    }

    SEXP ans;
    ans=PROTECT(allocVector(REALSXP,1));
    REAL(ans)[0]=logProb;
    UNPROTECT(1);
    return(ans);


}
