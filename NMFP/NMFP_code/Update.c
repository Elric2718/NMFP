#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>



SEXP Update(SEXP RW, SEXP RH, SEXP RV,SEXP RCset, SEXP Ralpha){
  const int alpha= asReal(Ralpha);		//character	-->char

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



    //some transition parameters
    double SumH=0;
    double SumW=0;
    double SumProduct=0;
    double newH [nCol];
    double newW [nRow];
    double Cpar=0;



    //update now
    for(int i=0;i<rank;i++){
        //update H
        for(int j=0;j<nCol;j++){
            SumH=0;
            SumW=0;

            for(int k=0;k<nRow;k++){
            SumW=SumW+W[k][i];
            }

            for(int k=0;k<nRow;k++){
            SumProduct=0;
                for(int n=0;n<rank;n++){
                SumProduct=SumProduct+W[k][n]*H[n][j];
                }
            SumH=SumH+V[k][j]*W[k][i]/SumProduct;
            }

            newH[j]=SumH*(H[i][j]/SumW);
        }


        for(int j=0;j<nCol;j++){
            H[i][j]=newH[j];
        }




        //update W
        for(int j=0;j<nRow;j++){
            SumH=0;
            SumW=0;
            Cpar=0;

            for(int k=0;k<nRow;k++){
                if(Cset[j][k]==1){
                    Cpar=Cpar+W[k][i];
                }
            }

            for(int k=0;k<nCol;k++){
            SumH=SumH+H[i][k];
            }

            for(int k=0;k<nCol;k++){
                SumProduct=0;
                for(int n=0;n<rank;n++){
                    SumProduct=SumProduct+W[j][n]*H[n][k];
                }
            SumW=SumW+V[j][k]*H[i][k]/SumProduct;
            }

            newW[j]=SumW*W[j][i]/(SumH+2*alpha*Cpar);
        }

       // SumW=0;
       // for(int j=0;j<nRow;j++){
       //     SumW=SumW+newW[j];
       // }
        for(int j=0;j<nRow;j++){
            W[j][i]=newW[j];  // /SumW;
        }
       // for(int j=0;j<nCol;j++){
       //    H[i][j]=H[i][j]*SumW;
       // }
    }




    //from c type back to R type
//	SEXP num;
//	PROTECT(num = allocVector(REALSXP, 1));
//	REAL(num)[0] = d;
//	UNPROTECT(1);//<-1,or 2... is determined by how many PROTECT you are using
//	return num;
    //create a matrix for R

    SEXP ans;
    PROTECT(ans = allocMatrix(REALSXP, nRow+nCol, rank));
    for(int m=0; m<nRow+nCol; m++)
        for(int n=0; n<rank; n++){
            if(m<nRow){
            REAL(ans)[m + n*(nRow+nCol)] = W[m][n];
            }else{
            REAL(ans)[m + n*(nRow+nCol)] = H[n][m-nRow];
            }
        }
    UNPROTECT(1);
    return(ans);
}
