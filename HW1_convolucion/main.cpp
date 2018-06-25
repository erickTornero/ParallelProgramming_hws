#include <iostream>
#include <vector>
#include <ctime>
#include <chrono>
#include <thread>
#include <array> 
#include <sstream>
#include "omp.h"
#define MATRIX_SIZE 2000
#define PATERN_SIZE 5 //9
void convolution(int * mat[MATRIX_SIZE],int pat[PATERN_SIZE][PATERN_SIZE], int *answer[MATRIX_SIZE], int cores);

void print_matrix(int *v[MATRIX_SIZE]){
    for(int i = 0; i < MATRIX_SIZE;i++){
        for(int j = 0; j < MATRIX_SIZE;j++){
            std::cout<<v[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}
int main(int argc, char *argv[]){

    int cores = 1;
    if (argc > 1){
        std::istringstream iss( argv[1]);  
        if (iss >> cores){
            // Conversion successful
        }
    }
    std::cout <<"Cores: "<<cores<<std::endl;
    std::chrono::high_resolution_clock::time_point t1,t2;
    int * matrix[MATRIX_SIZE];
    auto thread_count = std::thread::hardware_concurrency();
    /*srand(time(NULL));
    for(unsigned int i = 0; i < matrix.size();i++){
        for(unsigned int j = 0; j < matrix.at(0).size();j++){
            matrix.at(i).at(j) = rand()%255;
        }
    }*/
    for(unsigned int i = 0; i < MATRIX_SIZE;i++){
        matrix[i] = new int[MATRIX_SIZE];
        for(unsigned int j = 0; j < MATRIX_SIZE;j++){
            matrix[i][j] = j;
        }
    }
    int patern[PATERN_SIZE][PATERN_SIZE] = {{1,1,1,1,1},
                                        {1,1,1,1,1},
                                        {1,1,1,1,1},
                                        {1,1,1,1,1},
                                        {1,1,1,1,1}};
    int * ans[MATRIX_SIZE];
    for(int a = 0; a < MATRIX_SIZE; a++){
        ans[a] = new int[MATRIX_SIZE];
    }
    t1 = std::chrono::high_resolution_clock::now();
    convolution(matrix,patern,ans,cores);
    t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
    std::cout<<"Duration: "<<duration<<std::endl;
    //print_matrix(matrix);
    //std::cout<<std::endl<<std::endl<<std::endl;
    //print_matrix(ans);
    int xx = 21;
}
void convolution(int *mat[MATRIX_SIZE],int pat[PATERN_SIZE][PATERN_SIZE], int *answer[MATRIX_SIZE],int cores){
        //Variables para recorrer donde la matriz patron puede multiplicar todos los indices.
    int fromRow, toRow, fromCol, toCol;

    const int rowCenter = (PATERN_SIZE-1)/2;
    const int colCenter = (PATERN_SIZE-1)/2;
    fromRow = rowCenter;
    fromCol = colCenter;
    toRow = MATRIX_SIZE-rowCenter;
    toCol = MATRIX_SIZE-colCenter;

    int rowSizePatern = PATERN_SIZE;
    int colSizePatern = PATERN_SIZE;
    const int numberOfElements = rowSizePatern*colSizePatern;
    #pragma omp parallel for num_threads(cores) shared(answer)  schedule(static,9)
    //Iterate over all elements center and this is multiplied by all pattern matrix-avoid corners
    for(int i = fromRow; i < toRow;i++){
        for(int j = fromCol; j < toCol;j++){
            //Convolucionar pixel
            int acumm = 0;
            int rowMult = i-rowCenter;
            int colMult = j-colCenter;
           // #pragma omp simd reduction(+:acumm)
            for(int iPatern = 0; iPatern < PATERN_SIZE;iPatern++){
                for(int jPatern = 0; jPatern < PATERN_SIZE;jPatern++){
                    acumm+=pat[iPatern][jPatern]*mat[rowMult+iPatern][colMult+jPatern];
                }
            }
            answer[i][j] = acumm/numberOfElements;
        }
    }   
    //#pragma omp parallel for num_threads(cores) shared(answer)
    //This iterate over all top-up of matrix
    for(int i = 0; i < fromRow;i++){
        for(int j = 0; j < MATRIX_SIZE;j++){
            //Convolucionar pixel
            int acumm = 0;
            int rowMult = i-rowCenter;
            int colMult = j-colCenter;
            int numberOperations = 0;
            for(int iPatern = 0; iPatern < rowSizePatern;iPatern++){
                for(int jPatern = 0; jPatern < colSizePatern;jPatern++){
                    int whatrow = rowMult+iPatern;
                    int whatcol = colMult+jPatern;
                    if(( whatrow >= 0)&&(whatcol >= 0)&&(whatcol<MATRIX_SIZE)){
                        acumm+=pat[iPatern][jPatern]*mat[whatrow][whatcol];
                        numberOperations++;
                    }
                }
            }
            answer[i][j] = acumm/numberOperations;
        }
    }
    //This part is the part of above:
    // #pragma omp parallel for num_threads(cores) shared(answer)
    for(int i = toRow; i < MATRIX_SIZE;i++){
        for(int j = 0; j < MATRIX_SIZE;j++){
            //Convolucionar pixel
            int acumm = 0;
            int rowMult = i-rowCenter;
            int colMult = j-colCenter;
            int numberOperations = 0;
            for(int iPatern = 0; iPatern < rowSizePatern;iPatern++){
                for(int jPatern = 0; jPatern < colSizePatern;jPatern++){
                    int whatrow = rowMult+iPatern;
                    int whatcol = colMult+jPatern;
                    if((whatcol >= 0)&&(whatrow<MATRIX_SIZE)&&(whatcol<MATRIX_SIZE)){
                        acumm+=pat[iPatern][jPatern]*mat[whatrow][whatcol];
                        numberOperations++;
                    }
                }
            }
            answer[i][j] = acumm/numberOperations;
        }
    }

    //Lateral left part:
    // #pragma omp parallel for num_threads(cores) shared(answer)
    for(int i = fromRow; i < toRow;i++){
        for(int j = 0; j < fromCol;j++){
            //Convolucionar pixel
            int acumm = 0;
            int rowMult = i-rowCenter;
            int colMult = j-colCenter;
            int numberOperations = 0;
            for(int iPatern = 0; iPatern < rowSizePatern;iPatern++){
                for(int jPatern = 0; jPatern < colSizePatern;jPatern++){
                    int whatrow = rowMult+iPatern;
                    int whatcol = colMult+jPatern;
                    if(( whatrow >= 0)&&(whatcol >= 0)&&(whatrow<MATRIX_SIZE)){
                        acumm+=pat[iPatern][jPatern]*mat[whatrow][whatcol];
                        numberOperations++;
                    }
                }
            }
            answer[i][j] = acumm/numberOperations;
        }
    }
    //Right lateral part
    // #pragma omp parallel for num_threads(cores) shared(answer)
    for(int i = fromRow; i < toRow;i++){
        for(int j = toCol; j < MATRIX_SIZE;j++){
            //Convolucionar pixel
            int acumm = 0;
            int rowMult = i-rowCenter;
            int colMult = j-colCenter;
            int numberOperations = 0;
            for(int iPatern = 0; iPatern < rowSizePatern;iPatern++){
                for(int jPatern = 0; jPatern < colSizePatern;jPatern++){
                    int whatrow = rowMult+iPatern;
                    int whatcol = colMult+jPatern;
                    if(( whatrow >= 0)&&(whatrow<MATRIX_SIZE)&&(whatcol<MATRIX_SIZE)){
                        acumm+=pat[iPatern][jPatern]*mat[whatrow][whatcol];
                        numberOperations++;
                    }
                }
            }
            answer[i][j] = acumm/numberOperations;
        }
    }

}