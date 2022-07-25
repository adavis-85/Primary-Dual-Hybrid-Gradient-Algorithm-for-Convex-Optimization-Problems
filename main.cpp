/******************************************************************************
    Primary Dual Hybrid Gradient Algorithm.

  Part of convex optimization where there are two objectives.  This example lends
  itself to game theory where there are two players and accoring to a given
  score matrix, the optimal path to maximize each players total scores are chosen.
*******************************************************************************/

#include <iostream>
#include <bits/stdc++.h>
#include <vector>
#include <iterator>
#include <numeric>
#include <cmath>

using namespace std;
using std::vector;

vector<double> proj(vector<double> y,int ysize){

    //Vector "u" to contain the vector "y" sorted in descending order.
    vector<double> u(ysize);
    //Copying "y" into "u".
    u=y;
    //Initializing "K" to track largest entry in "y" that
    int K = 0;
    double tau=0.0;
    vector<double> x_n(ysize);
    vector<double> check_sum(ysize);
    sort(u.begin(),u.end(),greater<double>());
    double sum_amnt=0;

    //Iterating sums to find spot most affected by change.
    for (int i = 0; i < ysize; i++) {
      sum_amnt=0;
          for (int j=0;j < i+1; j++){
            sum_amnt=sum_amnt+u[j];
          }
      check_sum[i]=(sum_amnt-1)/(i+1);
    }

    //"K" is chosen.
    for (int i=0;i<ysize;i++){
        if (check_sum[i]<u[i]){
            K=i+1;
        }else{
            break;
        }
    }

    //Tau is initialized.
    for (int i=0;i<K;i++){
        tau=tau+u[i];
    }

    tau=(tau-1)/K;

    //Vector x_n is the probability vector for each projection.
    //Sum of entries must add up to 1 i.e. 100%.
    for (int i=0;i<ysize;i++){
        x_n[i]=max((y[i]-tau),0.0);
    }

    return x_n;
}



vector<double> proj(vector<double> y,int ysize);

int main()
{
    //Initial score matrix.
    double x[2][2] ={{10,3},{2,9}};

   int rows = sizeof(x)/sizeof(x[0]);
   int cols = sizeof(x[0])/sizeof(x[0][0]);

    double new_x[rows][rows];
    double transpose[cols][rows];

    //Transpose of matrix to compute the norm for rates of sigma and tau.
    for (int i=0;i<rows;i++){
        for (int j=0;j<cols;j++){
            transpose[j][i]=x[i][j];
        }
    }

    //Computing matrix norm.  norm(A*A') .
    double mult_sum;
    int c2 = sizeof(transpose[0])/sizeof(transpose[0][0]);
    for(int i=0; i<rows; ++i)
      for(int j=0; j<c2; ++j)
      for(int k=0; k<cols; ++k) {
         new_x[i][j]=x[i][k]*transpose[k][j];
      }

    int sumSq = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            sumSq += sqrt(new_x[i][j]);
        }
    }

    // Return the square root of
    // the sum of squares
    double x_L = sumSq;

    ////////////////////////////////end of norm of matrix

    ///////Initializing variables for probabilities/simplex
    vector<double> X_n(rows,.1);
    vector<double> Y_n(rows,.1);
    vector<double> X_bar_zero(rows,.1);
    vector<double> prob_x(rows,.1);
    vector<double> prob_y(rows,.1);
    vector<double> vec_mat_y(rows,.1);
    vector<double> vec_mat_x(rows,.1);
    vector<double> prev_y(rows,10);
    vector<double> prev_x(rows,10);

    int count=0;
    double sigma=1;
    double theta=1/sqrt(x_L);
    double tow=1.0;
    float tol=1e-8;
    double check_tol=1;
    double prev_tol=10000;
    double for_x;
    double for_y;

    while (prev_tol>tol){

        //Vector matrix multiplication per iteration for probabilities of column player.
        for(int i=0;i<rows;i++){
            for (int j=0;j<cols;j++){
                vec_mat_y[i]=prev_y[i]+(X_bar_zero[i]*x[i][j]*sigma);
            }
        }
        //Projecting onto simplex.
        prob_y=proj(vec_mat_y,rows);

        ///Vector matrix multiplication per iteration for probabilities of row player.
        for (int i=0; i<rows;i++){
            for (int j=0;j<cols;j++){
               vec_mat_x[i]=prev_x[i]-(prob_y[i]*transpose[j][i]*tow);
         }
        }

        //Projecting onto simplex.
        prob_x=proj(vec_mat_x,rows);

        for (int i=0;i<rows;i++){
            X_bar_zero[i]=prob_x[i]+(theta*(prob_x[i]-prev_x[i])*tow);
        }

       //Computing tolerance each iteration as a condition of convergence.
       for_x=0;
       for_y=0;

       for (int i=0;i<rows;i++){
            for_x=for_x+pow(prob_x[i]-prev_x[i],2);
            for_y=for_y+pow(prob_y[i]-prev_y[i],2);
        }

        if (for_x==0 && for_y==0){
           check_tol==0;
        }else{
                check_tol=(sqrt(for_x)/(2*sigma))+(sqrt(for_y)/(2*tow));
        }
        //Adjusting tolerance.//////////////////////////////////////
        prev_tol=prev_tol-check_tol;

        //Initializing for previous iteration for projecting at each iteration.
        prev_x=prob_x;
        prev_y=prob_y;

        //Visualizing progress each 1000 iterations.
        int mod_statement;
        mod_statement=count%1000;
        if (mod_statement==0){
            cout << "Probabilities for y" << endl;
            for (int i=0;i<rows;i++){
                cout << prob_y[i] << " ";
            }
            cout << endl;
            cout << "Probabilities for x" << endl;
            for (int i=0;i<rows;i++){
                cout << prob_x[i] << " ";
            }
            cout << endl;
        }

        //If tolerance is met.
         if (prev_tol<tol){
             cout << "Converged successfully after " << count << " iterations." << endl;
             break;
         }

        count++;
}
    return 0;
}
