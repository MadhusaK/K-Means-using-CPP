#include<iostream>
#include<armadillo>
#include<string>

using arma::uword; // 64 bit unsigned integer

/// PGM to N-dimensional Matrix conversion
//  Converts PGM's to N dimensional training data
//
arma::Mat<double> PicToMatrix(const std::string path, const int seg_size_rows, const int seg_size_cols)
{
/*

path = path of the image
seg_size_*   = Respective dimensions of each segment
*/

    arma::Mat<double> imgMatrix; // Image Matrix
    imgMatrix.load(path, arma::pgm_binary);

    const int seg_dim = seg_size_cols*seg_size_rows;
    const int imgMatrixSize = imgMatrix.n_rows*imgMatrix.n_cols;

    arma::Mat<double> convMatrix(seg_dim, imgMatrixSize/seg_dim, arma::fill::zeros); // Conversion Matrix with N-dimensional vectors


    int col_counter = 0; //Column Counter for loop

    for(int i = 0; i <  (imgMatrix.n_rows / seg_size_rows); i++)
    {
        for(int j = 0 ; j < (imgMatrix.n_cols / seg_size_cols); j++)
        {
        //Extract matrix segments and convert them into vectorts
        convMatrix.col(col_counter) = arma::vectorise(imgMatrix.submat(i*seg_size_rows, j*seg_size_cols, i*seg_size_rows + seg_size_rows - 1, j*seg_size_cols + seg_size_cols - 1));
        col_counter++;

        }
    }

    convMatrix.save("Cat Pictures/clustData", arma::pgm_binary);

    return convMatrix;

}

/// Rebuild Image
//  Converts data matrix to the image
//

arma::Mat<double> rebMatrix(const arma::Mat<double>& clusterData, const int seg_size_rows,const int seg_size_cols)
{
    arma::Mat<double> temp;
    arma::Mat<double> zombie;
    arma::Mat<double> zombie2;

    int const imgRow = 13;
    int const imgCol = 14;


    std::cout << clusterData.n_cols;
    for(int i = 0; i < clusterData.n_cols; i++)
    {
        temp = arma::reshape(clusterData.col(i), seg_size_rows, seg_size_cols);
        // temp += randu(size(temp)) * 128.0;
        zombie = arma::join_rows(zombie,temp);

    }

    for(int i = 0; i < imgRow; i++)
    {
        zombie2 = arma::join_cols(zombie2,zombie.cols(i*imgCol*seg_size_cols, i*imgCol*seg_size_cols + imgCol*seg_size_cols -1));
    }

    zombie2 = arma::flipud(zombie2);
    zombie2.save("Cat Pictures/rebuild", arma::pgm_binary);

}


int main()
{
    /////////////
    /// Variables
    //  Initialise variables and constants
    //
    //Settings
    arma::arma_rng::set_seed_random();                                                 // Randomise the RNG seed

    //Constants
    const int seg_size_rows = 10; // Size of the individual segments extracted. Dimension = seg_size_rows * seg_size_cols
    const int seg_size_cols = 10;
    const int seg_dim = seg_size_cols * seg_size_rows;
    const int total_pictures = 3;

    //Variables
    arma::Mat<double> trainingData;

    for(int i = 1; i <= total_pictures; i++)
    {
        std::string path = std::string("Cat Pictures/cat_test") + std::to_string(i) + ".pgm";
        trainingData = arma::join_horiz(trainingData, PicToMatrix(path, seg_size_rows, seg_size_cols));
    }

    trainingData.save("Cluster Data/clusterData", arma::pgm_binary);
    //rebMatrix(trainingData, seg_size_rows, seg_size_cols);




}
