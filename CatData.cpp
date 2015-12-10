#include<iostream>
#include<armadillo>
#include<string>
#include "K-Means V2.cpp"

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

    zombie2.save("Cat Pictures/rebuild", arma::pgm_binary);

}


/// Compress Image
//  Converts data matrix to the image
//

arma::Mat<double> compImage(const arma::Mat<double> dict, const std::string image, const int seg_size_rows, const int seg_size_cols)
{
/*
dict    = K-Means dictionary
image   = Path to full image
*/

    //Variables
    arma::Mat<double> comp = PicToMatrix(image, seg_size_rows, seg_size_cols); // Begin by matricies out vectors from the image

    double temp_distance {0}; // Tmep variable for calculating the distance a vector is from a centroid
    double min_distance {0};  // The norm of the closest centroid found so far
    arma::uword min_position {0};     // The position of the minimum vector


    // Compress the image by replacing each vector with the closest centroid

    for(int i = 0; i < comp.n_cols; i++)
    {
        // Initialise the centroid as the closest centroid
        temp_distance = arma::norm(comp.col(i) - dict.col(0));
        min_distance = temp_distance;
        min_position = 0;

        // Iterate through the rest of the centroids and find which ever is closer
        for(int j = 1; j < dict.n_cols; j++)
        {
            temp_distance = arma::norm(comp.col(i) - dict.col(j));

            // If a closer centroid is found, set it the current minimum
            if(temp_distance < min_distance)
            {
                min_distance = temp_distance;
                min_position = j;
            }
        }

        comp.col(i) = dict.col(min_position); // Compress the vector by replacing it with the closest centroid
    }

    comp.save("Cat Pictures/compImage", arma::pgm_binary);
    return comp;
}

//void distMeasure

int main()
{
    //Settings
    arma::arma_rng::set_seed_random();                                                 // Randomise the RNG seed

    //Constants
    const int seg_size_rows = 5; // Size of the individual segments extracted. Dimension = seg_size_rows * seg_size_cols
    const int seg_size_cols = 5;
    const int seg_dim = seg_size_cols * seg_size_rows;
    const int total_pictures = 4;

    const arma::uword iter = 20; // Number of Iterations
    const arma::uword nCentroids = 20; // Must be a common factor of points                     // Number of Centroids to look for

    const std::string imgPath = "Cat Pictures/cat.pgm"; // Path of the full image
    //Variables
    arma::Mat<double> compImg; // Matrix to store image to be compressed
    arma::Mat<double> trainingData; // Matrix to store training data
    arma::Mat<double> kmDict; //Dictionary buitl from training data

    arma::Mat<double> test = PicToMatrix(imgPath, 5, 5);
    rebMatrix(test,  5, 5);

    /*
    // Create Training data
    for(int i = 1; i <= total_pictures; i++)
    {
        std::string path = std::string("Cat Pictures/cat_test") + std::to_string(i) + ".pgm";
        trainingData = arma::join_horiz(trainingData, PicToMatrix(path, seg_size_rows, seg_size_cols));
    }

    trainingData.save("Cluster Data/clusterData", arma::pgm_binary);

    /// Compress Image
    kmDict = createDict("Cluster Data/clusterData", nCentroids, iter); // Create dictionary from training data

    compImg = compImage(kmDict, imgPath, seg_size_rows, seg_size_cols);
    rebMatrix(compImg, seg_size_rows, seg_size_cols);
    //rebMatrix(trainingData, seg_size_rows, seg_size_cols)

*/
}
