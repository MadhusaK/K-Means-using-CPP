#include<iostream>
#include<armadillo>

/////////////////////////
/// Generate Cluster Data
//  Generate Clusters for testing
//
arma::Mat<double> genCluster(const int points, const int dim, const int nCluster)
{
/*
points          = Number of points
dim             = The number of dimenions of each point
nCluster        = Size of Clusters
*/

    int numClusters = points/nCluster;
    arma::Mat<double> dataPoints(dim, points, arma::fill::randn);
    arma::Mat<double> clustMean(dim, numClusters, arma::fill::randu);
    clustMean = clustMean*50 -25; // Change Range from R[0 to 1,...] to R[-25 to 25,...]


    if(points % nCluster !=0) { std::cerr << "The cluster size must be a common factor of the total points";}

    for(int i = 0; i < numClusters; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            dataPoints(j, arma::span(i*nCluster, i*nCluster + nCluster - 1)) +=  clustMean(j,i);
        }
    }

    return dataPoints;
}

/////////////////////////
/// Initialise Centroids
//  Generate initial centroids
//
arma::Mat<double> initCent(const int nCentroids, const arma::Mat<double>& dataPoints)
{
/*
nCentroids              = Number of centroids
dataPoints              = Point data
*/
    arma::Mat<double> centroids(dataPoints.n_rows, nCentroids, arma::fill::zeros); //Centroid matrix
    arma::Col<int> tempVec = arma::randi<arma::Col<int>>(nCentroids, arma::distr_param(0, dataPoints.n_cols-1)); // Select three random indexes from the data points

    for(int i = 0; i < nCentroids; i++)
    {
        // Add the data point at each index to the list of centroids
        centroids.col(i) = dataPoints.col(tempVec(i));
    }
    tempVec.print();
    centroids.print();

    return centroids;

}

/////////////////////////
/// Recalculate Centroids
//  Recalculate centroids according to nearest points
//
arma::Mat<double> recalCentroid(const arma::Mat<double>& dataPoints, arma::Mat<double>& centroids)
{
/*
centroids               = Centroid data
dataPoints              = Point data
*/
    arma::Mat<double> tempCentroids(centroids.n_rows, centroids.n_cols, arma::fill::zeros);
    arma::Col<int>    tempCount(centroids.n_cols, arma::fill::zeros);


    double temp_distance {0};
    double min_distance {0};
    int min_position {0};

    for(int i = 0; i < dataPoints.n_cols; i++)
    {
        // Set the first point as the initial minimum
        temp_distance = arma::norm(dataPoints.col(i) - centroids.col(0));
        min_distance = temp_distance;
        min_position = 0;

        for(int j=1; j < centroids.n_cols; j++)
        {
            // Iterate through the rest of the points to find the new minimum
            temp_distance = arma::norm(dataPoints.col(i) - centroids.col(j));
            if (temp_distance < min_distance)
            {
                min_distance = temp_distance;
                min_position = j;
            }
        }
        tempCentroids.col(min_position) += dataPoints.col(i);
        tempCount(min_position) += 1;
    }



    for(int i = 0; i < centroids.n_cols; i++)
    {
        if( tempCount(i) == 0) { tempCentroids.col(i) =  dataPoints.col(rand() % (dataPoints.n_cols-1));}
        else{ tempCentroids.col(i) /= tempCount(i);}
    }

    return tempCentroids;
}

//////////////////////////////
/// Create Centroid Dictionary
//  Creates a dictionary of centroids from a given set of cluster data
//
void createDict(std::string path, const arma::uword nCentroids, const arma:uword iter)
{
/*
path                = Path to cluster data
nCentroids          = Number of centroids
iter                = Number of iterations until centroids are expected to converge
*/

    arma::Mat<double> dataPoints; // matrix to hold cluster data
    dataPoints.load(path, arma::pgm_binary); // load PGM into the matrix
    arma::Mat<double> centroids =  initCent(nCentroids, dataPoints);  // Store the initial centroids

    for(int i = 0; i < itt; i++)
    {
        centroids = recalCentroid(dataPoints,centroids);
    }
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
    const int points = 100;                                                            // Number of points
    const int dim = 2;                                                                 // Dimensions of the random variable
    const int itt = 20;                                                                 // Number of Iterations
    const int nCentroids = 10; // Must be a common factor of points                     // Number of Centroids to look for
    const int nCluster = points/nCentroids;                                            // Cluster size

    //Variables
    //arma::Mat<double> dataPoints = genCluster(points,dim, nCluster);             // Generate data
    arma::Mat<double> dataPoints;
    dataPoints.load("Cluster Data/clusterData", arma::pgm_binary); //Load clusterdata

    arma::Mat<double> centroids =  initCent(nCentroids, dataPoints);  // Store the initial centroids

    for(int i = 0; i < itt; i++)
    {
        centroids = recalCentroid(dataPoints,centroids);
    }

    dataPoints.save("Cluster Data/data", arma::csv_ascii);
    centroids.save("Cluster Data/centroids", arma::csv_ascii);

}





