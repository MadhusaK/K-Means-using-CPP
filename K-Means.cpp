#include<iostream>
#include<armadillo>


/////////////////////////
/// Generate Testing Data
//  Data used for testing
//
arma::Mat<double> genData(const int points, const int dim)
{
/*
points          = Number of points
dim             = The number of dimenions of each point
*/
    //Generate data
    arma::Mat<double> dataPoints(dim, points, arma::fill::randu);
    //Change range
    dataPoints = dataPoints*20-10; //Change range from R[[0,1]x[0,1]]-> R[[-10,10]x-[-10,10]]

    return dataPoints;
}

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
    clustMean = clustMean*50 -25;


    if(points % nCluster !=0) { std::cerr << "The cluster size must be a common factor of the total points";}

    for(int i = 0; i < numClusters; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            dataPoints(j, arma::span(i*nCluster, i*nCluster + nCluster - 1)) +=  clustMean(j,i);
        }
    }
    //Generate data
    //Change range
    return dataPoints;
}


////////////////
/// Group Points
//  Group points by closest centroid
//
void groupPoints(int iteration, const int points, const int nCentroids, const arma::Mat<double>& dataPoints, const arma::Mat<double>& centroids, arma::Mat<int>& groupMatrix)
{
/*
iteration       = Current iteration of the grouping function -- 0 to itt
points          = Number of points
nCentroids      = Number of centroids
dataPoints      = Data
centroids       = Centroids in matrix form -- Each Column represents the coordinates of one centroid
groupMatrix     = matrix to store the group of each point
*/
    double temp_distance {0};
    double min_distance {0};
    int min_position {0};

    for(int i = 0; i < points; i++)
    {
        // Set the first point as the initial minimum
        temp_distance = arma::norm(dataPoints.col(i) - centroids.col(0));
        min_distance = temp_distance;
        min_position = 0;
        for(int j = 1; j < nCentroids; j++)
        {
            // Iterate through the rest of the points to find the new minimum
            temp_distance = arma::norm(dataPoints.col(i) - centroids.col(j));
            if (temp_distance < min_distance)
            {
                min_distance = temp_distance;
                min_position = j;
            }
        }
        groupMatrix(i,iteration) = min_position;
    }
}

///////////////////////
/// Calculate Centroids
//  Group points by closest centroid
//
void calCentroid(int iteration, const int points, const int dim, const int nCentroids, const arma::Mat<double>& dataPoints, const arma::Mat<int>& groupMatrix, arma::Mat<double>& centroids)
{
/*
iteration       = Current iteration of the grouping function -- 0 to itt
points          = Number of points
dim             = The number of dimenions of each point
nCentroids      = Number of centroids
dataPoints      = Data
groupMatrix     = matrix to store the group of each point
centroids       = Centroids in matrix form -- Each Column represents the coordinates of one centroid
*/
    arma::Mat<double> tempCentroids(dim, nCentroids , arma::fill::zeros);      // Temporary matrix to store new centroids
    arma::Col<int>  tempCount(nCentroids , arma::fill::zeros);                    // Temporary vector to store the amount of points allocated to each centroid

    //Calculate new centroids
    for(int i = 0; i < points; i++)
    {
        tempCentroids.col(groupMatrix(i, iteration)) += dataPoints.col(i);
        tempCount(groupMatrix(i,iteration)) += 1;
    }


    //Normalise the centroids
    for(int i =0; i < nCentroids; i++)
    {
        tempCentroids.col(i) /= tempCount(i);
    }


    centroids = tempCentroids;

}

int main()
{
    /////////////
    /// Variables
    //  Initialise variables and constants
    //
    //Settings
    arma::arma_rng::set_seed_random();                              // Randomise the RNG seed

    //Constants
    const int points = 1000;                                          // Number of points
    const int dim = 3;                                            // Dimensions of the random variable
    const int itt = 5;                                              // Number of Iterations
    const int nCluster = 100; // Must be a common factor of points
    const int nCentroids = points/nCluster;                                       // Number of Centroids to look for

    //Variables
    arma::Mat<double> dataPoints = genData(points,dim);             // Generate data
    arma::Mat<int> groupMatrix(points, itt, arma::fill::zeros);     // Store the grouping for each point at each iteration
    arma::Mat<double> centroids = dataPoints.cols(0,nCentroids-1);  // Store the initial centroids


    std::cout << '\n' << "Centroids";
    centroids.print();
    centroids.save("Centroids", arma::csv_ascii);

    //////////////////////////
    /// Group Data by Centroid
    //  Group the points by the closest centroid
    groupPoints(0, points, nCentroids, dataPoints, centroids, groupMatrix); // Initial grouping


    arma::Mat<double> testData = genCluster(points,dim, nCluster);
    testData.save("Test", arma::csv_ascii);
    ////////////////////////////////////////////
    /// Calculate New Centroids and Regroup data
    //  Calculate the new centroids and regroup the data

    for(int i = 1 ; i < itt; i ++)
    {
        calCentroid(i-1, points,dim, nCentroids, dataPoints, groupMatrix, centroids);
        groupPoints(i, points, nCentroids, dataPoints, centroids, groupMatrix); // Initial grouping
    }

    std::cout << '\n' << "Data";
    dataPoints.print();
    dataPoints.save("data", arma::csv_ascii);

    std::cout << '\n' << "Allocatioms" << '\n';
    groupMatrix.print();
    groupMatrix.save("tags", arma::csv_ascii);

    std::cout << '\n' << "New Centroids";
    centroids.print();

}
