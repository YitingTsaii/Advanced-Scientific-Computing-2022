#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

static void readCsv(const std::string& filename,
                    std::vector<double>& x,
                    std::vector<double>& y)
{
    std::ifstream fin(filename);

    // skip header
    {
        std::string s;
        std::getline(fin, s);
    }

    double a, b;
    char delim;

    while(fin >> a >> delim >> b)
    {
        x.push_back(a);
        y.push_back(b);
    }
}

class PieceWiseLinearCDF
{
public:
    PieceWiseLinearCDF(std::vector<double> x, std::vector<double> y)
        : x_(std::move(x))
        , y_(std::move(y))
    {}

    PieceWiseLinearCDF(const std::string& filename)
    {
        readCsv(filename, x_, y_);
    }

    double generateSample(std::mt19937& gen) const
    {
        // TODO write your code here
        // get a sample from uniform(0,1) and perform inversion method
        std::uniform_real_distribution<double> dist(0.0,1.0);
        auto y_sample = dist(gen);
        double x_sample = 0;
        if (y_sample == 0){
            x_sample = this->getMinX();
        } else if (y_sample == 1){
            x_sample = this->getMaxX();
        } else {
            for(int i = 0; i < this->x_.size(); i++){
                if (y_sample >= this->y_[i] && y_sample < this->y_[i+1]){
                    x_sample = ((y_sample - this->y_[i])/(this->y_[i+1] - this->y_[i])) * (this->x_[i+1] - this->x_[i]) + this->x_[i];
                    break;
                }
            }
        }
        return x_sample;
    }

    double getMinX() const {return x_.front();}
    double getMaxX() const {return x_.back();}

//private:
    std::vector<double> x_;
    std::vector<double> y_;
};


int main(int argc, char **argv)
{
    // Load the data
    PieceWiseLinearCDF cdf("cdf.csv");

    // histogram quantities
    const int nbins = 100;

    std::vector<double> histogramLoc(nbins); // center of each bin
    std::vector<int> histogramCounts(nbins); // number of samples per bin

    const double a = cdf.getMinX(); // lowest bound of the histogram
    const double b = cdf.getMaxX(); // highest bound of the histogram
    const double h = (b-a) / nbins; // bin size

    for (int i = 0; i < nbins; ++i)
        histogramLoc[i] = a + (i+0.5) * h;


    // collect samples
    const int nsamples = 1'000'000;
    std::vector<double> samples(nsamples);
    std::mt19937 gen(3456789);

    // after this loop, 'samples' will be a vector storing all the x samples
    for (auto& x : samples) // loop x over each element of samples, and modify them
    {
        x = cdf.generateSample(gen); // each one of the x sample

        // TODO fill in the histogram counts
        double bin_lower = 0;
        double bin_upper = 0;
        if (x == a) {
                histogramCounts[0] += 1;
        } else {
            for (int i = 0; i < nbins; ++i){
                bin_lower = a + i * h;
                bin_upper = a + (i+1) * h;
                if (x > bin_lower && x <= bin_upper) {
                    histogramCounts[i] += 1;
                    break;
                }
            }
        }
    }

    // TODO compute and printout the empirical median and the analytical median
    // empirical median (obtained from the samples)
    sort(samples.begin(), samples.end()); // sort the samples in place
    double empirical_median = 0;
    if (samples.size() % 2 == 0){ // size is even
        empirical_median = (samples[samples.size() / 2 - 1] + samples[samples.size() / 2]) / 2;
    } else {
        empirical_median = samples[samples.size() / 2];
    }
    std::cout << "empirical median from samples = " << empirical_median << "\n";

    // analytical median (median of the CDF)
    // perform inversion method in main (not in the calss)
    double y_sample = 0.5;
    double x_sample = 0;
    if (y_sample == 0){
        x_sample = a;
    } else if (y_sample == 1){
        x_sample = b;
    } else {
        for(int i = 0; i < cdf.x_.size(); i++){
            if (y_sample >= cdf.y_[i] && y_sample < cdf.y_[i+1]){
                x_sample = ((y_sample - cdf.y_[i])/(cdf.y_[i+1] - cdf.y_[i])) * (cdf.x_[i+1] - cdf.x_[i]) + cdf.x_[i];
                break;
            }
        }
    }
    double analytical_median = x_sample;
    std::cout << "analytical median from the CDF = " << analytical_median << "\n";
    

    // write the histogram to a csv file
    std::ofstream fout("histogram.csv");
    fout << "x,count" << std::endl;
    for (int i = 0; i < nbins; ++i)
        fout << histogramLoc[i] << ',' << histogramCounts[i] << std::endl;

    return 0;
}