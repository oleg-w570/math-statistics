#include "triangle_experiment.h"

TriangleExperiment::TriangleExperiment(double sideLength) {
    triangle = EquilateralTriangle(sideLength);
}

void TriangleExperiment::RunExperiment(int n) {
    this->n = n;
    const double a = triangle.getSide();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis1(0.0, a + 1e-8);
    std::uniform_real_distribution<double> dis2(0.0, a * sqrt(3) * 0.5 + 1e-8);
    variationalSeries = std::vector<double>(n);
    int i = 0;

    while (i < n) {
        Point randPoint = Point(dis1(gen), dis2(gen));
        if (triangle.isInside(randPoint)) {
            const double distance = triangle.distanceToCenter(randPoint);
            variationalSeries[i] = distance;
            i++;
        }
    }

    std::sort(variationalSeries.begin(), variationalSeries.end());
}

void TriangleExperiment::calcSampleMean() {
    EmpiricalData.sampleMean = std::accumulate(variationalSeries.begin(), variationalSeries.end(), 0.0);
    EmpiricalData.sampleMean /= n;
}

void TriangleExperiment::calcSampleVariance() {
    EmpiricalData.sampleVariance = 0.0;
    for (const double &x: variationalSeries) {
        EmpiricalData.sampleVariance += (x - EmpiricalData.sampleMean) * (x - EmpiricalData.sampleMean);
    }
    EmpiricalData.sampleVariance /= n;
}

void TriangleExperiment::calcRange() {
    EmpiricalData.range = variationalSeries.back() - variationalSeries.front();
}

void TriangleExperiment::calcMedian() {
    int k = n / 2;
    if (n % 2 == 0)
        EmpiricalData.median = (variationalSeries[k - 1] + variationalSeries[k]) / 2.0;
    else
        EmpiricalData.median = variationalSeries[k];
}

void TriangleExperiment::calcDivMeasure()
{
    EmpiricalData.divmeasure = 0.0;
    for (int i = 0; i < n; i++) {
        const double& x = variationalSeries[i];
        const double F  = DistributionFunction(x);
        const double expr1 = (double)(i+1) / (double)n - F;
        const double expr2 = F - (double)i / (double)n;
        const double loc_max = std::max(expr1, expr2);
        EmpiricalData.divmeasure = std::max(EmpiricalData.divmeasure, loc_max);
    }
}

void TriangleExperiment::CalculateEmpiricalData() {
    calcSampleMean();
    calcSampleVariance();
    calcRange();
    calcMedian();
    calcDivMeasure();
}

void TriangleExperiment::CalculateHistogramData(const std::vector<double> &borders) {
    size_t size = borders.size();
    HistogramData.h = std::vector<double>(size - 1);
    HistogramData.z = std::vector<double>(size - 1);
    HistogramData.f = std::vector<double>(size - 1);

    int j = 0;
    for (int i = 0; i < size - 1; i++) {
        const double& left = borders[i];
        const double& right = borders[i+1];
        const double mid = (right + left) * 0.5;
        const double delta = right - left;
        int count = 0;
        while (j < n && variationalSeries[j] < right) {
            count++;
            j++;
        }
        HistogramData.h[i] = count / (n * delta);
        HistogramData.z[i] = mid;
        HistogramData.f[i] = ProbabilityDensity(mid);
    }
}

void TriangleExperiment::calcMean() {
    const double a = triangle.getSide();
    TheoreticData.mean = a * (12.0 * sqrt(3.0) + log(1351.0 + 780.0 * sqrt(3.0))) / 108.0;
}

void TriangleExperiment::calcVariance() {
    const double a = triangle.getSide();
    TheoreticData.variance = a * a / 12.0 - TheoreticData.mean * TheoreticData.mean;
}

void TriangleExperiment::CalculateTheoreticData() {
    calcMean();
    calcVariance();
}

double TriangleExperiment::DistributionFunction(double x) const {
    const double a = triangle.getSide();
    if (x > 0.0) {
        if (x > a * 0.5 / sqrt(3)) {
            if (x > a / sqrt(3)) {
                return 1.0;
            } else {
                return x * (a * sqrt(12.0 - a * a / x / x) + M_PI * x - 6.0 * x * asin(1.0 - a * a / (6.0 * x * x)))
                       / (sqrt(3.0) * a * a);
            }
        } else {
            return 4.0 * M_PI * x * x / (a * a * sqrt(3.0));
        }
    }
    return 0.0;
}

double TriangleExperiment::EmpiricalDistributionFunction(double x)
{
    auto xi = variationalSeries.begin();
    int count = 0;
    while (xi != variationalSeries.end() && *xi < x) {
        count++;
        xi++;
    }
    return (double)count / (double)n;
}

double TriangleExperiment::ProbabilityDensity(double x) const {
    const double a = triangle.getSide();
    if (x > 0.0) {
        if (x > a * 0.5 / sqrt(3.0)) {
            if (x > a / sqrt(3.0)) {
                return 0.0;
            } else {
                return 2.0 * x * (M_PI - 6.0 * asin(1.0 - a * a / (6.0 * x * x))) / (sqrt(3.0) * a * a);
            }
        } else {
            return 8.0 * M_PI * x / (sqrt(3.0) * a * a);
        }
    }
    return 0.0;
}

double TriangleExperiment::ProbabilityDensityX2(double x, int r) {
    if (x > 0.0) {
        const double rhf = r * 0.5;
        const double gamma = std::tgamma(rhf);

        return std::pow(2.0, -rhf) * std::pow(x, rhf - 1.0) * std::exp(-x * 0.5) / gamma;
    }

    return 0.0;
}

void TriangleExperiment::calcFR0() {
    const int k = HypothesisData.q.size();
    const int N = 1000;
    double sum = 0.0;
    for (int i = 1; i < N; i++) {
        const double expr1 = ProbabilityDensityX2(HypothesisData.R0 * (double)(i - 1) / N, k - 1);
        const double expr2 = ProbabilityDensityX2(HypothesisData.R0 * (double)i / N, k - 1);
        sum += (expr1 + expr2);
    }
    sum *= HypothesisData.R0 / (2.0 * N);
    HypothesisData.FR0 = 1.0 - sum;
}

void TriangleExperiment::CalculateHypothesisData(const std::vector<double> &brdrs) {
    const int k = brdrs.size() - 1;
    HypothesisData.R0 = 0.0;
    HypothesisData.q = std::vector<double>(k);

    int j = 0;
    for (int i = 0; i < k; i++) {
        const double &left = brdrs[i];
        const double &right = brdrs[i+1];
        const double qi = DistributionFunction(right) - DistributionFunction(left);
        int ni = 0;
        while (j < n && variationalSeries[j] < right) {
            ni++;
            j++;
        }

        HypothesisData.R0 += (ni - n * qi) * (ni - n * qi) / (n * qi);
        HypothesisData.q[i] = qi;
    }
    calcFR0();
}

bool TriangleExperiment::AcceptHypothesis(double alpha) const {
    return HypothesisData.FR0 > alpha;
}

