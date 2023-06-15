#pragma once
#define _USE_MATH_DEFINES
#include "triangle.h"
#include <cmath>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>


class TriangleExperiment {
public:
    EquilateralTriangle triangle;
    std::vector<double> variationalSeries;
    int n;
    struct {
        double sampleMean;
        double sampleVariance;
        double range;
        double median;
        double divmeasure;
    } EmpiricalData;
    struct {
        double mean;
        double variance;
    } TheoreticData;
    struct {
        std::vector<double> h;
        std::vector<double> z;
        std::vector<double> f;
    } HistogramData;
    struct {
        std::vector<double> q;
        double R0;
        double FR0;
    } HypothesisData;

protected:
    void calcSampleMean();
    void calcSampleVariance();
    void calcRange();
    void calcMedian();
    void calcDivMeasure();
    void calcMean();
    void calcVariance();
    void calcFR0();
    static double ProbabilityDensityX2(double x, int r);
public:
    explicit TriangleExperiment(double sideLength = 10.0);
    void RunExperiment(int n);
    double DistributionFunction(double x) const;
    double EmpiricalDistributionFunction(double x);
    double ProbabilityDensity(double x) const;
    void CalculateEmpiricalData();
    void CalculateTheoreticData();
    void CalculateHistogramData(const std::vector<double>& brdrs);
    void CalculateHypothesisData(const std::vector<double>& brdrs);
    bool AcceptHypothesis(double alpha) const;
};
