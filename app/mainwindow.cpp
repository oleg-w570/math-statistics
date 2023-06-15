#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
//    ui->tableOfValues->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->tableONumCharacteristics->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->tableHistogram->verticalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->graphF->addGraph();
    ui->graphF->graph(0)->setName(tr("Теоретическая функция распределения"));
    ui->graphF->addGraph();
    ui->graphF->graph(1)->setName(tr("Выборочная функция распределения"));
    ui->graphF->graph(1)->setPen(QPen(Qt::black));
    ui->graphF->xAxis->setLabel(tr("x"));
    ui->graphF->yAxis->setLabel(tr("y"));
    ui->graphF->legend->setVisible(true);

    ui->graphH->addGraph();
    ui->graphH->graph(0)->setName(tr("Гистограмма"));
    ui->graphH->xAxis->setLabel(tr("x"));
    ui->graphH->yAxis->setRange(0.0, 1.0);
    ui->graphH->yAxis->setLabel(tr("y"));
    ui->graphH->legend->setVisible(true);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButtonRun_clicked()
{
    const int n = ui->lineN->text().toInt();
    E.RunExperiment(n);
    const double& side = E.triangle.getSide();
    const std::vector<double>& values = E.variationalSeries;
    E.CalculateEmpiricalData();
    const double& mean = E.TheoreticData.mean;
    const double& sampleMean = E.EmpiricalData.sampleMean;
    const double& variance = E.TheoreticData.variance;
    const double& sampleVariance = E.EmpiricalData.sampleVariance;
    const double& median = E.EmpiricalData.median;
    const double& range = E.EmpiricalData.range;
    const double& divmeasure = E.EmpiricalData.divmeasure;

    ui->tableOfValues->setColumnCount(n);
    for (int i = 0; i < n; i++) {
        ui->tableOfValues->setItem(0, i, new QTableWidgetItem(QString::number(values[i])));
    }
    ui->tableOfValues->resizeColumnsToContents();

    ui->tableONumCharacteristics->setItem(0, 0, new QTableWidgetItem(QString::number(mean)));
    ui->tableONumCharacteristics->setItem(0, 1, new QTableWidgetItem(QString::number(sampleMean)));
    ui->tableONumCharacteristics->setItem(0, 2, new QTableWidgetItem(QString::number(std::abs(mean - sampleMean))));
    ui->tableONumCharacteristics->setItem(0, 3, new QTableWidgetItem(QString::number(variance)));
    ui->tableONumCharacteristics->setItem(0, 4, new QTableWidgetItem(QString::number(sampleVariance)));
    ui->tableONumCharacteristics->setItem(0, 5, new QTableWidgetItem(QString::number(std::abs(variance - sampleVariance))));
    ui->tableONumCharacteristics->setItem(0, 6, new QTableWidgetItem(QString::number(median)));
    ui->tableONumCharacteristics->setItem(0, 7, new QTableWidgetItem(QString::number(range)));
    ui->tableONumCharacteristics->setItem(0, 8, new QTableWidgetItem(QString::number(divmeasure)));

    QVector<double> x = {-1.0, -0.5, 0.0};
    QVector<double> F = {0.0, 0.0, 0.0};
    const double maxDist = side / sqrt(3.0);
    const double step = 0.01;
    double xi = step;

    while (xi < maxDist) {
        x.push_back(xi);
        F.push_back(E.EmpiricalDistributionFunction(xi));
        xi += step;
    }
    x.push_back(maxDist);
    x.push_back(maxDist+0.5);
    x.push_back(maxDist+1.0);
    F.push_back(1.0);
    F.push_back(1.0);
    F.push_back(1.0);

    ui->graphF->graph(1)->setData(x, F);
    ui->graphF->graph(1)->rescaleAxes(true);
    ui->graphF->replot();
}


void MainWindow::on_lineSide_textEdited(const QString &arg1)
{
    ui->graphF->graph(1)->data().data()->clear();
    const double side = arg1.toDouble();
    const double maxDist = side / sqrt(3.0);
    E = TriangleExperiment(side);
    E.CalculateTheoreticData();
    QVector<double> x = {-1.0, -0.5, 0.0};
    QVector<double> F = {0.0, 0.0, 0.0};
    const double step = 0.01;
    double xi = step;

    while (xi < maxDist) {
        x.push_back(xi);
        F.push_back(E.DistributionFunction(xi));
        xi += step;
    }
    x.push_back(maxDist);
    x.push_back(maxDist+0.5);
    x.push_back(maxDist+1.0);
    F.push_back(1.0);
    F.push_back(1.0);
    F.push_back(1.0);

    ui->graphF->graph(0)->setData(x, F);
    ui->graphF->graph(0)->rescaleAxes();    
    ui->graphF->replot();
}

void MainWindow::on_pushButtonBorders_clicked()
{
    const std::vector<double> borders = InputBorders();
    if (borders.empty())
        return;
    E.CalculateHistogramData(borders);
    const std::vector<double>& z = E.HistogramData.z;
    const std::vector<double>& f = E.HistogramData.f;
    const std::vector<double>& h = E.HistogramData.h;

    size_t size = borders.size()-1;
    ui->tableHistogram->setColumnCount(size);
    if (size <= 10)
        ui->tableHistogram->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    for (size_t i = 0; i < size; i++) {
        ui->tableHistogram->setItem(0, i, new QTableWidgetItem(QString::number(z[i])));
        ui->tableHistogram->setItem(1, i, new QTableWidgetItem(QString::number(f[i])));
        ui->tableHistogram->setItem(2, i, new QTableWidgetItem(QString::number(h[i])));
    }

    QVector<double> x;
    QVector<double> y;
    for (size_t i = 0; i < size; i++) {
        const double& left = borders[i];
        const double& right = borders[i+1];
        const double& height = h[i];
        x.push_back(left);
        x.push_back(right);
        y.push_back(height);
        y.push_back(height);
    }

    ui->graphH->graph(0)->setData(x, y);
//    ui->graphH->rescaleAxes(true);
    ui->graphH->replot();
}

std::vector<double> MainWindow::InputBorders()
{
    bool ok = false;
    const int k = QInputDialog::getInt(this, tr("Задайте промежутки"), tr("Количество промежутков"), 1, 1, 25, 1, &ok);
    if (ok) {
        const double maxDist = E.triangle.getSide() / sqrt(3.0);
        std::vector<double> borders(k+1);
        for (int i = 1; i < k; i++) {
            const double border = QInputDialog::getDouble(this, tr("Задайте промежутки"),
                                                          tr("Правая граница ") + QString::number(i) + tr(" промежутка"),
                                                          1, 0.001, maxDist, 3, &ok);
            if (ok)
                borders[i] = border;
            else
                return {};
        }
        borders[k] = maxDist;
        return borders;
    }

    return {};
}

void MainWindow::on_pushButtonHypothesis_clicked()
{
    const std::vector<double> borders = InputBorders();
    bordersHypothesis = borders;
    if (borders.empty())
        return;

    E.CalculateHypothesisData(borders);
    const std::vector<double>& q = E.HypothesisData.q;
    const size_t size = q.size();
    if (size <= 10)
        ui->tableHypothesis->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->tableHypothesis->setColumnCount(size);
    for (size_t i = 0; i < size; i++) {
        ui->tableHypothesis->setHorizontalHeaderItem(i, new QTableWidgetItem("q" + QString::number(i+1)));
        ui->tableHypothesis->setItem(0, i, new QTableWidgetItem(QString::number(q[i])));
    }
    ui->tableHypothesis->resizeColumnsToContents();
    QString FR0 = "F(R0) = " + QString::number(E.HypothesisData.FR0);
    ui->labelFR0->setText(FR0);
}


void MainWindow::on_pushButtonCheck_clicked()
{
    bool ok = false;
    const double alpha = QInputDialog::getDouble(this,
                                                 tr("Проверка гипотезы"),
                                                 tr("Введите уровень значимости"),
                                                 0.0, 0.0, 1.0, 5, &ok);
    if (ok) {
        const bool check = E.AcceptHypothesis(alpha);
        QString message = check ? "Гипотеза принята" : "Гипотеза отвергнута";
        QMessageBox::information(this, tr("Проверка гипотезы"), message);
    }

}

