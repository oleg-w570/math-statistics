#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <vector>
#include "triangle_experiment.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_pushButtonRun_clicked();

    void on_lineSide_textEdited(const QString &arg1);

    void on_pushButtonBorders_clicked();

    void on_pushButtonHypothesis_clicked();

    void on_pushButtonCheck_clicked();

private:
    Ui::MainWindow *ui;
    TriangleExperiment E;
    std::vector<double> bordersHypothesis;
    std::vector<double> InputBorders();
};
#endif // MAINWINDOW_H
