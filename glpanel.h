#ifndef GLPANEL_H
#define GLPANEL_H

#include <QGLWidget>

class Controller;

class GLPanel : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLPanel(QWidget *parent = 0);
    void setController(Controller *cont);

signals:

public slots:
    virtual void resizeGL(int w, int h);
    virtual void paintGL();
    virtual void mousePressEvent(QMouseEvent *me);

private:
    Controller *cont_;
};

#endif // GLPANEL_H
