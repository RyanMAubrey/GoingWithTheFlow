#include <QApplication>
#include <QSurfaceFormat>

#include "bodydebugwidget.h"
#include "fluiddebugwidget.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QSurfaceFormat fmt;
    fmt.setVersion(2, 1);
    fmt.setProfile(QSurfaceFormat::NoProfile);
    QSurfaceFormat::setDefaultFormat(fmt);

    //BodyDebugWidget w;
    //FluidDebugWidget w;
    //w.resize(800, 600);
    //w.show();

    return app.exec();
}
