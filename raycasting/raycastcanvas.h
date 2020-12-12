#pragma once

#include <functional>
#include <vector>

#include <QtMath>

#include <QOpenGLWidget>
#include <QOpenGLExtraFunctions>
#include <QOpenGLShaderProgram>

#include "mesh.h"
#include "raycastvolume.h"
#include "trackball.h"
//#include "vtkvolume.h"
#include "../data_importer.h"

/*!
 * \brief Class for a raycasting canvas widget.
 */
class RayCastCanvas : public QOpenGLWidget, protected QOpenGLExtraFunctions
{
    Q_OBJECT
public:
    explicit RayCastCanvas(QWidget *parent = nullptr);
    ~RayCastCanvas();

    void setStepLength(const GLfloat step_length) {
        m_stepLength = step_length;
        update();
    }
    void initVolume(DataImporter *_data_importer) {
        data_importer = _data_importer;
        setVolume();
    }
    void setVolume(long frame4display = 0) {
        if (!data_importer)
        {
            throw std::runtime_error("data_importer has not been initialized.");
        }
        if (!data_importer->p_vmin){// if max min value not defined
            data_importer->updateminmaxvalues();
        }
        double p_min = data_importer->p_vmin[0];
        double p_max = data_importer->p_vmax[0];
        long sx = data_importer->image4d->getXDim();
        long sy = data_importer->image4d->getYDim();
        long sz = data_importer->image4d->getZDim();
        long sc = data_importer->image4d->getCDim(); // for gray image stacked in channel, sc is the time indeed

        if (frame4display>=data_importer->image4d->getTDim()){
            throw std::runtime_error("data to show is not gray.");
        }
        long offsets = frame4display*sx*sy*sz;
        unsigned char *datahead = (unsigned char *)(data_importer->image4d->getRawDataAtChannel(0));

        if (!m_raycasting_volume)
            m_raycasting_volume = new RayCastVolume(this);
        m_raycasting_volume->transfer_volume(datahead + offsets, p_min, p_max, sx,
                                             sy, sz, sc);
        update();
    }

    void setThreshold(const double threshold) {
        auto range = m_raycasting_volume ? getRange() : std::pair<double, double>{0.0, 1.0};
        m_threshold = threshold / (range.second - range.first);
        update();
    }

    void setMode(const QString& mode) {
        m_active_mode = mode;
        update();
    }

    void setBackground(const QColor& colour) {
        m_background = colour;
        update();
    }

    std::vector<QString> getModes(void) {
        std::vector<QString> modes;
        for (const auto& [key, val] : m_modes) {
            modes.push_back(key);
        }
        return modes;
    }

    QColor getBackground(void) {
        return m_background;
    }

    std::pair<double, double> getRange(void) {
        return m_raycasting_volume->range();
    }

    DataImporter* getDataImporter(){return data_importer;}

    void handleKeyPressEvent(QKeyEvent * event); //for hook to MainWindow
    void handleKeyReleaseEvent(QKeyEvent * event); //for hook to MainWindow

signals:
    void changeVolumeTimePoint(int);
public slots:
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void wheelEvent(QWheelEvent * event);
    virtual void setVolumeTimePoint(int t);
    virtual void setLightPositionZero();
    virtual void setContrast(int relative_contrast/*[-100:100]*/);
    virtual void keyPressEvent(QKeyEvent *e){handleKeyPressEvent(e);}
    virtual void keyReleaseEvent(QKeyEvent *e){handleKeyReleaseEvent(e);}
protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

private:
    DataImporter *data_importer;
    QMatrix4x4 m_viewMatrix;
    QMatrix4x4 m_modelViewProjectionMatrix;
    QMatrix3x3 m_normalMatrix;
    // m_fov is the maximum vertical angle of cammera, it defines how large the fov will be
    const GLfloat m_fov = 30.0f;                                          /*!< Vertical field of view. */
    const GLfloat m_focalLength = 1.0 / qTan(M_PI / 180.0 * m_fov / 2.0); /*!< Focal length. */
    GLfloat m_aspectRatio;                                                /*!< width / height */

    QVector2D m_viewportSize;
    QVector3D m_rayOrigin; /*!< Camera position in model space coordinates. */

    QVector3D m_lightPosition {3.0, 0.0, 3.0};    /*!< In camera coordinates. */
    QVector3D m_diffuseMaterial {1.0, 1.0, 1.0};  /*!< Material colour. */
    GLfloat m_stepLength;                         /*!< Step length for ray march. */
    // no use
    GLfloat m_threshold = 50;                     /*!< Isosurface intensity threshold. */

    QColor m_background;                          /*!< Viewport background colour. */

    GLfloat m_gamma = 1.0f;//2.2 /*!< Gamma correction parameter. */

    RayCastVolume *m_raycasting_volume;
    //QPainter *canvas_painter;
    std::map<QString, QOpenGLShaderProgram*> m_shaders;
    std::map<QString, std::function<void(void)>> m_modes;
    QString m_active_mode;

    TrackBall m_trackBall {};       /*!< Trackball holding the model rotation. */
    TrackBall m_scene_trackBall {}; /*!< Trackball holding the scene rotation. */

    //center shift
    QPointF *centerShift = new QPointF(0.0, 0.0);
    GLint m_distExp = -200;

    GLuint scaled_width();
    GLuint scaled_height();

    void raycasting(const QString& shader);

    QPointF pixel_pos_to_view_pos(const QPointF& p);
    void create_noise(void);
    void add_shader(const QString& name, const QString& vector, const QString& fragment);
public:
    // rendering text
    void renderText(double x, double y, double z, QString text);
    inline GLint project(GLdouble objx, GLdouble objy, GLdouble objz,
                        const GLdouble model[16], const GLdouble proj[16],
                        const GLint viewport[4],
                        GLdouble * winx, GLdouble * winy, GLdouble * winz);
    inline void transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4]);
};
