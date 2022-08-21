/**
  \file App.h

  The G3D 10.00 default starter app is configured for OpenGL 4.1 and
  relatively recent GPUs.
 */
#pragma once
#include <G3D/G3D.h>

/** Application framework. */
class App : public GApp {
public:
    
    App(const GApp::Settings& settings = GApp::Settings());

    virtual void onInit() override;

    virtual void onGraphics3D(RenderDevice* rd, Array<shared_ptr<Surface> >& surfaceArray) override;

protected:

    shared_ptr<Texture>         m_environmentMap;
    Array<shared_ptr<Texture> > m_veniceTexture;

    const String m_dirCurrent;
    const String m_dirTextures;
    const String m_dirShaders;
};
