/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012  Daniel Robinson
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#ifndef SRVF_RENDER_H
#define SRVF_RENDER_H 1


namespace srvf
{


struct Color
{
  Color(float r, float g, float b)
   : red(r), green(g), blue(b)
  { }

  float red;
  float green;
  float blue;
};


enum DrawingMode
{
  POINTS=0,
  LINES=1,
  TUBES=2
};
  

/**
 * Abstract base class for renderers.
 */
class Renderer
{
public:

  virtual void begin(DrawingMode mode) = 0;
  virtual void vertex(double x, double y) = 0;
  virtual void vertex(double x, double y, double z) = 0;
  virtual void end() = 0;
  
  virtual void set_color(Color c) = 0;
  virtual void set_thickness(double t) = 0;

  virtual void translate(double x, double y) = 0;
  virtual void translate(double x, double y, double z) = 0;
  virtual void scale(double sfx, double sfy) = 0;
  virtual void scale(double sfx, double sfy, double sfz) = 0;
  virtual void rotate(double angle) = 0;
  virtual void rotate(double ax, double ay, double az, double angle) = 0;

protected:
  
  Renderer() { }
};

/**
 * OpenGL renderer.
 */
class OpenGlRenderer : public Renderer
{
public:
  
  OpenGlRenderer();

  virtual void begin(DrawingMode mode);
  virtual void vertex(double x, double y);
  virtual void vertex(double x, double y, double z);
  virtual void end();
  
  virtual void set_color(Color c);
  virtual void set_thickness(double t);

  virtual void translate(double x, double y);
  virtual void translate(double x, double y, double z);
  virtual void scale(double sfx, double sfy);
  virtual void scale(double sfx, double sfy, double sfz);
  virtual void rotate(double angle);
  virtual void rotate(double ax, double ay, double az, double angle);

private:
  DrawingMode mode_;
  double thickness_;
};


} // namespace srvf

#endif // SRVF_RENDER_H
