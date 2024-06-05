// include necessary libraries
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

class Vector {
public:
  explicit Vector(double x = 0, double y = 0, double z = 0) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }
  double norm2() const {
    return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
  }
  double norm() const {
    return sqrt(norm2());
  }
  void normalize() {
    double n = norm();
    data[0] /= n;
    data[1] /= n;
    data[2] /= n;
  }
  double operator[](int i) const { return data[i]; };
  double &operator[](int i) { return data[i]; };
  double data[3];
};

Vector operator+(const Vector &a, const Vector &b) {
  return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b) {
  return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector &a) {
  return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const double a, const Vector &b) {
  return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &a, const double b) {
  return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator/(const Vector &a, const double b) {
  return Vector(a[0] / b, a[1] / b, a[2] / b);
}
bool operator==(const Vector &a, const Vector &b) {
  return (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]);
}
std::ostream &operator<<(std::ostream &os, const Vector &obj) {
  os << "(" << obj[0] << ", " << obj[1] << ", " << obj[2] << ")";
  return os;
}

double dot(const Vector &a, const Vector &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector &a, const Vector &b) {
  return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

// if the Polygon class name conflicts with a class in wingdi.h on Windows, use a namespace or change the name
class Polygon {
public:
  // Needs to be described in the trigonometric order
  std::vector<Vector> vertices;
};

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
  FILE *f = fopen(filename.c_str(), "w+");
  fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
  for (int i = 0; i < polygons.size(); i++) {
    fprintf(f, "<g>\n");
    fprintf(f, "<polygon points = \"");
    for (int j = 0; j < polygons[i].vertices.size(); j++) {
      fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
    }
    fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
    fprintf(f, "</g>\n");
  }
  fprintf(f, "</svg>\n");
  fclose(f);
}

// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Polygon> &polygons, std::string filename, int frameid, int nbframes) {
  FILE *f;
  if (frameid == 0) {
    f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    fprintf(f, "<g>\n");
  } else {
    f = fopen(filename.c_str(), "a+");
  }
  fprintf(f, "<g>\n");
  for (int i = 0; i < polygons.size(); i++) {
    fprintf(f, "<polygon points = \"");
    for (int j = 0; j < polygons[i].vertices.size(); j++) {
      fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
    }
    fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
  }
  fprintf(f, "<animate\n");
  fprintf(f, "    id = \"frame%u\"\n", frameid);
  fprintf(f, "    attributeName = \"display\"\n");
  fprintf(f, "    values = \"");
  for (int j = 0; j < nbframes; j++) {
    if (frameid == j) {
      fprintf(f, "inline");
    } else {
      fprintf(f, "none");
    }
    fprintf(f, ";");
  }
  fprintf(f, "none\"\n    keyTimes = \"");
  for (int j = 0; j < nbframes; j++) {
    fprintf(f, "%2.3f", j / (double)(nbframes));
    fprintf(f, ";");
  }
  fprintf(f, "1\"\n   dur = \"5s\"\n");
  fprintf(f, "    begin = \"0s\"\n");
  fprintf(f, "    repeatCount = \"indefinite\"/>\n");
  fprintf(f, "</g>\n");
  if (frameid == nbframes - 1) {
    fprintf(f, "</g>\n");
    fprintf(f, "</svg>\n");
  }
  fclose(f);
}

struct Intersection {
  Vector position;
  double t;
  bool intersected;
  Intersection(Vector position, double t, bool intersected) : position(position), t(t), intersected(intersected) {}
};

Intersection intersect_line(Vector segment_a, Vector segment_b, Vector line_u, Vector line_v) {
  Vector outwards_normal = Vector(line_v[1] - line_u[1], line_u[0] - line_v[0]); // TODO verify
  double t = dot(line_u - segment_a, outwards_normal) / dot(segment_b - segment_a, outwards_normal);
  Vector P = segment_a + t * (segment_b - segment_a);
  if (t >= 0 && t <= 1) {
    return Intersection(P, t, true);
  } else {
    return Intersection(P, t, false);
  }
}

bool towards_inside(Vector point, Vector line_u, Vector line_v) {
  Vector outwards_normal = Vector(line_v[1] - line_u[1], line_u[0] - line_v[0]); // TODO verify
  return dot(point - line_u, outwards_normal) <= 0;
}

Polygon clip_by_polygon(Polygon to_clip, Polygon clipper) {
  // Sutherland-Hodgman algorithm
  for (int i = 0; i < clipper.vertices.size(); i++) {
    // We clip on this line
    Polygon clipped;
    Vector curClipVertex = clipper.vertices[i];
    Vector prevClipVertex = clipper.vertices[(i - 1) >= 0 ? (i-1):  clipper.vertices.size()-1];
    // We clip what's on the left of the line
    for (int j = 0; j < to_clip.vertices.size(); j++) {
      Vector curVertex = to_clip.vertices[j];
      Vector prevVertex = to_clip.vertices[(j - 1) >= 0 ? (j-1):  to_clip.vertices.size()-1];
      Intersection intersect = intersect_line(prevVertex, curVertex, prevClipVertex, curClipVertex);
      if (towards_inside(curVertex, prevClipVertex, curClipVertex)) {
        if (!towards_inside(prevVertex, prevClipVertex, curClipVertex)) {
          clipped.vertices.push_back(intersect.position);
        }
        clipped.vertices.push_back(curVertex);
      } else if (towards_inside(prevVertex, prevClipVertex, curClipVertex)) {
        clipped.vertices.push_back(intersect.position);
      }
    }
    to_clip = clipped;
  }
  return to_clip;
}

int main() {
  // create a polygon in trigonometric order
  Polygon poly;
  poly.vertices.push_back(Vector(0.5, 0.1));
  poly.vertices.push_back(Vector(0.9, 0.9));
  poly.vertices.push_back(Vector(0.1, 0.9));
  poly.vertices.push_back(Vector(0.1, 0.5));
  Polygon clipping_hexagon;
  // define a centered hexagon
  for (int i = 0; i < 6; i++) {
    clipping_hexagon.vertices.push_back(Vector(0.6 + 0.4 * cos(i * 2 * M_PI / 6), 0.5 + 0.4 * sin(i * 2 * M_PI / 6)));
  }
  // save the polygon to an svg file
  Polygon polynice = clip_by_polygon(poly, clipping_hexagon);
  save_svg({polynice}, "polygon.svg");
  return 0;
}
