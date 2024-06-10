// include necessary libraries
#include "lbfgs.c"
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <ctime>

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
  Polygon(std::vector<Vector> vertices = {}) : vertices(vertices) {}
};

// New saver for PNG
#include <sstream>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
void save_frame(const std::vector<Polygon> &cells, std::string filename, int N, int frameid = 0) {
  int W = 1000, H = 1000;
  std::vector<unsigned char> image(W * H * 3, 255);
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < cells.size(); i++) {

    double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
    for (int j = 0; j < cells[i].vertices.size(); j++) {
      bminx = std::min(bminx, cells[i].vertices[j][0]);
      bminy = std::min(bminy, cells[i].vertices[j][1]);
      bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
      bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
    }
    bminx = std::min(W - 1., std::max(0., W * bminx));
    bminy = std::min(H - 1., std::max(0., H * bminy));
    bmaxx = std::max(W - 1., std::max(0., W * bmaxx));
    bmaxy = std::max(H - 1., std::max(0., H * bmaxy));

    for (int y = bminy; y < bmaxy; y++) {
      for (int x = bminx; x < bmaxx; x++) {
        int prevSign = 0;
        bool isInside = true;
        double mindistEdge = 1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
          double x0 = cells[i].vertices[j][0] * W;
          double y0 = cells[i].vertices[j][1] * H;
          double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
          double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
          double det = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0);
          int sign = det / std::abs(det);
          if (prevSign == 0)
            prevSign = sign;
          else if (sign == 0)
            sign = prevSign;
          else if (sign != prevSign) {
            isInside = false;
            break;
          }
          prevSign = sign;
          double edgeLen = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
          double distEdge = std::abs(det) / edgeLen;
          double dotp = (x - x0) * (x1 - x0) + (y - y0) * (y1 - y0);
          if (dotp < 0 || dotp > edgeLen * edgeLen)
            distEdge = 1E9;
          mindistEdge = std::min(mindistEdge, distEdge);
        }
        if (isInside) {
          if (i < N) { // the N first particles may represent fluid, displayed in blue
            image[((H - y - 1) * W + x) * 3] = 0;
            image[((H - y - 1) * W + x) * 3 + 1] = 0;
            image[((H - y - 1) * W + x) * 3 + 2] = 255;
          }
          if (mindistEdge <= 2) {
            image[((H - y - 1) * W + x) * 3] = 0;
            image[((H - y - 1) * W + x) * 3 + 1] = 0;
            image[((H - y - 1) * W + x) * 3 + 2] = 0;
          }
        }
      }
    }
  }
  std::ostringstream os;
  os << filename << frameid << ".png";
  stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
  FILE *f = fopen(filename.c_str(), "w+");
  fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
  for (size_t i = 0; i < polygons.size(); i++) {
    fprintf(f, "<g>\n");
    fprintf(f, "<polygon points = \"");
    for (size_t j = 0; j < polygons[i].vertices.size(); j++) {
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
  for (size_t i = 0; i < polygons.size(); i++) {
    fprintf(f, "<polygon points = \"");
    for (size_t j = 0; j < polygons[i].vertices.size(); j++) {
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
  Vector outwards_normal = Vector(line_v[1] - line_u[1], line_u[0] - line_v[0]);
  double t = dot(line_u - segment_a, outwards_normal) / dot(segment_b - segment_a, outwards_normal);
  Vector P = segment_a + t * (segment_b - segment_a);
  if (t >= 0 && t <= 1) {
    return Intersection(P, t, true);
  } else {
    return Intersection(P, t, false);
  }
}

bool towards_inside(Vector point, Vector line_u, Vector line_v) {
  Vector outwards_normal = Vector(line_v[1] - line_u[1], line_u[0] - line_v[0]);
  return dot(point - line_u, outwards_normal) <= 0;
}

Polygon clip_by_polygon(Polygon to_clip, Polygon clipper) {
  // Sutherland-Hodgman algorithm
  for (size_t i = 0; i < clipper.vertices.size(); i++) {
    // We clip on this line
    Polygon clipped;
    Vector curClipVertex = clipper.vertices[i];
    Vector prevClipVertex = clipper.vertices[i >= 1 ? (i - 1) : clipper.vertices.size() - 1];
    // We clip what's on the left of the line
    for (size_t j = 0; j < to_clip.vertices.size(); j++) {
      Vector curVertex = to_clip.vertices[j];
      Vector prevVertex = to_clip.vertices[j >= 1 ? (j - 1) : to_clip.vertices.size() - 1];
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

Intersection intersect_bissec_voronoi(Vector segment_a, Vector segment_b, Vector voronoi_center, Vector other_point, Vector M) {
  double t = dot(M - segment_a, voronoi_center - other_point) / dot(segment_b - segment_a, voronoi_center - other_point);
  Vector P = segment_a + t * (segment_b - segment_a);
  if (t >= 0 && t <= 1) {
    return Intersection(P, t, true);
  } else {
    return Intersection(P, t, false);
  }
}

bool towards_inside_voronoi(Vector point, Vector voronoi_center, Vector other_point, Vector M) {
  return dot(point - M, other_point - voronoi_center) < 0;
}

Polygon clip_by_bissec_voronoi(Polygon to_clip, Vector voronoi_center, Vector other_point, double weight_center = 1, double weight_other = 1) {
  // Modified Sutherland-Hodgman algorithm
  Polygon clipped;
  Vector M = (voronoi_center + other_point) / 2 + (((weight_center - weight_other) / (2 * (voronoi_center - other_point).norm2())) * (other_point - voronoi_center));
  for (size_t j = 0; j < to_clip.vertices.size(); j++) {
    Vector B = to_clip.vertices[j];
    Vector A = to_clip.vertices[j >= 1 ? (j - 1) : to_clip.vertices.size() - 1];
    Intersection intersect = intersect_bissec_voronoi(A, B, voronoi_center, other_point, M);
    if (towards_inside_voronoi(B, voronoi_center, other_point, M)) {
      if (!towards_inside_voronoi(A, voronoi_center, other_point, M)) {
        clipped.vertices.push_back(intersect.position);
      }
      clipped.vertices.push_back(B);
    } else if (towards_inside_voronoi(A, voronoi_center, other_point, M)) {
      clipped.vertices.push_back(intersect.position);
    }
  }
  return clipped;
}

std::vector<Polygon> triangulate(Polygon polygon) {
  std::vector<Polygon> triangles;
  for (size_t i = 1; i + 1 < polygon.vertices.size(); i++) {
    triangles.push_back(Polygon({polygon.vertices[0], polygon.vertices[i], polygon.vertices[i + 1]}));
  }
  return triangles;
}
double triangle_area(Polygon tri) {
  if (tri.vertices.size() != 3) {
    std::cout << "Error: triangle has more than 3 vertices" << std::endl;
    throw "Error: triangle has more than 3 vertices";
  }
  return 0.5 * cross(tri.vertices[1] - tri.vertices[0], tri.vertices[2] - tri.vertices[0]).norm();
}

class PointCloud {
public:
  std::vector<Vector> points;
  std::vector<double> weights;
  std::vector<double> lambda; // Wanted area of the voronoi cell
  std::vector<Vector> velocities;
  int N_particles;
  int nb_liquid = 0;
  double liquid_proportion = 0.1;

  PointCloud(std::vector<Vector> points) : points(points) {
    for (size_t i = 0; i < points.size(); i++) {
      weights.push_back(1);
    }
  }
  PointCloud(std::vector<Vector> points, std::vector<double> lambda_) : points(points) {
    for (size_t i = 0; i < points.size(); i++) {
      weights.push_back(1);
    }
    optimize_for_lambda(lambda_);
  }
  PointCloud(double liquid_proportion, int N_part) :  N_particles(N_part), liquid_proportion(liquid_proportion) {
    if (liquid_proportion <= 0 || liquid_proportion >= 1) {
      std::cout << "Error: proportion must be between 0 and 1" << std::endl;
      throw "Error: proportion must be between 0 and 1";
    }
    // Generate points randomly
    nb_liquid = std::floor(N_part * liquid_proportion);
    std::mt19937 generator(std::clock());
    std::uniform_real_distribution<double> pos(0, 1);
    for (int i = 0; i < N_particles; i++) {
      points.push_back(Vector(pos(generator), pos(generator)));
      weights.push_back(1);
    }
    // Lloyd iterations
    int n_lloyd = 10;
    if (points.size() > 1000) {
      n_lloyd = 5;
      if (points.size() > 10000) {
        n_lloyd = 3;
      }
    }
    if (points.size() <= 400) {
      n_lloyd = 30;
    }
    for (int i = 0; i < n_lloyd; ++i){
      lloyd();
      if (i % (n_lloyd/10) == 0) {
        std::cout << "Lloyd iteration " << (double)100.0*i/n_lloyd << "%" << std::endl;
      }
    }
    std::cout << "Lloyd done" << std::endl;

    for (size_t i = 0; i < nb_liquid; i++) {
      velocities.push_back(Vector(0, 0));
    }
  }

  void lloyd(){
    std::vector<Polygon> voron = generate_voronoi();
    std::vector<Vector> new_points;
    for (size_t i = 0; i < points.size(); i++) {
      Polygon voronoi_cell = voron[i];
      double area = 0;
      Vector centroid = Vector(0, 0);
      for (int j = 0; j < voronoi_cell.vertices.size()-1; j++) {
        Polygon tri = Polygon({points[i], voronoi_cell.vertices[j], voronoi_cell.vertices[j+1]});
        double T = triangle_area(tri);
        area += T;
        centroid = centroid + T * (tri.vertices[0] + tri.vertices[1] + tri.vertices[2]) / 3;
      }
      centroid = centroid / area;
      new_points.push_back(centroid + (centroid-points[i])*0.3); // Over relax
    }
    points = new_points;
  } 

  void select(bool (*selection) (Vector)){
    std::vector<Vector> new_points;
    nb_liquid = 0;
    for (size_t i = 0; i < points.size(); i++) {
      if (selection(points[i])){
        new_points.push_back(points[i]);
        ++nb_liquid;
      }
    }
    for (size_t i = 0; i < points.size(); i++) {
      if (!selection(points[i])){
        new_points.push_back(points[i]);
      }
    }
    points = new_points;
  }

  void optimize_for_lambda(std::vector<double> lambda_) {
    set_lambda(lambda_);
    std::vector<double> new_weights = optimize();
    for (size_t i = 0; i < points.size(); i++) {
      weights[i] = new_weights[i];
    }
  }

  void set_lambda(std::vector<double> lambda_) {
    double tot = 0;
    for (double l : lambda_) {
      tot += l;
    }
    if (lambda.size() == 0) {
      for (size_t i = 0; i < points.size(); i++) {
        lambda.push_back(lambda_[i] / tot);
      }
    } else if (lambda.size() == points.size()) {
      for (size_t i = 0; i < points.size(); i++) {
        lambda[i] = lambda_[i] / tot;
      }
    } else {
      std::cout << "Error: lambda size does not match the number of points" << std::endl;
      throw "Error: lambda size does not match the number of points";
    }
  }

  std::vector<double> optimize() {
    // We want to optimize the weights
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(points.size());
    lbfgs_parameter_t param;

    // Initialize the weights
    for (size_t i = 0; i < points.size(); i++) {
      x[i] = 0.1;
    }
    // Do the parameter thing
    lbfgs_parameter_init(&param);

    // Call lbfgs
    lbfgs(points.size(), x, &fx, evaluate, progress, this, &param);

    // Copy the optimized weights
    std::vector<double> new_weights;
    for (size_t i = 0; i < points.size(); i++) {
      new_weights.push_back(x[i]);
    }

    // Don't forget to free the memory
    lbfgs_free(x);

    // Return the optimized weights
    return new_weights;
  }

  void optimize_for_liquid() {
    // We want to optimize the weights
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(nb_liquid+1); // 1 for air
    lbfgs_parameter_t param;

    // Initialize the weights
    for (size_t i = 0; i < nb_liquid + 1; i++) {
      x[i] = 0.1;
    }
    // Do the parameter thing
    lbfgs_parameter_init(&param);

    // Call lbfgs
    lbfgs(nb_liquid, x, &fx, evaluate_liquid, progress, this, &param);

    // Copy the optimized weights
    std::vector<double> new_weights;
    for (size_t i = 0; i < nb_liquid; i++) {
      new_weights.push_back(x[i]);
    }
    for (size_t i = nb_liquid; i < points.size(); i++) {
      new_weights.push_back(x[nb_liquid]);
    }

    // Don't forget to free the memory
    lbfgs_free(x);

    // Write the optimized weights
    weights = new_weights;
  }

  std::vector<Polygon> generate_voronoi() {
    std::vector<Polygon> voronoi;
    for (size_t i = 0; i < points.size(); i++) {
      Polygon voronoi_cell = Polygon({Vector(0, 0), Vector(1, 0), Vector(1, 1), Vector(0, 1)});
      for (size_t j = 0; j < points.size(); j++) {
        if (i != j) {
          voronoi_cell = clip_by_bissec_voronoi(voronoi_cell, points[i], points[j], weights[i], weights[j]);
        }
      }
      voronoi.push_back(voronoi_cell);
    }
    return voronoi;
  }

  static lbfgsfloatval_t evaluate( // contains extra parameters
      void *instance,
      const lbfgsfloatval_t *x,
      lbfgsfloatval_t *g,
      const int n,
      const lbfgsfloatval_t step) {
    // x is W, g is where i put gradient, n is the dimension? step i don't think i care
    // return g(W) ?
    (void)step;
    PointCloud pc = *(PointCloud *)(instance);
    if ((size_t)n != pc.lambda.size()) {
      std::cout << "Error: n and lambda size do not match" << std::endl;
      return -1;
    }
    for (int i = 0; i < n; i++) {
      pc.weights[i] = x[i];
    }

    // We have the point cloud, now we can calculate the voronoi diagram
    std::vector<Polygon> voronoi = pc.generate_voronoi();
    // say the weights
    lbfgsfloatval_t fx = 0.0;
    // We make the actual computations
    for (int i = 0; i < n; i++) {
      double v_area = 0;
      // We add integral of ||x - yi||^2 for each i
      for (Polygon tri : triangulate(voronoi[i])) {
        double T = triangle_area(tri);
        v_area += T;
        for (int k = 0; k < 3; k++) {
          for (int l = k; l < 3; l++) {
            fx += (T / 6) * dot(tri.vertices[k] - pc.points[i], tri.vertices[l] - pc.points[i]);
          }
        }
      }
      // Finally
      g[i] = -(-v_area + pc.lambda[i]);
      fx += -v_area * x[i] + pc.lambda[i] * x[i];
    }
    return -fx;
  }

  static Polygon disk(Vector origin, double radius, int n = 50){
    // Take 30 sides
    std::vector<Vector> vertices;
    for (int i = 0; i < n; i++){
      double angle = 2*3.141592653289793236*i/n;
      vertices.push_back(Vector(origin[0] + radius*cos(angle), origin[1] + radius*sin(angle)));
    }
    return Polygon(vertices);
  }

  static lbfgsfloatval_t evaluate_liquid( // contains extra parameters
      void *instance,
      const lbfgsfloatval_t *x,
      lbfgsfloatval_t *g,
      const int n,
      const lbfgsfloatval_t step) {
    // x is W, g is where i put gradient, n is the dimension? step i don't think i care
    // return g(W) ?
    (void)step;
    PointCloud pc = *(PointCloud *)(instance);
    for (int i = 0; i < pc.nb_liquid; i++) {
      pc.weights[i] = x[i];
    }
    for (int i = pc.nb_liquid; i < pc.points.size(); i++) {
      pc.weights[i] = x[pc.nb_liquid];
    }

    // We have the point cloud, now we can calculate the voronoi diagram
    std::vector<Polygon> voronoi = pc.generate_voronoi();
    // say the weights
    lbfgsfloatval_t fx = 0.0;
    // We make the actual computations
    double total_liquid_area = 0;
    for (int i = 0; i < pc.nb_liquid; i++) {
      Polygon intersected_cell = clip_by_polygon(voronoi[i], disk(pc.points[i], sqrt(x[i] - x[pc.nb_liquid])));
      double v_area = 0;
      // We add integral of ||x - yi||^2 for each i
      for (Polygon tri : triangulate(intersected_cell)) {
        double T = triangle_area(tri);
        v_area += T;
        for (int k = 0; k < 3; k++) {
          for (int l = k; l < 3; l++) {
            fx += (T / 6) * dot(tri.vertices[k] - pc.points[i], tri.vertices[l] - pc.points[i]);
          }
        }
      }
      total_liquid_area += v_area;
      // Finally
      g[i] = -(pc.liquid_proportion/pc.nb_liquid -v_area);
      fx += -v_area * x[i] + pc.liquid_proportion/pc.nb_liquid * x[i];
    }
    double air_area = 1 - total_liquid_area;
    fx += x[pc.nb_liquid] * (1 - pc.liquid_proportion - air_area);
    g[pc.nb_liquid] = pc.liquid_proportion/pc.nb_liquid - air_area;
    return -fx;
  }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
  static int progress(
      void *instance,
      const lbfgsfloatval_t *x,
      const lbfgsfloatval_t *g,
      const lbfgsfloatval_t fx,
      const lbfgsfloatval_t xnorm,
      const lbfgsfloatval_t gnorm,
      const lbfgsfloatval_t step,
      int n,
      int k,
      int ls) {

    /*
    printf("Iteration %d:\n", k);
    printf("  fx = %f", fx);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");*/
    return 0;
  }
#pragma GCC diagnostic pop
};

bool in_circle(Vector point) {
  return (point - Vector(0.5, 0.5)).norm() < 0.3;
}

int main() {
  // Create particles cloud
  std::vector<Vector> particles = {Vector(0.3, 0.5), Vector(0.2, 0.2), Vector(0.8, 0.2), Vector(0.8, 0.8), Vector(0.2, 0.8)};
  PointCloud pc(0.4, 400);
  pc.select(in_circle);
  // Create voronoi diagram
  std::vector<Polygon> voronoi = pc.generate_voronoi();
  save_frame(voronoi, "voronoi", pc.nb_liquid, 0);
  std::cout << "Finished" << std::endl;

  return 0;
}
