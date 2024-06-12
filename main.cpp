// include necessary libraries
#include "lbfgs.c"
#include <cmath>
#include <ctime>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>
#define N_THREADS 16
#define N_THREADS_EVALUATE 4

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

void save_line(size_t i, std::vector<unsigned char> &image, const std::vector<Polygon> &cells, int W, int H, int N) {
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

void save_lines(size_t i0, size_t size, std::vector<unsigned char> &image, const std::vector<Polygon> &cells, int W, int H, int N) {
  for (size_t i = i0; i < i0 + size; i++) {
    save_line(i, image, cells, W, H, N);
  }
}

void save_frame(const std::vector<Polygon> &cells, std::string filename, int N, int frameid = 0) {
  int W = 1000, H = 1000;
  std::vector<unsigned char> image(W * H * 3, 255);
  std::vector<std::thread> threads(N_THREADS - 1);
  size_t block_size = N / N_THREADS;
  for (int n_thread = 0; n_thread < N_THREADS - 1; n_thread++) {
    threads[n_thread] = std::thread(&save_lines, n_thread * block_size, block_size, std::ref(image), cells, W, H, N);
  }
  save_lines((N_THREADS - 1) * block_size, N - (N_THREADS - 1) * block_size, std::ref(image), cells, W, H, N);
  for (int n_thread = 0; n_thread < N_THREADS - 1; n_thread++) {
    threads[n_thread].join();
  }
  std::ostringstream os;
  os << filename << frameid << ".png";
  stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
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
  std::vector<Vector> velocities;
  int nb_liquid = 0;
  double liquid_proportion = 0.1;
  double air_weight = 0.1;

  PointCloud(double liquid_proportion, int N_part, std::vector<Vector> (*repartition)(int) = nullptr) : nb_liquid(N_part), liquid_proportion(liquid_proportion) {
    if (liquid_proportion <= 0 || liquid_proportion >= 1) {
      std::cout << "Error: proportion must be between 0 and 1" << std::endl;
      throw "Error: proportion must be between 0 and 1";
    }

    if (repartition == nullptr) {
      for (int i = 0; i < nb_liquid; i++) {
        weights.push_back(1);
        velocities.push_back(Vector(0, 0));
      }
      // Generate points randomly
      std::mt19937 generator(std::clock());
      std::uniform_real_distribution<double> pos(0, 1);
      for (int i = 0; i < nb_liquid; i++) {
        points.push_back(Vector(pos(generator), pos(generator)));
        weights.push_back(1);
      }
      // Lloyd iterations to make the points more regular
      int n_lloyd = 15;
      if (points.size() <= 400) {
        n_lloyd = 30;
      }
      for (int i = 0; i < n_lloyd; ++i) {
        lloyd();
        if (n_lloyd < 10 || i % (n_lloyd / 10) == 0) {
          std::cout << "Lloyd iteration " << (double)100.0 * i / n_lloyd << "%" << std::endl;
        }
      }
      std::cout << "Lloyd done" << std::endl;
    } else {
      points = repartition(nb_liquid);
      nb_liquid = points.size();
      for (int i = 0; i < nb_liquid; i++) {
        weights.push_back(1);
        velocities.push_back(Vector(0, 0));
      }
    }
  }

  void lloyd(size_t start_moving = 0) {
    std::vector<Polygon> voron = generate_voronoi();
    std::vector<Vector> new_points;
    for (size_t i = 0; i < points.size(); i++) {
      Polygon voronoi_cell = voron[i];
      double area = 0;
      Vector centroid = Vector(0, 0);
      for (int j = 0; j < voronoi_cell.vertices.size() - 1; j++) {
        Polygon tri = Polygon({points[i], voronoi_cell.vertices[j], voronoi_cell.vertices[j + 1]});
        double T = triangle_area(tri);
        area += T;
        centroid = centroid + T * (tri.vertices[0] + tri.vertices[1] + tri.vertices[2]) / 3;
      }
      centroid = centroid / area;
      if (i < start_moving) {
        new_points.push_back(points[i]);
      } else {
        new_points.push_back(centroid + (centroid - points[i]) * 0.3); // Over relax
      }
    }
    points = new_points;
  }

  void animate(int nb_frames, double epsilon, double dt, double mass) {
    for (int i = 0; i < nb_frames; ++i) {
      // Optimization and save time is the long part
      auto start = std::chrono::high_resolution_clock::now();
      optimize_for_liquid();
      auto end = std::chrono::high_resolution_clock::now();
      std::vector<Polygon> voronoi = generate_voronoi_liquid();
      auto start2 = std::chrono::high_resolution_clock::now();
      save_frame(voronoi, "animation", nb_liquid, i);
      auto end2 = std::chrono::high_resolution_clock::now();
      std::cout << "Frame " << i << std::endl;
      time_step_post_optimization(epsilon, dt, mass, voronoi);
      std::cout << "Optimization time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n" << "Save time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count() << "ms\n" << std::endl;
      // TODO parallelize more
    }
  }

  void time_step_post_optimization(double epsilon, double dt, double mass, std::vector<Polygon> voronoi) {
    for (int i = 0; i < nb_liquid; ++i) {
      Polygon drop = voronoi[i];
      Vector centroid = Vector(0, 0);
      double area = 0;
      for (int j = 0; j < drop.vertices.size() - 1; j++) {
        Polygon tri = Polygon({points[i], drop.vertices[j], drop.vertices[j + 1]});
        double T = triangle_area(tri);
        area += T;
        centroid = centroid + T * (tri.vertices[0] + tri.vertices[1] + tri.vertices[2]) / 3;
      }
      centroid = centroid / area;
      // Actual computations
      Vector fi_spring = (centroid - points[i]) / pow(epsilon, 2);
      Vector fi = fi_spring + mass * Vector(0, -9.81); // Gravity
      velocities[i] = velocities[i] + dt * fi / mass;
      points[i] = points[i] + dt * velocities[i];
      // Bounce with dampened velocity by 2
      double velocity_dampening[4] = {0.25, 0.25, 0.1, 0.2}; // Bounces less strongly on ground
      for (int j = 0; j < 2; j++) {
        if (points[i][j] < 0) {
          points[i][j] = -points[i][j] * velocity_dampening[2*j];
          velocities[i][j] = -velocities[i][j] * velocity_dampening[2*j];
        } else if (points[i][j] > 1) {
          points[i][j] = 1 - (points[i][j] - 1) * velocity_dampening[2*j + 1];
          velocities[i][j] = -velocities[i][j] * velocity_dampening[2*j + 1];
        }
      }
    }
  }

  void optimize_for_liquid() {
    // We want to optimize the weights
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(nb_liquid + 1); // 1 for air
    lbfgs_parameter_t param;

    // Initialize the weights
    for (size_t i = 0; i < nb_liquid; i++) {
      x[i] = weights[i];
    }
    x[nb_liquid] = air_weight;
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
    std::vector<Polygon> voronoi(points.size());
    std::vector<std::thread> threads(N_THREADS - 1);
    size_t block_size = points.size() / N_THREADS;
    for (int n_thread = 0; n_thread < N_THREADS - 1; n_thread++) {
      threads[n_thread] = std::thread([this, n_thread, block_size, &voronoi]() {
        for (size_t i = n_thread * block_size; i < (n_thread + 1) * block_size; i++) {
          voronoi[i] = Polygon({Vector(0, 0), Vector(1, 0), Vector(1, 1), Vector(0, 1)});
          for (size_t j = 0; j < points.size(); j++) {
            if (i != j) {
              voronoi[i] = clip_by_bissec_voronoi(voronoi[i], points[i], points[j], weights[i], weights[j]);
            }
          }
        }
      });
    }
    for (size_t i = (N_THREADS - 1) * block_size; i < points.size(); i++) {
      voronoi[i] = Polygon({Vector(0, 0), Vector(1, 0), Vector(1, 1), Vector(0, 1)});
      for (size_t j = 0; j < points.size(); j++) {
        if (i != j) {
          voronoi[i] = clip_by_bissec_voronoi(voronoi[i], points[i], points[j], weights[i], weights[j]);
        }
      }
    }
    for (int n_thread = 0; n_thread < N_THREADS - 1; n_thread++) {
      threads[n_thread].join();
    }
    return voronoi;
  }

  std::vector<Polygon> generate_voronoi_liquid() {
    std::vector<Polygon> voronoi(points.size());
    std::vector<std::thread> threads(N_THREADS - 1);
    size_t block_size = points.size() / N_THREADS;
    for (int n_thread = 0; n_thread < N_THREADS - 1; n_thread++) {
      threads[n_thread] = std::thread([this, n_thread, block_size, &voronoi]() {
        for (size_t i = n_thread * block_size; i < (n_thread + 1) * block_size; i++) {
          voronoi[i] = clip_by_polygon(Polygon({Vector(0, 0), Vector(1, 0), Vector(1, 1), Vector(0, 1)}), disk(points[i], sqrt(weights[i] - air_weight)));
          for (size_t j = 0; j < points.size(); j++) {
            if (i != j) {
              voronoi[i] = clip_by_bissec_voronoi(voronoi[i], points[i], points[j], weights[i], weights[j]);
            }
          }
        }
      });
    }
    for (size_t i = (N_THREADS - 1) * block_size; i < points.size(); i++) {
      voronoi[i] = clip_by_polygon(Polygon({Vector(0, 0), Vector(1, 0), Vector(1, 1), Vector(0, 1)}), disk(points[i], sqrt(weights[i] - air_weight)));
      for (size_t j = 0; j < points.size(); j++) {
        if (i != j) {
          voronoi[i] = clip_by_bissec_voronoi(voronoi[i], points[i], points[j], weights[i], weights[j]);
        }
      }
    }
    for (int n_thread = 0; n_thread < N_THREADS - 1; n_thread++) {
      threads[n_thread].join();
    }
    return voronoi;
  }

  static Polygon disk(Vector origin, double radius, int n = 30) {
    // Take 30 sides
    std::vector<Vector> vertices;
    for (int i = 0; i < n; i++) {
      double angle = 2 * 3.141592653289793236 * i / n;
      vertices.push_back(Vector(origin[0] + radius * cos(angle), origin[1] + radius * sin(angle)));
    }
    return Polygon(vertices);
  }

  static void evaluate_liquid_thread(const std::vector<Polygon> &voronoi, std::vector<lbfgsfloatval_t> &fx, std::vector<double> &liquid_area, lbfgsfloatval_t *g, const PointCloud &pc, int i_start, int i_end, int n_thread) {
    // Store TOTAL result in n_thread for fx and liquid_area and PARTIAL result in g at i
    fx[n_thread] = 0;
    liquid_area[n_thread] = 0;
    for (int i = i_start; i < i_end; i++) {
      double v_area = 0;
      // We add integral of ||x - yi||^2 for each i
      for (Polygon tri : triangulate(voronoi[i])) {
        double T = triangle_area(tri);
        v_area += T;
        for (int k = 0; k < 3; k++) {
          for (int l = k; l < 3; l++) {
            fx[n_thread] += (T / 6) * dot(tri.vertices[k] - pc.points[i], tri.vertices[l] - pc.points[i]);
          }
        }
      }
      // Finally
      liquid_area[n_thread] += v_area;
      g[i] = -(pc.liquid_proportion / pc.nb_liquid - v_area);
      fx[n_thread] += -v_area * pc.weights[i] + pc.liquid_proportion / pc.nb_liquid * pc.weights[i];
    }
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
    pc.air_weight = x[pc.nb_liquid];
    // We have the point cloud, now we can calculate the voronoi diagram
    std::vector<Polygon> voronoi = pc.generate_voronoi_liquid();
    // Initialize the vectors
    std::vector<std::thread> threads(N_THREADS_EVALUATE-1);
    int block_size = pc.nb_liquid / N_THREADS_EVALUATE;
    std::vector<lbfgsfloatval_t> fx_vec(N_THREADS_EVALUATE);
    std::vector<double> liquid_area_vec(N_THREADS_EVALUATE);
    // We make the actual computations
    for (int i = 0; i < N_THREADS_EVALUATE-1; i++) {
      threads[i] = std::thread(&evaluate_liquid_thread, voronoi, std::ref(fx_vec), std::ref(liquid_area_vec), std::ref(g), pc, i * block_size, (i + 1) * block_size, i);
    }
    evaluate_liquid_thread(voronoi, fx_vec, liquid_area_vec, g, pc, (N_THREADS_EVALUATE-1) * block_size, pc.nb_liquid, N_THREADS_EVALUATE-1);
    // Join the threads and merge the results
    lbfgsfloatval_t fx = fx_vec[N_THREADS_EVALUATE-1];
    double total_liquid_area = liquid_area_vec[N_THREADS_EVALUATE-1];
    for (int i = 0; i < N_THREADS_EVALUATE-1; i++) {
      threads[i].join();
      fx += fx_vec[i];
      total_liquid_area += liquid_area_vec[i];
    }
    double air_area = 1 - total_liquid_area;
    fx += x[pc.nb_liquid] * (1 - pc.liquid_proportion - air_area);
    g[pc.nb_liquid] = -(pc.liquid_proportion / pc.nb_liquid - air_area);
    return -fx;
  }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
  static int progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls) {
    return 0;
  }
#pragma GCC diagnostic pop
};

std::vector<Vector> semi_donut(int N) {
  std::vector<Vector> points;
  // 0.2 air, 0.15 liquid, 0.3 air, 0.15 liquid, 0.2 air
  // semicircumference: 0.5 to 1 (a bit less), say 0.75, width: 0.15
  // n_per_radius*n_per_width = N && n_per_radius / 0.75 = n_per_width / 0.15
  // n_per_radius = 0.75 * n_per_width / 0.15 = 5 * n_per_width
  int n_per_width = std::floor(sqrt(N / 5));
  int n_per_radius = std::floor(sqrt(5*N));
  while (n_per_width*n_per_radius > N){n_per_radius--;}
  for (int i = 0; i < n_per_width; i++) {
    double radius = 0.15 + 0.25 * i / (n_per_width - 1);
    int per_radius = n_per_radius;
    if (i == n_per_width - 1) {
      per_radius = N - n_per_radius * (n_per_width - 1);
    }
    for (int j = 0; j < per_radius; j++) {
      double angle = 3.141592653589793238 * (((double) j / (double) per_radius));
      points.push_back(Vector(0.5 - radius * cos(angle), 0.5 - radius * sin(angle)));
    }
  }
  return points;
}

int main() {
  // Create particles cloud
  auto start = std::chrono::high_resolution_clock::now();
  PointCloud pc(0.2, 1000);

  // Create animation
  pc.animate(100, 0.004, 0.01, 250);
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Finished in " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << "s" << std::endl;

  return 0;
}
