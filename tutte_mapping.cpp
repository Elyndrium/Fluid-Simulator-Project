#include <iostream>
#include <vector>
#include <thread>
#include <cmath>
#define N_THREADS 8

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
  size_t block_size = cells.size() / N_THREADS;
  for (int n_thread = 0; n_thread < N_THREADS - 1; n_thread++) {
    threads[n_thread] = std::thread(&save_lines, n_thread * block_size, block_size, std::ref(image), cells, W, H, N);
  }
  save_lines((N_THREADS - 1) * block_size, cells.size() - (N_THREADS - 1) * block_size, std::ref(image), cells, W, H, N);
  for (int n_thread = 0; n_thread < N_THREADS - 1; n_thread++) {
    threads[n_thread].join();
  }
  std::ostringstream os;
  os << filename << frameid << ".png";
  stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

struct Graph {
  std::vector<Vector> exterior_vertices; // Ordered to make the boundary
  std::vector<Vector> interior_vertices;
  //std::vector<std::vector<Vector>> exterior_edges; // Edges from exterior vertices
  std::vector<std::vector<size_t>> interior_edges; // Indices of vertices in "ev+iv" next to interior vertices
  Graph(std::vector<Vector> ev, std::vector<Vector> iv, std::vector<std::vector<size_t>> ie) : exterior_vertices(ev), interior_vertices(iv), interior_edges(ie) {}
};

Graph generate_mesh() {
  std::vector<Vector> ev = {Vector(0, 0, 0), Vector(1, 0, 0), Vector(1, 0.5, 0), Vector(0.5, 1, 0), Vector(0, 1, 0)};
  std::vector<Vector> iv = {Vector(0.5, 0.5, 1)};
  //std::vector<std::vector<Vector>> ee = {{ev[1], ev[4], iv[0]}, {ev[0], ev[2], iv[0]}, {ev[1], ev[3], iv[0]}, {ev[2], ev[4], iv[0]}, {ev[3], ev[0], iv[0]}};
  std::vector<std::vector<size_t>> ie = {{0,1,2,3,4}}; // Indices in "ev + iv"
  return Graph(ev, iv, ie);
}

std::vector<Polygon> Polygonize(const std::vector<Vector> &vecs) {
  std::vector<Polygon> polys;
  for (Vector v: vecs) {
    // Make a little square polygon around each block of side 0.02
    v = v / 3 + Vector(0.5, 0.5);
    polys.push_back(Polygon({v + Vector(-0.01, -0.01), v + Vector(0.01, -0.01), v + Vector(0.01, 0.01), v + Vector(-0.01, 0.01)}));
  }
  return polys;
}

std::vector<Vector> tutti(Graph g, int nb_iter = 3) {
  double s = (g.exterior_vertices[g.exterior_vertices.size() - 1] - g.exterior_vertices[0]).norm(); // Boundary length
  for (int i = 0; i < g.exterior_vertices.size() - 1; i++) {
    s += (g.exterior_vertices[i] - g.exterior_vertices[i + 1]).norm();
  }
  double cs = 0;
  std::vector<Vector> v_curr(g.exterior_vertices.size() + g.interior_vertices.size());
  // Layout external vertices
  for (int i = 0; i < g.exterior_vertices.size(); i++) {
    double angle = 2 * 3.141592653589793236 * cs / s;
    v_curr[i] = Vector(cos(angle), sin(angle));
    if (i < g.exterior_vertices.size() - 1) {
      cs += (g.exterior_vertices[i] - g.exterior_vertices[i + 1]).norm();
    }
  }
  // Layout internal vertices (simple projection)
  int offset = g.exterior_vertices.size();
  for (int i = 0; i < g.interior_vertices.size(); i++) {
    v_curr[i + offset] = Vector(g.interior_vertices[i][0], g.interior_vertices[i][1]);
  }

  // Iterate on internal vertices
  std::vector<Vector> v_next(g.exterior_vertices.size() + g.interior_vertices.size());
  for (int it = 0; it < nb_iter; it++) {
    for (int i = 0; i < g.interior_vertices.size(); i++) {
      Vector sum_linked = Vector(0, 0);
      for (size_t neighb_index : g.interior_edges[i]) {
        sum_linked = sum_linked + v_curr[neighb_index];
      }
      v_next[i + offset] = sum_linked / g.interior_edges[i].size();
    }
    // Save in file"
    std::cout << "Saving frame " << it << std::endl;
    save_frame(Polygonize(v_curr), "tutte", g.exterior_vertices.size(), it);
    // Copy new internal position to current position
    for (int i = 0; i < g.interior_vertices.size(); i++) {
      v_curr[i + offset] = v_next[i + offset];
    }
  }
  return v_curr;
}

int main() {
  tutti(generate_mesh());
  // Generates a sort of pyramid with cut corner and tutti maps it
  // Since we only have one interior vertex (displayed non-blue), it converges in one iteration, as seen in the saved frames
  return 0;
}
