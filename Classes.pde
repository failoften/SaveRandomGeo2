public class Delaunay {

  List<PVector> vertices;    
  ArrayList<Tetrahedron> tetras;  

  public List<Line> edges;

  public List<Line> surfaceEdges;
  public List<Triangle> triangles;


  public Delaunay() {
    vertices = new ArrayList<PVector>();
    tetras = new ArrayList<Tetrahedron>();
    edges = new ArrayList<Line>();
    surfaceEdges = new ArrayList<Line>();
    triangles = new ArrayList<Triangle>();
  }

  public void performTrinagularization(List<PVector> seq) {

    tetras.clear();
    edges.clear();



    PVector vMax = new PVector(-999, -999, -999);
    PVector vMin = new PVector( 999, 999, 999);
    for (PVector v : seq) {
      if (vMax.x < v.x) vMax.x = v.x;
      if (vMax.y < v.y) vMax.y = v.y;
      if (vMax.z < v.z) vMax.z = v.z;
      if (vMin.x > v.x) vMin.x = v.x;
      if (vMin.y > v.y) vMin.y = v.y;
      if (vMin.z > v.z) vMin.z = v.z;
    }

    PVector center = new PVector();    
    center.x = 0.5f * (vMax.x - vMin.x);
    center.y = 0.5f * (vMax.y - vMin.y);
    center.z = 0.5f * (vMax.z - vMin.z);
    float r = -1;                      
    for (PVector v : seq) {
      if (r < PVector.dist(center, v)) r = PVector.dist(center, v);
    }
    r += 0.1f;                         


    PVector v1 = new PVector();
    v1.x = center.x;
    v1.y = center.y + 3.0f*r;
    v1.z = center.z;

    PVector v2 = new PVector();
    v2.x = center.x - 2.0f*(float)Math.sqrt(2)*r;
    v2.y = center.y - r;
    v2.z = center.z;

    PVector v3 = new PVector();
    v3.x = center.x + (float)Math.sqrt(2)*r;
    v3.y = center.y - r;
    v3.z = center.z + (float)Math.sqrt(6)*r;

    PVector v4 = new PVector();
    v4.x = center.x + (float)Math.sqrt(2)*r;
    v4.y = center.y - r;
    v4.z = center.z - (float)Math.sqrt(6)*r;

    PVector[] outer = {v1, v2, v3, v4};
    tetras.add(new Tetrahedron(v1, v2, v3, v4));


    ArrayList<Tetrahedron> tmpTList = new ArrayList<Tetrahedron>();
    ArrayList<Tetrahedron> newTList = new ArrayList<Tetrahedron>();
    ArrayList<Tetrahedron> removeTList = new ArrayList<Tetrahedron>();
    for (PVector v : seq) {
      tmpTList.clear();
      newTList.clear();
      removeTList.clear();
      for (Tetrahedron t : tetras) {
        if ((t.o != null) && (t.r > PVector.dist(v, t.o))) {
          tmpTList.add(t);
        }
      }

      for (Tetrahedron t1 : tmpTList) {

        tetras.remove(t1);

        v1 = t1.vertices[0];
        v2 = t1.vertices[1];
        v3 = t1.vertices[2];
        v4 = t1.vertices[3];
        newTList.add(new Tetrahedron(v1, v2, v3, v));
        newTList.add(new Tetrahedron(v1, v2, v4, v));
        newTList.add(new Tetrahedron(v1, v3, v4, v));
        newTList.add(new Tetrahedron(v2, v3, v4, v));
      }

      boolean[] isRedundancy = new boolean[newTList.size()];
      for (int i = 0; i < isRedundancy.length; i++) isRedundancy[i] = false;
      for (int i = 0; i < newTList.size()-1; i++) {
        for (int j = i+1; j < newTList.size(); j++) {
          if (newTList.get(i).equals(newTList.get(j))) {
            isRedundancy[i] = isRedundancy[j] = true;
          }
        }
      }
      for (int i = 0; i < isRedundancy.length; i++) {
        if (!isRedundancy[i]) {
          tetras.add(newTList.get(i));
        }
      }
    }


    boolean isOuter = false;
    ArrayList<Tetrahedron> tetrasClone = new ArrayList<Tetrahedron>(tetras);
    for (Tetrahedron t4 : tetrasClone) {
      isOuter = false;
      for (PVector p1 : t4.vertices) {
        for (PVector p2 : outer) {
          if (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z) {
            isOuter = true;
          }
        }
      }
      if (isOuter) {
        tetras.remove(t4);
      }
    }

    triangles.clear();
    boolean isSame = false;
    for (Tetrahedron t : tetras) {
      for (Line l1 : t.getLines()) {
        isSame = false;
        for (Line l2 : edges) {
          if (l2.equals(l1)) {
            isSame = true;
            break;
          }
        }
        if (!isSame) {
          edges.add(l1);
        }
      }
    }




    ArrayList<Triangle> triList = new ArrayList<Triangle>();
    for (Tetrahedron t : tetras) {
      v1 = t.vertices[0];
      v2 = t.vertices[1];
      v3 = t.vertices[2];
      v4 = t.vertices[3];

      Triangle tri1 = new Triangle(v1, v2, v3);
      Triangle tri2 = new Triangle(v1, v3, v4);
      Triangle tri3 = new Triangle(v1, v4, v2);
      Triangle tri4 = new Triangle(v4, v3, v2);

      PVector n;

      n = tri1.getNormal();
      if (n.dot(v1) > n.dot(v4)) tri1.turnBack();

      n = tri2.getNormal();
      if (n.dot(v1) > n.dot(v2)) tri2.turnBack();

      n = tri3.getNormal();
      if (n.dot(v1) > n.dot(v3)) tri3.turnBack();

      n = tri4.getNormal();
      if (n.dot(v2) > n.dot(v1)) tri4.turnBack();

      triList.add(tri1);
      triList.add(tri2);
      triList.add(tri3);
      triList.add(tri4);
    }
    boolean[] isSameTriangle = new boolean[triList.size()];
    for (int i = 0; i < triList.size()-1; i++) {
      for (int j = i+1; j < triList.size(); j++) {

        if (triList.get(i).equals(triList.get(j))) isSameTriangle[i] = isSameTriangle[j] = true;
      }
    }
    for (int i = 0; i < isSameTriangle.length; i++) {
      if (!isSameTriangle[i]) triangles.add(triList.get(i));
    }

    surfaceEdges.clear();
    ArrayList<Line> surfaceEdgeList = new ArrayList<Line>();
    for (Triangle tri : triangles) {
      for (Line l : tri.getLines()) {
        surfaceEdgeList.add(l);
      }
      // surfaceEdgeList.addAll(Arrays.asList(tri.getLines()));
    }
    boolean[] isRedundancy = new boolean[surfaceEdgeList.size()];
    for (int i = 0; i < surfaceEdgeList.size()-1; i++) {
      for (int j = i+1; j < surfaceEdgeList.size(); j++) {
        if (surfaceEdgeList.get(i).equals(surfaceEdgeList.get(j))) isRedundancy[j] = true;
      }
    }

    for (int i = 0; i < isRedundancy.length; i++) {
      if (!isRedundancy[i]) surfaceEdges.add(surfaceEdgeList.get(i));
    }
  }
}

//**********************************************************************


class Line {
  public PVector start, end;
  public Line(PVector start, PVector end) {
    this.start = start;
    this.end = end;
  }


  public void reverse() {
    PVector tmp = this.start;
    this.start = this.end;
    this.end = tmp;
  }


  public boolean equals(Line l) {
    if ((this.start == l.start && this.end == l.end)
      || (this.start == l.end && this.end == l.start))
      return true;
    return false;
  }
}


//**********************************************************************


class Tetrahedron {

  PVector[] vertices;
  PVector o;     
  float   r;     

  public Tetrahedron(PVector[] v) {
    this.vertices = v;
    getCenterCircumcircle();
  }

  public Tetrahedron(PVector v1, PVector v2, PVector v3, PVector v4) {
    this.vertices = new PVector[4];
    vertices[0] = v1;
    vertices[1] = v2;
    vertices[2] = v3;
    vertices[3] = v4;
    getCenterCircumcircle();
  }

  public boolean equals(Tetrahedron t) {
    int count = 0;
    for (PVector p1 : this.vertices) {
      for (PVector p2 : t.vertices) {
        if (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z) {
          count++;
        }
      }
    }
    if (count == 4) return true;
    return false;
  }

  public Line[] getLines() {
    PVector v1 = vertices[0];
    PVector v2 = vertices[1];
    PVector v3 = vertices[2];
    PVector v4 = vertices[3];

    Line[] lines = new Line[6];

    lines[0] = new Line(v1, v2);
    lines[1] = new Line(v1, v3);
    lines[2] = new Line(v1, v4);
    lines[3] = new Line(v2, v3);
    lines[4] = new Line(v2, v4);
    lines[5] = new Line(v3, v4);
    return lines;
  }


  private void getCenterCircumcircle() {
    PVector v1 = vertices[0];
    PVector v2 = vertices[1];
    PVector v3 = vertices[2];
    PVector v4 = vertices[3];

    double[][] A = {
      {v2.x - v1.x, v2.y-v1.y, v2.z-v1.z}, 
      {v3.x - v1.x, v3.y-v1.y, v3.z-v1.z}, 
      {v4.x - v1.x, v4.y-v1.y, v4.z-v1.z}
    };
    double[] b = {
      0.5 * (v2.x*v2.x - v1.x*v1.x + v2.y*v2.y - v1.y*v1.y + v2.z*v2.z - v1.z*v1.z), 
      0.5 * (v3.x*v3.x - v1.x*v1.x + v3.y*v3.y - v1.y*v1.y + v3.z*v3.z - v1.z*v1.z), 
      0.5 * (v4.x*v4.x - v1.x*v1.x + v4.y*v4.y - v1.y*v1.y + v4.z*v4.z - v1.z*v1.z)
    };
    double[] x = new double[3];
    if (gauss(A, b, x) == 0) {
      o = null;
      r = -1;
    } else {
      o = new PVector((float)x[0], (float)x[1], (float)x[2]);
      r = PVector.dist(o, v1);
    }
  }

  /** LU分解による方程式の解法 **/
  private double lu(double[][] a, int[] ip) {
    int n = a.length;
    double[] weight = new double[n];

    for (int k = 0; k < n; k++) {
      ip[k] = k;
      double u = 0;
      for (int j = 0; j < n; j++) {
        double t = Math.abs(a[k][j]);
        if (t > u) u = t;
      }
      if (u == 0) return 0;
      weight[k] = 1/u;
    }
    double det = 1;
    for (int k = 0; k < n; k++) {
      double u = -1;
      int m = 0;
      for (int i = k; i < n; i++) {
        int ii = ip[i];
        double t = Math.abs(a[ii][k]) * weight[ii];
        if (t>u) { 
          u = t; 
          m = i;
        }
      }
      int ik = ip[m];
      if (m != k) {
        ip[m] = ip[k]; 
        ip[k] = ik;
        det = -det;
      }
      u = a[ik][k]; 
      det *= u;
      if (u == 0) return 0;
      for (int i = k+1; i < n; i++) {
        int ii = ip[i]; 
        double t = (a[ii][k] /= u);
        for (int j = k+1; j < n; j++) a[ii][j] -= t * a[ik][j];
      }
    }
    return det;
  }
  private void solve(double[][] a, double[] b, int[] ip, double[] x) {
    int n = a.length;
    for (int i = 0; i < n; i++) {
      int ii = ip[i]; 
      double t = b[ii];
      for (int j = 0; j < i; j++) t -= a[ii][j] * x[j];
      x[i] = t;
    }
    for (int i = n-1; i >= 0; i--) {
      double t = x[i]; 
      int ii = ip[i];
      for (int j = i+1; j < n; j++) t -= a[ii][j] * x[j];
      x[i] = t / a[ii][i];
    }
  }
  private double gauss(double[][] a, double[] b, double[] x) {
    int n = a.length;
    int[] ip = new int[n];
    double det = lu(a, ip);

    if (det != 0) { 
      solve(a, b, ip, x);
    }
    return det;
  }
}

//**********************************************************************


class Triangle {
  private PVector v1, v2, v3;
  public PVector[] points;
  public Triangle(PVector v1, PVector v2, PVector v3) {
    this.v1 = v1;
    this.v2 = v2;
    this.v3 = v3;
    this.points = new PVector[3];
    this.points[0] = v1;
    this.points[1] = v2;
    this.points[2] = v3;
  }



  public PVector getNormal() {
    PVector edge1 = new PVector(v2.x-v1.x, v2.y-v1.y, v2.z-v1.z);
    PVector edge2 = new PVector(v3.x-v1.x, v3.y-v1.y, v3.z-v1.z);


    PVector normal = edge1.cross(edge2);
    normal.normalize();
    return normal;
  }


  public void turnBack() {
    PVector tmp = this.v3;
    this.v3 = this.v1;
    this.v1 = tmp;
  }


  public Line[] getLines() {
    Line[] l = {
      new Line(v1, v2), 
      new Line(v2, v3), 
      new Line(v3, v1)
    };
    return l;
  }


  public boolean equals(Triangle t) {
    PVector[] points1 = this.points;
    PVector[] points2 = t.points;

    int cnt = 0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (points1[i] == points2[j])
          cnt++;
      }
    }
    if (cnt == 3) return true;
    else return false;
  }
}