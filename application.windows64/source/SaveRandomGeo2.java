import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import java.util.List; 
import java.util.Arrays; 
import java.util.Date; 
import java.text.DateFormat; 
import java.text.SimpleDateFormat; 
import de.bezier.guido.*; 
import java.util.concurrent.*; 
import nervoussystem.obj.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class SaveRandomGeo2 extends PApplet {











boolean record = false;

float jitter;

HScrollbar hs1, hs2, hs3, hs4, hs5, hs6;

int mode = 1;

// variables describing figure
List<PVector> vertexes;
int[][] triList;

int triNum;
float zoom=1.0f;
float e_sum=20.45f;
float mx=0, my=2;
float startX=0, startY=0;
int vertexNum;
float range;
float averageRadius;

SimpleButton newShapeButton;
SimpleButton saveShapeButton;
SimpleButton exitShapeButton;
SimpleButton modeButton;
SimpleButton PDFButton;

boolean mouseFree = true;


public void setup() {

  Interactive.make( this );

  
  // fullScreen(P3D);

  int x = 0;
  int y = 80;

  newShapeButton = new SimpleButton( 20 + x, y-60, 76, 20, "New Shape" );
  newShapeButton.setEvent("newShape");

  modeButton = new SimpleButton( 105 + x, y-60, 91, 20, "Change Mode" );
  modeButton.setEvent("mode");

  saveShapeButton = new SimpleButton( 204 + x, y-60, 82, 20, "Save Object" );
  saveShapeButton.setEvent("saveShape");

  exitShapeButton = new SimpleButton( width - x-75, y-60, 41, 20, "Exit" );
  exitShapeButton.setEvent("exit");

  hs1 = new HScrollbar(20, 85, 180, 16, 16, "Minimum Points");
  hs2 = new HScrollbar(20, 135, 180, 16, 16, "Maximum Points");
  hs3 = new HScrollbar(20, 185, 77, 16, 16, "Minimum Distance");
  hs4 = new HScrollbar(20, 235, 77, 16, 16, "Maximum Distance");
  hs5 = new HScrollbar(20, 285, 180, 16, 16, "Circularity");
  hs6 = new HScrollbar(20, 335, 180, 16, 16, "Jitter");
  
  createNewShape();
  
}
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
      0.5f * (v2.x*v2.x - v1.x*v1.x + v2.y*v2.y - v1.y*v1.y + v2.z*v2.z - v1.z*v1.z), 
      0.5f * (v3.x*v3.x - v1.x*v1.x + v3.y*v3.y - v1.y*v1.y + v3.z*v3.z - v1.z*v1.z), 
      0.5f * (v4.x*v4.x - v1.x*v1.x + v4.y*v4.y - v1.y*v1.y + v4.z*v4.z - v1.z*v1.z)
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

  /** LU\u5206\u89e3\u306b\u3088\u308b\u65b9\u7a0b\u5f0f\u306e\u89e3\u6cd5 **/
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
public void draw() {
  jitter = hs6.spos/20;
  // jitter = 0;
  if(jitter <= 1.1f){
    jitter = 0;
  }
  
  // setup the space
  background(75);

  //println(hs1.getPos() + " : " + hs2.getPos());

  if (record && mode != 4) {
    beginRecord("nervoussystem.obj.OBJExport", "obj/filename-###.obj");
  } 

  // save state of the matrix
  pushMatrix();

  translate(width/2+70, height/2, width/3);


  //gives a float between -1 and 1 depending on mouseX
  float ax= rotationDisplayment(mouseX);

  //as above, but for mouseY
  float ay= rotationDisplayment(mouseY);

  if (mousePressed== true && mouseX > 350) {
    my = ay - startY;
    mx = ax - startX;
  }

  rotateX(my);
  rotateY(mx);

  //scale(zoom);

  // add black frame around each triangle
  stroke(0);

  // fill each triangle with white color
  fill(180);

  if (mode==1) {
    beginShape(TRIANGLES);
  } else if (mode==2) {
    beginShape(TRIANGLE_STRIP);
  } else if (mode==3) {
    beginShape(TRIANGLE_FAN);
  } else if (mode==4) {
    beginShape(LINES);
  } else if (mode==5) {
    beginShape(QUADS);
  } else if (mode==6) {
    beginShape(QUAD_STRIP);
  } else {
    beginShape(TRIANGLES);
  }


  for (int[] tri : triList) {
    PVector v1 = vertexes.get(tri[0]);
    PVector v2 = vertexes.get(tri[1]);
    PVector v3 = vertexes.get(tri[2]);

    //strokeWeight(1.5);

    fill(180);
    vertex(v1.x, v1.y, v1.z + random(jitter));
    fill(170);
    vertex(v2.x+ random(jitter), v2.y, v2.z);
    fill(180);
    vertex(v3.x, v3.y+ random(jitter), v3.z);

  }
  endShape();

  // return the last saved state of the matrix (this is actually not necessary in this case)
  popMatrix();

  if (record) {
    endRecord();
    record = false;
  } 


  // scrollbars
  hs1.update();
  hs2.update();
  hs3.update();
  hs4.update();
  hs5.update();
  hs6.update();

  hs1.display();
  hs2.display();
  hs3.display();
  hs4.display();
  hs5.display();
  hs6.display();

  stroke(1.5f);
  // line(275,0,275,height);
}
//--------------------_GUI_--------------------------

public void buttonEvent(String event) {
  if (event == "newShape") {
    createNewShape();
  } else if (event == "saveShape") {
    record = true;
  } else if (event == "mode") {
    if (mode<6) {
      mode++;
    } else {
      mode = 1;
    }
  } else if (event == "exit") {
    exit();
  } 
}


public class SimpleButton {
  float x, y, width, height;
  boolean on;
  String button_text;
  String event = "none";

  SimpleButton ( float xx, float yy, float w, float h, String txt)
  {
    x = xx; 
    y = yy; 
    width = w; 
    height = h;
    button_text = txt;

    Interactive.add( this ); // register it with the manager
  }

  // called by manager

  public void mousePressed ()
  {
    buttonEvent(event);
  }

  public void setEvent(String newEvent) {
    event = newEvent;
  }

  public void mouseEntered () {
    mouseFree = false;
  }

  public void mouseExited () {
    mouseFree = true;
  }


  public void draw ()
  {
    //        fill( 100 );
    noFill();
    rect(x, y, width, height);
    fill( 255 );
    textAlign(CENTER, CENTER);
    text(button_text, x + width/2, y+height/2-2);
  }
}




class HScrollbar {
  int swidth, sheight;    // width and height of bar
  float xpos, ypos;       // x and y position of bar
  float spos, newspos;    // x position of slider
  float sposMin, sposMax; // max and min values of slider
  int loose;              // how loose/heavy
  boolean over;           // is the mouse over the slider?
  boolean locked;
  float ratio;
  String stitle;

  HScrollbar (float xp, float yp, int sw, int sh, int l, String newTitle) {
    swidth = sw;
    sheight = sh;
    int widthtoheight = sw - sh;
    ratio = (float)sw / (float)widthtoheight;
    xpos = xp;
    ypos = yp-sheight/2;
    spos = xpos + swidth/2 - sheight/2;
    newspos = spos;
    sposMin = xpos;
    sposMax = xpos + swidth - sheight;
    loose = l;
    stitle = newTitle;
  }

  public void update() {
    if (overEvent()) {
      over = true;
    } else {
      over = false;
    }

    if (mousePressed && over && mouseFree == true) {
      locked = true;
      mouseFree = false;
    }
    if (!mousePressed && mouseFree == false) {
      locked = false;
      mouseFree = true;
    }
    if (locked && mouseFree == false) {
      newspos = constrain(mouseX-sheight/2, sposMin, sposMax);
      mouseFree = true;
    }
    if (abs(newspos - spos) > 1) {
      spos = spos + (newspos-spos)/loose;
    }
  }

  public float constrain(float val, float minv, float maxv) {
    return min(max(val, minv), maxv);
  }

  public boolean overEvent() {
    if (mouseX > xpos && mouseX < xpos+swidth &&
      mouseY > ypos && mouseY < ypos+sheight) {
      return true;
    } else {
      return false;
    }
  }

  public void display() {
    noStroke();
    fill(204);
    rect(xpos, ypos, swidth, sheight);
    if (over || locked && mouseFree == false) {
      fill(0, 0, 0);
    } else {
      fill(102, 102, 102);
    }
    rect(spos, ypos, sheight, sheight);
    fill(255);
    textAlign(LEFT);
    text(stitle, xpos, ypos - 10);
    textAlign(CENTER);
    text(PApplet.parseInt(spos), swidth + xpos + 25, ypos + 12);
  }

  public float getPos() {
    // Convert spos to be values between
    // 0 and the total width of the scrollbar
    return spos * ratio;
  }
}
public void keyPressed() {
  if (key == 'c' || key == 'C') {
    createNewShape();
  } else if (key == 'r' || key == 'R') {
    record = true;
  } else if (key == 'q' || key == 'Q') {
    exit();
  } else if (key == 'f' || key == 'F') {
    if (mode<6) {
      mode++;
    } else {
      mode = 1;
    }
  }
}

//**********************************************************************


public void mousePressed() {
  if (mouseFree) {
    startX = rotationDisplayment(mouseX) - mx;
    startY = rotationDisplayment(mouseY) - my;
  }
}

public float rotationDisplayment(int mouseCoord) {
  return ((float)mouseCoord-(width/2.0f))/(width/2.0f)*4;
}

public void mouseWheel(MouseEvent event) {
  // for zooming the figure
  float e = event.getCount();

  e_sum+=e;
  if (e_sum>100) {
    e_sum=100;
  };
  if (e_sum<0) {
    e_sum=0;
  }
  zoom = map(e_sum, 0, 100, 0.1f, 4.5f);
}
public void createNewShape() {
  // I decided to use Delaunay triangularization for creating triangles from cloud of random points on a sphere - this ensures us that resulting figure makes sense
  zoom = 1;

  changeFigureType();
  Delaunay d = new Delaunay();
  //background(20);
  vertexes = new ArrayList<PVector>();
  float[][] vertexAngles = new float[vertexNum][2];

  
  // here are created random points. They are laying on sphere.
  for (int i = 0; i < vertexNum; i++) {
    // we add some small random number to radius in order to avoid computational bugs
    float r = 90 + random(0.01f);
    float phi = radians(random(-90, 90));
    float theta = radians(random(hs5.spos*2));
    vertexAngles[i][0] = phi;
    vertexAngles[i][1] = theta;
    
    //turn into sliders
      //vertexAngles[vAngles][0] = 179;
    
    vertexes.add(new PVector(r*cos(phi)*cos(theta), r*sin(phi), r*cos(phi)*sin(theta)));
  }

  // println("start triangularization");
  // here Delaunay triangularization take place
  d.performTrinagularization(vertexes);
  // println("triangularization finished");

  // after creating list of points we map achieved triangles to saved points so each triangle vertex is saved by reference to each vertex rather than absolute position
  // this enables us quick movement of vertex without changes in triangles lists
  mapTriangles(d);
  // than we randomize vertexes changing only the radius of each vertex. That way it's sure that resulting figure makes sense.
  randomizeVertexes(vertexAngles);
}

public void mapTriangles(Delaunay d) {
  // map triangles to vertexes
  triNum = d.triangles.size();
  // println("triangles: " + triNum);
  triList = new int[triNum][3];
  // for each triangle search for each vertex in vertex list and save reference to this vertex
  for (int i = 0; i < triNum; ++i) {
    PVector v1 = d.triangles.get(i).v1;
    PVector v2 = d.triangles.get(i).v2;
    PVector v3 = d.triangles.get(i).v3;
    for (int n = 0; n < vertexNum; ++n) {
      if (v1==vertexes.get(n)) {
        // println("success: tri=" + i + " ver:" + n );
        triList[i][0] = n;
      } else if (v2==vertexes.get(n)) {
        // println("success: tri=" + i + " ver:" + n );
        triList[i][1] = n;
      } else if (v3==vertexes.get(n)) {
        // println("success: tri=" + i + " ver:" + n );
        triList[i][2] = n;
      }
    }
  }
}

public void randomizeVertexes(float[][] vertexAngles) {
  // println("range: " + range);
  float sumRadius = 0;
  for (int i = 0; i < vertexNum; ++i) {
    // for each vertex we don't change it's angles only the radius is changing
    float phi = vertexAngles[i][0];
    float theta = vertexAngles[i][1];
    
    
    
    
    // the range variable sets the range for whole figure so it will be more like sphere or hedgehog
      //turn to sliders
    float r = random(70, 70+range);
    
    
    
    vertexes.add(i, new PVector(r*cos(phi)*cos(theta), r*sin(phi), r*cos(phi)*sin(theta)));
    sumRadius +=r;
  }
  averageRadius = sumRadius/vertexNum;
}

public void changeFigureType() {
  // number of vertexes:
  vertexNum = PApplet.parseInt(random(hs1.spos*2, hs2.spos*2));
  //vertexNum = int(300);

  // how far points can be from original sphere
  range = random(hs3.spos, hs4.spos);
}
  public void settings() {  size(1000, 800, P3D);  noSmooth(); }
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "SaveRandomGeo2" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
