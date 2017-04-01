void createNewShape() {
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

void mapTriangles(Delaunay d) {
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

void randomizeVertexes(float[][] vertexAngles) {
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

void changeFigureType() {
  // number of vertexes:
  vertexNum = int(random(hs1.spos*2, hs2.spos*2));
  //vertexNum = int(300);

  // how far points can be from original sphere
  range = random(hs3.spos, hs4.spos);
}