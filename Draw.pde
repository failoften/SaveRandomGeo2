public void draw() {
  jitter = hs6.spos/20;
  // jitter = 0;
  if(jitter <= 1.1){
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

  stroke(1.5);
  // line(275,0,275,height);
}