void keyPressed() {
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


void mousePressed() {
  if (mouseFree) {
    startX = rotationDisplayment(mouseX) - mx;
    startY = rotationDisplayment(mouseY) - my;
  }
}

float rotationDisplayment(int mouseCoord) {
  return ((float)mouseCoord-(width/2.0))/(width/2.0)*4;
}

void mouseWheel(MouseEvent event) {
  // for zooming the figure
  float e = event.getCount();

  e_sum+=e;
  if (e_sum>100) {
    e_sum=100;
  };
  if (e_sum<0) {
    e_sum=0;
  }
  zoom = map(e_sum, 0, 100, 0.1, 4.5);
}