//--------------------_GUI_--------------------------

void buttonEvent(String event) {
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

  void mousePressed ()
  {
    buttonEvent(event);
  }

  void setEvent(String newEvent) {
    event = newEvent;
  }

  void mouseEntered () {
    mouseFree = false;
  }

  void mouseExited () {
    mouseFree = true;
  }


  void draw ()
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

  void update() {
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

  float constrain(float val, float minv, float maxv) {
    return min(max(val, minv), maxv);
  }

  boolean overEvent() {
    if (mouseX > xpos && mouseX < xpos+swidth &&
      mouseY > ypos && mouseY < ypos+sheight) {
      return true;
    } else {
      return false;
    }
  }

  void display() {
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
    text(int(spos), swidth + xpos + 25, ypos + 12);
  }

  float getPos() {
    // Convert spos to be values between
    // 0 and the total width of the scrollbar
    return spos * ratio;
  }
}