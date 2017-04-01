import java.util.List;
import java.util.Arrays;
import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import de.bezier.guido.*;
import java.util.concurrent.*;
import nervoussystem.obj.*;


boolean record = false;

float jitter;

HScrollbar hs1, hs2, hs3, hs4, hs5, hs6;

int mode = 1;

// variables describing figure
List<PVector> vertexes;
int[][] triList;

int triNum;
float zoom=1.0;
float e_sum=20.45;
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

  size(1000, 800, P3D);
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
  noSmooth();
}