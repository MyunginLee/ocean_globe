#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

#include "al/app/al_App.hpp"
#include "al/graphics/al_Image.hpp"
#include "al/math/al_Random.hpp"

using namespace al;
using namespace std;

class MyApp : public App {
 public:
  Mesh pic, rgb, hsv, extra, target;
  Vec3f img_pos, rgb_pos, hsv_pos, extra_pos;
  int key, prevkey; 
  
  void onCreate() override {
    const char *filename = "../../texture/sst/sst_2003.jpeg";
    auto imageData = Image(filename);
    if (imageData.array().size() == 0) {
      cout << "failed to load image" << endl;
      exit(1);
    }
    cout << "loaded image size: " << imageData.width() << ", "
         << imageData.height() << endl;

    int W = imageData.width();
    int H = imageData.height();
    cout << "sibal" << endl;    
    pic.primitive(Mesh::POINTS);
    rgb.primitive(Mesh::POINTS);
    hsv.primitive(Mesh::POINTS);
    extra.primitive(Mesh::LINES);


    // iterate through all the pixel, scanning each row
    for (int row = 0; row < H; row++) {
      for (int column = 0; column < W; column++) {
        auto pixel = imageData.at(column, H - row - 1);
        
        img_pos = {1.0 * column / W - 0.5, 1.0 * row / H - 0.5, 0.0}; // position pixels in the middle of the space by writing -0.5
        pic.vertex(img_pos);
        pic.color(pixel.r / 255.0, pixel.g / 255.0, pixel.b / 255.0);


////////////RGB/////////////////////

    rgb_pos = {pixel.r /255.0, 
               pixel.g /255.0, 
               pixel.b /255.0};

    Color color(pixel.r/255.0, pixel.g/255.0, pixel.b/255.0);

    rgb.vertex(rgb_pos - 0.5); // place in the middle
    rgb.color(color);


/////////////HSV///////////////////
// S × cos(H*2*pi) , S × sin(H*2*pi), V

        HSV c = color;
        hsv_pos = {c.s * 0.5 * cos(c.h * M_2PI), 
                   c.v - 0.5,                           // place z in the middle
                   c.s * 0.5 * sin(c.h * M_2PI)
                  };

        hsv.vertex(hsv_pos);
        hsv.color(color);
    
      }
    }
    // set the camera position back some (z=3) and center on (x=0.5, y=0.5)
    nav().pos(0, 0, 5);

    cout << "sibal" << endl;    

    extra = pic;
    target = pic;
  }


  
  double t = 0.0;
  double phase = 0;

  void onAnimate(double dt) override {
    t += dt;
    phase += dt;

    
    if (key == 4){
      for (int i = 0; i < extra.vertices().size(); i++) {
          target.vertices()[i].lerp(extra.vertices()[i]*((t + dt)* 0.5), 0.1);
      }
    }else{
      for (int i = 0; i < extra.vertices().size(); i++) {
          target.vertices()[i].lerp(extra.vertices()[i], 0.1);
      }
    }
    // cout << extra.vertices().size() << endl;

    
  }

  bool onKeyDown(const Keyboard &k) override {
   
    switch (k.key()) {
      case '1':
        prevkey = key;
        key = 1;
        extra = pic;
        if (prevkey == 1){
          target = pic;  
        } else if (prevkey == 2){
          target = rgb;
        } else if (prevkey == 3){
          target = hsv;
        }
        break;

      case '2':
      // previous = next;
      // next = rgb;

        prevkey = key;
        key = 2;
        extra = rgb;
        if (prevkey == 1){
          target = pic;  
        } else if (prevkey == 2){
          target = rgb;
        } else if (prevkey == 3){
          target = hsv;
        }
        break;

      case '3':
        prevkey = key;
        key = 3;
        extra = hsv;
        if (prevkey == 1){
          target = pic;  
        } else if (prevkey == 2){
          target = rgb;
        } else if (prevkey == 3){
          target = hsv;
        }
        break;

      case '4':
        key = 4;
        if (prevkey == 1){
          target = pic;  
        } else if (prevkey == 2){
          target = rgb;
        } else if (prevkey == 3){
          target = hsv;
        }
        extra = target;
        break;    

      default:
        break;
    }
    return true;
  }

  void onDraw(Graphics &g) override {
    g.clear(0.0f);
    g.pushMatrix();
    g.meshColor();
    g.pointSize(10);
    g.rotate(phase * 3.5, 0.0, 1.0, 0.0);
    g.draw(target);
    g.popMatrix();
    
  }
};

int main() {
  MyApp app;
  app.configureAudio(48000, 512, 2, 0);
  app.start();
}